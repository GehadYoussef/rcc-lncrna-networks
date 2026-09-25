# 02_import_dataset.R: dataset adapters and schema audit.
# Each adapter reads the deposited raw count files and returns one
# SingleCellExperiment: assay "counts" (raw integer UMI), rowData (gene id,
# symbol) and colData (cell_id = dataset|sample|barcode, donor, tissue, QC
# fields, TISCH2 and author labels). TISCH2 log-normalised expression is not
# stored. Outcome and severity columns are removed at import and logged.
# A dataset that fails a schema check is recorded with the reason and skipped.
# Outputs: data/derived/singlecell/<dataset>.sce.rds, results/singlecell/
#          02_schema_audit.tsv (per dataset), 02_schema_audit_<dataset>.tsv
#          (per sample), 02_label_join_audit.tsv, 02_removed_outcome_fields.tsv

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
suppressPackageStartupMessages(library(SingleCellExperiment))
banner("02  import and schema audit")
start_log("02_import_dataset")

m <- read_manifest(); man <- m$manifest; smp <- m$samples
acq <- read_tsv("01_acquisition_log.tsv")
ONLY <- sub("^--dataset=", "", grep("^--dataset=", commandArgs(TRUE), value = TRUE))

norm_bc <- function(x) sub("-[0-9]+$", "", x)
raw_path <- function(ds, sub, fn) file.path(RAW_DIR, ds, sub, fn)

# ---- per-dataset label join rules --------------------------------------------
# How a TISCH2 "Cell" string maps to (sample, barcode) in the raw files. The
# manifest notes record per-table quirks (e.g. the GSE159115 Sample/Tissue swap).
tisch_key <- function(ds, tt, samples) {
  cell <- tt$Cell
  if (ds %in% c("KIRC_GSE159115", "KIRC_GSE171306")) {
    gsm <- sub("@.*", "", cell); bc <- sub(".*@", "", cell)
    sid <- samples$sample_id[match(gsm, samples$gsm)]
  } else if (ds %in% c("KIRC_GSE111360", "KIRC_GSE139555")) {
    sid <- sub("@.*", "", cell); bc <- sub(".*@", "", cell)
  } else if (ds == "KIRC_GSE121636") {
    sid <- samples$sample_id[match(tt$Sample, samples$gsm)]; bc <- sub(".*@", "", cell)
  } else if (ds == "KIRC_GSE145281_aPDL1") {
    sid <- tt$Sample; bc <- sub(".*_", "", cell)
  } else {
    # default: "<sample or GSM>@<barcode>"
    pre <- sub("@.*", "", cell); bc <- sub(".*@", "", cell)
    sid <- ifelse(pre %in% samples$gsm, samples$sample_id[match(pre, samples$gsm)], pre)
    if (!all(grepl("@", cell))) stop("no TISCH2 join rule for ", ds)
  }
  paste(sid, norm_bc(bc), sep = "|")
}

author_labels <- function(ds, samples) {
  if (ds == "KIRC_GSE159115") {
    fs <- raw_path(ds, "author", paste0("GSE159115_", c("ccRCC", "normal", "chRCC"), "_anno.csv.gz"))
    a <- rbindlist(lapply(fs[file.exists(fs)], fread), fill = TRUE)
    a[, key := paste(sample, norm_bc(sub("^SI_[0-9]+_", "", cell)), sep = "|")]
    return(a[, .(key, author_label = anno, author_patient = patient, author_doublet = doublet,
                 author_cellfilt = cellfilt, author_pct_mt = pct_MT)])
  }
  if (ds == "KIRC_GSE222703") {
    f <- raw_path(ds, "extracted", "GSE222703_SeuratObj_ALL_integrated.rds")
    if (!file.exists(f)) stop("GSE222703 author Seurat object not extracted: run 01")
    suppressPackageStartupMessages(requireNamespace("SeuratObject"))
    md <- as.data.table(readRDS(f)@meta.data, keep.rownames = "cell")
    # cell names are "<barcode>_<k>", where k indexes the sample. The k -> orig.ident
    # map is derived from the table itself and must be one-to-one.
    md[, k := sub(".*_", "", cell)]
    km <- unique(md[, .(k, orig.ident)])
    if (anyDuplicated(km$k) || anyDuplicated(km$orig.ident))
      stop("GSE222703: cell-name suffix does not map one-to-one to orig.ident")
    ident_to_sample <- c(p022_Tumor = "p022_Tumoral", p027_Tumor = "p027_Tumoral", p029_Tumor = "p029_Tumoral",
                         p022_Juxta = "p022_Juxta", p027_Juxta = "p027_Juxta", p029_Juxta = "p029_Juxta")
    if (!all(md$orig.ident %in% names(ident_to_sample))) stop("GSE222703: unexpected orig.ident values")
    md[, key := paste(ident_to_sample[orig.ident], norm_bc(sub("_[0-9]+$", "", cell)), sep = "|")]
    return(md[, .(key, author_label = Celltype_Harmony, author_patient = Patient, author_tissue = Tissue,
                  author_doublet = ManDblt | doublets_consensus.class)])
  }
  if (ds == "KIRC_GSE139555") {
    f <- raw_path(ds, "author", "GSE139555_all_metadata.txt.gz")
    if (!file.exists(f)) return(NULL)
    a <- fread(f)
    msg("GSE139555 author metadata columns: ", paste(names(a), collapse = ", "))
    return(NULL)   # columns are logged, no author labels are joined
  }
  NULL
}

# ---- readers per sample format ---------------------------------------------------
read_sample <- function(ds, r) {
  dir <- file.path(RAW_DIR, ds, if (grepl("_extracted$", r$format)) "extracted" else "geo")
  x <- switch(r$format,
    "10x_h5"     = read_10x_h5(file.path(dir, paste0(r$file_prefix, ".h5"))),
    "10x_mtx_v2" = read_10x_mtx(dir, r$file_prefix),
    "10x_mtx_v3" = read_10x_mtx(dir, r$file_prefix),
    "10x_mtx_v3_extracted" = read_10x_mtx(dir, r$file_prefix),
    "dense_txt"  = read_dense_txt(file.path(dir, paste0(r$file_prefix, ".txt.gz"))),
    stop("no reader for format ", r$format))
  x
}

qc_fields <- function(counts, symbols) {
  mito <- grepl("^MT-", symbols, ignore.case = TRUE)
  n_umi <- Matrix::colSums(counts)
  data.table(n_umi = n_umi, n_genes = Matrix::colSums(counts > 0),
             pct_mito = if (any(mito)) 100 * Matrix::colSums(counts[mito, , drop = FALSE]) / pmax(n_umi, 1) else NA_real_)
}

import_geo_dataset <- function(ds) {
  S <- smp[dataset_name == ds & include == TRUE]
  if (!nrow(S)) stop("no included samples")
  tt_f <- raw_path(ds, "tisch2", paste0(ds, "_CellMetainfo_table.tsv"))
  tt <- if (file.exists(tt_f)) fread(tt_f) else NULL
  if (!is.null(tt)) {
    tt[, key := tisch_key(ds, tt, smp[dataset_name == ds])]
    setnames(tt, names(tt), paste0("tisch_", make.names(names(tt))))
    setnames(tt, "tisch_key", "key")
  }
  au <- author_labels(ds, S)
  unfiltered <- grepl("unfiltered", man[dataset_name == ds, expression_type])

  mats <- list(); cds <- list(); samp_audit <- list(); gene_ref <- NULL
  for (i in seq_len(nrow(S))) {
    r <- S[i]
    msg("  ", ds, " / ", r$sample_id)
    x <- read_sample(ds, r)
    cnt <- x$counts
    int_ok <- is_integer_counts(cnt)
    bc <- colnames(cnt)
    if (unfiltered) {
      # The deposited matrix holds every barcode, so the curated cell calls
      # define the cell set. Cell calling is not re-run.
      calls <- if (!is.null(au)) au$key else if (!is.null(tt)) tt$key else
        stop("unfiltered matrix and no cell calls for ", r$sample_id)
      keep_keys <- calls[startsWith(calls, paste0(r$sample_id, "|"))]
      keep <- paste(r$sample_id, norm_bc(bc), sep = "|") %in% keep_keys
      cnt <- cnt[, keep, drop = FALSE]; bc <- bc[keep]
    }
    if (is.null(gene_ref)) gene_ref <- x$genes
    if (length(mats) && !identical(rownames(cnt), rownames(mats[[1]]))) {
      # Dense per-sample text files (GSE145281) list different genes per file.
      # Absent genes cannot be assumed zero, so only genes reported in every
      # sample are kept and the loss is logged.
      if (r$format != "dense_txt") stop("gene rows differ between samples of ", ds, " (", r$sample_id, ")")
      common <- intersect(rownames(mats[[1]]), rownames(cnt))
      msg("  gene sets differ across dense files; keeping ", length(common), " genes common so far (dropped ",
          length(union(rownames(mats[[1]]), rownames(cnt))) - length(common), ")")
      mats <- lapply(mats, function(z) z[common, , drop = FALSE])
      cnt <- cnt[common, , drop = FALSE]
      gene_ref <- gene_ref[match(common, gene_ref$gene_symbol)]
    }
    cd <- data.table(dataset = ds, sample_id = r$sample_id, gsm = r$gsm, donor_id = r$donor_id,
                     tissue = r$tissue, sample_type = r$sample_type, barcode = bc)
    cd[, key := paste(sample_id, norm_bc(barcode), sep = "|")]
    cd <- cbind(cd, qc_fields(cnt, x$genes$gene_symbol[match(rownames(cnt), rownames(x$counts))]))
    samp_audit[[i]] <- data.table(dataset = ds, sample_id = r$sample_id, donor_id = r$donor_id,
                                  tissue = r$tissue, n_genes_matrix = nrow(cnt), n_cells = ncol(cnt),
                                  integer_counts = int_ok, unfiltered_input = unfiltered,
                                  n_cells_tisch2 = if (is.null(tt)) NA_integer_ else sum(startsWith(tt$key, paste0(r$sample_id, "|"))),
                                  n_cells_author = if (is.null(au)) NA_integer_ else sum(startsWith(au$key, paste0(r$sample_id, "|"))))
    mats[[i]] <- cnt; cds[[i]] <- cd
  }
  if (length(unique(lapply(mats, rownames))) > 1) {
    common <- Reduce(intersect, lapply(mats, rownames))
    mats <- lapply(mats, function(z) z[common, , drop = FALSE])
    gene_ref <- gene_ref[match(common, gene_ref$gene_symbol)]
  }
  counts <- do.call(cbind, mats)
  cd <- rbindlist(cds)
  cd[, cell_id := make_cell_ids(dataset, sample_id, barcode)]
  colnames(counts) <- cd$cell_id
  if (!is.null(tt)) cd <- merge(cd, tt, by = "key", all.x = TRUE, sort = FALSE)
  if (!is.null(au)) cd <- merge(cd, au, by = "key", all.x = TRUE, sort = FALSE)
  cd <- cd[match(colnames(counts), cell_id)]

  # join audit: how many annotated cells are found in the raw matrices
  join <- data.table(dataset = ds,
    tisch2_cells = if (is.null(tt)) NA_integer_ else nrow(tt[sub("\\|.*", "", key) %in% S$sample_id]),
    tisch2_matched = if (is.null(tt)) NA_integer_ else sum(tt$key %in% cd$key),
    author_cells = if (is.null(au)) NA_integer_ else nrow(au[sub("\\|.*", "", key) %in% S$sample_id]),
    author_matched = if (is.null(au)) NA_integer_ else sum(au$key %in% cd$key))
  if (!is.null(tt) && join$tisch2_cells > 0)
    assert_cells_match(cd$key, tt$key[sub("\\|.*", "", tt$key) %in% S$sample_id],
                       min_frac = 0.9, what = paste(ds, "TISCH2 cells"))
  if (!is.null(au) && join$author_cells > 0)
    assert_cells_match(cd$key, au$key[sub("\\|.*", "", au$key) %in% S$sample_id],
                       min_frac = 0.9, what = paste(ds, "author cells"))

  assert_donors(cd$donor_id, ds)
  cd <- drop_outcome_fields(cd)
  removed <- attr(cd, "dropped_outcome_fields")
  list(counts = counts, genes = gene_ref, coldata = cd, samples = rbindlist(samp_audit),
       join = join, removed = removed)
}

# ---- Li et al. 2022 (AnnData) -------------------------------------------------------
# Expected obs columns: patient (12 donors, PD*), orig.ident (Sanger channel),
# summaryDescription (Tumour, Tumour-normal, Normal kidney, Blood, Metastasis,
# Thrombus, Fat, Normal adrenal), region (a-i, n, t), broad_type (RCC, Epi_PT,
# Epi_non-PT, EC, Fibro, immune types) and annotation (fine labels).
# X holds integer UMI counts (CSC over genes). Genes are GRCh37-era symbols.
LI_TISSUE <- c("Tumour" = "tumour", "Tumour-normal" = "other", "Thrombus" = "other",
               "Metastasis" = "metastasis", "Normal kidney" = "normal_kidney", "Blood" = "blood",
               "Fat" = "other", "Normal adrenal" = "other")
import_li2022 <- function(ds) {
  f <- raw_path(ds, "author", "RCC_upload_final_raw_counts.h5ad")
  obs <- h5ad_read_obs(f)
  need <- c("patient", "orig.ident", "summaryDescription", "region", "broad_type", "annotation")
  if (length(setdiff(need, names(obs))))
    stop("Li2022 adapter: expected obs columns missing: ", paste(setdiff(need, names(obs)), collapse = ", "),
         " (found: ", paste(names(obs), collapse = ", "), ")")
  unknown <- setdiff(unique(obs$summaryDescription), names(LI_TISSUE))
  if (length(unknown)) stop("Li2022 adapter: unmapped summaryDescription values: ", paste(unknown, collapse = ", "))
  counts <- h5ad_read_counts(f)
  cd <- data.table(dataset = ds, sample_id = paste(obs$orig.ident, obs$summaryDescription, sep = "_"), gsm = "",
                   donor_id = obs$patient, tissue = unname(LI_TISSUE[obs$summaryDescription]),
                   sample_type = obs$summaryDescription, region = obs$region, barcode = obs$obs_names,
                   author_label = obs$broad_type, author_fine_label = obs$annotation)
  cd[author_fine_label %in% c("Low quality", "Unknown"), author_label := "ua"]
  cd[, key := paste(sample_id, barcode, sep = "|")]
  cd[, cell_id := make_cell_ids(dataset, sample_id, barcode)]
  colnames(counts) <- cd$cell_id
  genes <- data.table(gene_id = NA_character_, gene_symbol = rownames(counts))
  cd <- cbind(cd, qc_fields(counts, genes$gene_symbol))
  assert_donors(cd$donor_id, ds)
  cd <- drop_outcome_fields(cd)
  sa <- cd[, .(n_cells = .N, n_rcc_labelled = sum(author_label == "RCC")),
           by = .(dataset, sample_id, donor_id, tissue, sample_type)]
  sa[, `:=`(n_genes_matrix = nrow(counts), integer_counts = is_integer_counts(counts), unfiltered_input = FALSE)]
  save_tsv(unique(cd[, .(dataset, sample_id, donor_id, tissue, sample_type)]), paste0("02_derived_sample_sheet_", ds, ".tsv"))
  list(counts = counts, genes = genes, coldata = cd, samples = sa,
       join = data.table(dataset = ds, tisch2_cells = NA_integer_, tisch2_matched = NA_integer_,
                         author_cells = nrow(cd), author_matched = nrow(cd)),
       removed = attr(cd, "dropped_outcome_fields"))
}

# ---- run ----------------------------------------------------------------------------
res_man <- read_tsv("01_manifest_resolved.tsv")
todo <- res_man[acquisition_status == "complete", dataset_name]
if (length(ONLY)) todo <- intersect(todo, ONLY)
skipped <- setdiff(res_man$dataset_name, todo)

audit <- list(); joins <- list(); removed <- list()
for (ds in todo) {
  msg("importing ", ds)
  out <- tryCatch({
    x <- if (ds == "KIRC_Li2022") import_li2022(ds) else import_geo_dataset(ds)
    ints <- is_integer_counts(x$counts)
    genes <- x$genes
    ensembl <- mean(grepl("^ENSG", genes$gene_id), na.rm = TRUE)
    sce <- SingleCellExperiment(assays = list(counts = x$counts),
                                rowData = S4Vectors::DataFrame(genes),
                                colData = S4Vectors::DataFrame(x$coldata))
    metadata(sce) <- list(dataset = ds, created = format(Sys.time(), tz = "UTC", usetz = TRUE),
                          config_sha256 = CFG_HASH, assay_scale = "raw_integer_umi",
                          removed_outcome_fields = x$removed)
    saveRDS(sce, file.path(DERIVED_DIR, paste0(ds, ".sce.rds")))
    save_tsv(x$samples, paste0("02_schema_audit_", ds, ".tsv"))
    joins[[ds]] <- x$join
    if (length(x$removed)) removed[[ds]] <- data.table(dataset = ds, field = x$removed)
    data.table(dataset = ds, status = "imported", reason = "",
               n_samples = uniqueN(x$coldata$sample_id), n_donors = uniqueN(x$coldata$donor_id),
               n_cells = ncol(x$counts), n_genes = nrow(x$counts),
               integer_counts = ints, gene_id_system = if (ensembl > 0.9) "ensembl" else "symbol",
               duplicated_cell_ids = anyDuplicated(colnames(x$counts)) > 0,
               duplicated_gene_ids = sum(duplicated(na.omit(genes$gene_id))),
               duplicated_gene_symbols = sum(duplicated(na.omit(genes$gene_symbol))))
  }, error = function(e) {
    msg("FAILED ", ds, ": ", conditionMessage(e))
    data.table(dataset = ds, status = "failed", reason = conditionMessage(e))
  })
  audit[[ds]] <- out
  gc(verbose = FALSE)
}
for (ds in skipped)
  audit[[ds]] <- data.table(dataset = ds, status = "not_imported",
                            reason = paste("acquisition", res_man[dataset_name == ds, acquisition_status]))

if (!length(ONLY)) {
  save_tsv(rbindlist(audit, fill = TRUE), "02_schema_audit.tsv")
  save_tsv(rbindlist(joins, fill = TRUE), "02_label_join_audit.tsv")
  save_tsv(if (length(removed)) rbindlist(removed) else data.table(dataset = character(), field = character()),
           "02_removed_outcome_fields.tsv")
}
print(rbindlist(audit, fill = TRUE))
write_session_info("02_import_dataset")
