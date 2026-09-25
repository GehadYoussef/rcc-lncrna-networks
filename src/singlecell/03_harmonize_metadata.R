# 03_harmonize_metadata.R: compartments, label concordance and sample flags.
# Per-cell compartment priority: (1) author annotation (SingleR labels for datasets
# deposited without one), (2) TISCH2 malignancy "Malignant cells", (3) TISCH2 major
# lineage via map_compartment(). TISCH2 labels normal renal epithelium "Stromal cells"
# in the malignancy column, hence step 3.
# Author-TISCH2 agreement is tabulated. Samples are flagged for multi-patient
# donor fields, replicate libraries of one donor and libraries reused across
# datasets. No outcome field is used.
# Outputs: data/derived/singlecell/<dataset>.coldata.rds, results/singlecell/
#          03_label_concordance.tsv, 03_compartment_mapping.tsv,
#          03_sample_flags.tsv, 03_duplicate_library_check.tsv,
#          03_reference_annotation_summary.tsv

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
suppressPackageStartupMessages(library(SingleCellExperiment))
banner("03  metadata harmonisation")
start_log("03_harmonize_metadata")

audit <- read_tsv("02_schema_audit.tsv")
datasets <- audit[status == "imported", dataset]

# ---- reference-based annotation (datasets deposited without labels) -----------
# SingleR against the reference dataset's author labels, on protein-coding genes
# only (reference_annotation in singlecell.yml). Log-normalisation is used for
# label transfer only. Stored counts are not modified.
RA <- CFG$reference_annotation %||% list()
ref_models <- list()
protein_coding_symbols <- function() {
  f <- file.path(DERIVED_DIR, "gencode_genes.rds")
  g <- if (file.exists(f)) readRDS(f) else read_gencode_genes(GENCODE_GTF)
  list(ids = strip_ensembl_version(g$gene_id[g$gene_type == "protein_coding"]))
}
train_reference <- function(spec, query_genes) {
  key <- paste(spec$reference_dataset, spec$reference_label, digest::digest(sort(query_genes)))
  if (!is.null(ref_models[[key]])) return(ref_models[[key]])
  msg("training SingleR reference from ", spec$reference_dataset, " (", spec$reference_label, ")")
  rs <- readRDS(file.path(DERIVED_DIR, paste0(spec$reference_dataset, ".sce.rds")))
  lab <- as.character(colData(rs)[[spec$reference_label]])
  ok <- !is.na(lab) & !tolower(lab) %in% c("ua", "unknown", "")
  if ("author_doublet" %in% names(colData(rs))) ok <- ok & !as.logical(colData(rs)$author_doublet) %in% TRUE
  rs <- rs[, ok]; lab <- lab[ok]
  keep_lab <- names(which(table(lab) >= spec$min_reference_cells_per_label))
  rs <- rs[, lab %in% keep_lab]; lab <- lab[lab %in% keep_lab]
  pc <- protein_coding_symbols()$ids
  # protein-coding genes present in both reference and query
  rs <- rs[strip_ensembl_version(rownames(rs)) %in% intersect(pc, query_genes), ]
  rownames(rs) <- strip_ensembl_version(rownames(rs))
  rs <- scuttle::logNormCounts(rs)
  model <- SingleR::trainSingleR(rs, labels = lab, aggr.ref = TRUE)
  ref_models[[key]] <<- list(model = model, n_genes = nrow(rs), labels = table(lab))
  ref_models[[key]]
}
ref_ann_summary <- list()

conc <- list(); mapping <- list(); flags <- list(); lib_sig <- list()
for (ds in datasets) {
  msg("harmonising ", ds)
  sce <- readRDS(file.path(DERIVED_DIR, paste0(ds, ".sce.rds")))
  cd <- as.data.table(as.data.frame(colData(sce)))

  if (!is.null(RA[[ds]])) {
    spec <- RA[[ds]]
    if (!all(unlist(spec$feature_biotypes) == "protein_coding"))
      stop("reference annotation for ", ds, " must use protein-coding features only")
    tr <- train_reference(spec, strip_ensembl_version(rownames(sce)))
    q <- sce[, cd$n_genes >= CFG$cell_qc$min_genes]
    rownames(q) <- strip_ensembl_version(rownames(q))
    q <- scuttle::logNormCounts(q)
    pred <- SingleR::classifySingleR(q, tr$model, BPPARAM = BiocParallel::SerialParam())
    lab <- rep(NA_character_, nrow(cd))
    lab[cd$n_genes >= CFG$cell_qc$min_genes] <- pred$pruned.labels
    cd[, author_label := lab]
    cd[, author_label_method := paste0("SingleR(ref=", spec$reference_dataset, ", protein_coding, n_genes=", tr$n_genes, ")")]
    ref_ann_summary[[ds]] <- data.table(dataset = ds, reference = spec$reference_dataset,
                                        cells_scored = sum(cd$n_genes >= CFG$cell_qc$min_genes),
                                        cells_pruned = sum(is.na(pred$pruned.labels)),
                                        label = names(table(lab)), n = as.integer(table(lab)))
    rm(q, pred); gc(verbose = FALSE)
  }

  has_author <- "author_label" %in% names(cd)
  tm <- "tisch_Celltype..malignancy."; tl <- "tisch_Celltype..major.lineage."
  has_tisch <- tm %in% names(cd)

  comp_author <- if (has_author) map_compartment(cd$author_label) else rep(NA_character_, nrow(cd))
  comp_tisch <- if (has_tisch) {
    x <- map_compartment(cd[[tl]])
    x[!is.na(cd[[tm]]) & cd[[tm]] == "Malignant cells"] <- "malignant"
    x[is.na(cd[[tm]])] <- NA_character_
    x
  } else rep(NA_character_, nrow(cd))

  cd[, compartment_author := comp_author]
  cd[, compartment_tisch2 := comp_tisch]
  cd[, compartment := fifelse(!is.na(compartment_author), compartment_author, compartment_tisch2)]
  cd[, label_source := fifelse(!is.na(compartment_author), "author",
                        fifelse(!is.na(compartment_tisch2), "tisch2", "none"))]
  cd[is.na(compartment), compartment := "unlabelled"]

  # normal epithelium is only meaningful from non-tumour tissue. Epithelial-looking
  # cells inside a tumour sample are kept but flagged separately
  cd[compartment == "normal_epithelial" & tissue == "tumour", compartment := "epithelial_in_tumour"]
  cd[compartment == "normal_epithelial" & tissue != "normal_kidney", compartment := "epithelial_other_tissue"]
  # malignant calls outside tumour tissue are implausible and set aside
  cd[compartment == "malignant" & tissue != "tumour", compartment := "malignant_label_in_nontumour"]

  if (has_author && has_tisch) {
    both <- cd[!is.na(compartment_author) & !is.na(compartment_tisch2)]
    conc[[ds]] <- both[, .N, by = .(dataset, compartment_author, compartment_tisch2)][
      , agree := compartment_author == compartment_tisch2][]
  }
  lab_cols <- intersect(c("author_label", tl), names(cd))
  mapping[[ds]] <- rbindlist(lapply(lab_cols, function(lc)
    cd[, .N, by = c(lc, "compartment")][, .(dataset = ds, label_column = lc, original_label = get(lc),
                                            compartment, n_cells = N)]))

  # ---- sample-level flags --------------------------------------------------------
  f <- cd[, .(n_cells = .N), by = .(dataset, sample_id, donor_id, tissue)]
  f[, flags := ""]
  if ("tisch_Patient" %in% names(cd)) {
    multi <- unique(cd[grepl(",", tisch_Patient), sample_id])
    f[sample_id %in% multi, flags := paste0(flags, "tisch2_patient_field_lists_several_patients;")]
  }
  # replicate libraries: two libraries of one donor whose cell barcodes overlap
  # heavily are re-sequencing of one cell suspension, not independent samples
  by_donor <- split(cd[, .(sample_id, bc = sub("-[0-9]+$", "", barcode))], cd$donor_id)
  for (dn in names(by_donor)) {
    s <- split(by_donor[[dn]]$bc, by_donor[[dn]]$sample_id)
    if (length(s) < 2) next
    for (a in names(s)) for (b in names(s)) if (a < b) {
      j <- length(intersect(s[[a]], s[[b]])) / min(length(s[[a]]), length(s[[b]]))
      if (j > 0.5) f[sample_id %in% c(a, b),
                     flags := paste0(flags, sprintf("barcode_overlap_%.0f%%_with_%s;", 100 * j, setdiff(c(a, b), sample_id)))]
    }
  }
  f[tissue == "blood", flags := paste0(flags, "blood_sample_no_kidney_cells;")]
  flags[[ds]] <- f

  # library signature: sorted (barcode, n_umi) digest per sample, compared
  # across datasets to catch the same deposited library appearing twice
  lib_sig[[ds]] <- cd[, .(signature = digest::digest(paste(sort(paste(barcode, n_umi)), collapse = ";"), algo = "sha1"),
                          n_cells = .N), by = .(dataset, sample_id, donor_id)]

  cd <- drop_outcome_fields(cd)
  saveRDS(cd, file.path(DERIVED_DIR, paste0(ds, ".coldata.rds")))
  rm(sce); gc(verbose = FALSE)
}

ls_all <- rbindlist(lib_sig)
if (!nrow(ls_all)) ls_all <- data.table(dataset = character(), sample_id = character(), signature = character())
dup <- ls_all[, .(n = .N, datasets = paste(unique(dataset), collapse = ";"),
                  samples = paste(paste(dataset, sample_id), collapse = ";")), by = signature][n > 1]
save_tsv(if (nrow(dup)) dup else data.table(signature = character(), n = integer(), datasets = character(),
                                            samples = character(), result = character()),
         "03_duplicate_library_check.tsv")
msg("identical libraries across samples: ", nrow(dup))

save_tsv(if (length(ref_ann_summary)) rbindlist(ref_ann_summary) else data.table(dataset = character()),
         "03_reference_annotation_summary.tsv")
save_tsv(rbindlist(conc, fill = TRUE), "03_label_concordance.tsv")
save_tsv(rbindlist(mapping, fill = TRUE), "03_compartment_mapping.tsv")
save_tsv(rbindlist(flags, fill = TRUE), "03_sample_flags.tsv")
write_session_info("03_harmonize_metadata")
