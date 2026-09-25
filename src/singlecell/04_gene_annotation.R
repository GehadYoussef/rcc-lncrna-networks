# 04_gene_annotation.R: reference gene mapping and matrix-level lncRNA audit.
# Reference is GENCODE v36 (GRCh38), the gene model of GDC data release 46 used
# for TCGA-KIRC and CPTAC-3. Ensembl IDs are matched unversioned, and symbol-only
# matrices map through symbols unique in the reference (map_genes in lib_sc.R).
# Per dataset, counts lncRNAs present, detected and shared with the bulk
# cohorts. Bulk expression uses the bulk pipeline's filter (FPKM >= 0.30 in
# >= 20% of samples). No bulk clinical field is read.
# Outputs: data/derived/singlecell/<dataset>.genes.rds, results/singlecell/
#          04_reference_provenance.tsv, 04_lncRNA_matrix_audit.tsv,
#          04_gene_mapping_summary.tsv, 04_tisch2_gene_retention.tsv,
#          04_bulk_lncRNA_features.tsv

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
suppressPackageStartupMessages(library(SingleCellExperiment))
banner("04  gene annotation and lncRNA matrix audit")
start_log("04_gene_annotation")

BIOTYPES <- CFG$annotation$lncrna_biotypes

# ---- reference -----------------------------------------------------------------
if (!file.exists(GENCODE_GTF)) {
  msg("GENCODE GTF not found locally; downloading ", CFG$annotation$gencode_url)
  st <- download_resumable(CFG$annotation$gencode_url, GENCODE_GTF, CFG$download$user_agent)
  if (st$status == "failed") stop("could not obtain GENCODE v", CFG$annotation$gencode_release)
}
ref_cache <- file.path(DERIVED_DIR, "gencode_genes.rds")
ref_sha <- sha256_file(GENCODE_GTF)
ref <- if (file.exists(ref_cache) && identical(attr(readRDS(ref_cache), "sha256"), ref_sha)) readRDS(ref_cache) else {
  msg("parsing ", basename(GENCODE_GTF))
  r <- read_gencode_genes(GENCODE_GTF); attr(r, "sha256") <- ref_sha; saveRDS(r, ref_cache); r }
header <- readLines(gzfile(GENCODE_GTF), n = 5)
save_tsv(data.table(field = c("genome_build", "gencode_release", "source_url", "local_path", "sha256",
                              "file_date", "gtf_header", "n_genes", "n_lncRNA_genes"),
                    value = c(CFG$annotation$genome_build, CFG$annotation$gencode_release, CFG$annotation$gencode_url,
                              sub(paste0("^", dirname(SUB_ROOT), "/"), "", GENCODE_GTF), ref_sha,
                              format(file.mtime(GENCODE_GTF), "%Y-%m-%d"), gsub("	", " ", paste(header[startsWith(header, "#")], collapse = " | ")),
                              nrow(ref), sum(is_lncrna(ref$gene_type, BIOTYPES)))),
         "04_reference_provenance.tsv")
legacy <- NULL
if (!is.na(LEGACY_GTF)) {
  if (!file.exists(LEGACY_GTF)) {
    msg("downloading legacy GENCODE v", CFG$annotation$legacy_gencode_release, " for the symbol bridge")
    st <- download_resumable(CFG$annotation$legacy_gencode_url, LEGACY_GTF, CFG$download$user_agent)
    if (st$status == "failed") stop("could not obtain legacy GENCODE GTF")
  }
  leg_cache <- file.path(DERIVED_DIR, "gencode_legacy_genes.rds")
  leg_sha <- sha256_file(LEGACY_GTF)
  legacy <- if (file.exists(leg_cache) && identical(attr(readRDS(leg_cache), "sha256"), leg_sha)) readRDS(leg_cache) else {
    msg("parsing ", basename(LEGACY_GTF)); r <- read_gencode_genes(LEGACY_GTF); attr(r, "sha256") <- leg_sha; saveRDS(r, leg_cache); r }
  prov <- read_tsv("04_reference_provenance.tsv")
  save_tsv(rbind(prov, data.table(field = c("legacy_gencode_release", "legacy_source_url", "legacy_sha256", "legacy_n_genes"),
                                  value = c(CFG$annotation$legacy_gencode_release, CFG$annotation$legacy_gencode_url,
                                            leg_sha, nrow(legacy)))), "04_reference_provenance.tsv")
}
n_ref_lnc <- sum(is_lncrna(ref$gene_type, BIOTYPES))
msg("reference: ", nrow(ref), " genes, ", n_ref_lnc, " lncRNA")

# ---- bulk lncRNA features ----------------------------------------------------------
bulk_features <- function(path, label) {
  if (!file.exists(path)) return(NULL)
  x <- readRDS(path)
  fp <- x$fpkm; ga <- as.data.table(x$gene_ann)
  lnc <- ga$gene_type %in% BIOTYPES
  pass <- rowMeans(fp[lnc, , drop = FALSE] >= 0.30) >= 0.20
  out <- data.table(cohort = label, ens_key = strip_ensembl_version(ga$gene_id[lnc]), gene_id = ga$gene_id[lnc],
                    gene_name = ga$gene_name[lnc], bulk_expressed = pass)
  setnames(out, "ens_key", "key")
}
bulk <- rbindlist(list(bulk_features(BULK_TCGA, "TCGA-KIRC"), bulk_features(BULK_CPTAC, "CPTAC-3")))
if (nrow(bulk)) {
  save_tsv(bulk, "04_bulk_lncRNA_features.tsv")
  tcga_ok <- bulk[cohort == "TCGA-KIRC" & bulk_expressed == TRUE, key]
  cptac_ok <- bulk[cohort == "CPTAC-3" & bulk_expressed == TRUE, key]
  msg("bulk lncRNAs passing expression filter: TCGA ", length(tcga_ok), ", CPTAC ", length(cptac_ok))
} else { tcga_ok <- cptac_ok <- character() }

# ---- per dataset ------------------------------------------------------------------
audit <- read_tsv("02_schema_audit.tsv")
datasets <- audit[status == "imported", dataset]
out <- list(); maps <- list()
for (ds in datasets) {
  msg("auditing ", ds)
  sce <- readRDS(file.path(DERIVED_DIR, paste0(ds, ".sce.rds")))
  g <- map_genes(as.data.table(as.data.frame(rowData(sce)))[, .(gene_id, gene_symbol)], ref, legacy)
  cnt <- counts(sce)
  g[, n_cells_detected := as.numeric(Matrix::rowSums(cnt > 0))]
  g[, total_umi := as.numeric(Matrix::rowSums(cnt))]
  g[, is_lnc := is_lncrna(gene_type, BIOTYPES)]
  coll <- ensembl_collisions(g$gene_id[grepl("^ENSG", g$gene_id)])
  saveRDS(g, file.path(DERIVED_DIR, paste0(ds, ".genes.rds")))

  lnc <- g[is_lnc == TRUE]
  lnc_keys <- unique(lnc$key)
  det_keys <- unique(lnc[n_cells_detected > 0, key])
  out[[ds]] <- data.table(
    dataset = ds, id_system = if (mean(grepl("^ENSG", g$gene_id), na.rm = TRUE) > 0.9) "ensembl" else "symbol",
    n_matrix_genes = nrow(g), n_mapped = sum(!is.na(g$ref_gene_id)),
    n_unmapped = sum(g$map_status == "unmapped"), n_symbol_ambiguous = sum(g$map_status == "symbol_ambiguous"),
    n_legacy_bridged = sum(g$map_status == "symbol_legacy_bridge"),
    n_ensembl_version_collisions = nrow(coll),
    lnc_in_reference = n_ref_lnc, lnc_in_matrix = length(lnc_keys),
    lnc_detected_any_cell = length(det_keys),
    lnc_share_of_total_umi = round(sum(lnc$total_umi) / sum(g$total_umi), 5),
    pc_detected_any_cell = g[gene_type == "protein_coding" & n_cells_detected > 0, uniqueN(key)],
    lnc_detected_in_TCGA_expressed = length(intersect(det_keys, tcga_ok)),
    lnc_detected_in_CPTAC_expressed = length(intersect(det_keys, cptac_ok)),
    lnc_detected_in_both_bulk = length(Reduce(intersect, list(det_keys, tcga_ok, cptac_ok))))
  maps[[ds]] <- g[, .N, by = map_status][, dataset := ds][]
  rm(sce, cnt); gc(verbose = FALSE)
}
save_tsv(rbindlist(out), "04_lncRNA_matrix_audit.tsv")
save_tsv(rbindlist(maps), "04_gene_mapping_summary.tsv")
print(rbindlist(out))

# ---- TISCH2 gene retention -------------------------------------------------------
# The TISCH2 "Expression" download holds genes x cell-group tables of mean
# log-normalised expression. Only its gene list is read, to count lncRNAs that
# TISCH2 retains.
ret <- list()
for (ds in datasets) {
  z <- file.path(RAW_DIR, ds, "tisch2", paste0(ds, "_Expression.zip"))
  if (!file.exists(z)) next
  lst <- utils::unzip(z, list = TRUE)$Name
  tab <- grep("Celltype_malignancy[.]txt$", lst, value = TRUE)[1]
  if (is.na(tab)) tab <- grep("[.]txt$", lst, value = TRUE)[1]
  if (is.na(tab)) { ret[[ds]] <- data.table(dataset = ds, note = paste("no table in zip:", paste(lst, collapse = ","))); next }
  con <- unz(z, tab); gn <- fread(text = readLines(con), sep = "	", header = FALSE, skip = 1, select = 1)[[1]]
  tg <- map_genes(data.table(gene_id = NA_character_, gene_symbol = gn), ref, legacy)
  raw <- readRDS(file.path(DERIVED_DIR, paste0(ds, ".genes.rds")))
  raw_det <- raw[is_lnc == TRUE & n_cells_detected > 0, unique(key)]
  t_lnc <- tg[is_lncrna(gene_type, BIOTYPES), unique(key)]
  ret[[ds]] <- data.table(dataset = ds, tisch2_table = basename(tab), tisch2_content = "mean expression per cell group (not cell-level)",
                          tisch2_genes = length(gn), tisch2_symbols_unmapped_to_v36 = sum(tg$map_status == "unmapped"),
                          tisch2_symbols_ambiguous = sum(tg$map_status == "symbol_ambiguous"),
                          tisch2_lnc_mapped = length(t_lnc), raw_lnc_detected = length(raw_det),
                          raw_lnc_detected_absent_from_tisch2 = length(setdiff(raw_det, t_lnc)))
}
save_tsv(rbindlist(ret, fill = TRUE), "04_tisch2_gene_retention.tsv")
write_session_info("04_gene_annotation")
