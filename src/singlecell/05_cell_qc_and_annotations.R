# 05_cell_qc_and_annotations.R: cell QC, label checks and donor-level coverage.
# Applies the configured per-cell malignant rescue, then cell QC (thresholds in
# singlecell.yml, author doublet and cell-filter calls), counting exclusions by
# reason. Checks labels against marker panels without reassigning them, counts
# cells per donor x compartment, and sums raw counts per donor to test lncRNA
# coverage in malignant cells (lnc_min_counts_pseudobulk in at least
# lnc_min_donor_frac of evaluable donors). Donors, not cells, are the unit.
# Outputs: data/derived/singlecell/<dataset>.coldata_qc.rds, results/singlecell/
#          05_*.tsv (exclusions, rescue, markers, donor cells, coverage,
#          detection), figures/singlecell/05_composition, 05_lncRNA_coverage

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
suppressPackageStartupMessages({ library(SingleCellExperiment); library(ggplot2) })
banner("05  cell QC, label checks and donor-level lncRNA coverage")
start_log("05_cell_qc_and_annotations")

Q <- CFG$cell_qc; PB <- CFG$pseudobulk; FE <- CFG$feasibility
MIN_CELLS <- PB$min_cells_per_donor_compartment
audit <- read_tsv("02_schema_audit.tsv")
datasets <- audit[status == "imported", dataset]
man <- read_manifest()$manifest

excl <- list(); markers <- list(); dcc <- list(); cov <- list(); det <- list(); rescues <- list()
for (ds in datasets) {
  msg("QC ", ds)
  sce <- readRDS(file.path(DERIVED_DIR, paste0(ds, ".sce.rds")))
  cd <- readRDS(file.path(DERIVED_DIR, paste0(ds, ".coldata.rds")))
  g <- readRDS(file.path(DERIVED_DIR, paste0(ds, ".genes.rds")))
  stopifnot(identical(cd$cell_id, colnames(sce)))
  cnt <- counts(sce); rm(sce)
  require_raw_counts(cnt, paste(ds, "counts"))
  assert_no_outcome_fields(cd, paste(ds, "cell metadata"))

  # ---- malignant-cell rescue (configured compartments only, per-cell rule) --------
  sym_all <- fifelse(is.na(g$ref_symbol), g$gene_symbol, g$ref_symbol)
  cd[, compartment_before_rescue := compartment]
  RM <- CFG$rescue_malignant
  flagged <- unlist(RM$compartments[[ds]])
  if (length(flagged)) {
    resc <- rescue_malignant(cnt, sym_all, cd$compartment, flagged, unlist(CFG$markers$malignant_ccrcc),
                             RM$min_ccrcc_markers, unlist(RM$exclude_if_detected))
    cd[resc, compartment := "malignant"]
    rescues[[ds]] <- cd[compartment_before_rescue %in% flagged,
                        .(n_candidates = .N, n_rescued = sum(compartment == "malignant")),
                        by = .(dataset, donor_id, from_compartment = compartment_before_rescue)]
    msg("  rescued ", sum(resc), " malignant cells from ", paste(flagged, collapse = ", "))
  }

  # ---- QC ------------------------------------------------------------------------
  cd[, qc_reason := ""]
  add <- function(cond, why) cd[cond & qc_reason == "", qc_reason := why]
  excl_donors <- names(CFG$exclude_donors[[ds]] %||% list())
  add(cd$donor_id %in% excl_donors, "donor_excluded_not_ccRCC")
  if (isTRUE(Q$use_author_doublet_calls) && "author_doublet" %in% names(cd))
    add(!is.na(cd$author_doublet) & as.logical(cd$author_doublet), "author_doublet")
  if ("author_cellfilt" %in% names(cd))
    add(!is.na(cd$author_cellfilt) & as.logical(cd$author_cellfilt), "author_cell_filter")
  add(cd$n_genes < Q$min_genes, "low_genes")
  add(cd$n_umi < Q$min_umi, "low_umi")
  add(!is.na(cd$pct_mito) & cd$pct_mito > Q$max_pct_mito, "high_mito")
  add(cd$compartment %in% c("unlabelled", "other_uncertain", "malignant_label_in_nontumour"), "uncertain_or_unlabelled")
  add(cd$tissue == "blood", "blood_sample")
  excl[[ds]] <- cd[, .N, by = .(dataset, sample_id, donor_id, tissue, reason = fifelse(qc_reason == "", "retained", qc_reason))]
  keep <- cd$qc_reason == ""
  saveRDS(cd, file.path(DERIVED_DIR, paste0(ds, ".coldata_qc.rds")))
  cdk <- cd[keep]; ck <- cnt[, keep, drop = FALSE]
  msg("  retained ", sum(keep), " of ", nrow(cd), " cells")
  if (!sum(keep)) next

  # ---- marker check --------------------------------------------------------------
  sym <- g$ref_symbol %||% g$gene_symbol
  sym[is.na(sym)] <- g$gene_symbol[is.na(sym)]
  for (panel in names(CFG$markers)) {
    idx <- which(sym %in% CFG$markers[[panel]])
    if (!length(idx)) next
    detm <- ck[idx, , drop = FALSE] > 0
    frac <- sapply(split(seq_len(ncol(ck)), cdk$compartment), function(j)
      Matrix::rowMeans(detm[, j, drop = FALSE]))
    frac <- matrix(frac, nrow = length(idx), dimnames = list(sym[idx], sort(unique(cdk$compartment))))
    markers[[paste(ds, panel)]] <- data.table(dataset = ds, panel = panel, marker = rownames(frac),
                                              as.data.table(frac))
  }

  # ---- donor x compartment ------------------------------------------------------------
  d <- donor_compartment_table(cdk, MIN_CELLS)
  pooled <- cdk[tissue == "tumour" & compartment %in% c("immune", "endothelial", "fibroblast_stromal"),
                .(n_cells = .N), by = .(dataset, donor_id, tissue)][, `:=`(compartment = "pooled_nonmalignant_tme",
                                                                           evaluable = n_cells >= MIN_CELLS)]
  d <- rbind(d, pooled, use.names = TRUE)
  d[, analysis_role := man[dataset_name == ds, analysis_role]]
  dcc[[ds]] <- d

  # ---- donor-level lncRNA coverage (malignant compartment) ---------------------------
  lnc_rows <- which(g$is_lnc)
  mal <- cdk$compartment == "malignant"
  ev_mal <- d[compartment == "malignant" & evaluable == TRUE, donor_id]
  n_ev <- length(ev_mal)
  if (n_ev > 0) {
    sel <- mal & cdk$donor_id %in% ev_mal
    pb <- aggregate_pseudobulk(ck[lnc_rows, sel, drop = FALSE], cdk$donor_id[sel])
    pass_n <- Matrix::rowSums(pb >= FE$lnc_min_counts_pseudobulk)
    passing <- pass_n >= ceiling(FE$lnc_min_donor_frac * n_ev)
    pb_keys <- g$key[lnc_rows]
    cov[[ds]] <- data.table(dataset = ds, analysis_role = man[dataset_name == ds, analysis_role],
                            evaluable_malignant_donors = n_ev,
                            lnc_detected_in_malignant_cells = uniqueN(pb_keys[Matrix::rowSums(pb) > 0]),
                            lnc_passing_donor_filter = uniqueN(pb_keys[passing]))
    saveRDS(data.table(gene_key = pb_keys, symbol = g$ref_symbol[lnc_rows], donors_passing = pass_n, passing = passing),
            file.path(DERIVED_DIR, paste0(ds, ".lnc_malignant_passing.rds")))
  } else {
    cov[[ds]] <- data.table(dataset = ds, analysis_role = man[dataset_name == ds, analysis_role],
                            evaluable_malignant_donors = 0L, lnc_detected_in_malignant_cells = 0L,
                            lnc_passing_donor_filter = 0L)
  }

  # ---- detection prevalence by compartment -------------------------------------------
  if (length(lnc_rows)) {
    det_by <- lapply(split(seq_len(ncol(ck)), cdk$compartment), function(j) {
      fr <- Matrix::rowMeans(ck[lnc_rows, j, drop = FALSE] > 0)
      data.table(n_cells = length(j), lnc_detected_ge_1pct_cells = sum(fr >= 0.01),
                 lnc_detected_ge_5pct_cells = sum(fr >= 0.05), lnc_detected_ge_10pct_cells = sum(fr >= 0.10),
                 median_lnc_umi_per_cell = stats::median(Matrix::colSums(ck[lnc_rows, j, drop = FALSE])),
                 median_total_umi_per_cell = stats::median(Matrix::colSums(ck[, j, drop = FALSE])))
    })
    det[[ds]] <- rbindlist(det_by, idcol = "compartment")[, dataset := ds]
  }
  rm(cnt, ck); gc(verbose = FALSE)
}

save_tsv(rbindlist(excl), "05_cell_exclusions.tsv")
save_tsv(if (length(rescues)) rbindlist(rescues) else data.table(dataset = character()), "05_malignant_rescue.tsv")
MK <- rbindlist(markers, fill = TRUE)
save_tsv(MK, "05_marker_check.tsv")
# marker discordance: mean ccRCC-panel detection per assigned compartment
if (nrow(MK)) {
  comp_cols <- setdiff(names(MK), c("dataset", "panel", "marker"))
  disc <- melt(MK[panel == "malignant_ccrcc"], id.vars = c("dataset", "panel", "marker"), measure.vars = comp_cols,
               variable.name = "compartment", value.name = "detection")[!is.na(detection)]
  disc <- disc[, .(ccrcc_panel_mean_detection = mean(detection), n_markers = .N), by = .(dataset, compartment)]
  disc[, flag := fifelse(compartment != "malignant" & ccrcc_panel_mean_detection >= CFG$marker_discordance_threshold,
                         "probable_malignant_cells_under_nonmalignant_label", "")]
  save_tsv(disc, "05_marker_discordance.tsv")
}
DCC <- rbindlist(dcc, fill = TRUE)
save_tsv(DCC, "05_donor_compartment_cells.tsv")
save_tsv(rbindlist(cov, fill = TRUE), "05_lncRNA_donor_coverage.tsv")
save_tsv(rbindlist(det, fill = TRUE), "05_lncRNA_detection_by_compartment.tsv")

# ---- figures (donors, not cells, are the units shown) -----------------------------------
if (nrow(DCC)) {
  pd <- DCC[compartment != "pooled_nonmalignant_tme"]
  p1 <- ggplot(pd, aes(x = donor_id, y = n_cells, fill = compartment)) +
    geom_col() + facet_grid(tissue ~ dataset, scales = "free", space = "free_x") +
    labs(x = "Donor", y = "Cells after QC", fill = "Compartment",
         title = "Cell composition per donor after QC",
         subtitle = sprintf("Inference unit is the donor; a donor-compartment is evaluable with >= %d cells", MIN_CELLS)) +
    theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  save_fig(p1, "05_composition", width = 12, height = 6)
}
cv <- rbindlist(cov, fill = TRUE)
if (nrow(cv)) {
  p2 <- ggplot(melt(cv, id.vars = c("dataset", "analysis_role", "evaluable_malignant_donors")),
               aes(x = dataset, y = value, fill = variable)) +
    geom_col(position = "dodge") + coord_flip() +
    geom_hline(yintercept = FE$min_lnc_passing_discovery, linetype = 2) +
    labs(x = NULL, y = "lncRNA genes", fill = NULL,
         title = "lncRNA coverage in malignant cells (donor-level pseudobulk)",
         subtitle = sprintf("Dashed line: go threshold for the discovery dataset (%d)", FE$min_lnc_passing_discovery)) +
    theme_bw(base_size = 9)
  save_fig(p2, "05_lncRNA_coverage", width = 8, height = 4.5)
}
write_session_info("05_cell_qc_and_annotations")
