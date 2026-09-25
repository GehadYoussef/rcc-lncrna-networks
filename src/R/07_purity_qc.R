# 07_purity_qc.R: tumour composition, library quality and the leading lncRNA axis.
#
# The leading lncRNA axis is PC1 of the observed lncRNA matrix. It is tested
# against tumour composition (ESTIMATE stromal and immune scores) and STAR
# library quality, with protein-coding PC1 as a negative control.
# Inputs: caches dataset.rds, networks.rds, expr_raw.rds.
# Outputs: caches estimate_scores.rds, star_qc.rds and lnc_global_axis.rds,
#   results 07_* (PC variance, axis correlations, per-gene quality correlation,
#   nested Cox models, non-feature fraction as exposure) and figure
#   07_axis_vs_purity_qc.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(estimate); library(ggplot2)
})
banner("07 | Purity, library quality and the leading lncRNA axis")

ds     <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets   <- readRDS(file.path(CACHE_DIR, "networks.rds"))
cohort <- as.data.table(ds$cohort_full)
stopifnot(all(c("pct_noFeature", "pct_multimapping", "pct_unmapped",
                "pct_ambiguous", "libsize", "tss") %in% names(cohort)))

# -----------------------------------------------------------------------------
# 1. ESTIMATE stromal / immune / purity scores for every sample of cohort_full
# -----------------------------------------------------------------------------
# Computed on the full protein-coding matrix so no signature genes are lost to
# the WGCNA variance filter. Recomputed if the cache misses any sample.
est_rds <- file.path(CACHE_DIR, "estimate_scores.rds")
est <- if (file.exists(est_rds)) readRDS(est_rds) else NULL
if (!is.null(est) && !all(cohort$sample_barcode %in% est$sample_barcode)) {
  msg("Cached ESTIMATE scores cover ", nrow(est), " of ", nrow(cohort),
      " samples; recomputing")
  est <- NULL
}
if (is.null(est)) {
  raw  <- readRDS(file.path(CACHE_DIR, "expr_raw.rds"))
  pc   <- which(raw$gene_ann$gene_type == "protein_coding")
  # raw$fpkm has one column per file and some barcodes have two files, so
  # columns are selected by the file_id kept in 01.
  .ci  <- match(cohort$file_id, raw$sample_info$file_id)
  stopifnot(!anyNA(.ci), !anyDuplicated(.ci))
  m    <- raw$fpkm[pc, .ci, drop = FALSE]
  colnames(m) <- cohort$sample_barcode
  sym  <- raw$gene_ann$gene_name[pc]
  keep <- !is.na(sym) & nzchar(sym) & !duplicated(sym)
  m    <- log2(m[keep, , drop = FALSE] + 1)
  rownames(m) <- sym[keep]
  rm(raw)
  msg("ESTIMATE input: ", nrow(m), " protein-coding genes x ", ncol(m), " samples")

  in_f  <- file.path(CACHE_DIR, "estimate_input.txt")
  fil_f <- file.path(CACHE_DIR, "estimate_filtered.gct")
  sc_f  <- file.path(CACHE_DIR, "estimate_scores.gct")
  write.table(data.frame(GeneSymbol = rownames(m), m, check.names = FALSE),
              in_f, sep = "\t", quote = FALSE, row.names = FALSE)
  filterCommonGenes(input.f = in_f, output.f = fil_f, id = "GeneSymbol")
  # The "illumina" platform gives no tumour-purity value (its calibration is
  # Affymetrix-only), so ESTIMATEScore serves as an inverse purity proxy.
  estimateScore(fil_f, sc_f, platform = "illumina")

  # GCT turns barcode hyphens into dots but keeps column order, so barcodes
  # are restored by position from the input matrix.
  g   <- fread(sc_f, skip = 2, header = TRUE)
  mat <- as.matrix(g[, -c(1, 2), with = FALSE])
  rownames(mat) <- g[[1]]
  stopifnot(ncol(mat) == ncol(m),
            all(c("StromalScore", "ImmuneScore", "ESTIMATEScore") %in% rownames(mat)))
  est <- data.table(sample_barcode = colnames(m),
                    StromalScore  = as.numeric(mat["StromalScore", ]),
                    ImmuneScore   = as.numeric(mat["ImmuneScore", ]),
                    ESTIMATEScore = as.numeric(mat["ESTIMATEScore", ]))
  saveRDS(est, est_rds)
}
est <- as.data.table(est)[sample_barcode %in% cohort$sample_barcode]
msg("ESTIMATE scores for ", nrow(est), " samples")

# -----------------------------------------------------------------------------
# 2. STAR alignment summary for every sample of cohort_full
# -----------------------------------------------------------------------------
# Taken from the cohort table, where 01 read them from the expression file.
qc <-cohort[, .(sample_barcode, pct_unmapped, pct_multimapping, pct_noFeature,
                 pct_ambiguous, assigned_reads = libsize)]
saveRDS(qc, file.path(CACHE_DIR, "star_qc.rds"))
msg("STAR QC metrics cached for ", nrow(qc), " samples")

# Shared covariates. Composition and quality terms are standardised and
# clinical terms keep their natural scales.
cov <- merge(cohort[, .(sample_barcode, os_time, os_event, age, sex, T_stage,
                        N_pos, M1, grade_num, stage_group, grade_group, tss,
                        pct_unmapped, pct_multimapping, pct_noFeature,
                        pct_ambiguous, assigned_reads = libsize, libsize)],
             est, by = "sample_barcode")
cov[, `:=`(male    = as.numeric(sex == "male"),
           stromal = as.numeric(scale(StromalScore)),
           immune  = as.numeric(scale(ImmuneScore)),
           nf      = as.numeric(scale(pct_noFeature)),
           nf_log  = as.numeric(scale(log10(pct_noFeature))),
           mm      = as.numeric(scale(pct_multimapping)),
           dep     = as.numeric(scale(log10(assigned_reads))))]

# -----------------------------------------------------------------------------
# 3. Leading axes: PC1 of the observed lncRNA matrix, and the protein-coding
#    negative control
# -----------------------------------------------------------------------------
# Matrices are column-centred, unscaled. Each PC is signed to correlate
# positively with per-sample mean expression.
pc_summary <- function(X, label, k = 5) {
  s     <- svd(scale(X, center = TRUE, scale = FALSE), nu = k, nv = k)
  share <- s$d^2 / sum(s$d^2)
  q     <- cohort[match(rownames(X), sample_barcode)]
  mexp  <- rowMeans(X)
  pcs   <- s$u[, seq_len(k), drop = FALSE]
  for (j in seq_len(k)) if (cor(pcs[, j], mexp) < 0) pcs[, j] <- -pcs[, j]
  r <- function(a, b) cor(a, b, method = "spearman", use = "complete.obs")
  tbl <- rbindlist(lapply(seq_len(k), function(j) data.table(
    matrix = label, pc = j, n = nrow(X),
    var_explained  = round(share[j], 4),
    rho_noFeature  = round(r(pcs[, j], q$pct_noFeature), 3),
    rho_multimap   = round(r(pcs[, j], q$pct_multimapping), 3),
    rho_log_depth  = round(r(pcs[, j], log10(q$libsize)), 3),
    rho_mean_expr  = round(r(pcs[, j], mexp), 3))))
  list(tbl = tbl, pc1 = setNames(pcs[, 1], rownames(X)))
}
msg("PCA of the four expression matrices ...")
pc_lnc_obs   <- pc_summary(obs_expr(nets$lnc),  "lncRNA_observed")
pc_lnc_res   <- pc_summary(nets$lnc$expr,       "lncRNA_residualised")
pc_mrna_obs  <- pc_summary(obs_expr(nets$mrna), "protein_coding_observed")
pc_mrna_res  <- pc_summary(nets$mrna$expr,      "protein_coding_residualised")
pc_tbl <- rbindlist(list(pc_lnc_obs$tbl, pc_lnc_res$tbl,
                         pc_mrna_obs$tbl, pc_mrna_res$tbl))
save_tsv(pc_tbl, "07_pc_variance_explained.tsv")
print(pc_tbl[pc == 1])

axis <- data.table(sample_barcode = names(pc_lnc_obs$pc1), lnc_axis = unname(pc_lnc_obs$pc1))
saveRDS(axis, file.path(CACHE_DIR, "lnc_global_axis.rds"))
pc_axis <- data.table(sample_barcode = names(pc_mrna_obs$pc1), pc_axis = unname(pc_mrna_obs$pc1))
msg("lncRNA axis: ", nrow(axis), " samples; protein-coding PC1: ", nrow(pc_axis), " samples")

# -----------------------------------------------------------------------------
# 4. Correlation of each axis with composition and alignment quality
# -----------------------------------------------------------------------------
vars <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "pct_unmapped",
          "pct_multimapping", "pct_noFeature", "pct_ambiguous",
          "assigned_reads", "libsize")
axis_cor <- function(ax, label) {
  d <- merge(ax, cov, by = "sample_barcode")
  setnames(d, names(ax)[2], "axis")
  d[, axis_z := as.numeric(scale(axis))]
  out <- rbindlist(lapply(vars, function(v) {
    ct <- cor.test(d$axis_z, d[[v]], method = "spearman", exact = FALSE)
    data.table(matrix = label, variable = v,
               n = sum(complete.cases(d[, c("axis_z", v), with = FALSE])),
               spearman_rho = round(unname(ct$estimate), 3), p = signif(ct$p.value, 3))
  }))
  out[, abs_rho := abs(spearman_rho)]
  setorder(out, -abs_rho); out[, abs_rho := NULL]
  out[]
}
cor_tbl <- rbind(axis_cor(axis, "lncRNA"), axis_cor(pc_axis, "protein_coding"))
save_tsv(cor_tbl, "07_axis_vs_purity_qc_correlations.tsv")
print(cor_tbl)

# -----------------------------------------------------------------------------
# 5. Non-feature fraction versus clinical variables
# -----------------------------------------------------------------------------
wilcox_row <- function(variable, grp) {
  ok <- !is.na(grp) & !is.na(cohort$pct_noFeature)
  x  <- cohort$pct_noFeature[ok]; g <- droplevels(grp[ok]); lv <- levels(g)
  w  <- wilcox.test(x ~ g, exact = FALSE)
  data.table(variable = variable, test = "Wilcoxon rank-sum", n = sum(ok),
             group1 = lv[1], n_group1 = sum(g == lv[1]),
             median_group1 = round(median(x[g == lv[1]]), 2),
             group2 = lv[2], n_group2 = sum(g == lv[2]),
             median_group2 = round(median(x[g == lv[2]]), 2),
             rho = NA_real_, p = signif(w$p.value, 3))
}
spearman_row <- function(variable, v) {
  ok <- !is.na(v) & !is.na(cohort$pct_noFeature)
  ct <- cor.test(cohort$pct_noFeature[ok], v[ok], method = "spearman", exact = FALSE)
  data.table(variable = variable, test = "Spearman", n = sum(ok),
             group1 = NA_character_, n_group1 = NA_integer_, median_group1 = NA_real_,
             group2 = NA_character_, n_group2 = NA_integer_, median_group2 = NA_real_,
             rho = round(unname(ct$estimate), 3), p = signif(ct$p.value, 3))
}
nf_clin <- rbindlist(list(
  wilcox_row("stage_group", cohort$stage_group),
  wilcox_row("grade_group", cohort$grade_group),
  wilcox_row("M1", factor(cohort$M1, levels = c(0, 1), labels = c("M0", "M1"))),
  spearman_row("grade_num", cohort$grade_num),
  spearman_row("T_stage",   cohort$T_stage),
  spearman_row("age",       cohort$age)))
save_tsv(nf_clin, "07_noFeature_vs_clinical.tsv")
print(nf_clin)

# -----------------------------------------------------------------------------
# 6. Per-gene correlation of observed expression with the non-feature fraction
# -----------------------------------------------------------------------------
# Spearman rho as the Pearson correlation of average ranks.
gene_quality_cor <- function(net, biotype) {
  X  <- obs_expr(net)
  nf <- cohort$pct_noFeature[match(rownames(X), cohort$sample_barcode)]
  ok <- !is.na(nf)
  Xr <- apply(X[ok, , drop = FALSE], 2, rank)
  rho <- as.numeric(cor(Xr, rank(nf[ok])))
  gt <- as.data.table(net$gene_tbl)
  i  <- match(colnames(X), gt$gene_id)
  data.table(biotype = biotype, gene_id = colnames(X), gene_name = gt$gene_name[i],
             module = gt$module[i], n = sum(ok), rho_noFeature = round(rho, 4))
}
msg("Per-gene correlation with pct_noFeature ...")
gene_cor <- rbind(gene_quality_cor(nets$lnc,  "lncRNA"),
                  gene_quality_cor(nets$mrna, "protein_coding"))
save_tsv(gene_cor, "07_per_gene_quality_correlation.tsv")
gene_sum <- gene_cor[, .(n_genes = .N,
                         median_abs_rho = round(median(abs(rho_noFeature)), 3),
                         q25 = round(quantile(abs(rho_noFeature), 0.25), 3),
                         q75 = round(quantile(abs(rho_noFeature), 0.75), 3),
                         frac_abs_rho_gt_0.3 = round(mean(abs(rho_noFeature) > 0.3), 3),
                         frac_abs_rho_gt_0.5 = round(mean(abs(rho_noFeature) > 0.5), 3)),
                     by = biotype]
save_tsv(gene_sum, "07_per_gene_quality_summary.tsv")
print(gene_sum)

# -----------------------------------------------------------------------------
# 7. Nested Cox sequence for each axis, on one complete-case set
# -----------------------------------------------------------------------------
NESTED <- list(
  `Axis alone`                               = "axis_z",
  `+ age, sex, T, N, M1, ordinal grade`      = "axis_z + age + male + T_stage + N_pos + M1 + grade_num",
  `+ ESTIMATE stromal and immune scores`     = "axis_z + age + male + T_stage + N_pos + M1 + grade_num + stromal + immune",
  `+ non-feature, multimapping, log10 depth` = "axis_z + age + male + T_stage + N_pos + M1 + grade_num + stromal + immune + nf + mm + dep")
NEED <- c("os_time", "os_event", "axis_z", "age", "male", "T_stage", "N_pos", "M1",
          "grade_num", "stromal", "immune", "nf", "mm", "dep")

# VIF of the axis term from a linear model with the same terms (car::vif), or
# 1 / (1 - R2) of the axis on the other covariates if car is absent.
vif_axis <- function(rhs, data) {
  others <- setdiff(trimws(strsplit(rhs, "+", fixed = TRUE)[[1]]), "axis_z")
  if (!length(others)) return(1)
  if (requireNamespace("car", quietly = TRUE)) {
    v <- tryCatch(car::vif(lm(as.formula(paste("os_time ~", rhs)), data = data))[["axis_z"]],
                  error = function(e) NA_real_)
    if (is.finite(v)) return(v)
  }
  r2 <- summary(lm(as.formula(paste("axis_z ~", paste(others, collapse = " + "))),
                   data = data))$r.squared
  1 / (1 - r2)
}
ph_term <- function(fit, term = "axis_z")
  tryCatch(cox.zph(fit)$table[term, "p"], error = function(e) NA_real_)

nested_axis <- function(ax, label) {
  d <- merge(ax, cov, by = "sample_barcode")
  setnames(d, names(ax)[2], "axis")
  d[, axis_z := as.numeric(scale(axis))]
  dc <- d[complete.cases(d[, NEED, with = FALSE])]
  msg(label, " nested models: complete cases n = ", nrow(dc), ", events = ", sum(dc$os_event))
  rbindlist(lapply(names(NESTED), function(nm) {
    fit <- coxph(as.formula(paste("Surv(os_time, os_event) ~", NESTED[[nm]])), data = dc)
    s <- summary(fit)
    data.table(model = nm, n = s$n, events = s$nevent,
               HR = round(s$conf.int["axis_z", 1], 3),
               lo = round(s$conf.int["axis_z", 3], 3),
               hi = round(s$conf.int["axis_z", 4], 3),
               p  = signif(s$coefficients["axis_z", 5], 3),
               C  = round(unname(s$concordance[1]), 3),
               ph_p = signif(ph_term(fit), 3),
               vif_axis = round(vif_axis(NESTED[[nm]], dc), 3))
  }))
}
axis_nested <- nested_axis(axis, "lncRNA axis")
save_tsv(axis_nested, "07_axis_nested_adjustment.tsv")
save_tsv(axis_nested, "22_axis_nested_adjustment.tsv")   # read by 13_figures.R
print(axis_nested)

pc_nested <- nested_axis(pc_axis, "protein-coding PC1")
save_tsv(pc_nested, "07_pc_axis_nested_adjustment.tsv")
print(pc_nested)

# -----------------------------------------------------------------------------
# 8. The non-feature fraction as the exposure
# -----------------------------------------------------------------------------
# HR per SD of pct_noFeature, on the complete-case set of the fullest model.
# Tissue source sites with fewer than 10 patients are pooled into one stratum.
NEED_NF <- c("os_time", "os_event", "nf", "nf_log", "mm", "dep", "age", "male", "T_stage",
             "N_pos", "M1", "grade_num", "stromal", "immune", "tss")
dn <- cov[complete.cases(cov[, NEED_NF, with = FALSE])]
site_n <- table(dn$tss)
dn[, site := ifelse(tss %in% names(site_n)[site_n >= 10], tss, "other")]
msg("pct_noFeature as exposure: n = ", nrow(dn), ", events = ", sum(dn$os_event),
    ", sites (>= 10 patients) = ", sum(site_n >= 10), ", pooled = ", sum(site_n < 10))

CLIN_RHS <- "age + male + T_stage + N_pos + M1 + grade_num"
NF_MODELS <- list(
  `(a) pct_noFeature alone`                       = "nf",
  `(b) + age, sex, T, N, M1, ordinal grade`       = paste("nf +", CLIN_RHS),
  `(c) + ESTIMATE stromal and immune scores`      = paste("nf +", CLIN_RHS, "+ stromal + immune"),
  `(c) + multimapping, log10 depth`               = paste("nf +", CLIN_RHS, "+ stromal + immune + mm + dep"),
  `(d) as (c), stratified by tissue source site`  = paste("nf +", CLIN_RHS, "+ stromal + immune + strata(site)"))
nf_fit <- function(rhs) coxph(as.formula(paste("Surv(os_time, os_event) ~", rhs)), data = dn)
# The fraction is right-skewed, so HRs are given per SD on linear and log10 scales.
fit_nf_models <- function(term, label) rbindlist(lapply(names(NF_MODELS), function(nm) {
  rhs <- sub("^nf", term, NF_MODELS[[nm]])
  fit <- nf_fit(rhs); s <- summary(fit)
  data.table(exposure_scale = label, model = nm, n = s$n, events = s$nevent,
             HR = round(s$conf.int[term, 1], 3),
             lo = round(s$conf.int[term, 3], 3),
             hi = round(s$conf.int[term, 4], 3),
             p  = signif(s$coefficients[term, 5], 3),
             ph_p = signif(ph_term(fit, term), 3),
             n_strata = if (grepl("strata", rhs)) length(unique(dn$site)) else 1L,
             lrt_chisq = NA_real_, lrt_df = NA_integer_, lrt_p = NA_real_)
}))
nf_tbl <- rbind(fit_nf_models("nf", "linear, per SD"),
                fit_nf_models("nf_log", "log10, per SD"))
# Joint likelihood-ratio test of the three STAR metrics (non-feature,
# multimapping, log10 depth) added to the clinical + ESTIMATE model.
f_red  <- nf_fit(paste(CLIN_RHS, "+ stromal + immune"))
f_full <- nf_fit(paste("nf +", CLIN_RHS, "+ stromal + immune + mm + dep"))
lrt <- 2 * (f_full$loglik[2] - f_red$loglik[2])
nf_tbl[model == "(c) + multimapping, log10 depth" & exposure_scale == "linear, per SD",
       `:=`(lrt_chisq = round(lrt, 3), lrt_df = 3L,
            lrt_p = signif(pchisq(lrt, df = 3, lower.tail = FALSE), 3))]
save_tsv(nf_tbl, "07_noFeature_as_exposure.tsv")
print(nf_tbl)

# -----------------------------------------------------------------------------
# 9. Figure
# -----------------------------------------------------------------------------
d <- merge(axis, cov, by = "sample_barcode")
d[, axis_z := as.numeric(scale(lnc_axis))]
pl <- melt(d[, .(axis_z, StromalScore, ImmuneScore, ESTIMATEScore,
                 pct_noFeature, pct_multimapping)],
           id.vars = "axis_z")
g <- ggplot(pl, aes(value, axis_z)) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, colour = "#B2182B", linewidth = 0.7) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(title = "Global lncRNA axis vs tumour composition and alignment quality",
       x = NULL, y = "Global lncRNA axis (z)") +
  theme_bw()
save_fig(g, "07_axis_vs_purity_qc", 9, 5.5)

write_session_info("07_purity_qc")
banner("07 | done")
