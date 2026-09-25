# 29_normalisation_check.R: tests whether the lncRNA quality axis is an FPKM normalisation artefact
# GDC FPKM divides counts by the protein-coding read total, so libraries with more
# non-feature, intronic or lncRNA reads get larger lncRNA FPKMs. From raw STAR counts
# (528 libraries): TMM and DESeq2 size factors, compartment read shares, PC1 of the
# 3,442 production lncRNAs under five normalisations, and the production module
# eigengenes refitted on TMM and compared with the FPKM arms of stage 12 by a paired
# patient bootstrap over a seed sweep.
# Inputs: caches dataset, networks, star_qc, estimate_scores, lnc_global_axis, gdc_file_map,
#   and results/12_technical_adjustment_sensitivity.tsv. Run after stage 12.
# Outputs: results/29_*.tsv, caches 29_counts_unstranded.rds and 29_bootstrap_seed_sweep.rds.
# Spearman throughout (the non-feature fraction is skewed). HRs are per SD over the 511
# network samples. BH within the family named in each table.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(edgeR); library(DESeq2)
})
banner("29 | Is the lncRNA quality axis an FPKM normalisation artefact?")
set.seed(SEED)
t_start <- Sys.time()

ds     <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets   <- readRDS(file.path(CACHE_DIR, "networks.rds"))
qc     <- as.data.table(readRDS(file.path(CACHE_DIR, "star_qc.rds")))
est    <- as.data.table(readRDS(file.path(CACHE_DIR, "estimate_scores.rds")))
axis   <- as.data.table(readRDS(file.path(CACHE_DIR, "lnc_global_axis.rds")))
fmap   <- as.data.table(readRDS(file.path(CACHE_DIR, "gdc_file_map.rds")))
cohort <- as.data.table(ds$cohort_full)
stopifnot(!anyDuplicated(cohort$sample_barcode), !anyDuplicated(cohort$file_id))
msg("Discovery cohort: ", nrow(cohort), " patients, ", sum(cohort$os_event), " deaths")

spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(c(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}

# ---- 1. raw unstranded counts, all genes, 528 libraries ----
banner("1 | Raw STAR counts for the 528 discovery libraries")
counts_rds <- file.path(CACHE_DIR, "29_counts_unstranded.rds")
cnt <- if (file.exists(counts_rds)) readRDS(counts_rds) else NULL
if (!is.null(cnt) && !identical(sort(colnames(cnt$counts)), sort(cohort$sample_barcode))) {
  msg("Cached counts do not cover the cohort; re-reading"); cnt <- NULL
}
if (is.null(cnt)) {
  fi <- cohort[, .(sample_barcode, file_id)]
  fi[, file_name := fmap$file_name[match(file_id, fmap$file_id)]]
  stopifnot(!anyNA(fi$file_name))
  fi[, path := winlong(file.path(EXPR_DIR, file_id, file_name))]
  stopifnot(all(file.exists(fi$path)))
  msg("Reading ", nrow(fi), " STAR count files (unstranded column) ...")
  first <- fread(fi$path[1], skip = 1, showProgress = FALSE,
                 select = c("gene_id", "gene_name", "gene_type", "unstranded"))
  is_gene  <- grepl("^ENSG", first$gene_id)
  gene_ann <- first[is_gene, .(gene_id, gene_name, gene_type)]
  qc_rows  <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
  counts   <- matrix(NA_integer_, nrow(gene_ann), nrow(fi),
                     dimnames = list(gene_ann$gene_id, fi$sample_barcode))
  qmat     <- matrix(NA_real_, nrow(fi), length(qc_rows),
                     dimnames = list(fi$sample_barcode, qc_rows))
  for (i in seq_len(nrow(fi))) {
    d <- fread(fi$path[i], skip = 1, showProgress = FALSE, select = c("gene_id", "unstranded"))
    v <- setNames(as.numeric(d$unstranded), d$gene_id)
    qmat[i, ] <- v[qc_rows]
    d <- d[grepl("^ENSG", gene_id)]
    if (!identical(d$gene_id, gene_ann$gene_id)) d <- d[match(gene_ann$gene_id, d$gene_id)]
    counts[, i] <- as.integer(d$unstranded)
    if (i %% 100 == 0) msg("  ", i, " / ", nrow(fi))
  }
  stopifnot(!anyNA(counts))
  cnt <- list(counts = counts, gene_ann = gene_ann,
              sample_info = cbind(fi[, .(sample_barcode, file_id, file_name)], as.data.table(qmat)),
              read_at = Sys.time())
  saveRDS(cnt, counts_rds)
  msg("Counts cached: ", nrow(counts), " genes x ", ncol(counts), " libraries")
} else msg("Counts from cache: ", nrow(cnt$counts), " genes x ", ncol(cnt$counts), " libraries")

counts   <- cnt$counts[, cohort$sample_barcode]
gene_ann <- cnt$gene_ann
si       <- cnt$sample_info[match(cohort$sample_barcode, sample_barcode)]
stopifnot(identical(colnames(counts), cohort$sample_barcode), identical(si$file_id, cohort$file_id))

# Check the counts against star_qc (stage 07), which derives from the same files.
lib <- colSums(counts)
q   <- qc[match(cohort$sample_barcode, sample_barcode)]
tot <- lib + si$N_unmapped + si$N_multimapping + si$N_noFeature + si$N_ambiguous
verif <- data.table(sample_barcode = cohort$sample_barcode, file_id = cohort$file_id,
                    assigned_reads_counts = lib, assigned_reads_star_qc = q$assigned_reads,
                    diff_assigned = lib - q$assigned_reads,
                    pct_noFeature_counts = 100 * si$N_noFeature / tot,
                    pct_noFeature_star_qc = q$pct_noFeature,
                    diff_pct_noFeature = 100 * si$N_noFeature / tot - q$pct_noFeature)
save_tsv(verif, "29_count_verification.tsv")
msg("Assigned reads: max |difference| vs star_qc = ", max(abs(verif$diff_assigned)),
    " reads; non-feature fraction: max |difference| = ", signif(max(abs(verif$diff_pct_noFeature)), 3), " points")
stopifnot(max(abs(verif$diff_assigned)) == 0, max(abs(verif$diff_pct_noFeature)) < 1e-6)

# ---- 2. size factors and compartment read shares versus quality ----
banner("2 | TMM and median-of-ratios size factors; lncRNA read share")
dge <- DGEList(counts)
dge <- calcNormFactors(dge, method = "TMM")
tmm_nf   <- dge$samples$norm.factors
sf_deseq <- estimateSizeFactorsForMatrix(counts)
# Composition components: the DESeq2 factor divided by relative library size,
# and the TMM factor itself (edgeR multiplies it by depth).
lib_rel      <- lib / exp(mean(log(lib)))
sf_deseq_rel <- sf_deseq / lib_rel

is_lnc  <- gene_ann$gene_type == "lncRNA"
is_pc   <- gene_ann$gene_type == "protein_coding"
prod_g  <- colnames(nets$lnc$expr)                    # the 3,442 production lncRNAs
stopifnot(all(prod_g %in% rownames(counts)))
lnc_share_all  <- colSums(counts[is_lnc, ]) / lib
lnc_share_prod <- colSums(counts[prod_g, ]) / lib
pc_share       <- colSums(counts[is_pc, ]) / lib
# GDC FPKM = count x 1e9 / (length x protein-coding reads). Relative to a CPM
# on all assigned reads, every gene's FPKM is multiplied by assigned / coding.
fpkm_vs_cpm_scale <- lib / colSums(counts[is_pc, ])

per_sample <- data.table(
  sample_barcode = cohort$sample_barcode, assigned_reads = lib,
  pct_noFeature = q$pct_noFeature, pct_multimapping = q$pct_multimapping,
  tmm_norm_factor = tmm_nf, tmm_effective_libsize = lib * tmm_nf,
  deseq_size_factor = sf_deseq, deseq_size_factor_rel_libsize = sf_deseq_rel,
  lnc_read_share_all = lnc_share_all, lnc_read_share_production = lnc_share_prod,
  protein_coding_read_share = pc_share, fpkm_vs_cpm_scale = fpkm_vs_cpm_scale,
  lnc_axis = axis$lnc_axis[match(cohort$sample_barcode, axis$sample_barcode)])
save_tsv(per_sample, "29_sample_normalisation_factors.tsv")

QUANT <- c(
  tmm_norm_factor               = "edgeR TMM normalisation factor (all genes)",
  tmm_effective_libsize         = "TMM effective library size (assigned reads x TMM factor)",
  deseq_size_factor             = "DESeq2 median-of-ratios size factor (all genes)",
  deseq_size_factor_rel_libsize = "DESeq2 size factor relative to library size (composition component)",
  assigned_reads                = "assigned reads (library size)",
  lnc_read_share_all            = "lncRNA read share, all 16,901 lncRNA genes / assigned reads",
  lnc_read_share_production     = "lncRNA read share, 3,442 production lncRNAs / assigned reads",
  protein_coding_read_share     = "protein-coding read share / assigned reads",
  fpkm_vs_cpm_scale             = "assigned reads / protein-coding reads (per-sample FPKM scaling relative to CPM)")
# fpkm_vs_cpm_scale is exactly 1 / protein_coding_read_share, so those two rows
# are one measurement (mirrored rho, identical p). The assigned_reads versus
# log depth cell is a self-correlation and is set to NA.
stopifnot(max(abs(per_sample$fpkm_vs_cpm_scale - 1 / per_sample$protein_coding_read_share)) < 1e-12)
NOTE <- c(
  protein_coding_read_share = "exactly 1 / fpkm_vs_cpm_scale; that row's correlations are the sign-flipped mirror of this one's, not an independent test",
  fpkm_vs_cpm_scale         = "exactly 1 / protein_coding_read_share; this row and that one are one measurement reported twice, with mirrored rho and identical p",
  assigned_reads            = "log_depth correlation not reported: this quantity IS the library size, so it would be a correlation with itself")
size_tbl <- rbindlist(lapply(names(QUANT), function(v) {
  x <- per_sample[[v]]
  s_nf <- spearman(x, per_sample$pct_noFeature)
  s_mm <- spearman(x, per_sample$pct_multimapping)
  s_dp <- if (v == "assigned_reads") c(n = NA_real_, rho = NA_real_, p = NA_real_)
          else spearman(x, log10(per_sample$assigned_reads))
  s_ax <- spearman(x, per_sample$lnc_axis)
  data.table(quantity = v, description = QUANT[[v]], n = s_nf[["n"]],
             median = signif(median(x), 4), iqr_lo = signif(quantile(x, 0.25), 4),
             iqr_hi = signif(quantile(x, 0.75), 4),
             rho_noFeature = round(s_nf[["rho"]], 3), p_noFeature = signif(s_nf[["p"]], 3),
             rho_multimap  = round(s_mm[["rho"]], 3), p_multimap  = signif(s_mm[["p"]], 3),
             rho_log_depth = round(s_dp[["rho"]], 3), p_log_depth = signif(s_dp[["p"]], 3),
             n_axis = s_ax[["n"]], rho_axis_fpkm = round(s_ax[["rho"]], 3), p_axis_fpkm = signif(s_ax[["p"]], 3),
             note = if (v %in% names(NOTE)) NOTE[[v]] else "")
}))
# BH over all Spearman tests in the table as one family. The four duplicated
# tests make it slightly conservative.
SIZE_FDR_FAMILY <- paste(
  "all Spearman tests in 29_size_factors_vs_quality.tsv: 9 quantities x 4 quality",
  "metrics less the assigned-reads/log-depth self-correlation; the",
  "protein_coding_read_share and fpkm_vs_cpm_scale rows are the same measurement",
  "twice, so the family contains 4 duplicated tests")
.p_size <- c(size_tbl$p_noFeature, size_tbl$p_multimap, size_tbl$p_log_depth,
             size_tbl$p_axis_fpkm)
.f_size <- rep(NA_real_, length(.p_size))
.f_size[is.finite(.p_size)] <- p.adjust(.p_size[is.finite(.p_size)], "BH")
.n_size <- nrow(size_tbl)
size_tbl[, `:=`(fdr_noFeature = signif(.f_size[1:.n_size], 3),
                fdr_multimap  = signif(.f_size[.n_size + 1:.n_size], 3),
                fdr_log_depth = signif(.f_size[2 * .n_size + 1:.n_size], 3),
                fdr_axis_fpkm = signif(.f_size[3 * .n_size + 1:.n_size], 3),
                fdr_family    = SIZE_FDR_FAMILY)]
setcolorder(size_tbl, c("quantity", "description", "n", "median", "iqr_lo", "iqr_hi",
                        "rho_noFeature", "p_noFeature", "fdr_noFeature",
                        "rho_multimap", "p_multimap", "fdr_multimap",
                        "rho_log_depth", "p_log_depth", "fdr_log_depth",
                        "n_axis", "rho_axis_fpkm", "p_axis_fpkm", "fdr_axis_fpkm",
                        "fdr_family", "note"))
msg("Size-factor correlation family: ", sum(is.finite(.p_size)), " Spearman tests, ",
    sum(.p_size < 0.05, na.rm = TRUE), " with raw p < 0.05, ",
    sum(.f_size < FDR_ALPHA, na.rm = TRUE), " at BH < ", FDR_ALPHA)
save_tsv(size_tbl, "29_size_factors_vs_quality.tsv")
print(size_tbl[, .(quantity, n, median, rho_noFeature, p_noFeature, rho_log_depth, rho_axis_fpkm)], row.names = FALSE)

# ---- 3. the production lncRNA matrix under five normalisations ----
banner("3 | PC1 of the lncRNA matrix under count-based and quantile normalisation")
S <- nets$lnc$samples                                  # 511 network samples
stopifnot(all(S %in% cohort$sample_barcode), all(S %in% axis$sample_barcode))
E_fpkm <- t(ds$lnc$expr[prod_g, S])
stopifnot(max(abs(E_fpkm - obs_expr(nets$lnc)[S, prod_g])) < 1e-8)

# log2 CPM by library size only, and with TMM factors from all genes of all
# 528 libraries. cpm(log = TRUE) adds a library-scaled prior count.
dge_lib <- DGEList(counts)                             # norm.factors = 1
E_cpm <- t(cpm(dge_lib, log = TRUE, prior.count = 1)[prod_g, S])
E_tmm <- t(cpm(dge,     log = TRUE, prior.count = 1)[prod_g, S])

# Blind DESeq2 VST with size factors and a parametric dispersion trend from the 3,442
# lncRNAs. The steps are written out so the log shows the dispersion fit used (DESeq2
# silently falls back to a local fit). Equivalent to varianceStabilizingTransformation(
# blind = TRUE, fitType = "parametric").
msg("DESeq2 variance-stabilising transformation of the lncRNA counts ...")
.t_vst <- Sys.time()
dds_lnc <- DESeqDataSetFromMatrix(counts[prod_g, S],
                                  S4Vectors::DataFrame(row.names = S), design = ~ 1)
dds_lnc <- estimateSizeFactors(dds_lnc)
dds_lnc <- estimateDispersionsGeneEst(dds_lnc, quiet = TRUE)
dds_lnc <- estimateDispersionsFit(dds_lnc, fitType = "parametric", quiet = TRUE)
E_vst   <- t(getVarianceStabilizedData(dds_lnc))
msg("  VST done in ", round(as.numeric(difftime(Sys.time(), .t_vst, units = "secs")), 1),
    " s; dispersion trend fitted by '", attr(dispersionFunction(dds_lnc), "fitType"),
    "'; lncRNA-only size factors span ",
    paste(signif(range(sizeFactors(dds_lnc)), 3), collapse = " to "))
stopifnot(identical(dimnames(E_vst), dimnames(E_fpkm)))

# Quantile normalisation gives every sample the same distribution, so no
# per-sample offset or scale remains.
quantile_normalise <- function(M) {                    # genes x samples
  if (requireNamespace("limma", quietly = TRUE)) return(limma::normalizeQuantiles(M))
  r  <- apply(M, 2, rank, ties.method = "average")
  s  <- apply(M, 2, sort)
  mu <- rowMeans(s)
  apply(r, 2, function(rk) approx(seq_along(mu), mu, xout = rk, rule = 2)$y)
}
E_qn <- t(quantile_normalise(ds$lnc$expr[prod_g, S]))
dimnames(E_qn) <- dimnames(E_fpkm)

qS   <- qc[match(S, sample_barcode)]
nf_S <- qS$pct_noFeature; mm_S <- qS$pct_multimapping; dp_S <- log10(qS$assigned_reads)
ax_S <- axis$lnc_axis[match(S, axis$sample_barcode)]

pc_stats <- function(E, label) {
  x  <- scale(E, center = TRUE, scale = FALSE)
  sv <- svd(x, nu = 2, nv = 2)
  share <- sv$d^2 / sum(sv$d^2)
  mexp <- rowMeans(E)
  # PC1 is signed against the cached FPKM axis, since quantile normalisation makes the
  # per-sample mean constant. Rho with the mean and its SD are reported as diagnostics of
  # that alternative anchor. PC2 has no anchor, so only |rho| is reported.
  pc1 <- sv$u[, 1] * sv$d[1]; if (stats::cor(pc1, ax_S) < 0) pc1 <- -pc1
  pc2 <- sv$u[, 2] * sv$d[2]
  s_nf <- spearman(pc1, nf_S); s_mm <- spearman(pc1, mm_S); s_dp <- spearman(pc1, dp_S)
  s_mu <- spearman(pc1, mexp); s_ax <- spearman(pc1, ax_S); s_nf2 <- spearman(pc2, nf_S)
  # per-gene Spearman rho with the non-feature fraction (Pearson on ranks)
  rho_g <- as.numeric(stats::cor(apply(E, 2, rank), rank(nf_S)))
  data.table(matrix = label, n = nrow(E), n_genes = ncol(E),
             pc1_sign_anchor = "cached FPKM axis (lnc_global_axis.rds)",
             pc1_var_share = round(share[1], 4), pc2_var_share = round(share[2], 4),
             rho_noFeature = round(s_nf[["rho"]], 3), p_noFeature = signif(s_nf[["p"]], 3),
             rho_multimap = round(s_mm[["rho"]], 3), rho_log_depth = round(s_dp[["rho"]], 3),
             pearson_r_axis_fpkm = round(stats::cor(pc1, ax_S), 3),
             rho_axis_fpkm = round(s_ax[["rho"]], 3), p_axis_fpkm = signif(s_ax[["p"]], 3),
             pearson_r_mean_expr = round(stats::cor(pc1, mexp), 3),
             rho_mean_expr = round(s_mu[["rho"]], 3),
             abs_pearson_r_mean_expr = round(abs(stats::cor(pc1, mexp)), 3),
             sd_mean_expr = signif(sd(mexp), 3),
             mean_expr_range = paste(signif(range(mexp), 6), collapse = " to "),
             mean_expr_anchor_determinate = abs(stats::cor(pc1, mexp)) > 0.5,
             pc2_abs_rho_noFeature = round(abs(s_nf2[["rho"]]), 3),
             gene_median_abs_rho_noFeature = round(median(abs(rho_g)), 3),
             gene_frac_abs_rho_gt_0.3 = round(mean(abs(rho_g) > 0.3), 3),
             gene_frac_rho_positive = round(mean(rho_g > 0), 3))
}
msg("PCA under each normalisation ...")
pc_tbl <- rbindlist(list(
  pc_stats(E_fpkm, "log2(FPKM+1) as used"),
  pc_stats(E_cpm,  "log2 CPM, library size only"),
  pc_stats(E_tmm,  "log2 CPM TMM (edgeR, all-gene factors)"),
  pc_stats(E_vst,  "DESeq2 VST on lncRNA counts (blind)"),
  pc_stats(E_qn,   "quantile-normalised log2(FPKM+1)")))
save_tsv(pc_tbl, "29_normalisation_pc1.tsv")
print(pc_tbl[, .(matrix, n, pc1_var_share, rho_noFeature, rho_axis_fpkm,
                 abs_pearson_r_mean_expr, sd_mean_expr, mean_expr_anchor_determinate,
                 gene_median_abs_rho_noFeature)], row.names = FALSE)
stopifnot(abs(pc_tbl$rho_axis_fpkm[1] - 1) < 1e-3)   # the FPKM arm is the cached axis
if (any(!pc_tbl$mean_expr_anchor_determinate))
  msg("NOTE: the per-sample mean is degenerate as a sign anchor in: ",
      paste(pc_tbl$matrix[!pc_tbl$mean_expr_anchor_determinate], collapse = "; "),
      " -- |r(PC1, per-sample mean)| = ",
      paste(pc_tbl$abs_pearson_r_mean_expr[!pc_tbl$mean_expr_anchor_determinate], collapse = ", "),
      "; sd(per-sample mean) = ",
      paste(pc_tbl$sd_mean_expr[!pc_tbl$mean_expr_anchor_determinate], collapse = ", "),
      ". PC1 there is signed against the FPKM axis, not against the per-sample mean.")

# ---- 4. fixed production modules on the TMM matrix versus the FPKM arms of 12 ----
banner("4 | Module eigengenes on the TMM matrix, with and without residualisation")
gt   <- as.data.table(nets$lnc$gene_tbl)
cov_S <- tech_covariates(qS)
E_tmm_res  <- apply_technical(E_tmm,  fit_technical(E_tmm,  cov_S), cov_S)
E_fpkm_res <- nets$lnc$expr[S, prod_g]                 # the production residualised matrix
.d_res <- max(abs(E_fpkm_res - apply_technical(E_fpkm, fit_technical(E_fpkm, cov_S), cov_S)))
msg("Production residualised matrix vs fit_technical/apply_technical here: max |difference| = ", signif(.d_res, 2))
if (.d_res > 1e-6) warning("production residualisation not reproduced exactly (see message above)")

LOAD_LNC  <- { L <- discovery_loadings(nets); L[grep("^lnc_", names(L))] }
scores_on <- function(E) score_modules(E, fit_module_loadings(E, gt, "lnc_ME"))
SC <- list(
  fpkm_observed     = scores_on(E_fpkm),
  fpkm_residualised = score_modules(E_fpkm_res, LOAD_LNC),
  tmm_observed      = scores_on(E_tmm),
  tmm_residualised  = scores_on(E_tmm_res))
SC$fpkm_observed_STAR <- SC$fpkm_observed
SC$tmm_observed_STAR  <- SC$tmm_observed
mods <- Reduce(intersect, lapply(SC, colnames))
stopifnot(length(mods) == nets$lnc$n_modules)

cl  <- cohort[match(S, sample_barcode)]
eS  <- est[match(S, sample_barcode)]
X   <- clinical_design(cl, "augmented", eS, qS)
STAR_TERMS <- c("noFeature", "multimap", "libsize_z")
keep <- complete.cases(X[, c(CLIN_TERMS, STAR_TERMS), drop = FALSE]) &
        is.finite(cl$os_time) & is.finite(cl$os_event)
y <- Surv(cl$os_time[keep], cl$os_event[keep])
msg("Cox analysis set: n = ", sum(keep), ", events = ", sum(cl$os_event[keep]),
    " (network samples with complete clinical covariates)")

ARM_LABEL <- c(
  fpkm_observed      = "log2(FPKM+1), observed eigengene",
  fpkm_observed_STAR = "log2(FPKM+1), observed eigengene + STAR covariates",
  fpkm_residualised  = "log2(FPKM+1), residualised eigengene (production)",
  tmm_observed       = "log2 CPM TMM, observed eigengene",
  tmm_observed_STAR  = "log2 CPM TMM, observed eigengene + STAR covariates",
  tmm_residualised   = "log2 CPM TMM, residualised eigengene")
arms <- rbindlist(lapply(names(ARM_LABEL), function(arm) {
  covs <- c(CLIN_TERMS, if (grepl("_STAR$", arm)) STAR_TERMS)
  r <- rbindlist(lapply(mods, function(m) {
    d <- data.frame(ME = SC[[arm]][keep, m], X[keep, covs, drop = FALSE])
    s <- summary(coxph(y ~ ., data = d))
    data.table(biotype = "lncRNA", module = sub("^lnc_ME", "", m), arm = arm,
               description = ARM_LABEL[[arm]], covariate_set = "clinical",
               n = s$n, events = s$nevent,
               HR = s$conf.int["ME", 1], lo = s$conf.int["ME", 3], hi = s$conf.int["ME", 4],
               p = s$coefficients["ME", 5])
  }))
  r[, fdr := p.adjust(p, "BH")][]
}))
arms[, n_genes := gt[module != "grey", .N, by = module][match(arms$module, module), N]]
save_tsv(arms[, .(biotype, module, n_genes, arm, description, covariate_set, n, events,
                  HR = round(HR, 3), lo = round(lo, 3), hi = round(hi, 3),
                  p = signif(p, 3), fdr = signif(fdr, 3))],
         "29_normalisation_module_arms.tsv")

# The FPKM arms refitted here should reproduce stage 12.
s12 <- fread(file.path(RESULTS_DIR, "12_technical_adjustment_sensitivity.tsv"))
s12 <- s12[biotype == "lncRNA" & covariate_set == "clinical"]
stopifnot(setequal(s12$module, arms[arm == "fpkm_observed", module]))
chk <- merge(s12[, .(module, HR_unadjusted, HR_observed_STAR, HR_adjusted, n, events)],
             dcast(arms[grepl("^fpkm", arm), .(module, arm, HR)], module ~ arm, value.var = "HR"),
             by = "module")
msg("Reproduction of 12 (max |HR difference|): observed ",
    signif(max(abs(chk$HR_unadjusted - chk$fpkm_observed)), 2), "; + STAR covariates ",
    signif(max(abs(chk$HR_observed_STAR - chk$fpkm_observed_STAR)), 2), "; residualised ",
    signif(max(abs(chk$HR_adjusted - chk$fpkm_residualised)), 2),
    "; n ", unique(chk$n), " vs ", sum(keep))
if (max(abs(chk$HR_unadjusted - chk$fpkm_observed), abs(chk$HR_adjusted - chk$fpkm_residualised)) > 0.005)
  warning("FPKM arms refitted here differ from 12_technical_adjustment_sensitivity.tsv by more than 0.005 in HR")

# Eigengene agreement across normalisations (same samples, same membership).
me_cor <- rbindlist(lapply(mods, function(m) data.table(
  module = sub("^lnc_ME", "", m),
  r_ME_observed     = round(stats::cor(SC$fpkm_observed[, m],     SC$tmm_observed[, m]), 3),
  r_ME_residualised = round(stats::cor(SC$fpkm_residualised[, m], SC$tmm_residualised[, m]), 3))))

# ---- 4b. paired patient bootstrap of log HR(TMM) - log HR(FPKM) ----
# Both eigengenes are on the same patients, so the difference is estimated by resampling
# patients (same draws for both arms and all modules) and refitting both Cox models.
# Predictors keep full-cohort standardisation, as in paired_boot_delta_c. coxph.fit() is
# called directly for speed, checked against coxph() first.
# boot_p = 2 * min(P(delta <= 0), P(delta >= 0)), and 0 means p < 1 / BOOT_B.
# BH over the 16 modules within each arm. The contrast is repeated over N_SEED_SWEEP
# seeds at full BOOT_B, because the count of intervals excluding zero is itself a
# Monte-Carlo quantity.
banner("4b | Paired patient bootstrap of the TMM-minus-FPKM difference in log HR")
tm_k <- cl$os_time[keep]; ev_k <- cl$os_event[keep]
Xk   <- X[keep, CLIN_TERMS, drop = FALSE]
n_k  <- sum(keep)
.ctl <- survival::coxph.control()
.all <- seq_len(n_k)
MIN_BOOT_EVENTS <- 5             # same floor as paired_boot_delta_c
N_SEED_SWEEP    <- 50L           # seeds in the stability sweep (SEED is the first)
SWEEP_SEEDS     <- SEED + seq_len(N_SEED_SWEEP) - 1L
BOOT_FDR_FAMILY <- paste(
  "paired patient bootstrap of the TMM-minus-FPKM difference in log hazard ratio;",
  "16 lncRNA modules, Benjamini-Hochberg within arm (observed; residualised)")

.beta1 <- function(v)
  unname(survival::coxph.fit(cbind(ME = v, Xk), Surv(tm_k, ev_k), strata = NULL,
                             offset = NULL, init = NULL, control = .ctl,
                             weights = NULL, method = "efron",
                             rownames = NULL)$coefficients[1])
.chk <- max(abs(vapply(mods, function(m) exp(.beta1(SC$tmm_observed[keep, m])), numeric(1)) -
                arms[arm == "tmm_observed"][match(sub("^lnc_ME", "", mods), module), HR]))
msg("coxph.fit vs coxph point estimates (tmm_observed HR): max |difference| = ", signif(.chk, 2))
stopifnot(.chk < 1e-8)

# Per-module score vectors on the analysis set, passed to the sweep workers.
VMOD <- lapply(mods, function(m) list(
  tmm_observed      = as.numeric(SC$tmm_observed[keep, m]),
  fpkm_observed     = as.numeric(SC$fpkm_observed[keep, m]),
  tmm_residualised  = as.numeric(SC$tmm_residualised[keep, m]),
  fpkm_residualised = as.numeric(SC$fpkm_residualised[keep, m])))
names(VMOD) <- sub("^lnc_ME", "", mods)

# One seed of the contrast. Returns one row per module and arm with the
# within-arm BH adjustment, plus the fit and warning counts.
boot_one_seed <- function(seed, B = BOOT_B) {
  set.seed(seed)
  idx  <- lapply(seq_len(B), function(k) sample.int(n_k, n_k, replace = TRUE))
  nfit <- 0L; nwarn <- 0L
  bfit <- function(v, i) {
    nfit <<- nfit + 1L
    withCallingHandlers(
      unname(survival::coxph.fit(cbind(ME = v[i], Xk[i, , drop = FALSE]),
                                 survival::Surv(tm_k[i], ev_k[i]), strata = NULL,
                                 offset = NULL, init = NULL, control = .ctl,
                                 weights = NULL, method = "efron",
                                 rownames = NULL)$coefficients[1]),
      warning = function(w) { nwarn <<- nwarn + 1L; invokeRestart("muffleWarning") })
  }
  pair <- function(v_tmm, v_fpkm) {
    d <- vapply(idx, function(i) {
      if (sum(ev_k[i]) < MIN_BOOT_EVENTS) return(NA_real_)
      tryCatch(bfit(v_tmm, i) - bfit(v_fpkm, i), error = function(e) NA_real_)
    }, numeric(1))
    d <- d[is.finite(d)]
    c(delta = bfit(v_tmm, .all) - bfit(v_fpkm, .all),
      lo = unname(quantile(d, 0.025)), hi = unname(quantile(d, 0.975)),
      p_boot = 2 * min(mean(d <= 0), mean(d >= 0)), n_boot = length(d))
  }
  r <- data.table::rbindlist(lapply(names(VMOD), function(m) {
    o <- pair(VMOD[[m]]$tmm_observed,     VMOD[[m]]$fpkm_observed)
    s <- pair(VMOD[[m]]$tmm_residualised, VMOD[[m]]$fpkm_residualised)
    data.table::data.table(
      module = m, arm = c("observed", "residualised"), seed = seed, boot_B = B,
      delta  = c(o[["delta"]],  s[["delta"]]),
      lo     = c(o[["lo"]],     s[["lo"]]),
      hi     = c(o[["hi"]],     s[["hi"]]),
      p_boot = c(o[["p_boot"]], s[["p_boot"]]),
      n_boot = c(o[["n_boot"]], s[["n_boot"]]))
  }))
  r[, ci_excludes_0 := lo > 0 | hi < 0]
  r[, fdr := stats::p.adjust(p_boot, "BH"), by = arm]
  r[, fdr_sig := fdr < FDR_ALPHA]
  list(tbl = r[], n_fit = nfit, n_warn = nwarn)
}

# The sweep is cached with a signature of its inputs, and any change invalidates it.
sweep_rds <- file.path(CACHE_DIR, "29_bootstrap_seed_sweep.rds")
sweep_sig <- list(seeds = SWEEP_SEEDS, B = BOOT_B, n = n_k, modules = names(VMOD),
                  min_events = MIN_BOOT_EVENTS,
                  checksum = c(sum(tm_k), sum(ev_k), sum(Xk),
                               vapply(VMOD, function(z) sum(unlist(z)), numeric(1))))
sw <- if (file.exists(sweep_rds)) readRDS(sweep_rds) else NULL
if (!is.null(sw) && !isTRUE(all.equal(sw$sig, sweep_sig, tolerance = 1e-8))) {
  msg("Cached seed sweep does not match the current inputs; recomputing"); sw <- NULL
}
if (is.null(sw)) {
  n_work <- max(1L, min(10L, N_THREADS, N_SEED_SWEEP))
  msg("Seed-stability sweep: ", N_SEED_SWEEP, " seeds x ", BOOT_B, " draws x ",
      length(VMOD), " modules x 2 contrasts, on ", n_work, " worker(s) ...")
  .t_sw <- Sys.time()
  res <- NULL
  if (n_work > 1L) {
    cls <- try(parallel::makePSOCKcluster(n_work), silent = TRUE)
    if (inherits(cls, "try-error")) {
      warning("could not start a PSOCK cluster; running the seed sweep serially")
    } else {
      parallel::clusterEvalQ(cls, suppressPackageStartupMessages({
        library(data.table); library(survival) }))
      parallel::clusterExport(cls,
        c("VMOD", "Xk", "tm_k", "ev_k", "n_k", ".ctl", ".all", "BOOT_B",
          "FDR_ALPHA", "MIN_BOOT_EVENTS", "boot_one_seed"),
        envir = environment())
      res <- try(parallel::parLapply(cls, SWEEP_SEEDS, boot_one_seed), silent = TRUE)
      try(parallel::stopCluster(cls), silent = TRUE)
      if (inherits(res, "try-error")) {
        warning("parallel seed sweep failed (", as.character(res), "); running serially")
        res <- NULL
      }
    }
  }
  if (is.null(res)) res <- lapply(SWEEP_SEEDS, boot_one_seed)
  sw <- list(sig = sweep_sig, res = res, n_workers = n_work,
             minutes = as.numeric(difftime(Sys.time(), .t_sw, units = "mins")))
  saveRDS(sw, sweep_rds)
}
msg("Seed sweep: ", length(sw$res), " seeds on ", sw$n_workers, " worker(s) in ",
    round(sw$minutes, 1), " min (", length(sw$res) * length(VMOD) * 2L * (2L * BOOT_B + 2L),
    " Cox fits in total)")

# Per-module intervals are reported at the primary seed, the first of the sweep.
stopifnot(identical(unique(sw$res[[1]]$tbl$seed), SEED))
boot_raw <- copy(sw$res[[1]]$tbl)
.nfit    <- sw$res[[1]]$n_fit
.nwarn   <- sw$res[[1]]$n_warn
bo <- boot_raw[arm == "observed"][match(names(VMOD), module)]
br <- boot_raw[arm == "residualised"][match(names(VMOD), module)]
boot_tbl <- data.table(
  module = bo$module,
  boot_delta_logHR_observed = round(bo$delta, 3),
  boot_lo_observed = round(bo$lo, 3), boot_hi_observed = round(bo$hi, 3),
  boot_p_observed = signif(bo$p_boot, 3), boot_fdr_observed = signif(bo$fdr, 3),
  boot_n_observed = bo$n_boot,
  boot_delta_logHR_residualised = round(br$delta, 3),
  boot_lo_residualised = round(br$lo, 3), boot_hi_residualised = round(br$hi, 3),
  boot_p_residualised = signif(br$p_boot, 3), boot_fdr_residualised = signif(br$fdr, 3),
  boot_n_residualised = br$n_boot,
  boot_ci_excludes_0_observed     = bo$ci_excludes_0,
  boot_ci_excludes_0_residualised = br$ci_excludes_0,
  boot_sig_fdr_observed           = bo$fdr_sig,
  boot_sig_fdr_residualised       = br$fdr_sig,
  boot_fdr_family                 = BOOT_FDR_FAMILY)
msg("Paired bootstrap at the primary seed ", SEED, ": ", BOOT_B, " draws x ",
    length(VMOD), " modules x 2 contrasts (", .nfit, " Cox fits)")
msg("Convergence warnings during the bootstrap: ", .nwarn, " of ", .nfit,
    " fits (monotone likelihood in a NUISANCE covariate -- N_pos separates in a few",
    " resamples; the module coefficient is still estimated, the draw is retained,",
    " and both arms of a draw share the same resample)")
msg("Observed arm: median difference in log HR (TMM - FPKM) = ",
    round(median(boot_tbl$boot_delta_logHR_observed), 3), "; ",
    sum(boot_tbl$boot_sig_fdr_observed), " of ", nrow(boot_tbl),
    " modules differ at BH < ", FDR_ALPHA, " within the arm (",
    sum(boot_tbl$boot_ci_excludes_0_observed),
    " before the multiplicity adjustment), ",
    sum(boot_tbl$boot_sig_fdr_observed & boot_tbl$boot_delta_logHR_observed < 0),
    " of them negative. Residualised arm: median ",
    round(median(boot_tbl$boot_delta_logHR_residualised), 3), "; ",
    sum(boot_tbl$boot_sig_fdr_residualised), " at BH < ", FDR_ALPHA, " (",
    sum(boot_tbl$boot_ci_excludes_0_residualised), " before adjustment).")
print(boot_tbl[order(boot_delta_logHR_observed),
               .(module, boot_delta_logHR_observed, boot_lo_observed, boot_hi_observed,
                 boot_p_observed, boot_fdr_observed, boot_sig_fdr_observed,
                 boot_delta_logHR_residualised, boot_p_residualised,
                 boot_fdr_residualised)], row.names = FALSE)

# ---- 4c. seed stability, per seed and per module ----
banner("4c | Stability of the bootstrap count across seeds")
all_sw <- rbindlist(lapply(sw$res, `[[`, "tbl"))
# The point estimate does not depend on the draws, so it must match across seeds.
stopifnot(all_sw[, .(u = uniqueN(round(delta, 10))), by = .(module, arm)][, all(u == 1L)])

seed_tbl <- all_sw[, .(boot_B = BOOT_B, n_modules = .N,
                       n_ci_excludes_0 = sum(ci_excludes_0),
                       n_fdr_sig = sum(fdr_sig),
                       modules_ci_excludes_0 = paste(sort(module[ci_excludes_0]), collapse = ";"),
                       modules_fdr_sig = paste(sort(module[fdr_sig]), collapse = ";")),
                   by = .(arm, seed)]
seed_tbl[, is_primary_seed := seed == SEED]
setorder(seed_tbl, arm, seed)
save_tsv(seed_tbl[, .(arm, seed, is_primary_seed, boot_B, n_modules,
                      n_ci_excludes_0, n_fdr_sig,
                      modules_ci_excludes_0, modules_fdr_sig)],
         "29_bootstrap_seed_stability.tsv")

mod_stab <- all_sw[, .(n_seeds = .N, boot_B = BOOT_B,
                       delta_logHR = round(delta[1], 3),
                       n_seeds_ci_excludes_0 = sum(ci_excludes_0),
                       frac_ci_excludes_0 = round(mean(ci_excludes_0), 3),
                       n_seeds_fdr_sig = sum(fdr_sig),
                       frac_fdr_sig = round(mean(fdr_sig), 3),
                       median_lo = round(median(lo), 3), median_hi = round(median(hi), 3),
                       min_lo = round(min(lo), 3), max_hi = round(max(hi), 3),
                       median_p_boot = signif(median(p_boot), 3),
                       min_p_boot = signif(min(p_boot), 3),
                       max_p_boot = signif(max(p_boot), 3),
                       median_fdr = signif(median(fdr), 3)),
                   by = .(module, arm)]
mod_stab <- merge(mod_stab,
                  boot_raw[, .(module, arm,
                               primary_seed_ci_excludes_0 = ci_excludes_0,
                               primary_seed_fdr_sig = fdr_sig)],
                  by = c("module", "arm"))
mod_stab[, `:=`(seed_sweep_seeds = paste(range(SWEEP_SEEDS), collapse = " to "),
                boot_fdr_family = BOOT_FDR_FAMILY)]
setorder(mod_stab, arm, -frac_ci_excludes_0, module)
save_tsv(mod_stab, "29_bootstrap_seed_stability_by_module.tsv")

sc_o <- seed_tbl[arm == "observed"]; sc_r <- seed_tbl[arm == "residualised"]
# "module 16%" for each module whose call changes across seeds, or "none" if all
# are stable. The guard stops paste0() rendering an empty input as "%".
unstable_list <- function(mod, frac) {
  i <- which(frac > 0 & frac < 1)
  if (!length(i)) return("none")
  paste(paste0(mod[i], " ", round(100 * frac[i]), "%"), collapse = "; ")
}
msg("Observed arm across ", N_SEED_SWEEP, " seeds: intervals excluding 0, median ",
    median(sc_o$n_ci_excludes_0), " of 16 (range ", min(sc_o$n_ci_excludes_0), " to ",
    max(sc_o$n_ci_excludes_0), "; primary seed ", sc_o[is_primary_seed == TRUE, n_ci_excludes_0],
    "); at BH < ", FDR_ALPHA, ", median ", median(sc_o$n_fdr_sig), " (range ",
    min(sc_o$n_fdr_sig), " to ", max(sc_o$n_fdr_sig), "; primary seed ",
    sc_o[is_primary_seed == TRUE, n_fdr_sig], ")")
msg("Residualised arm across ", N_SEED_SWEEP, " seeds: intervals excluding 0, median ",
    median(sc_r$n_ci_excludes_0), " (range ", min(sc_r$n_ci_excludes_0), " to ",
    max(sc_r$n_ci_excludes_0), "); at BH < ", FDR_ALPHA, ", median ",
    median(sc_r$n_fdr_sig), " (range ", min(sc_r$n_fdr_sig), " to ",
    max(sc_r$n_fdr_sig), ")")
print(mod_stab[arm == "observed",
               .(module, delta_logHR, frac_ci_excludes_0, frac_fdr_sig,
                 median_lo, median_hi, min_p_boot, max_p_boot)], row.names = FALSE)

tw <- dcast(arms, module ~ arm, value.var = c("HR", "lo", "hi", "p", "fdr"))
status_of <- function(f_obs, f_res)
  fifelse(f_res < FDR_ALPHA & f_obs < FDR_ALPHA, "robust",
  fifelse(f_res >= FDR_ALPHA & f_obs < FDR_ALPHA, "LOST on adjustment",
  fifelse(f_res < FDR_ALPHA & f_obs >= FDR_ALPHA, "gained on adjustment", "not significant")))
contrast <- merge(s12[, .(module, n_12 = n, events_12 = events,
                          HR_fpkm_observed = HR_unadjusted, fdr_fpkm_observed = fdr_unadjusted,
                          HR_fpkm_observed_STAR = HR_observed_STAR, fdr_fpkm_observed_STAR = fdr_observed_STAR,
                          HR_fpkm_residualised = HR_adjusted, fdr_fpkm_residualised = fdr_adjusted,
                          status_fpkm = status)],
                  tw[, .(module,
                         HR_tmm_observed = HR_tmm_observed, lo_tmm_observed, hi_tmm_observed,
                         p_tmm_observed, fdr_tmm_observed,
                         HR_tmm_observed_STAR, lo_tmm_observed_STAR, hi_tmm_observed_STAR,
                         p_tmm_observed_STAR, fdr_tmm_observed_STAR,
                         HR_tmm_residualised, lo_tmm_residualised, hi_tmm_residualised,
                         p_tmm_residualised, fdr_tmm_residualised)],
                  by = "module")
contrast <- merge(contrast, me_cor, by = "module")
contrast <- merge(contrast, boot_tbl, by = "module")
# Fraction of seeds in which each module's interval excludes zero.
stab_wide <- dcast(mod_stab, module ~ arm,
                   value.var = c("frac_ci_excludes_0", "frac_fdr_sig"))
setnames(stab_wide,
         c("frac_ci_excludes_0_observed", "frac_ci_excludes_0_residualised",
           "frac_fdr_sig_observed", "frac_fdr_sig_residualised"),
         c("seed_frac_ci_excludes_0_observed", "seed_frac_ci_excludes_0_residualised",
           "seed_frac_fdr_sig_observed", "seed_frac_fdr_sig_residualised"))
stab_wide[, `:=`(seed_sweep_n_seeds = N_SEED_SWEEP, seed_sweep_B = BOOT_B)]
contrast <- merge(contrast, stab_wide, by = "module")
contrast[, n_genes := gt[module != "grey", .N, by = module][match(contrast$module, module), N]]
contrast[, `:=`(n = sum(keep), events = sum(cl$os_event[keep]))]
contrast[, status_tmm := status_of(fdr_tmm_observed, fdr_tmm_residualised)]
contrast[, `:=`(
  sig_fpkm_observed         = fdr_fpkm_observed < FDR_ALPHA,
  sig_tmm_observed          = fdr_tmm_observed < FDR_ALPHA,
  sig_fpkm_observed_STAR    = fdr_fpkm_observed_STAR < FDR_ALPHA,
  sig_tmm_observed_STAR     = fdr_tmm_observed_STAR < FDR_ALPHA,
  sig_fpkm_residualised     = fdr_fpkm_residualised < FDR_ALPHA,
  sig_tmm_residualised      = fdr_tmm_residualised < FDR_ALPHA)]
contrast[, `:=`(
  agree_sig_observed      = sig_fpkm_observed == sig_tmm_observed,
  agree_sig_observed_STAR = sig_fpkm_observed_STAR == sig_tmm_observed_STAR,
  agree_sig_residualised  = sig_fpkm_residualised == sig_tmm_residualised,
  agree_direction_observed     = sign(log(HR_fpkm_observed)) == sign(log(HR_tmm_observed)),
  agree_direction_residualised = sign(log(HR_fpkm_residualised)) == sign(log(HR_tmm_residualised)),
  agree_status = status_fpkm == status_tmm,
  log_HR_ratio_observed     = round(log(HR_tmm_observed / HR_fpkm_observed), 3),
  log_HR_ratio_residualised = round(log(HR_tmm_residualised / HR_fpkm_residualised), 3))]
setorder(contrast, fdr_tmm_residualised)
out <- contrast[, .(biotype = "lncRNA", module, n_genes, covariate_set = "clinical", n, events,
                    HR_fpkm_observed, fdr_fpkm_observed, HR_fpkm_observed_STAR, fdr_fpkm_observed_STAR,
                    HR_fpkm_residualised, fdr_fpkm_residualised, status_fpkm,
                    HR_tmm_observed = round(HR_tmm_observed, 3), lo_tmm_observed = round(lo_tmm_observed, 3),
                    hi_tmm_observed = round(hi_tmm_observed, 3), p_tmm_observed = signif(p_tmm_observed, 3),
                    fdr_tmm_observed = signif(fdr_tmm_observed, 3),
                    HR_tmm_observed_STAR = round(HR_tmm_observed_STAR, 3),
                    lo_tmm_observed_STAR = round(lo_tmm_observed_STAR, 3),
                    hi_tmm_observed_STAR = round(hi_tmm_observed_STAR, 3),
                    p_tmm_observed_STAR = signif(p_tmm_observed_STAR, 3),
                    fdr_tmm_observed_STAR = signif(fdr_tmm_observed_STAR, 3),
                    HR_tmm_residualised = round(HR_tmm_residualised, 3),
                    lo_tmm_residualised = round(lo_tmm_residualised, 3),
                    hi_tmm_residualised = round(hi_tmm_residualised, 3),
                    p_tmm_residualised = signif(p_tmm_residualised, 3),
                    fdr_tmm_residualised = signif(fdr_tmm_residualised, 3),
                    status_tmm, r_ME_observed, r_ME_residualised,
                    agree_sig_observed, agree_sig_observed_STAR, agree_sig_residualised,
                    agree_direction_observed, agree_direction_residualised, agree_status,
                    log_HR_ratio_observed, log_HR_ratio_residualised,
                    boot_delta_logHR_observed, boot_lo_observed, boot_hi_observed,
                    boot_p_observed, boot_fdr_observed, boot_ci_excludes_0_observed,
                    boot_sig_fdr_observed, boot_n_observed,
                    boot_delta_logHR_residualised, boot_lo_residualised, boot_hi_residualised,
                    boot_p_residualised, boot_fdr_residualised, boot_ci_excludes_0_residualised,
                    boot_sig_fdr_residualised, boot_n_residualised, boot_fdr_family,
                    seed_sweep_n_seeds, seed_sweep_B,
                    seed_frac_ci_excludes_0_observed, seed_frac_fdr_sig_observed,
                    seed_frac_ci_excludes_0_residualised, seed_frac_fdr_sig_residualised)]
save_tsv(out, "29_normalisation_module_contrast.tsv")
print(out[, .(module, HR_fpkm_observed, HR_tmm_observed, HR_fpkm_residualised, HR_tmm_residualised,
              status_fpkm, status_tmm, r_ME_observed, r_ME_residualised,
              boot_delta_logHR_observed, boot_lo_observed, boot_hi_observed,
              boot_fdr_observed, seed_frac_ci_excludes_0_observed,
              seed_frac_fdr_sig_observed)], row.names = FALSE)

# One row per normalisation. Pairwise quantities go to the agreement table.
# Stage 12 reports FDRs to three significant figures and some lie near 0.05,
# so the full-precision FPKM refit is given as a separate row (see `source`).
sig_counts <- function(f_obs, f_star, f_res, stat) list(
  n_sig_observed = sum(f_obs < FDR_ALPHA), n_sig_observed_STAR = sum(f_star < FDR_ALPHA),
  n_sig_residualised = sum(f_res < FDR_ALPHA),
  n_lost = sum(stat == "LOST on adjustment"), n_gained = sum(stat == "gained on adjustment"),
  n_robust = sum(stat == "robust"))
fpkm_refit <- tw[match(contrast$module, module)]
status_fpkm_refit <- status_of(fpkm_refit$fdr_fpkm_observed, fpkm_refit$fdr_fpkm_residualised)
summary_tbl <- rbindlist(list(
  c(list(normalisation = "log2(FPKM+1)",
         source = "12_technical_adjustment_sensitivity.tsv, FDR as published (3 significant figures)",
         n_modules = nrow(contrast), n = unique(contrast$n_12), events = unique(contrast$events_12)),
    sig_counts(contrast$fdr_fpkm_observed, contrast$fdr_fpkm_observed_STAR,
               contrast$fdr_fpkm_residualised, contrast$status_fpkm)),
  c(list(normalisation = "log2(FPKM+1)",
         source = "refitted in this script, full FDR precision",
         n_modules = nrow(contrast), n = sum(keep), events = sum(cl$os_event[keep])),
    sig_counts(fpkm_refit$fdr_fpkm_observed, fpkm_refit$fdr_fpkm_observed_STAR,
               fpkm_refit$fdr_fpkm_residualised, status_fpkm_refit)),
  c(list(normalisation = "log2 CPM TMM",
         source = "this script, full FDR precision",
         n_modules = nrow(contrast), n = sum(keep), events = sum(cl$os_event[keep])),
    sig_counts(contrast$fdr_tmm_observed, contrast$fdr_tmm_observed_STAR,
               contrast$fdr_tmm_residualised, contrast$status_tmm))))
save_tsv(summary_tbl, "29_normalisation_module_summary.tsv")
print(t(summary_tbl))

# The FPKM-versus-TMM contrast, one row. Agreement of significance calls does not
# test whether effects differ. The paired bootstrap columns do.
agree_tbl <- data.table(
  contrast = "log2(FPKM+1) vs log2 CPM TMM; same 16 modules, same patients, same covariates",
  n_modules = nrow(contrast), n = sum(keep), events = sum(cl$os_event[keep]),
  n_modules_same_call_observed      = sum(contrast$agree_sig_observed),
  n_modules_same_call_observed_STAR = sum(contrast$agree_sig_observed_STAR),
  n_modules_same_call_residualised  = sum(contrast$agree_sig_residualised),
  n_modules_same_status             = sum(contrast$agree_status),
  n_modules_same_direction_observed     = sum(contrast$agree_direction_observed),
  n_modules_same_direction_residualised = sum(contrast$agree_direction_residualised),
  median_r_ME_observed     = round(median(contrast$r_ME_observed), 3),
  min_r_ME_observed        = round(min(contrast$r_ME_observed), 3),
  median_r_ME_residualised = round(median(contrast$r_ME_residualised), 3),
  min_r_ME_residualised    = round(min(contrast$r_ME_residualised), 3),
  # Overlap of the per-arm 95% CIs, both refitted here.
  n_modules_HR_ci_overlap_observed =
    sum(pmax(fpkm_refit$lo_fpkm_observed, contrast$lo_tmm_observed) <=
        pmin(fpkm_refit$hi_fpkm_observed, contrast$hi_tmm_observed)),
  n_modules_tmm_point_inside_fpkm_ci_observed =
    sum(fpkm_refit$lo_fpkm_observed <= contrast$HR_tmm_observed &
        contrast$HR_tmm_observed <= fpkm_refit$hi_fpkm_observed),
  boot_B = BOOT_B, boot_n_cox_fits = .nfit, boot_n_convergence_warnings = .nwarn,
  boot_primary_seed = SEED, boot_fdr_family = BOOT_FDR_FAMILY,
  median_boot_delta_logHR_observed = round(median(contrast$boot_delta_logHR_observed), 3),
  # BH-adjusted counts, with the unadjusted interval counts beside them.
  n_boot_fdr_sig_observed          = sum(contrast$boot_sig_fdr_observed),
  n_boot_fdr_sig_negative_observed =
    sum(contrast$boot_sig_fdr_observed & contrast$boot_delta_logHR_observed < 0),
  n_boot_ci_excludes_0_observed    = sum(contrast$boot_ci_excludes_0_observed),
  n_boot_delta_negative_ci_excludes_0_observed =
    sum(contrast$boot_ci_excludes_0_observed & contrast$boot_delta_logHR_observed < 0),
  median_boot_delta_logHR_residualised = round(median(contrast$boot_delta_logHR_residualised), 3),
  n_boot_fdr_sig_residualised          = sum(contrast$boot_sig_fdr_residualised),
  n_boot_fdr_sig_negative_residualised =
    sum(contrast$boot_sig_fdr_residualised & contrast$boot_delta_logHR_residualised < 0),
  n_boot_ci_excludes_0_residualised    = sum(contrast$boot_ci_excludes_0_residualised),
  n_boot_delta_negative_ci_excludes_0_residualised =
    sum(contrast$boot_ci_excludes_0_residualised & contrast$boot_delta_logHR_residualised < 0),
  # Distribution of those counts over the seed sweep.
  seed_sweep_n_seeds = N_SEED_SWEEP, seed_sweep_B = BOOT_B,
  seed_sweep_seeds   = paste(range(SWEEP_SEEDS), collapse = " to "),
  seed_count_ci_excludes_0_observed_at_primary_seed = sc_o[is_primary_seed == TRUE, n_ci_excludes_0],
  seed_count_ci_excludes_0_observed_median = as.numeric(median(sc_o$n_ci_excludes_0)),
  seed_count_ci_excludes_0_observed_min    = min(sc_o$n_ci_excludes_0),
  seed_count_ci_excludes_0_observed_max    = max(sc_o$n_ci_excludes_0),
  seed_count_ci_excludes_0_observed_p10    = unname(quantile(sc_o$n_ci_excludes_0, 0.10)),
  seed_count_ci_excludes_0_observed_p90    = unname(quantile(sc_o$n_ci_excludes_0, 0.90)),
  seed_count_fdr_sig_observed_at_primary_seed = sc_o[is_primary_seed == TRUE, n_fdr_sig],
  seed_count_fdr_sig_observed_median = as.numeric(median(sc_o$n_fdr_sig)),
  seed_count_fdr_sig_observed_min    = min(sc_o$n_fdr_sig),
  seed_count_fdr_sig_observed_max    = max(sc_o$n_fdr_sig),
  seed_count_fdr_sig_observed_p10    = unname(quantile(sc_o$n_fdr_sig, 0.10)),
  seed_count_fdr_sig_observed_p90    = unname(quantile(sc_o$n_fdr_sig, 0.90)),
  seed_count_ci_excludes_0_residualised_at_primary_seed = sc_r[is_primary_seed == TRUE, n_ci_excludes_0],
  seed_count_ci_excludes_0_residualised_median = as.numeric(median(sc_r$n_ci_excludes_0)),
  seed_count_ci_excludes_0_residualised_min    = min(sc_r$n_ci_excludes_0),
  seed_count_ci_excludes_0_residualised_max    = max(sc_r$n_ci_excludes_0),
  seed_count_fdr_sig_residualised_at_primary_seed = sc_r[is_primary_seed == TRUE, n_fdr_sig],
  seed_count_fdr_sig_residualised_median = as.numeric(median(sc_r$n_fdr_sig)),
  seed_count_fdr_sig_residualised_min    = min(sc_r$n_fdr_sig),
  seed_count_fdr_sig_residualised_max    = max(sc_r$n_fdr_sig),
  # Modules whose call is the same in every seed, and those whose call changes.
  seed_n_modules_ci_excludes_0_always_observed =
    mod_stab[arm == "observed", sum(frac_ci_excludes_0 == 1)],
  seed_n_modules_ci_excludes_0_never_observed =
    mod_stab[arm == "observed", sum(frac_ci_excludes_0 == 0)],
  seed_n_modules_ci_excludes_0_unstable_observed =
    mod_stab[arm == "observed", sum(frac_ci_excludes_0 > 0 & frac_ci_excludes_0 < 1)],
  seed_modules_ci_excludes_0_unstable_observed =
    mod_stab[arm == "observed", unstable_list(module, frac_ci_excludes_0)],
  seed_n_modules_fdr_sig_always_observed =
    mod_stab[arm == "observed", sum(frac_fdr_sig == 1)],
  seed_n_modules_fdr_sig_never_observed =
    mod_stab[arm == "observed", sum(frac_fdr_sig == 0)],
  seed_n_modules_fdr_sig_unstable_observed =
    mod_stab[arm == "observed", sum(frac_fdr_sig > 0 & frac_fdr_sig < 1)],
  seed_modules_fdr_sig_unstable_observed =
    mod_stab[arm == "observed", unstable_list(module, frac_fdr_sig)],
  seed_n_modules_ci_excludes_0_unstable_residualised =
    mod_stab[arm == "residualised", sum(frac_ci_excludes_0 > 0 & frac_ci_excludes_0 < 1)],
  seed_modules_ci_excludes_0_unstable_residualised =
    mod_stab[arm == "residualised", unstable_list(module, frac_ci_excludes_0)],
  seed_n_modules_fdr_sig_unstable_residualised =
    mod_stab[arm == "residualised", sum(frac_fdr_sig > 0 & frac_fdr_sig < 1)],
  seed_modules_fdr_sig_unstable_residualised =
    mod_stab[arm == "residualised", unstable_list(module, frac_fdr_sig)])
save_tsv(agree_tbl, "29_normalisation_agreement.tsv")
print(t(agree_tbl))

msg("Elapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
write_session_info("29_normalisation_check")
banner("29 | done")
