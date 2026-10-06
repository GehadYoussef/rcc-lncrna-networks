# 42_singlecell_truth_correction.R: tests whether library-quality correction removes artefact or biology
# The malignant-vs-normal-epithelium specificity of an lncRNA in single cells does
# not depend on bulk library quality, so it serves as ground truth. If the
# correction removes artefact and keeps biology, expect:
#   T1  bulk tumour-minus-normal differences agree with single-cell fold changes as well or better
#   T2  a score of single-cell malignant-specific lncRNAs loses its correlation with
#       the non-feature fraction but keeps it with a malignant-marker score
#   T3  candidates set aside in bulk as technically sensitive become single-cell concordant
#   T4  abundance-matched non-specific lncRNA sets show no marker coherence.
# Inputs: data/derived/cache/pan_TCGA-KIRC.rds (pan-cancer build) and results/singlecell.
# Outputs: 42_sc_truth_*.tsv tables and the single-cell correction supplementary figure.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(matrixStats); library(ggplot2); library(patchwork) })
banner("42 | Single-cell ground truth for the library-quality correction")
set.seed(SEED)
B_BOOT <- 2000; N_NULL <- 1000
MARKERS <- c("CA9", "NDUFA4L2", "ANGPTL4", "VEGFA", "EGLN3")

SC_DIR <- Filter(dir.exists, c(file.path(PROJECT_ROOT, "submission", "results", "singlecell"),
                               file.path(PROJECT_ROOT, "results", "singlecell")))[1]
if (is.na(SC_DIR)) stop("single-cell results directory not found")

# ---- 1. bulk: KIRC tumours and matched normals ------------------------------
obj <- readRDS(file.path(CACHE_DIR, "pan_TCGA-KIRC.rds"))
s <- obj$samples[keep == TRUE]
sT <- s[group == "tumour"]; sN <- s[group == "normal"]
lnc_ids <- obj$ann[gene_type == "lncRNA", gene_id]
fT <- obj$fpkm[lnc_ids, sT$file_id]; fN <- obj$fpkm[lnc_ids, sN$file_id]
expressed <- rowMeans(fT >= LNC_MIN_FPKM) >= LNC_MIN_FRAC | rowMeans(fN >= LNC_MIN_FPKM) >= LNC_MIN_FRAC
XT <- t(log2(fT[expressed, ] + 1)); XN <- t(log2(fN[expressed, ] + 1))
cT <- tech_covariates(sT); cN <- tech_covariates(sN)
msg(nrow(XT), " tumours, ", nrow(XN), " normals, ", ncol(XT), " expressed lncRNAs")

fitT <- fit_technical(XT, cT)
corr <- list(
  observed = list(T = XT, N = XN),
  tumour_fitted = list(T = apply_technical(XT, fitT, cT), N = apply_technical(XN, fitT, cN)))
# Main arm: per-gene OLS on the standardised STAR metrics, fitted in tumours and
# applied unchanged to normals. Joint arm: expr ~ tumour + metrics on both groups,
# removing only the metric terms so the tissue difference is not absorbed.
Xall <- rbind(XT, XN); call <- rbind(cT, cN)
mu <- colMeans(call); sdv <- apply(call, 2, sd)
Z <- scale(call, center = mu, scale = sdv)
D <- cbind(1, tumour = c(rep(1, nrow(XT)), rep(0, nrow(XN))), Z)
Bj <- qr.solve(crossprod(D), crossprod(D, Xall))
Xj <- Xall - Z %*% Bj[colnames(Z), , drop = FALSE]
corr$joint_fitted <- list(T = Xj[seq_len(nrow(XT)), ], N = Xj[nrow(XT) + seq_len(nrow(XN)), ])

pairs <- merge(sT[, .(patient, fT = file_id, plateT = plate)], sN[, .(patient, fN = file_id, plateN = plate)],
               by = "patient")
msg(nrow(pairs), " tumour-normal pairs, ", pairs[plateT == plateN, .N], " on a shared plate")
paired_diff <- function(arm, sel = rep(TRUE, nrow(pairs))) {
  p <- pairs[sel]
  colMeans(corr[[arm]]$T[p$fT, , drop = FALSE] - corr[[arm]]$N[p$fN, , drop = FALSE])
}

# ---- 2. single-cell truth ---------------------------------------------------
# The single-cell candidates were locked before any bulk step. To avoid
# circularity the bulk-retained set is never used as a selector.
meta <- fread(file.path(SC_DIR, "08_meta_lncRNA.tsv.gz"))
scN <- meta[contrast == "malignant_vs_normal_epithelial", .(gk = gene_key, sc_lfc = pooled_log2FC,
                                                           sc_fdr = pooled_FDR, sc_k = k)]
scT <- meta[contrast == "malignant_vs_pooled_tme", .(gk = gene_key, tme_lfc = pooled_log2FC, tme_fdr = pooled_FDR)]
lock <- fread(file.path(SC_DIR, "locked_lncRNA_candidates.tsv"))[, gk := sub("[.].*$", "", gene_key)]
dec  <- fread(file.path(SC_DIR, "10_bulk_candidate_decisions.tsv"))[, gk := sub("[.].*$", "", gene_key)]
key_of <- function(x) sub("[.].*$", "", x)
gene_key_bulk <- key_of(colnames(XT))

# ---- 3. T1: concordance of bulk tumour-normal with single-cell fold change ---
gl <- data.table(gk = gene_key_bulk, gene_id = colnames(XT))
for (arm in names(corr)) {
  gl[, (paste0("d_", arm)) := paired_diff(arm)]
  gl[, (paste0("dshared_", arm)) := paired_diff(arm, pairs$plateT == pairs$plateN)]
}
gl[, rho_metric_obs := as.numeric(cor(XT, cT[, "pct_noFeature"], method = "spearman"))]
gl[, rho_metric_corr := as.numeric(cor(corr$tumour_fitted$T, cT[, "pct_noFeature"], method = "spearman"))]
gl <- merge(gl, scN, by = "gk")
gl <- merge(gl, obj$ann[, .(gene_id, gene_name)], by = "gene_id")
save_tsv(gl, "42_sc_truth_gene_level.tsv")

conc <- function(d, sc, sig) c(rho = cor(d, sc, method = "spearman"),
                                agree = mean(sign(d[sig]) == sign(sc[sig])))
t1 <- rbindlist(lapply(c("all pairs", "shared-plate pairs"), function(pset) {
  pre <- if (pset == "all pairs") "d_" else "dshared_"
  rbindlist(lapply(c("tumour_fitted", "joint_fitted"), function(arm) {
    sig <- gl$sc_fdr < 0.05
    o <- conc(gl[[paste0(pre, "observed")]], gl$sc_lfc, sig)
    a <- conc(gl[[paste0(pre, arm)]], gl$sc_lfc, sig)
    bt <- replicate(B_BOOT, {
      i <- sample.int(nrow(gl), replace = TRUE); sg <- sig[i]
      conc(gl[[paste0(pre, arm)]][i], gl$sc_lfc[i], sg) - conc(gl[[paste0(pre, "observed")]][i], gl$sc_lfc[i], sg)
    })
    data.table(pairs = pset, correction = arm, n_genes = nrow(gl), n_sc_significant = sum(sig),
               rho_observed = o[["rho"]], rho_corrected = a[["rho"]], delta_rho = a[["rho"]] - o[["rho"]],
               delta_rho_lo = quantile(bt["rho", ], 0.025), delta_rho_hi = quantile(bt["rho", ], 0.975),
               agree_observed = o[["agree"]], agree_corrected = a[["agree"]],
               delta_agree_lo = quantile(bt["agree", ], 0.025, na.rm = TRUE),
               delta_agree_hi = quantile(bt["agree", ], 0.975, na.rm = TRUE))
  }))
}))
save_tsv(t1, "42_sc_truth_concordance.tsv"); print(t1)

# ---- 3b. T1u: the unpaired contrast -------------------------------------------
# Paired tumour and normal libraries differ little in the covariates, so the
# correction barely moves the paired difference. The unpaired contrast carries
# the between-library quality gap, so it is the contrast the correction can change.
cov_gap <- colMeans(scale(rbind(cT, cN))[seq_len(nrow(cT)), ]) - colMeans(scale(rbind(cT, cN))[nrow(cT) + seq_len(nrow(cN)), ])
msg("unpaired tumour-minus-normal difference in the covariates (pooled SD units): ",
    paste(names(cov_gap), round(cov_gap, 2), collapse = ", "))
pair_gap <- colMeans(scale(rbind(cT, cN))[match(pairs$fT, sT$file_id), ] - scale(rbind(cT, cN))[nrow(cT) + match(pairs$fN, sN$file_id), ])
unp <- function(arm) colMeans(corr[[arm]]$T) - colMeans(corr[[arm]]$N)
gu <- data.table(gk = gene_key_bulk)
for (arm in names(corr)) gu[, (paste0("u_", arm)) := unp(arm)]
gu <- merge(gu, scN, by = "gk")
sigu <- gu$sc_fdr < 0.05
t1u <- rbindlist(lapply(c("tumour_fitted", "joint_fitted"), function(arm) {
  o <- conc(gu$u_observed, gu$sc_lfc, sigu); a <- conc(gu[[paste0("u_", arm)]], gu$sc_lfc, sigu)
  bt <- replicate(B_BOOT, { i <- sample.int(nrow(gu), replace = TRUE); sg <- sigu[i]
    conc(gu[[paste0("u_", arm)]][i], gu$sc_lfc[i], sg) - conc(gu$u_observed[i], gu$sc_lfc[i], sg) })
  data.table(pairs = "unpaired, all tumours vs all normals", correction = arm,
             n_genes = nrow(gu), n_sc_significant = sum(sigu),
             rho_observed = o[["rho"]], rho_corrected = a[["rho"]], delta_rho = a[["rho"]] - o[["rho"]],
             delta_rho_lo = quantile(bt["rho", ], 0.025), delta_rho_hi = quantile(bt["rho", ], 0.975),
             agree_observed = o[["agree"]], agree_corrected = a[["agree"]],
             delta_agree_lo = quantile(bt["agree", ], 0.025, na.rm = TRUE),
             delta_agree_hi = quantile(bt["agree", ], 0.975, na.rm = TRUE))
}))
gaps <- data.table(covariate = names(cov_gap), unpaired_gap_sd = cov_gap, paired_mean_gap_sd = pair_gap[names(cov_gap)])
save_tsv(gaps, "42_sc_truth_covariate_gaps.tsv"); print(gaps)
t1 <- rbind(t1, t1u)
save_tsv(t1, "42_sc_truth_concordance.tsv"); print(t1u)

# ---- 4. T2 and T4: score coherence in tumours --------------------------------
pc_ids <- obj$ann[gene_type == "protein_coding" & gene_name %in% MARKERS, gene_id]
mk <- scale(t(log2(obj$fpkm[pc_ids, sT$file_id] + 1)))
marker_score <- rowMeans(mk)
nf <- cT[, "pct_noFeature"]
marker_vs_metric <- cor(marker_score, nf, method = "spearman")
msg("malignant-marker score vs non-feature fraction: Spearman ", round(marker_vs_metric, 3))
score_of <- function(X, ids) rowMeans(scale(X[, ids, drop = FALSE]))
sets <- list(
  `single-cell tier 1 (tumour-specific)` = lock[tier == "tier1_tumour_specific", gk],
  `single-cell tier 2 (malignant-compartment)` = lock[tier != "tier1_tumour_specific", gk],
  `single-cell malignant-up vs normal epithelium (meta FDR < 0.05)` = scN[sc_fdr < 0.05 & sc_lfc > 0, gk])
coh <- rbindlist(lapply(names(sets), function(nm) {
  ids <- colnames(XT)[gene_key_bulk %in% sets[[nm]]]
  if (length(ids) < 3) return(NULL)
  rbindlist(lapply(names(corr), function(arm) {
    sc_ <- score_of(corr[[arm]]$T, ids)
    bt <- replicate(B_BOOT, { i <- sample.int(length(sc_), replace = TRUE)
      c(cor(sc_[i], marker_score[i], method = "spearman"), cor(sc_[i], nf[i], method = "spearman")) })
    data.table(set = nm, n_genes = length(ids), correction = arm,
               rho_marker = cor(sc_, marker_score, method = "spearman"),
               rho_marker_lo = quantile(bt[1, ], 0.025), rho_marker_hi = quantile(bt[1, ], 0.975),
               rho_metric = cor(sc_, nf, method = "spearman"),
               rho_metric_lo = quantile(bt[2, ], 0.025), rho_metric_hi = quantile(bt[2, ], 0.975))
  }))
}))
# T4: abundance-matched sets of single-cell non-specific lncRNAs
nonspec <- scN[abs(sc_lfc) < 0.25 & sc_fdr > 0.5, gk]
pool_ids <- colnames(XT)[gene_key_bulk %in% nonspec]
tier1_ids <- colnames(XT)[gene_key_bulk %in% sets[[1]]]
abund <- colMeans(XT); bins <- cut(abund, quantile(abund, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
names(bins) <- colnames(XT)
need <- table(bins[tier1_ids])
null <- rbindlist(lapply(seq_len(N_NULL), function(k) {
  ids <- unlist(lapply(names(need), function(b) { cand <- pool_ids[bins[pool_ids] == as.integer(b)]
    if (length(cand)) sample(cand, min(length(cand), need[[b]])) }))
  data.table(draw = k, arm = names(corr),
             rho_marker = vapply(names(corr), function(a) cor(score_of(corr[[a]]$T, ids), marker_score, method = "spearman"), 1),
             rho_metric = vapply(names(corr), function(a) cor(score_of(corr[[a]]$T, ids), nf, method = "spearman"), 1))
}))
nsum <- null[, .(set = "abundance-matched single-cell non-specific (null, 1000 draws)",
                 n_genes = as.integer(sum(need)), rho_marker = median(rho_marker),
                 rho_marker_lo = quantile(rho_marker, 0.025), rho_marker_hi = quantile(rho_marker, 0.975),
                 rho_metric = median(rho_metric), rho_metric_lo = quantile(rho_metric, 0.025),
                 rho_metric_hi = quantile(rho_metric, 0.975)), by = .(correction = arm)]
coh <- rbind(coh, nsum, fill = TRUE)
coh[, marker_score_rho_metric := marker_vs_metric]
save_tsv(coh, "42_sc_truth_score_coherence.tsv"); print(coh)

# ---- 5. T3: the candidates set aside in bulk as technically sensitive -------
tech <- dec[decision_reason == "technically_sensitive", .(gk, symbol, tier, technical_r2_TCGA_KIRC)]
t3 <- merge(tech, gl[, .(gk, sc_lfc, sc_fdr, d_observed, d_tumour_fitted, d_joint_fitted)],
            by = "gk", all.x = TRUE)
t3 <- merge(t3, scT, by = "gk", all.x = TRUE)
mcor <- function(X, ids) vapply(ids, function(g) if (!is.na(g) && g %in% colnames(X)) cor(X[, g], marker_score, method = "spearman") else NA_real_, 1)
id_of <- setNames(colnames(XT), gene_key_bulk)
# metric correlation for every set-aside gene expressed in bulk, with or without
# a single-cell estimate
qcor <- function(X, ids) vapply(ids, function(g) if (!is.na(g) && g %in% colnames(X)) cor(X[, g], nf, method = "spearman") else NA_real_, 1)
t3[, in_bulk := gk %in% names(id_of)]
t3[, rho_metric_obs := qcor(corr$observed$T, id_of[gk])]
t3[, rho_metric_corr := qcor(corr$tumour_fitted$T, id_of[gk])]
t3[, rho_marker_obs := mcor(corr$observed$T, id_of[gk])]
t3[, rho_marker_corr := mcor(corr$tumour_fitted$T, id_of[gk])]
t3[, direction_agrees_obs := sign(d_observed) == sign(sc_lfc)]
t3[, direction_agrees_corr := sign(d_tumour_fitted) == sign(sc_lfc)]
save_tsv(t3[order(-technical_r2_TCGA_KIRC)], "42_sc_truth_technical_candidates.tsv"); print(t3)

# ---- 6. figure ---------------------------------------------------------------
# Each panel carries its letter and a short title; the statistics are in the
# Supplementary Fig. 6 legend.
BLUE <- "#0072B2"; VERM <- "#D55E00"
th <- theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 8.5),
                                      plot.title.position = "plot")
pa <- ggplot(gl, aes(sc_lfc, d_tumour_fitted)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(aes(colour = sc_fdr < 0.05), size = 0.8, alpha = 0.7) +
  scale_colour_manual(values = c(`TRUE` = VERM, `FALSE` = "grey60"), labels = c(`TRUE` = "single-cell FDR < 0.05", `FALSE` = "not significant"), name = NULL) +
  labs(title = "Tumour-normal agreement",
       x = "Single-cell log2 fold change, malignant vs normal epithelium", y = "Bulk paired tumour - normal, corrected (log2)") +
  th + theme(legend.position = "bottom")
cd <- coh[correction != "joint_fitted"]
cd[, set := factor(sub(" [(].*", "", set), levels = unique(sub(" [(].*", "", set)))]
cl <- melt(cd, id.vars = c("set", "correction"), measure.vars = list(c("rho_marker", "rho_metric"), c("rho_marker_lo", "rho_metric_lo"), c("rho_marker_hi", "rho_metric_hi")),
           variable.name = "with", value.name = c("rho", "lo", "hi"))
cl[, with := factor(ifelse(with == 1, "malignant-marker score", "non-feature fraction"))]
cl[, correction := factor(ifelse(correction == "observed", "observed", "corrected"), levels = c("observed", "corrected"))]
pb <- ggplot(cl, aes(rho, set, colour = correction)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 1.6, position = position_dodge(0.5)) + facet_wrap(~ with) +
  scale_colour_manual(values = c(observed = VERM, corrected = BLUE), name = NULL) +
  labs(title = "Score coherence",
       x = "Spearman rho", y = NULL) + th + theme(legend.position = "bottom")
t3p <- t3[in_bulk == TRUE]
t3p[, symbol := factor(symbol, levels = symbol[order(rho_metric_obs)])]
pc_ <- ggplot(t3p) +
  geom_segment(aes(x = rho_metric_obs, xend = rho_metric_corr, y = symbol, yend = symbol), colour = "grey70") +
  geom_point(aes(rho_metric_obs, symbol), colour = VERM, size = 1.6) +
  geom_point(aes(rho_metric_corr, symbol), colour = BLUE, size = 1.6) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
  labs(title = "Set-aside candidates",
       x = "Spearman rho with the non-feature fraction", y = NULL) + th
# Bold lower-case tags, as in the other supplementary figures.
save_fig((pa | pb) / (pc_ | plot_spacer()) +
           plot_annotation(tag_levels = "a") &
           theme(plot.tag = element_text(face = "bold", size = 9.5)),
         "SupplementaryFigureS6_sc_truth_correction", 7.2, 7.6)

write_session_info("42_singlecell_truth_correction")
msg("42 done")
