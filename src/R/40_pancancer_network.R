# 40_pancancer_network.R: does the quality axis reorganise the lncRNA network in every cancer type? (rule C)
# For each eligible tumour project (>= PAN_MIN_TUMOURS retained tumours,
# TCGA-LAML excluded), builds the lncRNA network on observed log2(FPKM + 1)
# and on expression residualised on the three STAR metrics, with the discovery
# lncRNA WGCNA parameters and sample-outlier rule (standardised connectivity
# < -2.5). Parameters are not tuned per project. Reports module counts,
# largest-module and grey fractions, the correlation of the largest observed
# eigengene with the non-feature fraction, and the adjusted Rand index between
# partitions, then applies rule C (00_config.R). Run after 38.
# Outputs: results/40_pancancer_network.tsv, 40_pancancer_network_prevalence.tsv,
#   cache/40_pan_net_<project>.rds (both partitions, so the stage resumes).

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
source(file.path(R_DIR, "pan_helpers.R"))
suppressPackageStartupMessages({ library(data.table); library(WGCNA); library(matrixStats) })
enableWGCNAThreads(N_THREADS)
banner("40 | Pan-cancer lncRNA network reorganisation (rule C)")
verify_decision_lock()

# PAN_TEST_ALL=1 builds networks for ineligible projects too, writing *_partial
# tables, so the mechanics can be tested before the eligible set is complete.
TEST <- nzchar(Sys.getenv("PAN_TEST_ALL")); SUFFIX <- if (TEST) "_partial" else ""

ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  s_ij <- sum(choose(tab, 2)); s_i <- sum(choose(rowSums(tab), 2)); s_j <- sum(choose(colSums(tab), 2))
  e <- s_i * s_j / choose(n, 2)
  (s_ij - e) / ((s_i + s_j) / 2 - e)
}
build <- function(E) {
  set.seed(SEED)
  bw <- blockwiseModules(E, power = LNC_POWER, networkType = NETWORK_TYPE, TOMType = TOM_TYPE,
                         minModuleSize = MIN_MODULE_SIZE, mergeCutHeight = LNC_MERGE,
                         deepSplit = LNC_DEEPSPLIT, minKMEtoStay = LNC_MINKME,
                         numericLabels = FALSE, pamRespectsDendro = FALSE,
                         maxBlockSize = MAX_BLOCK_SIZE, saveTOMs = FALSE, verbose = 0)
  cols <- setNames(bw$colors, colnames(E))
  list(colors = cols, MEs = moduleEigengenes(E, colors = cols)$eigengenes)
}

rows <- list()
for (proj in pan_projects()) {
  obj <- readRDS(file.path(CACHE_DIR, paste0("pan_", proj, ".rds")))
  n_ret <- obj$samples[keep == TRUE & group == "tumour", .N]
  eligible <- n_ret >= PAN_MIN_TUMOURS && !proj %in% PAN_EXCLUDE_INFERENCE
  if (!eligible && !TEST) { rm(obj); next }
  f <- file.path(CACHE_DIR, paste0("40_pan_net_", proj, ".rds"))
  if (file.exists(f) && !TEST) { nt <- readRDS(f); msg(proj, ": networks from cache") } else {
    t0 <- Sys.time()
    G <- pan_group(obj, "tumour"); E <- G$lnc; s <- G$samples
    gsg <- goodSamplesGenes(E, verbose = 0)
    E <- E[gsg$goodSamples, gsg$goodGenes]; s <- s[gsg$goodSamples]
    A <- adjacency(t(E), type = "distance"); k <- colSums(A) - 1; rm(A)
    outl <- which(scale(k)[, 1] < -2.5)
    if (length(outl)) { E <- E[-outl, , drop = FALSE]; s <- s[-outl] }
    cov <- tech_covariates(s)
    obs <- build(E); adj <- build(remove_technical(E, cov))
    nt <- list(project = proj, samples = s$file_id, noFeature = s$pct_noFeature,
               n_outliers = length(outl), obs = obs, adj = adj, eligible = eligible)
    if (!TEST) saveRDS(nt, f)
    msg(sprintf("%s: %d samples x %d lncRNAs, both networks in %.1f min", proj, nrow(E), ncol(E),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  summ <- function(cols) {
    tb <- table(cols[cols != "grey"])
    list(n_modules = length(tb), largest = names(tb)[which.max(tb)],
         largest_frac = if (length(tb)) max(tb) / length(cols) else 0,
         grey_frac = mean(cols == "grey"))
  }
  so <- summ(nt$obs$colors); sa <- summ(nt$adj$colors)
  me <- nt$obs$MEs[[paste0("ME", so$largest)]]
  r_me <- if (is.null(me)) NA_real_ else as.numeric(cor(as.numeric(me), nt$noFeature))
  both <- nt$obs$colors != "grey" & nt$adj$colors != "grey"
  rows[[length(rows) + 1]] <- data.table(
    project = proj, eligible = nt$eligible, n_samples = length(nt$samples),
    n_outliers_removed = nt$n_outliers, n_lncRNA = length(nt$obs$colors),
    n_modules_observed = so$n_modules, n_modules_adjusted = sa$n_modules,
    largest_frac_observed = as.numeric(so$largest_frac), largest_frac_adjusted = as.numeric(sa$largest_frac),
    grey_frac_observed = so$grey_frac, grey_frac_adjusted = sa$grey_frac,
    r_largest_observed_ME_noFeature = r_me,
    ari_all_genes = as.numeric(ari(nt$obs$colors, nt$adj$colors)),
    ari_assigned_both = if (sum(both) > 1) as.numeric(ari(nt$obs$colors[both], nt$adj$colors[both])) else NA_real_)
  rm(obj, nt); gc(verbose = FALSE)
}
nw <- rbindlist(rows)
stopifnot(all(c("r_largest_observed_ME_noFeature", "ari_all_genes") %in% names(nw)))
nw[, rule_C_pass := largest_frac_observed >= PAN_NET_LARGEST_MIN &
                    abs(r_largest_observed_ME_noFeature) >= PAN_NET_EIGEN_R_MIN &
                    ari_all_genes < PAN_NET_ARI_MAX]
setorder(nw, -largest_frac_observed)
save_tsv(nw, paste0("40_pancancer_network", SUFFIX, ".tsv"))
el <- nw[eligible == TRUE]
prev <- data.table(
  quantity = c("eligible tumour projects with networks", "rule C met",
               "largest observed module >= 40% of genes", "largest observed eigengene |r| >= 0.7",
               "ARI < 0.2", "more modules after adjustment than before"),
  n = c(nrow(el), sum(el$rule_C_pass), sum(el$largest_frac_observed >= PAN_NET_LARGEST_MIN),
        sum(abs(el$r_largest_observed_ME_noFeature) >= PAN_NET_EIGEN_R_MIN, na.rm = TRUE),
        sum(el$ari_all_genes < PAN_NET_ARI_MAX), sum(el$n_modules_adjusted > el$n_modules_observed)),
  of = c(NA, rep(nrow(el), 5)))
save_tsv(prev, paste0("40_pancancer_network_prevalence", SUFFIX, ".tsv")); print(prev)
write_session_info("40_pancancer_network")
msg("40 done")
