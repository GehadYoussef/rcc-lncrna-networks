# 22_locked_model_export.R: export the locked prediction model so it can be applied to new patients.
#
# Covers both comparators (clinical, and augmented with ESTIMATE and STAR metrics).
# Inputs: locked_model.rds (09), networks.rds (gene names only),
# estimate_scores.rds and star_qc.rds (07). Nothing is refitted.
# Outputs: model coefficients, standardisation constants, discovery baseline
# survival at 1, 3 and 5 years, gene-level module loadings with their centring
# and scaling constants, and score-level constants (22_*.tsv in results).

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(survival) })
banner("22 | Locked model tables and module loadings")

L <- readRDS(file.path(CACHE_DIR, "locked_model.rds"))
stopifnot(identical(L$version, "v9"))
yt  <- L$yt
tcl <- as.data.table(L$tcl)
stopifnot(identical(tcl$sample_barcode, L$tsamp))
msg("Discovery locked set: n = ", length(L$tsamp), ", events = ", sum(yt[, 2]))

lp_of <- function(X, b) as.numeric(X[, names(b), drop = FALSE] %*% b)
COMPARATORS <- c("clinical", "augmented")

cat("\n############ D. THE LOCKED PREDICTION MODEL ############\n")

# ---- coefficients -----------------------------------------------------------
coefs <- rbindlist(lapply(COMPARATORS, function(cmp) {
  M <- L$models[[cmp]]
  rbindlist(list(
    data.table(comparator = cmp, model = "comparator only", term = names(M$b_clin),
               coefficient = round(M$b_clin, 4), HR = round(exp(M$b_clin), 3),
               penalised = FALSE, lambda_min = NA_real_),
    data.table(comparator = cmp, model = "combined", term = names(M$b_full),
               coefficient = round(M$b_full, 4), HR = round(exp(M$b_full), 3),
               penalised = names(M$b_full) %in% colnames(L$Et),
               lambda_min = signif(M$lambda, 4))))
}))
print(coefs, row.names = FALSE)
save_tsv(coefs, "22_locked_model_coefficients.tsv")

# ---- standardisation constants ----------------------------------------------
# Age uses the discovery mean and SD. Augmented terms are standardised within
# the locked discovery set from the same sources as stage 09. Every constant is
# checked against the locked design matrix.
est_t <- as.data.table(readRDS(file.path(CACHE_DIR, "estimate_scores.rds")))[
           match(L$tsamp, sample_barcode)]
qc_t  <- as.data.table(readRDS(file.path(CACHE_DIR, "star_qc.rds")))[
           match(L$tsamp, sample_barcode)]
stopifnot(!anyNA(est_t$StromalScore), !anyNA(qc_t$pct_noFeature))
raw <- list(age = tcl$age, stromal = est_t$StromalScore, immune = est_t$ImmuneScore,
            noFeature = qc_t$pct_noFeature, multimap = qc_t$pct_multimapping,
            libsize_z = log10(qc_t$assigned_reads))
src <- c(age = "clinical", stromal = "ESTIMATE", immune = "ESTIMATE",
         noFeature = "STAR summary (07)", multimap = "STAR summary (07)",
         libsize_z = "STAR summary (07)")
consts <- rbindlist(lapply(COMPARATORS, function(cmp) {
  Xt <- L$X[[cmp]]$Xt
  vars <- if (cmp == "clinical") "age" else names(raw)
  rbindlist(lapply(vars, function(v) {
    x <- raw[[v]]
    if (v == "age") { m <- L$age_const[["mean"]]; s <- L$age_const[["sd"]] }
    else            { m <- mean(x);               s <- sd(x) }
    data.table(comparator = cmp, covariate = v, source = src[[v]],
               mean = round(m, 4), sd = round(s, 4),
               max_abs_dev = signif(max(abs((x - m) / s - Xt[, v])), 3))
  }))
}))
print(consts, row.names = FALSE)
stopifnot(all(consts$max_abs_dev < 1e-8))
msg("Standardisation constants reproduce the locked design matrices exactly.")
save_tsv(consts[, .(comparator, covariate, source, mean, sd)],
         "22_locked_model_standardisation.tsv")

# ---- discovery baseline survival --------------------------------------------
# Breslow baseline with the locked predictor as an offset, so 1 - risk at
# lp = 0 is the baseline survival of the locked model. Survival at the mean
# discovery predictor is given for reference.
base <- rbindlist(lapply(COMPARATORS, function(cmp) {
  M  <- L$models[[cmp]]
  Xt <- L$X[[cmp]]$Xt
  lps <- list(`comparator only` = lp_of(Xt, M$b_clin),
              combined          = lp_of(cbind(Xt, L$Et), M$b_full))
  rbindlist(lapply(names(lps), function(nm) {
    lp <- lps[[nm]]; rf <- baseline_risk_fun(yt, lp)
    rbindlist(lapply(EVAL_TIMES_YRS, function(yr) {
      t <- yr * 365.25
      data.table(comparator = cmp, model = nm, years = yr,
                 baseline_survival_lp0 = round(1 - rf(0, t), 4),
                 baseline_survival_at_mean_lp = round(1 - rf(mean(lp), t), 4),
                 mean_lp_discovery = round(mean(lp), 4),
                 apparent_C_discovery = round(cindex(yt, lp), 4),
                 n_train = length(lp), events_train = sum(yt[, 2]))
    }))
  }))
}))
print(base, row.names = FALSE)
save_tsv(base, "22_locked_model_baseline.tsv")

cat("\n############ G. MODULE LOADINGS ############\n")
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
gene_names <- rbind(as.data.table(nets$mrna$gene_tbl)[, .(gene_id, gene_name)],
                    as.data.table(nets$lnc$gene_tbl)[, .(gene_id, gene_name)])
gene_names <- gene_names[!duplicated(gene_id)]
rm(nets)

loadings <- L$loadings
load_tbl <- rbindlist(lapply(names(loadings), function(m) {
  Lm <- loadings[[m]]
  data.table(module = m, gene_id = Lm$genes,
             loading = signif(Lm$rot[Lm$genes], 6),
             gene_center = signif(Lm$gene_center[Lm$genes], 6),
             gene_scale  = signif(Lm$gene_scale[Lm$genes], 6))
}))
load_tbl <- merge(load_tbl, gene_names, by = "gene_id", all.x = TRUE, sort = FALSE)
setcolorder(load_tbl, c("module", "gene_id", "gene_name", "loading", "gene_center", "gene_scale"))
load_tbl <- load_tbl[order(module, -abs(loading))]
msg("Loadings written for ", uniqueN(load_tbl$module), " modules, ",
    nrow(load_tbl), " gene rows")
save_tsv(load_tbl, "22_module_loadings.tsv")

const_tbl <- rbindlist(lapply(names(loadings), function(m) {
  Lm <- loadings[[m]]
  data.table(module = m,
             score_center = signif(Lm$score_center, 6),
             score_scale  = signif(Lm$score_scale, 6),
             var_explained = round(Lm$var_explained, 4),
             n_genes = length(Lm$genes),
             retained_clinical  = m %in% L$models$clinical$sel_me,
             retained_augmented = m %in% L$models$augmented$sel_me)
}))
print(const_tbl, row.names = FALSE)
save_tsv(const_tbl, "22_module_score_constants.tsv")

write_session_info("22_locked_model_export")
banner("22 | done")
