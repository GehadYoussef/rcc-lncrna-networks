# 12_model_diagnostics.R: diagnostics of the locked prognostic model.
#
# Run after 09. Every table has a `comparator` column: clinical (primary: age,
# sex, T, N, M1, grade) or augmented (clinical plus ESTIMATE and the STAR
# quality metrics). Sections:
#   A  calibration and Brier score at 1, 3 and 5 years, discovery baseline hazard
#   B  module hazard ratios with and without technical adjustment
#   C  full multivariable tables
#   D  networks rebuilt inside each CV fold
#   E  proportional hazards of the risk scores
#   F  decision curves in CPTAC-3
# Inputs (cache): dataset.rds, networks.rds, locked_model.rds,
# estimate_scores.rds, star_qc.rds.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet); library(WGCNA)
  library(ggplot2)
})
banner("12 | Model diagnostics")
set.seed(SEED)
enableWGCNAThreads(N_THREADS)

ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
L    <- readRDS(file.path(CACHE_DIR, "locked_model.rds"))
EST  <- as.data.table(readRDS(file.path(CACHE_DIR, "estimate_scores.rds")))
QC   <- as.data.table(readRDS(file.path(CACHE_DIR, "star_qc.rds")))
cohort <- as.data.table(ds$cohort_full)
stopifnot(identical(L$version, "v9"))

COMPARATORS <- c("clinical", "augmented")
lpv <- function(X, b) as.numeric(X[, names(b), drop = FALSE] %*% b)

# One entry per cohort x comparator x standardisation, holding the linear
# predictors of the comparator ("comp") and comparator + modules ("full")
# models and their discovery baselines. Discovery-fixed standardisation is
# used for the clinical comparator only, because composition and STAR metric
# scales do not transport across library protocols.
arms <- list()
for (cmp in COMPARATORS) {
  Xt <- L$X[[cmp]]$Xt; Xv <- L$X[[cmp]]$Xv
  bC <- L$models[[cmp]]$b_clin; bF <- L$models[[cmp]]$b_full
  lp_t <- list(comp = lpv(Xt, bC), full = lpv(cbind(Xt, L$Et), bF))
  rf   <- list(comp = baseline_risk_fun(L$yt, lp_t$comp),
               full = baseline_risk_fun(L$yt, lp_t$full))
  arms[[length(arms) + 1]] <- list(cohort = "discovery", comparator = cmp,
    standardisation = "cohort", y = L$yt, lp = lp_t, rf = rf)
  arms[[length(arms) + 1]] <- list(cohort = "validation", comparator = cmp,
    standardisation = "cohort", y = L$yv,
    lp = list(comp = lpv(Xv, bC), full = lpv(cbind(Xv, L$Ev), bF)), rf = rf)
  if (cmp == "clinical") {
    Xd <- L$X$clinical$Xv_disc
    arms[[length(arms) + 1]] <- list(cohort = "validation", comparator = cmp,
      standardisation = "discovery", y = L$yv,
      lp = list(comp = lpv(Xd, bC), full = lpv(cbind(Xd, L$Ev_disc), bF)), rf = rf)
  }
}
arm_id <- function(a) data.table(cohort = a$cohort, comparator = a$comparator,
                                 standardisation = a$standardisation)
msg("Locked sets: discovery n = ", length(L$tsamp), " (events ", sum(L$yt[, 2]),
    "), CPTAC-3 n = ", length(L$vsamp), " (events ", sum(L$yv[, 2]), ")")

# ---- A. calibration and Brier score ----
banner("A | Calibration")
# brier_clinical is the Brier score of the comparator named in the row. The
# null model gives every patient the mean predicted risk of the combined model.
brier_rows <- list(); cal_rows <- list(); cs_rows <- list()
for (a in arms) {
  for (yr in EVAL_TIMES_YRS) {
    t  <- yr * 365.25
    pf <- a$rf$full(a$lp$full, t); pc <- a$rf$comp(a$lp$comp, t)
    brier_rows[[length(brier_rows) + 1]] <- cbind(arm_id(a), data.table(
      years = yr, primary = yr == PRIMARY_HORIZON_YR,
      n = length(pf), events = sum(a$y[, 2]),
      brier_clinical = round(brier_ipcw(a$y, pc, t), 4),
      brier_combined = round(brier_ipcw(a$y, pf, t), 4),
      brier_null     = round(brier_ipcw(a$y, rep(mean(pf), length(pf)), t), 4)))
    for (mdl in c("comparator", "comparator + modules")) {
      pr <- if (mdl == "comparator") pc else pf
      lp <- if (mdl == "comparator") a$lp$comp else a$lp$full
      cb <- calibration_bins(a$y, pr, t)
      cal_rows[[length(cal_rows) + 1]] <- cbind(arm_id(a)[rep(1L, nrow(cb))], data.table(
        model = mdl, years = yr, primary = yr == PRIMARY_HORIZON_YR,
        group = cb$bin, n = cb$n, events_by_t = cb$events_by_t,
        predicted = cb$predicted, observed = cb$observed,
        obs_lo = cb$obs_lo, obs_hi = cb$obs_hi))
      cs <- calibration_summary(a$y, lp, pr, t)
      cs_rows[[length(cs_rows) + 1]] <- cbind(arm_id(a), data.table(
        model = mdl, years = yr, primary = yr == PRIMARY_HORIZON_YR,
        slope = round(cs[["slope"]], 4), slope_lo = round(cs[["slope_lo"]], 4),
        slope_hi = round(cs[["slope_hi"]], 4),
        observed = round(cs[["observed"]], 4), expected = round(cs[["expected"]], 4),
        OE = round(cs[["OE"]], 4), n_at_risk = cs[["n_at_risk"]]))
    }
  }
}
brier <- rbindlist(brier_rows); cal <- rbindlist(cal_rows); calsum <- rbindlist(cs_rows)
save_tsv(brier, "12_brier_scores.tsv"); print(brier)
save_tsv(cal, "12_calibration_bins.tsv")
save_tsv(calsum, "12_calibration_summary.tsv")
print(calsum[primary == TRUE])

cal_fig <- cal[comparator == "clinical" & standardisation == "cohort" &
               model == "comparator + modules"]
gcal <- ggplot(cal_fig, aes(predicted, observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = obs_lo, ymax = obs_hi), width = 0.01,
                colour = "grey60") +
  geom_point(size = 1.8, colour = "#B2182B") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, colour = "#2166AC") +
  facet_grid(cohort ~ years, labeller = labeller(years = function(x)
    paste0(x, " year"))) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Calibration of the locked risk model (clinical comparator + modules)",
       subtitle = "Quintiles of predicted risk; observed from Kaplan-Meier. Baseline hazard fixed in discovery.",
       x = "Predicted risk of death", y = "Observed risk of death") +
  theme_bw()
save_fig(gcal, "12_calibration", 8, 5.5)

# ---- B. sensitivity to the technical adjustment ----
banner("B | With vs without library-quality adjustment")
# Same module membership in every arm:
#   observed                      eigengene on the observed matrix
#   observed_plus_STAR_covariates as observed, with the STAR metrics as covariates
#   residualised                  eigengene on the residualised matrix
# Each arm is fitted under the clinical and clinical + ESTIMATE covariate sets,
# on one complete-case set per covariate set.
ARMS     <- c("observed", "observed_plus_STAR_covariates", "residualised")
COV_SETS <- list(clinical = CLIN_TERMS,
                 clinical_ESTIMATE = c(CLIN_TERMS, "stromal", "immune"))
STAR_TERMS <- c("noFeature", "multimap", "libsize_z")

tech_sens <- function(net, tag) {
  E_obs <- obs_expr(net); samp <- rownames(E_obs)
  cl  <- cohort[match(samp, sample_barcode)]
  est <- EST[match(samp, sample_barcode)]
  qc  <- QC[match(samp, sample_barcode)]
  stopifnot(identical(cl$sample_barcode, samp), !anyNA(qc$pct_noFeature),
            !anyNA(est$StromalScore))
  E_res <- remove_technical(E_obs, tech_covariates(qc))
  S <- list(
    observed     = score_modules(E_obs, fit_module_loadings(E_obs, net$gene_tbl, "")),
    residualised = score_modules(E_res, fit_module_loadings(E_res, net$gene_tbl, "")))
  S$observed_plus_STAR_covariates <- S$observed
  mods <- intersect(colnames(S$observed), colnames(S$residualised))
  X <- clinical_design(cl, "augmented", est, qc)
  rbindlist(lapply(names(COV_SETS), function(cs) {
    keep <- complete.cases(X[, c(COV_SETS[[cs]], STAR_TERMS), drop = FALSE]) &
            !is.na(cl$os_time) & !is.na(cl$os_event)
    msg("  ", tag, " | ", cs, ": n = ", sum(keep), ", events = ",
        sum(cl$os_event[keep]))
    rbindlist(lapply(ARMS, function(arm) {
      cov <- c(COV_SETS[[cs]],
               if (arm == "observed_plus_STAR_covariates") STAR_TERMS)
      r <- rbindlist(lapply(mods, function(m) {
        d <- data.frame(os_time = cl$os_time[keep], os_event = cl$os_event[keep],
                        ME = S[[arm]][keep, m], X[keep, cov, drop = FALSE])
        s <- summary(coxph(Surv(os_time, os_event) ~ ., data = d))
        data.table(biotype = tag, module = m, covariate_set = cs, arm = arm,
                   n = s$n, events = s$nevent,
                   HR = s$conf.int["ME", "exp(coef)"],
                   lo = s$conf.int["ME", "lower .95"],
                   hi = s$conf.int["ME", "upper .95"],
                   p  = s$coefficients["ME", "Pr(>|z|)"])
      }))
      r[, fdr := p.adjust(p, "BH")][]
    }))
  }))
}
sens <- rbind(tech_sens(nets$mrna, "mRNA"), tech_sens(nets$lnc, "lncRNA"))
save_tsv(sens[, .(biotype, module, covariate_set, arm, n, events,
                  HR = round(HR, 3), lo = round(lo, 3), hi = round(hi, 3),
                  p = signif(p, 3), fdr = signif(fdr, 3))],
         "12_technical_adjustment_arms.tsv")

wide <- dcast(sens, biotype + module + covariate_set + n + events ~ arm,
              value.var = c("HR", "fdr"))
setnames(wide,
         c("HR_observed", "fdr_observed", "HR_residualised", "fdr_residualised",
           "HR_observed_plus_STAR_covariates", "fdr_observed_plus_STAR_covariates"),
         c("HR_unadjusted", "fdr_unadjusted", "HR_adjusted", "fdr_adjusted",
           "HR_observed_STAR", "fdr_observed_STAR"))
wide[, sig_adjusted      := fdr_adjusted      < FDR_ALPHA]
wide[, sig_unadjusted    := fdr_unadjusted    < FDR_ALPHA]
wide[, sig_observed_STAR := fdr_observed_STAR < FDR_ALPHA]
wide[, status := fifelse(sig_adjusted & sig_unadjusted, "robust",
                 fifelse(!sig_adjusted & sig_unadjusted, "LOST on adjustment",
                 fifelse(sig_adjusted & !sig_unadjusted, "gained on adjustment",
                         "not significant")))]
setorder(wide, covariate_set, biotype, fdr_adjusted)
save_tsv(wide[, .(biotype, module, covariate_set, n, events,
                  HR_unadjusted = round(HR_unadjusted, 3),
                  fdr_unadjusted = signif(fdr_unadjusted, 3),
                  HR_observed_STAR = round(HR_observed_STAR, 3),
                  fdr_observed_STAR = signif(fdr_observed_STAR, 3),
                  HR_adjusted = round(HR_adjusted, 3),
                  fdr_adjusted = signif(fdr_adjusted, 3),
                  status, sig_observed_STAR)],
         "12_technical_adjustment_sensitivity.tsv")
print(wide[, .(biotype, module, covariate_set, HR_unadjusted = round(HR_unadjusted, 2),
               HR_observed_STAR = round(HR_observed_STAR, 2),
               HR_adjusted = round(HR_adjusted, 2), status)])
for (cs in names(COV_SETS)) {
  w <- wide[covariate_set == cs]
  msg(cs, " set -- significant WITHOUT adjustment: ", sum(w$sig_unadjusted),
      "; with STAR covariates: ", sum(w$sig_observed_STAR),
      "; WITH residualisation: ", sum(w$sig_adjusted),
      "; lost: ", sum(w$status == "LOST on adjustment"),
      "; gained: ", sum(w$status == "gained on adjustment"))
}

gs <- ggplot(wide[covariate_set == "clinical"],
             aes(HR_unadjusted, HR_adjusted, colour = status)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_hline(yintercept = 1, linewidth = 0.2, colour = "grey80") +
  geom_vline(xintercept = 1, linewidth = 0.2, colour = "grey80") +
  geom_point(size = 2.2) + facet_wrap(~ biotype) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Effect of library-quality adjustment on module hazard ratios",
       subtitle = "Same module membership both axes; clinical covariate set; only the expression adjustment differs",
       x = "HR without adjustment", y = "HR with adjustment", colour = NULL) +
  theme_bw() + theme(legend.position = "bottom")
save_fig(gs, "12_technical_adjustment_sensitivity", 8.5, 4.8)

# ---- C. full multivariable tables ----
banner("C | Full multivariable models")
# Locked terms refitted unpenalised in each cohort. CPTAC-3 rows have a low
# events-per-parameter ratio (`epv`) and are descriptive.
full_tbl <- function(cohort_name, cmp, model_name, X, E, y, b) {
  d <- as.data.frame(cbind(X, E))
  keep <- names(b)[names(b) %in% colnames(d)]
  d <- data.frame(os_time = y[, 1], os_event = y[, 2], d[, keep, drop = FALSE])
  f <- coxph(Surv(os_time, os_event) ~ ., data = d)
  s <- summary(f)
  data.table(cohort = cohort_name, comparator = cmp, model = model_name,
             n = s$n, events = s$nevent, n_terms = length(keep),
             epv = round(s$nevent / length(keep), 1),
             term = rownames(s$conf.int),
             HR = round(s$conf.int[, "exp(coef)"], 3),
             lo = round(s$conf.int[, "lower .95"], 3),
             hi = round(s$conf.int[, "upper .95"], 3),
             p  = signif(s$coefficients[, "Pr(>|z|)"], 3))
}
mv <- rbindlist(lapply(COMPARATORS, function(cmp) rbind(
  full_tbl("TCGA-KIRC", cmp, "comparator", L$X[[cmp]]$Xt, NULL, L$yt,
           L$models[[cmp]]$b_clin),
  full_tbl("TCGA-KIRC", cmp, "comparator + modules", L$X[[cmp]]$Xt, L$Et, L$yt,
           L$models[[cmp]]$b_full),
  full_tbl("CPTAC-3", cmp, "comparator", L$X[[cmp]]$Xv, NULL, L$yv,
           L$models[[cmp]]$b_clin),
  full_tbl("CPTAC-3", cmp, "comparator + modules", L$X[[cmp]]$Xv, L$Ev, L$yv,
           L$models[[cmp]]$b_full))))
save_tsv(mv, "12_full_multivariable_models.tsv")
print(mv[model == "comparator + modules" & comparator == "clinical"])

# ---- D. fold-wise network rebuild ----
banner("D | Modules recomputed inside CV folds")
# Analysis set: the discovery locked set from 09. Within each training fold:
# residualisation, both networks (stage 02 parameters), module loadings and
# the elastic net (comparator unpenalised) are fitted on the training fold and
# applied to the held-out fold. Sample QC from 02 is not repeated.
tsamp <- L$tsamp
tcl   <- as.data.table(L$tcl)[match(tsamp, sample_barcode)]
stopifnot(identical(tcl$sample_barcode, tsamp))
est_t <- EST[match(tsamp, sample_barcode)]
qc_t  <- QC[match(tsamp, sample_barcode)]
# Standardising unpenalised covariates only shifts the linear predictor by a
# constant, so comparator designs are built once. Penalised module scores use
# training-fold constants.
XD <- list(clinical  = clinical_design(tcl, "clinical"),
           augmented = clinical_design(tcl, "augmented", est_t, qc_t))
ok_d <- complete.cases(XD$augmented) & !is.na(tcl$os_time) & !is.na(tcl$os_event)
if (!all(ok_d)) {
  msg("Dropping ", sum(!ok_d), " samples with an incomplete augmented design")
  tsamp <- tsamp[ok_d]; tcl <- tcl[ok_d]; qc_t <- qc_t[ok_d]
  XD <- lapply(XD, function(X) X[ok_d, , drop = FALSE])
}
y_d   <- Surv(tcl$os_time, tcl$os_event)
cov_d <- tech_covariates(qc_t)
stopifnot(all(tsamp %in% rownames(obs_expr(nets$mrna))),
          all(tsamp %in% rownames(obs_expr(nets$lnc))))
E_obs <- list(mrna = obs_expr(nets$mrna)[tsamp, , drop = FALSE],
              lnc  = obs_expr(nets$lnc)[tsamp, , drop = FALSE])
NET_PAR <- list(
  mrna = list(prefix = "mRNA_ME", power = nets$mrna$power, deep_split = DEEP_SPLIT,
              merge_cut = MERGE_CUT_HEIGHT, min_kme = MIN_KME_TO_STAY),
  lnc  = list(prefix = "lnc_ME", power = LNC_POWER,
              deep_split = if (is.na(LNC_DEEPSPLIT)) DEEP_SPLIT else LNC_DEEPSPLIT,
              merge_cut  = if (is.na(LNC_MERGE)) MERGE_CUT_HEIGHT else LNC_MERGE,
              min_kme    = if (is.na(LNC_MINKME)) MIN_KME_TO_STAY else LNC_MINKME))
msg("Fold-wise analysis set: n = ", length(tsamp), ", events = ", sum(y_d[, 2]),
    "; K = ", FOLDWISE_K, ", assignments = ", FOLDWISE_N_ASSIGN)

rebuild_modules <- function(E, par) {
  set.seed(SEED)
  net <- blockwiseModules(E, power = par$power, networkType = NETWORK_TYPE,
                          TOMType = TOM_TYPE, minModuleSize = MIN_MODULE_SIZE,
                          mergeCutHeight = par$merge_cut, deepSplit = par$deep_split,
                          minKMEtoStay = par$min_kme, numericLabels = FALSE,
                          pamRespectsDendro = FALSE, maxBlockSize = MAX_BLOCK_SIZE,
                          saveTOMs = FALSE, verbose = 0)
  data.table(gene_id = names(net$colors), module = unname(net$colors))
}

fold_rows <- list(); assign_rows <- list()
for (a in seq_len(FOLDWISE_N_ASSIGN)) {
  t_assign <- Sys.time()
  set.seed(SEED + a)
  folds <- sample(rep(seq_len(FOLDWISE_K), length.out = length(tsamp)))
  oof_comp <- lapply(COMPARATORS, function(z) rep(NA_real_, length(tsamp)))
  oof_comb <- oof_comp; names(oof_comp) <- names(oof_comb) <- COMPARATORS
  for (k in seq_len(FOLDWISE_K)) {
    t_fold <- Sys.time()
    tr <- folds != k; te <- !tr
    S_tr <- list(); S_te <- list(); n_mod <- integer(2); names(n_mod) <- names(E_obs)
    for (nm in names(E_obs)) {
      tf   <- fit_technical(E_obs[[nm]][tr, , drop = FALSE], cov_d[tr, , drop = FALSE])
      Etr  <- apply_technical(E_obs[[nm]][tr, , drop = FALSE], tf, cov_d[tr, , drop = FALSE])
      Ete  <- apply_technical(E_obs[[nm]][te, , drop = FALSE], tf, cov_d[te, , drop = FALSE])
      gt   <- rebuild_modules(Etr, NET_PAR[[nm]])
      n_mod[nm] <- length(setdiff(unique(gt$module), "grey"))
      Ltr  <- fit_module_loadings(Etr, gt, NET_PAR[[nm]]$prefix)
      S_tr[[nm]] <- score_modules(Etr, Ltr, gene_standardise = "discovery",
                                  score_standardise = "discovery")
      S_te[[nm]] <- score_modules(Ete, Ltr, gene_standardise = "discovery",
                                  score_standardise = "discovery")
    }
    S_tr <- do.call(cbind, S_tr); S_te <- do.call(cbind, S_te)
    stopifnot(identical(colnames(S_tr), colnames(S_te)))
    y_tr <- y_d[tr]
    n_ret   <- setNames(integer(length(COMPARATORS)),   COMPARATORS)
    dropped <- setNames(character(length(COMPARATORS)), COMPARATORS)
    for (cmp in COMPARATORS) {
      Xtr <- XD[[cmp]][tr, , drop = FALSE]; Xte <- XD[[cmp]][te, , drop = FALSE]
      fc <- coxph(y_tr ~ ., data = as.data.frame(Xtr))
      # A rare indicator (e.g. N_pos) can be constant in a training fold. It is
      # dropped from both models for that fold and recorded in `dropped_terms`.
      bad <- is.na(coef(fc))
      if (any(bad)) {
        dropped[cmp] <- paste(names(coef(fc))[bad], collapse = ";")
        warning("Comparator '", cmp, "': dropping inestimable term(s) ", dropped[cmp],
                " in assignment ", a, ", fold ", k, call. = FALSE, immediate. = TRUE)
        Xtr <- Xtr[, !bad, drop = FALSE]; Xte <- Xte[, !bad, drop = FALSE]
        fc  <- coxph(y_tr ~ ., data = as.data.frame(Xtr))
        stopifnot(!anyNA(coef(fc)))
      }
      oof_comp[[cmp]][te] <- lpv(Xte, coef(fc))
      if (ncol(S_tr) == 0) {
        oof_comb[[cmp]][te] <- oof_comp[[cmp]][te]; n_ret[cmp] <- 0L; next
      }
      set.seed(SEED)
      fit <- cv.glmnet(cbind(Xtr, S_tr), y_tr, family = "cox", alpha = ML_ALPHA,
                       nfolds = 10, nlambda = ML_NLAMBDA,
                       lambda.min.ratio = ML_LAMBDA_MIN_RATIO,
                       thresh = ML_THRESH, maxit = ML_MAXIT,
                       penalty.factor = c(rep(0, ncol(Xtr)), rep(1, ncol(S_tr))))
      newx <- cbind(Xte, S_te)
      stopifnot(identical(colnames(newx), colnames(cbind(Xtr, S_tr))))
      oof_comb[[cmp]][te] <- as.numeric(predict(fit, newx = newx, s = "lambda.min"))
      bf <- as.matrix(coef(fit, s = "lambda.min"))[, 1]
      n_ret[cmp] <- sum(bf[colnames(S_tr)] != 0)
    }
    el <- as.numeric(difftime(Sys.time(), t_fold, units = "mins"))
    # One row per fold and comparator. Both share the fold's network rebuild.
    for (cmp in COMPARATORS) fold_rows[[length(fold_rows) + 1]] <- data.table(
      assignment = a, fold = k, comparator = cmp, n_train = sum(tr), n_test = sum(te),
      events_train = sum(y_tr[, 2]),
      n_modules_mRNA = n_mod[["mrna"]], n_modules_lncRNA = n_mod[["lnc"]],
      n_retained = n_ret[[cmp]], dropped_terms = dropped[[cmp]],
      elapsed_min = round(el, 2))
    msg(sprintf("  assignment %d/%d, fold %d/%d: %d mRNA + %d lncRNA modules, %.1f min",
                a, FOLDWISE_N_ASSIGN, k, FOLDWISE_K, n_mod[["mrna"]], n_mod[["lnc"]], el))
  }
  el_a <- as.numeric(difftime(Sys.time(), t_assign, units = "mins"))
  for (cmp in COMPARATORS) {
    Cc <- cindex(y_d, oof_comp[[cmp]]); Cf <- cindex(y_d, oof_comb[[cmp]])
    assign_rows[[length(assign_rows) + 1]] <- data.table(
      assignment = a, seed = SEED + a, comparator = cmp,
      n = length(tsamp), events = sum(y_d[, 2]),
      C_comparator_oof = round(Cc, 4), C_combined_oof = round(Cf, 4),
      delta = round(Cf - Cc, 4), elapsed_min = round(el_a, 1))
  }
  msg(sprintf("  assignment %d done in %.1f min", a, el_a))
}
per_fold <- rbindlist(fold_rows); per_assign <- rbindlist(assign_rows)
save_tsv(per_fold, "12_foldwise_per_fold.tsv")
save_tsv(per_assign, "12_foldwise_per_assignment.tsv")
print(per_assign)

# Repeated-CV estimate from stage 05 (modules fixed on the full cohort).
read_cv05 <- function(cmp) {
  out <- list(delta = NA_real_, sd = NA_real_, lo = NA_real_, hi = NA_real_,
              n = NA_integer_, C_c = NA_real_, C_f = NA_real_)
  f <- file.path(RESULTS_DIR, "05_delta_cindex.tsv")
  if (!file.exists(f)) { msg("05_delta_cindex.tsv not found; row left empty"); return(out) }
  d <- fread(f)
  r <- d[gsub(" ", "", comparison) == paste0(cmp, "_eig-", cmp)]
  if (nrow(r) != 1) { msg("No ", cmp, "_eig - ", cmp, " row in 05_delta_cindex.tsv"); return(out) }
  out$delta <- r$delta_mean
  if ("pct2.5"  %in% names(r)) out$lo <- r[["pct2.5"]]
  if ("pct97.5" %in% names(r)) out$hi <- r[["pct97.5"]]
  if ("n_repeats" %in% names(r)) out$n <- r$n_repeats
  fr <- file.path(RESULTS_DIR, "05_cv_cindex_per_repeat.tsv")
  if (file.exists(fr)) {
    p <- fread(fr); eig <- paste0(cmp, "_eig")
    if (all(c(cmp, eig) %in% names(p))) {
      out$C_c <- mean(p[[cmp]]); out$C_f <- mean(p[[eig]]); out$sd <- sd(p[[eig]] - p[[cmp]])
    }
  }
  out
}
foldwise <- rbindlist(lapply(COMPARATORS, function(cmp) {
  ap <- arms[[which(vapply(arms, function(a) a$cohort == "discovery" &&
                                            a$comparator == cmp, logical(1)))]]
  bt <- paired_boot_delta_c(L$yt, ap$lp$full, ap$lp$comp)
  cv <- read_cv05(cmp)
  pa <- per_assign[comparator == cmp]
  rbind(
    data.table(analysis = "modules fixed on full cohort (apparent)", comparator = cmp,
               C_clinical = round(cindex(L$yt, ap$lp$comp), 4),
               C_combined = round(cindex(L$yt, ap$lp$full), 4),
               delta_C = round(bt[["delta"]], 4), delta_sd = NA_real_,
               delta_lo = round(bt[["lo"]], 4), delta_hi = round(bt[["hi"]], 4),
               n_assignments = NA_integer_,
               interval_type = "paired patient bootstrap 95% (in-sample)"),
    data.table(analysis = "modules fixed on full cohort (10x10 CV, from 05)", comparator = cmp,
               C_clinical = round(cv$C_c, 4), C_combined = round(cv$C_f, 4),
               delta_C = round(cv$delta, 4), delta_sd = round(cv$sd, 4),
               delta_lo = round(cv$lo, 4), delta_hi = round(cv$hi, 4),
               n_assignments = cv$n,
               interval_type = "2.5 and 97.5 percentiles across repeats"),
    data.table(analysis = sprintf("modules recomputed inside each fold (%d x %d assignments)",
                                  FOLDWISE_K, FOLDWISE_N_ASSIGN), comparator = cmp,
               C_clinical = round(mean(pa$C_comparator_oof), 4),
               C_combined = round(mean(pa$C_combined_oof), 4),
               delta_C = round(mean(pa$delta), 4), delta_sd = round(sd(pa$delta), 4),
               delta_lo = round(min(pa$delta), 4), delta_hi = round(max(pa$delta), 4),
               n_assignments = nrow(pa),
               interval_type = "range across fold assignments"))
}))
save_tsv(foldwise, "12_foldwise_module_sensitivity.tsv")
print(foldwise)

# ---- E. proportional hazards of the risk scores ----
banner("E | Proportional hazards, final risk scores")
COHORT_LAB <- c(discovery = "TCGA-KIRC", validation = "CPTAC-3")
ph <- rbindlist(lapply(arms, function(a) rbindlist(lapply(
  c("comparator", "comparator + modules"), function(mdl) {
    sc <- if (mdl == "comparator") a$lp$comp else a$lp$full
    f  <- coxph(a$y ~ score, data = data.frame(score = sc))
    tb <- cox.zph(f)$table
    data.table(cohort = COHORT_LAB[[a$cohort]], comparator = a$comparator,
               standardisation = a$standardisation, model = mdl,
               n = length(sc), events = sum(a$y[, 2]),
               term = rownames(tb), chisq = round(tb[, "chisq"], 3),
               df = tb[, "df"], p = signif(tb[, "p"], 3))
  }))))
save_tsv(ph, "12_proportional_hazards_riskscore.tsv")
print(ph)

# ---- F. decision curves in CPTAC-3 (discovery baseline) ----
banner("F | Decision curves, CPTAC-3")
# Summary: median net-benefit gain of the combined model over its comparator
# in two threshold windows, with a paired patient bootstrap.
thr  <- seq(0.02, 0.60, by = 0.01)
WIN  <- list(`0.10-0.40` = round(thr, 2) >= 0.10 & round(thr, 2) <= 0.40,
             `0.05-0.30` = round(thr, 2) >= 0.05 & round(thr, 2) <= 0.30)
DCA_B <- 500
dca_rows <- list(); dsum_rows <- list()
for (a in arms[vapply(arms, function(z) z$cohort == "validation", logical(1))]) {
  for (yr in c(PRIMARY_HORIZON_YR, SECONDARY_HORIZON_YR)) {
    t <- yr * 365.25
    rc <- a$rf$comp(a$lp$comp, t); rf_ <- a$rf$full(a$lp$full, t)
    nbc <- net_benefit(a$y, rc, t, thr); nbf <- net_benefit(a$y, rf_, t, thr)
    id <- data.table(baseline = "discovery", comparator = a$comparator,
                     standardisation = a$standardisation, years = yr)
    idt <- id[rep(1L, length(thr))]
    dca_rows[[length(dca_rows) + 1]] <- rbind(
      cbind(idt, data.table(threshold = thr, strategy = "comparator", net_benefit = nbc)),
      cbind(idt, data.table(threshold = thr, strategy = "comparator + modules", net_benefit = nbf)),
      cbind(idt, data.table(threshold = thr, strategy = "treat all",
                            net_benefit = net_benefit_treat_all(a$y, t, thr))),
      cbind(idt, data.table(threshold = thr, strategy = "treat none", net_benefit = 0)))
    gain <- nbf - nbc
    set.seed(SEED); n <- length(rc)
    bmed <- vapply(seq_len(DCA_B), function(b) {
      i  <- sample.int(n, n, replace = TRUE); yi <- a$y[i]
      gi <- net_benefit(yi, rf_[i], t, thr) - net_benefit(yi, rc[i], t, thr)
      vapply(WIN, function(s) median(gi[s]), numeric(1))
    }, numeric(length(WIN)))
    rownames(bmed) <- names(WIN)
    for (w in names(WIN)) {
      bw <- bmed[w, ]; bw <- bw[is.finite(bw)]
      dsum_rows[[length(dsum_rows) + 1]] <- cbind(id, data.table(
        window = w, n_thresholds = sum(WIN[[w]]),
        median_gain = round(median(gain[WIN[[w]]]), 4),
        n_positive = sum(gain[WIN[[w]]] > 0),
        boot_lo = round(unname(quantile(bw, 0.025)), 4),
        boot_hi = round(unname(quantile(bw, 0.975)), 4),
        p_boot = signif(2 * min(mean(bw <= 0), mean(bw >= 0)), 3),
        n_boot = length(bw)))
    }
  }
}
dca <- rbindlist(dca_rows); dsum <- rbindlist(dsum_rows)
save_tsv(dca, "22_validation_decision_curve.tsv")
save_tsv(dca, "12_validation_decision_curve.tsv")
save_tsv(dsum, "12_validation_decision_curve_summary.tsv")
print(dsum)

dca_fig <- dca[comparator == "clinical" & standardisation == "cohort"]
gd <- ggplot(dca_fig, aes(threshold, net_benefit, colour = strategy, linetype = strategy)) +
  geom_line(linewidth = 0.5) + facet_wrap(~ years, labeller = labeller(years = function(x)
    paste0(x, " years"))) +
  coord_cartesian(ylim = c(-0.02, max(dca_fig$net_benefit, na.rm = TRUE) * 1.1)) +
  scale_colour_manual(values = c(`comparator + modules` = "#E66101", comparator = "#4D4D4D",
                                 `treat all` = "grey45", `treat none` = "grey70"), name = NULL) +
  scale_linetype_manual(values = c(`comparator + modules` = 1, comparator = 1,
                                   `treat all` = 2, `treat none` = 3), name = NULL) +
  labs(title = "Decision curves in CPTAC-3 (clinical comparator, discovery baseline)",
       x = "Threshold probability", y = "Net benefit") +
  theme_bw() + theme(legend.position = "bottom")
save_fig(gd, "12_validation_decision_curve", 8, 4.5)

write_session_info("12_model_diagnostics")
banner("12 | done")
