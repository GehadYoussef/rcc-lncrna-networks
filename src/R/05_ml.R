# 05_ml.R: cross-validated incremental value of co-expression modules.
#
# Penalised Cox models in TCGA-KIRC, scored out of fold over repeated K-fold CV.
# Comparators (unpenalised): clinical (age, sex, T, N, M1, grade) and augmented
# (plus ESTIMATE and the STAR quality metrics). Module representations
# (elastic net): "eig", the discovery module scores on the residualised
# matrices, and "hub", lncRNA hub genes on observed expression with the module
# pool chosen inside each training fold.
# Inputs: dataset.rds, networks.rds, estimate_scores.rds, star_qc.rds.
# Outputs: C-index per repeat, paired-bootstrap increments, time-dependent AUC,
# decision curves, final coefficients and a likelihood-ratio test. External
# validation is in 09.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet)
  library(timeROC); library(ggplot2)
  library(foreach); library(doParallel)
})
banner("05 | Penalised Cox models: incremental value of module scores")
set.seed(SEED)

ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
est  <- as.data.table(readRDS(file.path(CACHE_DIR, "estimate_scores.rds")))
qc   <- as.data.table(readRDS(file.path(CACHE_DIR, "star_qc.rds")))
cohort <- as.data.table(ds$cohort_full)

# ---- analysis set ----
# Samples in both networks with complete clinical inputs and survival, so every
# arm uses the same patients.
CLIN_INPUTS <- c("os_time", "os_event", "age", "sex", "T_stage", "N_pos", "M1", "grade_num")
in_nets  <- intersect(rownames(nets$lnc$expr), rownames(nets$mrna$expr))
cand     <- cohort[sample_barcode %in% in_nets]
complete <- complete.cases(cand[, ..CLIN_INPUTS])
msg("Samples in both networks: ", nrow(cand),
    "; incomplete clinical inputs or survival: ", sum(!complete))
cl <- cand[complete]

aux <- merge(est, qc, by = "sample_barcode")
aux <- aux[match(cl$sample_barcode, sample_barcode)]
stopifnot(identical(aux$sample_barcode, cl$sample_barcode))

# Drop rows with incomplete comparator designs. The count is logged.
Xc <- clinical_design(cl, "clinical")
Xa <- clinical_design(cl, "augmented", est = aux, qc = aux)
ok <- complete.cases(Xc) & complete.cases(Xa)
msg("Rows with incomplete augmented inputs: ", sum(complete.cases(Xc) & !complete.cases(Xa)))
cl <- cl[ok]; Xc <- Xc[ok, , drop = FALSE]; Xa <- Xa[ok, , drop = FALSE]
samp <- cl$sample_barcode
rownames(Xc) <- samp; rownames(Xa) <- samp

y  <- Surv(cl$os_time, cl$os_event)
cl_event <- cl$os_event          # plain vector for the parallel workers
n  <- length(samp)
msg("Analysis set: n = ", n, ", events = ", sum(cl_event))

save_tsv(data.table(
  step = c("in_both_networks", "incomplete_clinical_inputs_or_survival",
           "incomplete_augmented_inputs", "analysis_set", "events"),
  n    = c(nrow(cand), sum(!complete), sum(complete) - n, n, sum(cl_event))),
  "05_analysis_set.tsv")

# ---- module scores (eig representation) ----
# Unit-variance discovery scores on the residualised matrices (00_config.R).
sc <- discovery_scores(nets)
Xe <- cbind(sc$mrna[match(samp, rownames(sc$mrna)), , drop = FALSE],
            sc$lnc [match(samp, rownames(sc$lnc)),  , drop = FALSE])
rownames(Xe) <- samp
stopifnot(!anyNA(Xe))
lnc_cols <- grep("^lnc_ME", colnames(Xe), value = TRUE)
msg("Module scores available as predictors: ", ncol(Xe),
    " (", ncol(Xe) - length(lnc_cols), " protein-coding, ", length(lnc_cols), " lncRNA)")

# ---- observed lncRNA expression (hub representation) ----
# Observed log2(FPKM + 1), so the predictor is computable for a single patient.
# Candidate genes per module are ordered by kME.
Xo <- obs_expr(nets$lnc)[samp, , drop = FALSE]
gt <- as.data.frame(nets$lnc$gene_tbl[module != "grey" & gene_id %in% colnames(Xo)])
gt <- gt[order(-gt$kME), ]
hub_pool <- split(gt, gt$module)
msg("Observed lncRNA matrix for the hub representation: ", ncol(Xo), " transcripts in ",
    length(hub_pool), " modules")

# ---- per-fold hub selection ----
# Training rows only: Cox model per lncRNA module (score + clinical terms), BH
# across modules, keep FDR < FDR_ALPHA (else the three smallest p), then the
# top ML_MAX_FEATURES genes by kME.
select_hub_genes <- function(tr) {
  p <- vapply(lnc_cols, function(m) {
    d <- data.frame(ME = Xe[tr, m], Xc[tr, , drop = FALSE])
    tryCatch(summary(coxph(y[tr] ~ ., data = d))$coefficients["ME", "Pr(>|z|)"],
             error = function(e) NA_real_)
  }, numeric(1))
  fdr <- p.adjust(p, "BH")
  sel <- names(p)[which(fdr < FDR_ALPHA)]
  fallback <- !length(sel)
  if (fallback) sel <- names(sort(p))[seq_len(min(3L, sum(is.finite(p))))]
  mods  <- intersect(sub("^lnc_ME", "", sel), names(hub_pool))
  if (!length(mods)) mods <- names(hub_pool)   # only if every per-module fit failed
  pool  <- do.call(rbind, hub_pool[mods])
  pool  <- pool[order(-pool$kME), ]
  list(modules = mods, genes = head(pool$gene_id, ML_MAX_FEATURES), fallback = fallback)
}

# Truncated lambda path: the Cox solver converges poorly on near-collinear
# features at small lambda. Comparator columns have penalty factor 0.
fit_cv <- function(x, yy, pfac) {
  cv.glmnet(x, yy, family = "cox", alpha = ML_ALPHA, nfolds = 5,
            nlambda = ML_NLAMBDA, lambda.min.ratio = ML_LAMBDA_MIN_RATIO,
            thresh = ML_THRESH, maxit = ML_MAXIT, penalty.factor = pfac)
}
pf_for <- function(Xcomp, Xmod) c(rep(0, ncol(Xcomp)), rep(1, ncol(Xmod)))

# Comparator linear predictor on held-out rows. A coefficient is NA only when
# its column is constant in the training fold, so it contributes zero.
lp_cox <- function(fit, Xte) {
  b <- coef(fit); b[is.na(b)] <- 0
  as.numeric(Xte[, names(b), drop = FALSE] %*% b)
}
lp_net <- function(fit, Xte) as.numeric(predict(fit, newx = Xte, s = "lambda.min"))

# ---- repeated cross-validation ----
ARMS <- c("clinical", "augmented", "clinical_eig", "augmented_eig",
          "clinical_hub", "augmented_hub")
msg("Running ", ML_N_REPEATS, " x ", ML_N_FOLDS,
    "-fold cross-validation (hub selection and lambda nested inside each fold) ...")

t0 <- Sys.time()
par_cl <- makeCluster(min(10L, N_THREADS))
registerDoParallel(par_cl)
clusterExport(par_cl, c("y", "Xc", "Xa", "Xe", "Xo", "hub_pool", "lnc_cols",
                        "select_hub_genes", "fit_cv", "pf_for", "lp_cox", "lp_net",
                        "FDR_ALPHA", "ML_MAX_FEATURES", "ML_ALPHA", "ML_NLAMBDA",
                        "ML_LAMBDA_MIN_RATIO", "ML_THRESH", "ML_MAXIT"))
CV_list <- foreach(r = seq_len(ML_N_REPEATS),
                   .packages = c("survival", "glmnet")) %dopar% {
  set.seed(SEED + r)
  out   <- matrix(NA_real_, n, length(ARMS), dimnames = list(samp, ARMS))
  hub   <- vector("list", ML_N_FOLDS)
  folds <- sample(rep(seq_len(ML_N_FOLDS), length.out = n))
  Xce <- cbind(Xc, Xe); Xae <- cbind(Xa, Xe)
  for (k in seq_len(ML_N_FOLDS)) {
    tr <- folds != k; te <- !tr
    if (sum(cl_event[tr]) < 5) next

    out[te, "clinical"]  <- lp_cox(coxph(y[tr] ~ ., data = as.data.frame(Xc[tr, , drop = FALSE])),
                                   Xc[te, , drop = FALSE])
    out[te, "augmented"] <- lp_cox(coxph(y[tr] ~ ., data = as.data.frame(Xa[tr, , drop = FALSE])),
                                   Xa[te, , drop = FALSE])

    out[te, "clinical_eig"]  <- lp_net(fit_cv(Xce[tr, ], y[tr], pf_for(Xc, Xe)), Xce[te, , drop = FALSE])
    out[te, "augmented_eig"] <- lp_net(fit_cv(Xae[tr, ], y[tr], pf_for(Xa, Xe)), Xae[te, , drop = FALSE])

    hs  <- select_hub_genes(tr)
    Xg  <- Xo[, hs$genes, drop = FALSE]
    Xch <- cbind(Xc, Xg); Xah <- cbind(Xa, Xg)
    out[te, "clinical_hub"]  <- lp_net(fit_cv(Xch[tr, ], y[tr], pf_for(Xc, Xg)), Xch[te, , drop = FALSE])
    out[te, "augmented_hub"] <- lp_net(fit_cv(Xah[tr, ], y[tr], pf_for(Xa, Xg)), Xah[te, , drop = FALSE])

    hub[[k]] <- data.frame(repeat_id = r, fold = k, n_train = sum(tr),
                           events_train = sum(cl_event[tr]),
                           n_modules_selected = length(hs$modules),
                           modules = paste(hs$modules, collapse = ";"),
                           n_genes = length(hs$genes), fallback = hs$fallback)
  }
  list(lp = out, hub = do.call(rbind, hub))
}
stopCluster(par_cl)
msg("Cross-validation finished in ",
    round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

LP <- array(NA_real_, dim = c(n, length(ARMS), ML_N_REPEATS),
            dimnames = list(samp, ARMS, NULL))
for (r in seq_len(ML_N_REPEATS)) LP[, , r] <- CV_list[[r]]$lp

hub_sel <- rbindlist(lapply(CV_list, `[[`, "hub"))
save_tsv(hub_sel, "05_hub_selection_per_fold.tsv")
hub_freq <- hub_sel[, .(module = unlist(strsplit(modules, ";"))), by = .(repeat_id, fold)][
  , .(n_folds_selected = .N, frac_folds_selected = .N / nrow(hub_sel)), by = module]
setorder(hub_freq, -n_folds_selected)
save_tsv(hub_freq, "05_hub_module_selection_frequency.tsv")
msg("Hub pool: ", sum(hub_sel$fallback), " of ", nrow(hub_sel),
    " folds fell back to the three smallest p; median genes per fold ",
    median(hub_sel$n_genes))
print(hub_freq)

# ---- C-index per repeat and on the repeat-averaged predictor ----
DELTAS <- list(c("clinical_eig",  "clinical"),
               c("augmented_eig", "augmented"),
               c("clinical_hub",  "clinical"),
               c("augmented_hub", "augmented"),
               c("augmented",     "clinical"))
delta_name <- function(z) paste0("delta_", z[1], "_vs_", z[2])

per_rep <- rbindlist(lapply(seq_len(ML_N_REPEATS), function(r) {
  x <- as.list(setNames(lapply(ARMS, function(m) cindex(y, LP[, m, r])), ARMS))
  cbind(data.table(repeat_id = r), as.data.table(x))
}))
for (z in DELTAS) set(per_rep, j = delta_name(z), value = per_rep[[z[1]]] - per_rep[[z[2]]])
save_tsv(per_rep, "05_cv_cindex_per_repeat.tsv")

lp_avg <- apply(LP, c(1, 2), mean, na.rm = TRUE)
ci_tbl <- rbindlist(lapply(ARMS, function(m) {
  cc <- cindex_ci(y, lp_avg[, m])
  data.table(model = m, n = n, events = sum(cl_event),
             C_cv_mean    = mean(per_rep[[m]]),
             C_cv_pct2.5  = unname(quantile(per_rep[[m]], 0.025)),
             C_cv_pct97.5 = unname(quantile(per_rep[[m]], 0.975)),
             C_cv_min     = min(per_rep[[m]]),
             C_cv_max     = max(per_rep[[m]]),
             n_repeats    = nrow(per_rep),
             C_avgscore   = cc[["C"]], C_avg_se = cc[["se"]],
             C_avg_lo     = cc[["lo"]], C_avg_hi = cc[["hi"]])
}))
save_tsv(ci_tbl, "05_cv_cindex_summary.tsv")
print(ci_tbl)

# Across-repeat spread of the increment describes CV noise. The paired patient
# bootstrap on the repeat-averaged predictors gives the sampling interval.
msg("Paired bootstrap of the C-index increments (", BOOT_B, " resamples) ...")
delta_tbl <- rbindlist(lapply(DELTAS, function(z) {
  d  <- per_rep[[delta_name(z)]]
  bt <- paired_boot_delta_c(y, lp_avg[, z[1]], lp_avg[, z[2]])
  data.table(comparison = paste(z[1], "-", z[2]),
             delta_mean = mean(d),
             range_lo = min(d), range_hi = max(d),
             pct2.5 = unname(quantile(d, 0.025)), pct97.5 = unname(quantile(d, 0.975)),
             n_repeats_positive = sum(d > 0), n_repeats = length(d),
             boot_delta = bt[["delta"]], boot_lo = bt[["lo"]], boot_hi = bt[["hi"]],
             boot_p = bt[["p_boot"]], boot_prop_positive = bt[["prop_positive"]],
             n_boot = bt[["n_boot"]], n = n, events = sum(cl_event))
}))
save_tsv(delta_tbl, "05_delta_cindex.tsv")
banner("Incremental value of module information over each comparator")
print(delta_tbl)

# ---- time-dependent AUC ----
times <- EVAL_TIMES_YRS * 365.25
auc <- rbindlist(lapply(ARMS, function(m) {
  tr <- timeROC(T = cl$os_time, delta = cl$os_event, marker = lp_avg[, m],
                cause = 1, times = times, iid = TRUE)
  data.table(model = m, years = EVAL_TIMES_YRS, n = n, events = sum(cl_event),
             AUC = as.numeric(tr$AUC),
             se  = as.numeric(tr$inference$vect_sd_1))
}))
auc[, `:=`(lo = AUC - 1.96 * se, hi = AUC + 1.96 * se)]
auc[, model := factor(model, levels = ARMS)]
save_tsv(auc, "05_time_dependent_auc.tsv")
print(auc)

g <- ggplot(auc, aes(factor(years), AUC, fill = model)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.8),
                width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(title = "Cross-validated time-dependent AUC",
       subtitle = "Repeat-averaged out-of-fold risk scores, TCGA-KIRC (internal only)",
       x = "Years from diagnosis", y = "AUC (95% CI)") +
  theme_bw()
save_fig(g, "05_time_dependent_auc", 7.5, 4)

# ---- decision curves ----
# Absolute risk from a Breslow baseline with the out-of-fold predictor at
# slope 1. The baseline is fitted on the same patients, so the curves are
# internal.
DCA_ARMS <- c(clinical     = "Clinical (age, sex, TNM, grade)",
              clinical_eig = "Clinical + module scores")
thr <- seq(0.02, 0.60, by = 0.01)
dca <- rbindlist(lapply(c(PRIMARY_HORIZON_YR, SECONDARY_HORIZON_YR), function(yr) {
  t_eval <- yr * 365.25
  rbind(
    rbindlist(lapply(names(DCA_ARMS), function(m) {
      rf <- baseline_risk_fun(y, lp_avg[, m])
      data.table(years = yr, threshold = thr, strategy = DCA_ARMS[[m]],
                 net_benefit = net_benefit(y, rf(lp_avg[, m], t_eval), t_eval, thr))
    })),
    data.table(years = yr, threshold = thr, strategy = "Treat all",
               net_benefit = net_benefit_treat_all(y, t_eval, thr)),
    data.table(years = yr, threshold = thr, strategy = "Treat none", net_benefit = 0))
}))
dca[, `:=`(n = n, events = sum(cl_event))]
save_tsv(dca, "05_decision_curve.tsv")

g <- ggplot(dca, aes(threshold, net_benefit, colour = strategy, linetype = strategy)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ years, labeller = labeller(years = function(x) paste0(x, "-year overall survival"))) +
  coord_cartesian(ylim = c(-0.05, max(dca$net_benefit, na.rm = TRUE) * 1.1)) +
  labs(title = "Decision curve analysis",
       subtitle = "Cross-validated risk scores, TCGA-KIRC",
       x = "Threshold probability of death by the horizon",
       y = "Net benefit", colour = NULL, linetype = NULL) +
  theme_bw() + theme(legend.position = "bottom")
save_fig(g, "05_decision_curve", 9, 4.5)

# ---- final models on the full analysis set (apparent) ----
# Comparator (unpenalised) plus all module scores (penalised), per comparator.
final <- lapply(c(clinical = "clinical", augmented = "augmented"), function(cmp) {
  Xcomp <- if (cmp == "clinical") Xc else Xa
  X <- cbind(Xcomp, Xe)
  set.seed(SEED)
  gf <- fit_cv(X, y, pf_for(Xcomp, Xe))
  co <- as.matrix(coef(gf, s = "lambda.min"))
  tab <- data.table(comparator = cmp, feature = rownames(co), coefficient = co[, 1])[coefficient != 0]
  tab[, HR := exp(coefficient)]
  tab[, type := ifelse(feature %in% colnames(Xcomp), "comparator", "module")]
  tab[, abs_coef := abs(coefficient)]      # setorder() takes columns, not calls
  setorder(tab, -abs_coef)
  tab[, abs_coef := NULL]
  tab[, `:=`(lambda_min = gf$lambda.min, n = n, events = sum(cl_event))]
  b_mod <- co[colnames(Xe), 1]
  list(table = tab, lp_modules = as.numeric(Xe %*% b_mod), Xcomp = Xcomp)
})
coef_tbl <- rbindlist(lapply(final, `[[`, "table"))
save_tsv(coef_tbl, "05_final_model_coefficients.tsv")
for (cmp in names(final)) {
  msg("Final ", cmp, " + modules model retains ",
      sum(final[[cmp]]$table$type == "module"), " module scores")
}
print(coef_tbl)

# ---- likelihood-ratio test for incremental value (apparent) ----
# The module part of each final linear predictor enters the comparator model
# as one term (1 df). Fitted and tested on the same data, so optimistic.
lrt_tbl <- rbindlist(lapply(names(final), function(cmp) {
  Xcomp <- final[[cmp]]$Xcomp
  m0 <- coxph(y ~ ., data = as.data.frame(Xcomp))
  m1 <- coxph(y ~ ., data = data.frame(Xcomp, module_score = final[[cmp]]$lp_modules))
  a  <- as.data.frame(anova(m0, m1))
  data.table(comparator = cmp,
             model  = c(cmp, paste(cmp, "+ module score")),
             n = n, events = sum(cl_event),
             loglik = a$loglik, Chisq = a$Chisq, Df = a$Df, p = a[["Pr(>|Chi|)"]])
}))
save_tsv(lrt_tbl, "05_likelihood_ratio_test.tsv")
print(lrt_tbl)

banner("CAVEAT")
cat(
"All performance estimates above are internal to TCGA-KIRC. Cross-validation\n",
"controls optimism from hub selection and penalty tuning but not from the\n",
"module definitions themselves (see 12) and cannot detect cohort-specific\n",
"artefacts. External validation of the locked model is in 09.\n", sep = "")

write_session_info("05_ml")
banner("05 | done")
