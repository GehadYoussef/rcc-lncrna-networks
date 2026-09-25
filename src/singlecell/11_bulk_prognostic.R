# 11_bulk_prognostic.R: prognostic development in TCGA-KIRC, validation in CPTAC-3.
# Run after 10. The only single-cell stage that reads survival.
# A. TCGA repeated nested CV, all preprocessing inside training folds. Arms:
#    clinical (CLIN_TERMS), clinical_sc_score (+ unweighted tier-1 score,
#    primary) and clinical_enet (elastic net over retained candidates).
# B. Final TCGA fits locked to bulk_model_lock.json before CPTAC-3 survival
#    is read. C. CPTAC-3 validation with locked coefficients: discrimination,
#    calibration, IPCW Brier and decision curves. D. Score coefficient
#    re-estimated in CPTAC-3 (association replication).
# STAR residualisation and z-scores are fitted within each cohort. Helpers come
# from src/R/00_config.R. Outputs: results/singlecell/11_*.tsv,
# figures/singlecell/11_*, reports/singlecell_bulk_validation.md

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
suppressPackageStartupMessages({ library(survival); library(glmnet); library(ggplot2) })
banner("11  prognostic development (TCGA-KIRC) and validation (CPTAC-3)")
start_log("11_bulk_prognostic")

BH <- CFG$bulk_handoff; PG <- CFG$prognostic
pth <- function(x) normalizePath(file.path(SUB_ROOT, x), winslash = "/", mustWork = TRUE)
H <- new.env(); H$R_DIR <- dirname(bulk_path(BH$helpers))
sys.source(bulk_path(BH$helpers), envir = H)
SEED <- H$SEED
thr <- seq(PG$nb_thresholds[[1]], PG$nb_thresholds[[2]], by = PG$nb_thresholds[[3]])
HZ <- unlist(PG$horizons_years); HZ_PRIMARY <- PG$primary_horizon_years

# ---- locks and candidate expression ---------------------------------------------------
bl <- jsonlite::read_json(file.path(RESULTS_DIR, "bulk_candidate_lock.json"))
if (!identical(sha256_file(file.path(RESULTS_DIR, "bulk_retained_candidates.tsv")), bl$retained_sha256))
  stop("bulk_retained_candidates.tsv does not match bulk_candidate_lock.json")
kept <- read_tsv("bulk_retained_candidates.tsv")
BX <- readRDS(file.path(DERIVED_DIR, "bulk_candidate_expression.rds"))
tier1 <- kept[tier == "tier1_tumour_specific", gene_key]
all_keys <- kept$gene_key
if (length(tier1) < BH$min_retained_tier1) stop("fewer than ", BH$min_retained_tier1, " retained tier-1 candidates")
msg("bulk lock ", substr(bl$bulk_lock_id, 1, 16), ": ", length(tier1), " tier-1, ", length(all_keys), " retained candidates")

inner_foldid <- function(n, seed, k = PG$inner_folds) {   # as src/R/30_endpoint_sensitivity.R
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed); f <- sample(rep(seq_len(k), length.out = n))
  if (is.null(old)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old, envir = .GlobalEnv)
  f
}
lp_lin <- function(b, X) { b[is.na(b)] <- 0; as.numeric(X[, names(b), drop = FALSE] %*% b) }
fit_enet <- function(X, y, pf, fid) cv.glmnet(X, y, family = "cox", alpha = H$ML_ALPHA, foldid = fid, nlambda = H$ML_NLAMBDA,
                                              lambda.min.ratio = H$ML_LAMBDA_MIN_RATIO, thresh = H$ML_THRESH,
                                              maxit = H$ML_MAXIT, penalty.factor = pf)
enet_coef <- function(fit) { b <- as.matrix(coef(fit, s = "lambda.min"))[, 1]; b[b != 0] }

# prepare a cohort's feature blocks given STAR and z-score fits
features <- function(E, S, tfit, zfit1, zfit_all) {
  R <- H$apply_technical(E, tfit, S)
  list(score = apply_zscore_mean(R[, zfit1$genes, drop = FALSE], zfit1),
       Z = sweep(sweep(R[, zfit_all$genes, drop = FALSE], 2, zfit_all$center, "-"), 2, zfit_all$scale, "/"))
}
fit_all_arms <- function(Xc, E, S, y, fid_seed) {
  tfit <- H$fit_technical(E, S)
  R <- H$apply_technical(E, tfit, S)
  z1 <- fit_zscore(R[, tier1, drop = FALSE]); za <- fit_zscore(R[, all_keys, drop = FALSE])
  f <- features(E, S, tfit, z1, za)
  b_clin <- coef(coxph(y ~ ., data = as.data.frame(Xc)))
  b_sc <- coef(coxph(y ~ ., data = as.data.frame(cbind(Xc, sc_score = f$score))))
  Xe <- cbind(Xc, f$Z); pf <- c(rep(0, ncol(Xc)), rep(1, ncol(f$Z)))
  en <- fit_enet(Xe, y, pf, inner_foldid(nrow(Xe), fid_seed))
  list(tfit = tfit, z1 = z1, za = za, b_clin = b_clin, b_sc = b_sc, b_enet = enet_coef(en), lambda = en$lambda.min,
       train_features = f)
}
predict_arms <- function(m, Xc, E, S, tfit = m$tfit, z1 = m$z1, za = m$za) {
  f <- features(E, S, tfit, z1, za)
  c(list(clinical = lp_lin(m$b_clin, Xc),
         clinical_sc_score = lp_lin(m$b_sc, cbind(Xc, sc_score = f$score)),
         clinical_enet = lp_lin(m$b_enet, cbind(Xc, f$Z))), list(score = f$score))
}
ARMS <- c("clinical", "clinical_sc_score", "clinical_enet")

# ---- A. TCGA-KIRC development ----
d <- readRDS(bulk_path(BH$tcga_dataset))
tcl <- as.data.table(d$cohort_adj); rm(d)
tx <- BX[["TCGA-KIRC"]]
tcl <- tcl[match(tx$samples, sample_barcode)]
Xc_t <- H$clinical_design(tcl)
ok_t <- complete.cases(Xc_t) & is.finite(tcl$os_time) & tcl$os_time > 0 & !is.na(tcl$os_event)
tcl <- tcl[ok_t]; Xc_t <- Xc_t[ok_t, , drop = FALSE]
E_t <- tx$log2fpkm[ok_t, , drop = FALSE]; S_t <- tx$star[ok_t, , drop = FALSE]
y_t <- Surv(tcl$os_time, tcl$os_event)
msg("TCGA-KIRC development: n=", nrow(tcl), ", deaths=", sum(tcl$os_event))

cv_rows <- list(); oof <- list()
for (r in seq_len(PG$cv_repeats)) {
  set.seed(SEED + r)
  folds <- sample(rep(seq_len(PG$cv_folds), length.out = nrow(tcl)))
  lp <- matrix(NA_real_, nrow(tcl), length(ARMS), dimnames = list(NULL, ARMS))
  for (k in seq_len(PG$cv_folds)) {
    tr <- folds != k; te <- !tr
    if (sum(tcl$os_event[tr]) < H$MIN_EPV) next
    m <- fit_all_arms(Xc_t[tr, , drop = FALSE], E_t[tr, , drop = FALSE], S_t[tr, , drop = FALSE], y_t[tr],
                      SEED + 1000L * r + k)
    p <- predict_arms(m, Xc_t[te, , drop = FALSE], E_t[te, , drop = FALSE], S_t[te, , drop = FALSE])
    for (a in ARMS) lp[te, a] <- p[[a]]
  }
  oof[[r]] <- lp
  cv_rows[[r]] <- data.table(repeat_id = r, arm = ARMS, C = vapply(ARMS, function(a) H$cindex(y_t, lp[, a]), numeric(1)))
  msg("  CV repeat ", r, ": ", paste(sprintf("%s %.3f", ARMS, cv_rows[[r]]$C), collapse = "; "))
}
CVR <- rbindlist(cv_rows)
CVR[, delta_vs_clinical := C - C[arm == "clinical"], by = repeat_id]
save_tsv(CVR, "11_tcga_cv_repeats.tsv")
avg_lp <- Reduce(`+`, oof) / length(oof)
cv_sum <- rbindlist(lapply(ARMS, function(a) {
  x <- CVR[arm == a]
  ci <- H$cindex_ci(y_t, avg_lp[, a])
  bd <- if (a == "clinical") rep(NA_real_, 6) else H$paired_boot_delta_c(y_t, avg_lp[, a], avg_lp[, "clinical"], B = H$BOOT_B, seed = SEED)
  data.table(cohort = "TCGA-KIRC", analysis = "repeated_nested_cv", arm = a, n = nrow(tcl), events = sum(tcl$os_event),
             C_median_over_repeats = median(x$C), C_p2.5 = quantile(x$C, 0.025), C_p97.5 = quantile(x$C, 0.975),
             deltaC_median_over_repeats = median(x$delta_vs_clinical),
             deltaC_p2.5 = quantile(x$delta_vs_clinical, 0.025), deltaC_p97.5 = quantile(x$delta_vs_clinical, 0.975),
             C_repeat_averaged_lp = ci[["C"]], deltaC_repeat_averaged = bd[1], deltaC_boot_lo = bd[2],
             deltaC_boot_hi = bd[3], deltaC_p_boot = bd[4])
}))
save_tsv(cv_sum, "11_tcga_cv_summary.tsv")
print(cv_sum[, .(arm, C_median_over_repeats, deltaC_median_over_repeats, deltaC_p2.5, deltaC_p97.5)])

# ---- B. final TCGA fits and model lock (before any CPTAC-3 outcome is read) ----
final <- fit_all_arms(Xc_t, E_t, S_t, y_t, SEED)
lp_t <- predict_arms(final, Xc_t, E_t, S_t)
risk_fun <- lapply(setNames(ARMS, ARMS), function(a) H$baseline_risk_fun(y_t, lp_t[[a]]))
sc_fit_t <- coxph(y_t ~ ., data = as.data.frame(cbind(Xc_t, sc_score = lp_t$score)))
zph <- cox.zph(sc_fit_t)
coefs <- rbindlist(lapply(c("clinical" = "b_clin", "clinical_sc_score" = "b_sc", "clinical_enet" = "b_enet"), function(nm)
  data.table(term = names(final[[nm]]), coefficient = unname(final[[nm]]))), idcol = "arm")
coefs <- merge(coefs, kept[, .(term = gene_key, symbol)], by = "term", all.x = TRUE)[order(arm, term)]
save_tsv(coefs, "11_locked_model_coefficients.tsv")
save_tsv(data.table(term = rownames(zph$table), chisq = zph$table[, "chisq"], df = zph$table[, "df"], p = zph$table[, "p"]),
         "11_tcga_proportional_hazards.tsv")
tcga_hr <- summary(sc_fit_t)$conf.int["sc_score", ]
msg("TCGA apparent: HR per unit sc_score ", sprintf("%.2f (%.2f-%.2f)", tcga_hr[1], tcga_hr[3], tcga_hr[4]),
    "; enet kept ", sum(names(final$b_enet) %in% all_keys), " candidate lncRNAs")

model <- list(bulk_lock_id = bl$bulk_lock_id, tier1 = tier1, all_keys = all_keys, clin_terms = colnames(Xc_t),
              b_clin = final$b_clin, b_sc = final$b_sc, b_enet = final$b_enet, lambda_min = final$lambda,
              tcga_z_tier1 = final$z1, tcga_z_all = final$za, tcga_technical_fit = final$tfit,
              tcga_baseline = list(y = y_t, lp = lp_t[ARMS]), sc_score_sd_tcga = sd(lp_t$score),
              cptac_rules = list(technical = "fit_technical within CPTAC-3 (own STAR metrics)", zscore = "within CPTAC-3",
                                 clinical_design = "within-cohort age z (primary)"),
              seed = SEED, cv = list(repeats = PG$cv_repeats, folds = PG$cv_folds, inner = PG$inner_folds))
model_id <- digest::digest(paste(bl$bulk_lock_id, table_digest(coefs), signif(final$lambda, 10), sep = "|"), algo = "sha256", serialize = FALSE)
mj <- file.path(RESULTS_DIR, "bulk_model_lock.json"); mr <- file.path(DERIVED_DIR, "bulk_locked_model.rds")
if (file.exists(mj)) {
  old <- jsonlite::read_json(mj)
  if (!identical(old$model_lock_id, model_id)) {
    if (!"--new-lock" %in% commandArgs(TRUE)) stop("a different locked bulk model exists; rerun with --new-lock to replace it")
    arch <- file.path(RESULTS_DIR, "lock_archive", paste0("model_", old$model_lock_id)); dir.create(arch, recursive = TRUE, showWarnings = FALSE)
    Sys.chmod(mj, "0644"); file.copy(mj, arch); unlink(mj)
  } else msg("locked bulk model reproduced exactly (", substr(model_id, 1, 16), ")")
}
if (!file.exists(mj)) {
  saveRDS(model, mr)
  jsonlite::write_json(list(model_lock_id = model_id, parent_bulk_lock_id = bl$bulk_lock_id,
                            parent_single_cell_lock_id = bl$parent_single_cell_lock_id,
                            created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                            written_before_validation_outcomes_read = TRUE,
                            development = list(cohort = "TCGA-KIRC", n = nrow(tcl), events = sum(tcl$os_event)),
                            primary_arm = "clinical_sc_score", secondary_arm = "clinical_enet",
                            coefficients = coefs, lambda_min = final$lambda, model_rds_sha256 = sha256_file(mr),
                            rules = PG, helpers_sha256 = sha256_file(bulk_path(BH$helpers))),
                       mj, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
  Sys.chmod(mj, "0444")
  msg("bulk model lock written: ", substr(model_id, 1, 16))
}

# ---- C. CPTAC-3 external validation with the locked model ----
v <- readRDS(bulk_path(BH$cptac_dataset))
vcl <- as.data.table(v$cohort); rm(v)
vx <- BX[["CPTAC-3"]]
vcl <- vcl[match(vx$samples, sample_barcode)]
Xc_v <- H$clinical_design(vcl)
ok_v <- complete.cases(Xc_v) & is.finite(vcl$os_time) & vcl$os_time > 0 & !is.na(vcl$os_event)
vcl <- vcl[ok_v]; Xc_v <- Xc_v[ok_v, , drop = FALSE]
E_v <- vx$log2fpkm[ok_v, , drop = FALSE]; S_v <- vx$star[ok_v, , drop = FALSE]
y_v <- Surv(vcl$os_time, vcl$os_event)
msg("CPTAC-3 validation: n=", nrow(vcl), ", deaths=", sum(vcl$os_event), ", max follow-up ", round(max(vcl$os_time) / 365.25, 1), " y")

tfit_v <- H$fit_technical(E_v, S_v)
Rv <- H$apply_technical(E_v, tfit_v, S_v)
z1_v <- fit_zscore(Rv[, tier1, drop = FALSE]); za_v <- fit_zscore(Rv[, all_keys, drop = FALSE])
missing_enet <- setdiff(intersect(names(final$b_enet), all_keys), za_v$genes)
if (length(missing_enet)) stop("elastic-net genes without variance in CPTAC-3: ", paste(missing_enet, collapse = ", "))
lp_v <- predict_arms(final, Xc_v, E_v, S_v, tfit = tfit_v, z1 = z1_v, za = za_v)

disc <- rbindlist(lapply(ARMS, function(a) {
  ci <- H$cindex_ci(y_v, lp_v[[a]])
  bd <- if (a == "clinical") rep(NA_real_, 6) else H$paired_boot_delta_c(y_v, lp_v[[a]], lp_v$clinical, B = H$BOOT_B, seed = SEED)
  data.table(cohort = "CPTAC-3", analysis = "locked_model_external_validation", arm = a, n = nrow(vcl), events = sum(vcl$os_event),
             C = ci[["C"]], C_lo = ci[["lo"]], C_hi = ci[["hi"]], deltaC_vs_clinical = bd[1], deltaC_lo = bd[2],
             deltaC_hi = bd[3], deltaC_p_boot = bd[4], deltaC_prop_positive = bd[5])
}))

hz_rows <- list(); bins <- list(); nb <- list(); auc_rows <- list()
set.seed(SEED)
boot_idx <- lapply(seq_len(H$BOOT_B), function(i) sample.int(nrow(vcl), nrow(vcl), replace = TRUE))
for (yr in HZ) {
  t_d <- yr * 365.25
  km <- summary(survfit(y_v ~ 1), times = t_d, extend = TRUE)
  at_risk <- km$n.risk
  preds <- lapply(setNames(ARMS, ARMS), function(a) vapply(lp_v[[a]], function(l) risk_fun[[a]](l, t_d), numeric(1)))
  for (a in ARMS) {
    ta <- tryCatch(timeROC::timeROC(T = vcl$os_time, delta = vcl$os_event, marker = lp_v[[a]], cause = 1, times = t_d, iid = TRUE),
                   error = function(e) NULL)
    cs <- H$calibration_summary(y_v, lp_v[[a]], preds[[a]], t_d)
    br <- H$brier_ipcw(y_v, preds[[a]], t_d)
    dbr <- if (a == "clinical") c(NA, NA, NA) else {
      b <- vapply(boot_idx, function(i) H$brier_ipcw(y_v[i], preds[[a]][i], t_d) - H$brier_ipcw(y_v[i], preds$clinical[i], t_d), numeric(1))
      c(br - H$brier_ipcw(y_v, preds$clinical, t_d), quantile(b, c(0.025, 0.975), na.rm = TRUE))
    }
    hz_rows[[paste(yr, a)]] <- data.table(horizon_years = yr, primary = yr == HZ_PRIMARY, arm = a, n_at_risk = at_risk,
      low_at_risk_flag = at_risk < PG$min_at_risk_report,
      AUC = if (!is.null(ta)) unname(ta$AUC[2]) else NA_real_,
      AUC_se = if (!is.null(ta)) unname(ta$inference$vect_sd_1[2]) else NA_real_,
      calibration_slope = cs[["slope"]], slope_lo = cs[["slope_lo"]], slope_hi = cs[["slope_hi"]],
      observed_risk = cs[["observed"]], expected_risk = cs[["expected"]], OE = cs[["OE"]],
      brier = br, delta_brier_vs_clinical = dbr[1], delta_brier_lo = dbr[2], delta_brier_hi = dbr[3])
    if (yr == HZ_PRIMARY) bins[[a]] <- H$calibration_bins(y_v, preds[[a]], t_d)[, arm := a]
  }
  if (yr == HZ_PRIMARY) {
    nb_all <- H$net_benefit_treat_all(y_v, t_d, thr)
    nbs <- lapply(setNames(ARMS, ARMS), function(a) H$net_benefit(y_v, preds[[a]], t_d, thr))
    nb <- rbind(rbindlist(lapply(ARMS, function(a) data.table(arm = a, threshold = thr, net_benefit = nbs[[a]]))),
                data.table(arm = "treat_all", threshold = thr, net_benefit = nb_all),
                data.table(arm = "treat_none", threshold = thr, net_benefit = 0))
    win <- thr >= PG$nb_summary_window[[1]] & thr <= PG$nb_summary_window[[2]]
    nb_sum <- rbindlist(lapply(setdiff(ARMS, "clinical"), function(a) {
      obs <- median(nbs[[a]][win] - nbs$clinical[win])
      b <- vapply(boot_idx[seq_len(min(500, H$BOOT_B))], function(i) median(H$net_benefit(y_v[i], preds[[a]][i], t_d, thr[win]) -
                                                                                H$net_benefit(y_v[i], preds$clinical[i], t_d, thr[win])), numeric(1))
      data.table(horizon_years = yr, arm = a, window = paste(PG$nb_summary_window, collapse = "-"),
                 median_delta_net_benefit = obs, boot_lo = quantile(b, 0.025), boot_hi = quantile(b, 0.975), n_boot = length(b))
    }))
  }
}
HM <- rbindlist(hz_rows)
save_tsv(disc, "11_cptac_discrimination.tsv")
save_tsv(HM, "11_cptac_horizon_metrics.tsv")
save_tsv(rbindlist(bins), "11_cptac_calibration_bins.tsv")
save_tsv(nb, "11_cptac_net_benefit.tsv")
save_tsv(nb_sum, "11_cptac_net_benefit_summary.tsv")
print(disc[, .(arm, C, C_lo, C_hi, deltaC_vs_clinical, deltaC_lo, deltaC_hi)])
print(HM[primary == TRUE, .(arm, n_at_risk, AUC, calibration_slope, OE, brier, delta_brier_vs_clinical, delta_brier_lo, delta_brier_hi)])

# ---- D. association replication (coefficient re-estimated in CPTAC-3) ----
fit_v <- coxph(y_v ~ ., data = as.data.frame(cbind(Xc_v, sc_score = lp_v$score)))
epv <- sum(vcl$os_event) / (ncol(Xc_v) + 1)
hr_v <- summary(fit_v)$conf.int["sc_score", ]; p_v <- summary(fit_v)$coefficients["sc_score", "Pr(>|z|)"]
uni_t <- summary(coxph(y_t ~ lp_t$score))$conf.int[1, ]; uni_v <- summary(coxph(y_v ~ lp_v$score))$conf.int[1, ]
assoc <- data.table(
  cohort = c("TCGA-KIRC", "CPTAC-3", "TCGA-KIRC", "CPTAC-3"),
  model = c("clinical + sc_score", "clinical + sc_score", "sc_score alone", "sc_score alone"),
  analysis = c("development_apparent", "association_replication", "development_apparent", "association_replication"),
  n = c(nrow(tcl), nrow(vcl), nrow(tcl), nrow(vcl)), events = c(sum(tcl$os_event), sum(vcl$os_event), sum(tcl$os_event), sum(vcl$os_event)),
  HR_per_unit_score = c(tcga_hr[1], hr_v[1], uni_t[1], uni_v[1]), HR_lo = c(tcga_hr[3], hr_v[3], uni_t[3], uni_v[3]),
  HR_hi = c(tcga_hr[4], hr_v[4], uni_t[4], uni_v[4]),
  p = c(summary(sc_fit_t)$coefficients["sc_score", "Pr(>|z|)"], p_v,
        summary(coxph(y_t ~ lp_t$score))$coefficients[1, "Pr(>|z|)"], summary(coxph(y_v ~ lp_v$score))$coefficients[1, "Pr(>|z|)"]),
  events_per_parameter = c(sum(tcl$os_event) / (ncol(Xc_t) + 1), epv, sum(tcl$os_event), sum(vcl$os_event)),
  epv_below_min = c(sum(tcl$os_event) / (ncol(Xc_t) + 1), epv, sum(tcl$os_event), sum(vcl$os_event)) < H$MIN_EPV)
save_tsv(assoc, "11_association_replication.tsv")
print(assoc)

# ---- figures --------------------------------------------------------------------------
B3 <- rbindlist(bins)
p1 <- ggplot(B3, aes(x = predicted, y = observed, colour = arm)) + geom_abline(linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = obs_lo, ymax = obs_hi), position = position_dodge(width = 0.01)) + geom_line() +
  coord_equal(xlim = c(0, max(c(B3$predicted, B3$obs_hi), na.rm = TRUE)), ylim = c(0, max(c(B3$predicted, B3$obs_hi), na.rm = TRUE))) +
  labs(x = sprintf("Predicted %d-year risk (TCGA-locked model)", HZ_PRIMARY), y = "Observed risk (Kaplan-Meier, 95% CI)",
       title = "CPTAC-3 calibration of locked models", subtitle = "Quintiles of predicted risk") + theme_bw(base_size = 9)
save_fig(p1, "11_cptac_calibration", width = 6, height = 5.5)
p2 <- ggplot(nb, aes(x = threshold, y = net_benefit, colour = arm)) + geom_line() +
  coord_cartesian(ylim = c(-0.05, max(nb$net_benefit, na.rm = TRUE) + 0.02)) +
  labs(x = "Threshold probability", y = "Net benefit", title = sprintf("CPTAC-3 decision curves at %d years", HZ_PRIMARY)) +
  theme_bw(base_size = 9)
save_fig(p2, "11_cptac_decision_curves", width = 6.5, height = 4.5)

# ---- report -----------------------------------------------------------------------------
fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")
L <- c("# Single-cell lncRNA candidates in bulk RNA-seq: development and validation\n",
  sprintf("Single-cell lock `%s` -> bulk candidate lock `%s` -> model lock `%s`.\n",
          substr(bl$parent_single_cell_lock_id, 1, 16), substr(bl$bulk_lock_id, 1, 16), substr(model_id, 1, 16)),
  "Single-cell discovery validates cellular localisation only; the survival results below are an independent evaluation.\n",
  sprintf("## Candidates\n\n%d locked; %d retained after pre-specified bulk availability and library-quality rules (%d tier-1 used in the primary score).\n",
          bl$n_locked, bl$n_retained, length(tier1)),
  "## TCGA-KIRC repeated nested cross-validation (development)\n",
  paste0("| arm | median C | median dC vs clinical | dC 2.5-97.5 percentile over repeats |\n|---|---|---|---|\n",
         paste(cv_sum[, sprintf("| %s | %s | %s | %s to %s |", arm, fmt(C_median_over_repeats), fmt(deltaC_median_over_repeats, 4),
                                fmt(deltaC_p2.5, 4), fmt(deltaC_p97.5, 4))], collapse = "\n"), "\n"),
  "## CPTAC-3 external validation of the locked models (predictive validation)\n",
  paste0("| arm | C (95% CI) | dC vs clinical (bootstrap 95% CI) |\n|---|---|---|\n",
         paste(disc[, sprintf("| %s | %s (%s-%s) | %s |", arm, fmt(C), fmt(C_lo), fmt(C_hi),
                              ifelse(is.na(deltaC_vs_clinical), "-", sprintf("%s (%s to %s)", fmt(deltaC_vs_clinical, 4), fmt(deltaC_lo, 4), fmt(deltaC_hi, 4))))], collapse = "\n"), "\n"),
  sprintf("\n%d-year (primary) metrics:\n", HZ_PRIMARY),
  paste0("| arm | at risk | AUC | calibration slope | O/E | Brier | dBrier vs clinical (95% CI) |\n|---|---|---|---|---|---|---|\n",
         paste(HM[primary == TRUE, sprintf("| %s | %d | %s | %s | %s | %s | %s |", arm, as.integer(n_at_risk), fmt(AUC), fmt(calibration_slope, 2),
                                           fmt(OE, 2), fmt(brier, 4), ifelse(is.na(delta_brier_vs_clinical), "-",
                                           sprintf("%s (%s to %s)", fmt(delta_brier_vs_clinical, 4), fmt(delta_brier_lo, 4), fmt(delta_brier_hi, 4))))], collapse = "\n"), "\n"),
  "\nNet benefit (median difference vs clinical over the threshold window, bootstrap CI):\n",
  paste(nb_sum[, sprintf("- %s: %s (%s to %s)", arm, fmt(median_delta_net_benefit, 4), fmt(boot_lo, 4), fmt(boot_hi, 4))], collapse = "\n"),
  "\n\n## Association replication (coefficient re-estimated in CPTAC-3; not predictive validation)\n",
  paste0("| cohort | model | HR per unit score (95% CI) | p | EPV |\n|---|---|---|---|---|\n",
         paste(assoc[, sprintf("| %s | %s | %s (%s-%s) | %s | %s%s |", cohort, model, fmt(HR_per_unit_score, 2), fmt(HR_lo, 2), fmt(HR_hi, 2),
                               formatC(p, digits = 2, format = "g"), fmt(events_per_parameter, 1), ifelse(epv_below_min, " (below MIN_EPV)", ""))], collapse = "\n"), "\n"),
  sprintf("\nHorizon metrics with fewer than %d patients at risk are flagged in `11_cptac_horizon_metrics.tsv`.\n", PG$min_at_risk_report))
writeLines(L, file.path(REPORT_DIR, "singlecell_bulk_validation.md"))
msg("report written: reports/singlecell_bulk_validation.md")
write_session_info("11_bulk_prognostic")
