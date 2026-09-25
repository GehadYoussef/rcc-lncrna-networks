# 41_pancancer_metric_survival.R: is the non-feature fraction prognostic across cancer types?
# Per eligible tumour project, Cox models of TCGA-CDR overall survival
# (censored at OS_CENSOR_DAYS) on the non-feature fraction per project SD:
#   primary    metric + age + sex (+ ordinal stage) + strata(sequencing plate)
#   secondary  the same without the plate stratum
# Plates with fewer than PAN_PLATE_MIN_N patients are pooled into one stratum.
# Projects with >= PAN_SURV_MIN_EVENTS events enter a REML random-effects
# meta-analysis of the log hazard ratio, per model.
# Stage: AJCC pathologic stage collapsed to I-IV, else clinical stage (FIGO or
# AJCC), each used only if available for PAN_STAGE_MIN_FRAC of patients.
# Age is standardised within project. Sex is omitted for single-sex projects.
# The decision rules are checked against results/38_decision_rules_lock.json.
# Outputs: 41_pancancer_metric_survival.tsv (project x model),
#          41_pancancer_metric_meta.tsv (one row per model).

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
source(file.path(R_DIR, "pan_helpers.R"))
suppressPackageStartupMessages({ library(data.table); library(survival); library(readxl); library(metafor) })
banner("41 | Pan-cancer prognosis of the non-feature fraction")
verify_decision_lock()
TEST <- nzchar(Sys.getenv("PAN_TEST_ALL")); SUFFIX <- if (TEST) "_partial" else ""

cdr <- as.data.table(read_excel(file.path(CACHE_DIR, "30_TCGA-CDR.xlsx"), sheet = "TCGA-CDR",
                                guess_max = 20000))
cdr <- cdr[, .(patient = bcr_patient_barcode, type,
               age = suppressWarnings(as.numeric(age_at_initial_pathologic_diagnosis)),
               sex = tolower(gender), path_stage = ajcc_pathologic_tumor_stage,
               clin_stage = clinical_stage,
               os = suppressWarnings(as.numeric(OS)), os_time = suppressWarnings(as.numeric(OS.time)))]
roman <- function(x) {
  x <- toupper(trimws(as.character(x)))
  r <- sub("^STAGE\\s*", "", x); r <- sub("[A-C][0-9]?$", "", r)
  o <- rep(NA_integer_, length(x))
  o[r == "I"] <- 1L; o[r == "II"] <- 2L; o[r == "III"] <- 3L; o[r == "IV"] <- 4L
  o
}
cdr[, `:=`(path_stage_n = roman(path_stage), clin_stage_n = roman(clin_stage))]
cdr[!is.na(os_time) & os_time > OS_CENSOR_DAYS, os := 0]
cdr[!is.na(os_time), os_time := pmin(os_time, OS_CENSOR_DAYS)]

fit_one <- function(d, rhs, strat) {
  f <- as.formula(paste("Surv(os_time, os) ~", paste(c(rhs, if (strat) "strata(plate_s)"), collapse = " + ")))
  m <- tryCatch(coxph(f, data = d), error = function(e) NULL)
  if (is.null(m) || !"metric_z" %in% names(coef(m))) return(NULL)
  s <- summary(m)$coefficients["metric_z", ]
  data.table(log_hr = s[["coef"]], se = s[["se(coef)"]], hr = exp(s[["coef"]]),
             lo = exp(s[["coef"]] - 1.96 * s[["se(coef)"]]), hi = exp(s[["coef"]] + 1.96 * s[["se(coef)"]]),
             p = s[["Pr(>|z|)"]])
}

rows <- list()
for (proj in pan_projects()) {
  obj <- readRDS(file.path(CACHE_DIR, paste0("pan_", proj, ".rds")))
  s <- obj$samples[keep == TRUE & group == "tumour"]
  eligible <- nrow(s) >= PAN_MIN_TUMOURS && !proj %in% PAN_EXCLUDE_INFERENCE
  rm(obj)
  if (!eligible && !TEST) next
  d <- merge(s[, .(patient, pct_noFeature, plate)], cdr, by = "patient")
  d <- d[is.finite(os_time) & os_time > 0 & !is.na(os) & is.finite(age) & sex %in% c("female", "male")]
  stage_src <- if (mean(!is.na(d$path_stage_n)) >= PAN_STAGE_MIN_FRAC) "ajcc_pathologic" else
               if (mean(!is.na(d$clin_stage_n)) >= PAN_STAGE_MIN_FRAC) "clinical" else "none"
  d[, stage := switch(stage_src, ajcc_pathologic = path_stage_n, clinical = clin_stage_n, NA_integer_)]
  if (stage_src != "none") d <- d[!is.na(stage)]
  d[, `:=`(metric_z = as.numeric(scale(pct_noFeature)), age_z = as.numeric(scale(age)),
           male = as.integer(sex == "male"), plate_s = pool_levels(plate))]
  rhs <- c("metric_z", "age_z", if (uniqueN(d$male) > 1) "male", if (stage_src != "none") "stage")
  for (strat in c(TRUE, FALSE)) {
    r <- fit_one(d, rhs, strat)
    if (is.null(r)) next
    rows[[length(rows) + 1]] <- cbind(data.table(
      project = proj, eligible = eligible, model = if (strat) "primary (plate stratum)" else "secondary (no plate stratum)",
      n = nrow(d), events = sum(d$os), stage_source = stage_src, n_plate_strata = nlevels(d$plate_s),
      covariates = paste(rhs, collapse = " + ")), r)
  }
  msg(sprintf("%-10s n=%4d events=%3d stage=%s", proj, nrow(d), sum(d$os), stage_src))
}
sv <- rbindlist(rows)
sv[, in_meta := eligible & events >= PAN_SURV_MIN_EVENTS]
save_tsv(sv[order(model, project)], paste0("41_pancancer_metric_survival", SUFFIX, ".tsv"))

meta <- rbindlist(lapply(unique(sv$model), function(mdl) {
  x <- sv[model == mdl & in_meta == TRUE]
  if (nrow(x) < 3) return(data.table(model = mdl, k = nrow(x)))
  m <- rma(yi = x$log_hr, sei = x$se, method = "REML")
  pr <- predict(m, transf = exp)
  data.table(model = mdl, k = nrow(x), pooled_hr = exp(as.numeric(m$b)), lo = exp(m$ci.lb), hi = exp(m$ci.ub),
             p = m$pval, tau2 = m$tau2, i2 = m$I2, q_p = m$QEp, pi_lo = pr$pi.lb, pi_hi = pr$pi.ub,
             n_projects_hr_above_1 = sum(x$hr > 1), n_projects_p_below_0.05 = sum(x$p < 0.05))
}), fill = TRUE)
save_tsv(meta, paste0("41_pancancer_metric_meta", SUFFIX, ".tsv")); print(meta)
write_session_info("41_pancancer_metric_survival")
msg("41 done")
