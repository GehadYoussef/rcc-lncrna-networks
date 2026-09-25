# 33_sex_stratified.R: sex-disaggregated module effects and X/Y-encoded lncRNA module members (SAGER)
# For the five prognostic modules, a module x sex interaction model and separate fits in
# women and men in TCGA-KIRC (principal and clinical sets), CPTAC-3 (locked-model scores)
# and TCGA-KIRP. Also tests scores, STAR metrics and the lncRNA axis for association
# with sex in four cohorts, lists chrX/chrY lncRNA module members with Fisher enrichment
# against GENCODE v36, and runs two sensitivity analyses (X/Y members deleted, recorded
# sex checked against XIST and chrY lncRNAs).
# Inputs: caches dataset, networks, survival, locked_model, estimate_scores,
#   lnc_global_axis, subtype_*, validation_dataset, 28_gencode_v36.gtf.gz (if present),
#   and results of stages 02, 07, 09 and 11.
# Outputs: results/33_*.tsv.
# The sex-specific fits let the baseline hazard and every covariate differ by sex and
# carry their own Wald and LR difference tests. The common-covariate interaction model
# tests the module x sex term alone. The two meet only in 33_sex_contrast_reconciliation.tsv.
# Scores are per SD over the whole set. Quality metrics and the axis are standardised
# before complete-case restriction, as in stage 07.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(httr); library(jsonlite)
})
banner("33 | Sex-disaggregated module effects and X/Y-encoded lncRNA members")
set.seed(SEED)
t_start <- Sys.time()

ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
surv <- readRDS(file.path(CACHE_DIR, "survival.rds"))
L    <- readRDS(file.path(CACHE_DIR, "locked_model.rds"))
stopifnot(identical(L$version, "v9"))
EST  <- as.data.table(readRDS(file.path(CACHE_DIR, "estimate_scores.rds")))
AXIS <- as.data.table(readRDS(file.path(CACHE_DIR, "lnc_global_axis.rds")))
LOAD <- discovery_loadings(nets)
SC   <- discovery_scores(nets, LOAD)
cohort <- as.data.table(ds$cohort_full)
if (!"assigned_reads" %in% names(cohort)) cohort[, assigned_reads := libsize]
stopifnot(!anyNA(cohort$sex))
msg("TCGA-KIRC: ", nrow(cohort), " patients, ", sum(cohort$os_event), " deaths; ",
    sum(cohort$sex == "female"), " women (", sum(cohort$os_event[cohort$sex == "female"]),
    " deaths), ", sum(cohort$sex == "male"), " men (",
    sum(cohort$os_event[cohort$sex == "male"]), " deaths)")

# ---- prognostic modules (within-network FDR under the full model) ----
TARGET <- rbind(
  as.data.table(surv$mrna)[fdr_full < FDR_ALPHA, .(biotype = "mRNA",   module, HR_full)],
  as.data.table(surv$lnc)[fdr_full < FDR_ALPHA,  .(biotype = "lncRNA", module, HR_full)])
TARGET[, feature := paste0(ifelse(biotype == "mRNA", "mRNA_ME", "lnc_ME"), module)]
msg("Prognostic modules: ",
    paste(TARGET[, paste(biotype, module)], collapse = ", "))

NET_TAG <- c(mRNA = "mrna", lncRNA = "lnc")
Z975    <- qnorm(0.975)
per_sd  <- function(x) as.numeric((x - mean(x)) / sd(x))

# ---- Cox helpers ----
# Cox fit on complete cases. A covariate constant on the analysed set is
# dropped and named. With interaction = TRUE the exposure is crossed with the
# male indicator and `term` should be "ME:male".
cox_fit <- function(d, exposure, covs, term = exposure, interaction = FALSE) {
  vars <- c("os_time", "os_event", exposure, covs)
  dc   <- d[complete.cases(d[, vars, with = FALSE])]
  keep <- vapply(covs, function(v) length(unique(dc[[v]])) > 1, logical(1))
  dropped <- covs[!keep]; covs <- covs[keep]
  rhs <- if (interaction) c(paste0(exposure, " * male"), setdiff(covs, "male"))
         else c(exposure, covs)
  f  <- coxph(as.formula(paste("Surv(os_time, os_event) ~",
                               paste(rhs, collapse = " + "))), data = dc)
  s  <- summary(f); cf <- s$coefficients
  list(fit = f, data = dc, n = s$n, events = s$nevent, n_par = nrow(cf),
       HR = unname(s$conf.int[term, "exp(coef)"]),
       lo = unname(s$conf.int[term, "lower .95"]),
       hi = unname(s$conf.int[term, "upper .95"]),
       p  = unname(cf[term, "Pr(>|z|)"]),
       coef = unname(cf[term, "coef"]), se = unname(cf[term, "se(coef)"]),
       dropped = dropped)
}
implausible <- function(HR, se) !is.finite(HR) | HR > 20 | HR < 0.05 | !is.finite(se) | se > 3
wald_p <- function(b, se) 2 * pnorm(-abs(b / se))

# Model labels written into every row.
MODEL_STRATIFIED  <- "sex-specific fit: separate baseline hazard and separate covariate effects in women and in men"
MODEL_INTERACTION <- "common-covariate interaction model: one baseline hazard, covariate effects common to both sexes, module x sex term"

# Likelihood-ratio test for the sex-specific fits: a strata(male) model with
# sex-specific covariates, with and without a common exposure effect. The
# unconstrained fit reproduces the summed log-likelihood of the two
# sex-specific fits, and the gap is returned and checked.
strat_lrt <- function(d, covs, sex_fits) {
  vars <- c("os_time", "os_event", "ME", covs)
  dc   <- d[complete.cases(d[, vars, with = FALSE])]
  Z <- data.table(os_time = dc$os_time, os_event = dc$os_event, male = dc$male,
                  ME = as.numeric(dc$ME), ME_male = as.numeric(dc$ME) * dc$male)
  sex_terms <- character(0)
  for (v in setdiff(covs, "male")) for (s in c(0, 1)) {
    inside <- dc$male == s
    if (length(unique(dc[[v]][inside])) > 1) {
      nm <- paste0(v, "_sex", s)
      Z[, (nm) := fifelse(inside, as.numeric(dc[[v]]), 0)]
      sex_terms <- c(sex_terms, nm)
    }
  }
  mk <- function(extra) coxph(as.formula(paste(
    "Surv(os_time, os_event) ~",
    paste(c("strata(male)", "ME", extra, sex_terms), collapse = " + "))), data = Z)
  f0 <- mk(NULL); f1 <- mk("ME_male")
  ll_sep <- sum(vapply(sex_fits, function(f) f$fit$loglik[2], numeric(1)))
  stat   <- 2 * (f1$loglik[2] - f0$loglik[2])
  list(p = pchisq(stat, 1, lower.tail = FALSE), chisq = stat,
       loglik_gap_vs_sex_fits = abs(f1$loglik[2] - ll_sep))
}

# One exposure (column ME of d) under one covariate set. Returns three blocks:
# the interaction model (with the additive fit), the sex-specific fits, and
# the two side by side.
sex_models <- function(d, covs, cohort_name, spec, family, biotype, module,
                       is_transform = FALSE) {
  stopifnot("male" %in% covs, all(c("ME", covs) %in% names(d)))
  add <- cox_fit(d, "ME", covs)
  int <- cox_fit(d, "ME", covs, term = "ME:male", interaction = TRUE)
  V <- vcov(int$fit); b <- coef(int$fit)
  se_f <- sqrt(V["ME", "ME"])
  b_m  <- b[["ME"]] + b[["ME:male"]]
  se_m <- sqrt(V["ME", "ME"] + V["ME:male", "ME:male"] + 2 * V["ME", "ME:male"])
  lrt  <- anova(add$fit, int$fit)
  p_lrt <- lrt[2, ncol(lrt)]
  ev  <- int$data[, .(n = .N, events = sum(os_event)), by = male]
  gv  <- function(sx, col) { v <- ev[male == sx][[col]]; if (length(v)) v[1] else 0L }
  # ---- block 1: common-covariate interaction model ----
  # The female and male HRs implied by this model differ from those of the
  # sex-specific fits, hence the explicit column names.
  interaction <- data.table(
    cohort = cohort_name, spec = spec, family = family, biotype = biotype, module = module,
    model = MODEL_INTERACTION,
    is_transform = is_transform, covariates = paste(covs, collapse = ","),
    n = int$n, events = int$events,
    n_female = gv(0, "n"), events_female = gv(0, "events"),
    n_male   = gv(1, "n"), events_male   = gv(1, "events"),
    HR_additive = add$HR, lo_additive = add$lo, hi_additive = add$hi, p_additive = add$p,
    HR_female_interaction_model = exp(b[["ME"]]),
    lo_female_interaction_model = exp(b[["ME"]] - Z975 * se_f),
    hi_female_interaction_model = exp(b[["ME"]] + Z975 * se_f),
    p_female_interaction_model  = wald_p(b[["ME"]], se_f),
    HR_male_interaction_model = exp(b_m),
    lo_male_interaction_model = exp(b_m - Z975 * se_m),
    hi_male_interaction_model = exp(b_m + Z975 * se_m),
    p_male_interaction_model  = wald_p(b_m, se_m),
    HR_interaction = int$HR, lo_interaction = int$lo, hi_interaction = int$hi,
    p_interaction = int$p, p_interaction_LRT = p_lrt,
    n_parameters = int$n_par, epv = int$events / int$n_par,
    epv_below_min = int$events / int$n_par < MIN_EPV,
    covariates_dropped = paste(int$dropped, collapse = ";"))
  # ---- block 2: sex-specific fits ----
  sf <- lapply(c(female = "female", male = "male"), function(sx)
    cox_fit(d[male == as.numeric(sx == "male")], "ME", setdiff(covs, "male")))
  # Difference test for these fits: Wald contrast of two independent log HRs.
  d_b   <- sf$male$coef - sf$female$coef
  d_se  <- sqrt(sf$female$se^2 + sf$male$se^2)
  p_w   <- wald_p(d_b, d_se)
  ll    <- strat_lrt(d, covs, sf)
  stratified <- rbindlist(lapply(names(sf), function(sx) {
    f   <- sf[[sx]]
    sdv <- sd(f$data$ME)
    data.table(
      cohort = cohort_name, spec = spec, family = family, biotype = biotype, module = module,
      model = MODEL_STRATIFIED,
      is_transform = is_transform, sex = sx,
      covariates = paste(setdiff(covs, "male"), collapse = ","),
      n = f$n, events = f$events,
      HR = f$HR, lo = f$lo, hi = f$hi, p = f$p,
      sd_score_in_stratum = sdv,
      HR_per_stratum_SD = exp(f$coef * sdv),
      lo_per_stratum_SD = exp((f$coef - Z975 * f$se) * sdv),
      hi_per_stratum_SD = exp((f$coef + Z975 * f$se) * sdv),
      # pair-level contrast, repeated on both rows
      HR_ratio_male_over_female = exp(d_b),
      ratio_lo = exp(d_b - Z975 * d_se), ratio_hi = exp(d_b + Z975 * d_se),
      p_sex_difference = p_w, p_sex_difference_LRT = ll$p,
      n_parameters = f$n_par, epv = f$events / f$n_par,
      epv_below_min = f$events / f$n_par < MIN_EPV,
      implausible = implausible(f$HR, f$se),
      covariates_dropped = paste(f$dropped, collapse = ";"))
  }))
  # ---- block 3: both models side by side ----
  contrast <- data.table(
    cohort = cohort_name, spec = spec, family = family, biotype = biotype, module = module,
    is_transform = is_transform,
    n = int$n, events = int$events,
    n_female = gv(0, "n"), events_female = gv(0, "events"),
    n_male   = gv(1, "n"), events_male   = gv(1, "events"),
    strat_model = MODEL_STRATIFIED,
    strat_HR_female = sf$female$HR, strat_lo_female = sf$female$lo,
    strat_hi_female = sf$female$hi, strat_p_female = sf$female$p,
    strat_HR_male = sf$male$HR, strat_lo_male = sf$male$lo,
    strat_hi_male = sf$male$hi, strat_p_male = sf$male$p,
    strat_HR_ratio_male_over_female = exp(d_b),
    strat_ratio_lo = exp(d_b - Z975 * d_se), strat_ratio_hi = exp(d_b + Z975 * d_se),
    strat_p_sex_difference = p_w, strat_p_sex_difference_LRT = ll$p,
    strat_lrt_loglik_gap = ll$loglik_gap_vs_sex_fits,
    int_model = MODEL_INTERACTION,
    int_HR_female = exp(b[["ME"]]), int_lo_female = exp(b[["ME"]] - Z975 * se_f),
    int_hi_female = exp(b[["ME"]] + Z975 * se_f), int_p_female = wald_p(b[["ME"]], se_f),
    int_HR_male = exp(b_m), int_lo_male = exp(b_m - Z975 * se_m),
    int_hi_male = exp(b_m + Z975 * se_m), int_p_male = wald_p(b_m, se_m),
    int_HR_interaction = int$HR, int_lo_interaction = int$lo, int_hi_interaction = int$hi,
    int_p_interaction = int$p, int_p_interaction_LRT = p_lrt)
  list(interaction = interaction, stratified = stratified, contrast = contrast)
}

INT <- list(); STR <- list(); CON <- list()
collect <- function(r) {
  INT[[length(INT) + 1]] <<- r$interaction
  STR[[length(STR) + 1]] <<- r$stratified
  CON[[length(CON) + 1]] <<- r$contrast
  invisible(NULL)
}

# ---- 1. TCGA-KIRC: principal and clinical covariate sets ----
banner("1 | TCGA-KIRC: module x sex interaction and sex-stratified models")
# Covariate frame as in stage 03: STAR metrics and ESTIMATE scores are
# standardised over the samples given, before any complete-case restriction.
kirc_base <- function(samples) {
  cl <- cohort[match(samples, sample_barcode)]
  stopifnot(identical(cl$sample_barcode, samples))
  e  <- EST[match(cl$sample_barcode, sample_barcode)]
  data.table(
    sample_barcode = cl$sample_barcode,
    os_time  = cl$os_time, os_event = cl$os_event,
    age      = cl$age,
    male     = as.numeric(cl$sex == "male"),
    T_stage  = as.numeric(cl$T_stage), N_pos = as.numeric(cl$N_pos),
    M1       = as.numeric(cl$M1), grade = as.numeric(cl$grade_num),
    stromal  = as.numeric(scale(e$StromalScore)),
    immune   = as.numeric(scale(e$ImmuneScore)),
    noFeat   = as.numeric(scale(cl$pct_noFeature)),
    multimap = as.numeric(scale(cl$pct_multimapping)),
    depth    = as.numeric(scale(log10(cl$assigned_reads))))
}
KIRC_SPECS <- list(
  clinical  = c("age", "male", "T_stage", "N_pos", "M1", "grade"),
  principal = c("age", "male", "T_stage", "N_pos", "M1", "grade",
                "stromal", "immune", "noFeat", "multimap", "depth"))
# Keyed by biotype. NET_TAG maps it to the network names ("mrna", "lnc").
FR <- lapply(NET_TAG, function(k) list(scores = SC[[k]], base = kirc_base(rownames(SC[[k]]))))
stopifnot(identical(names(FR), names(NET_TAG)), all(TARGET$biotype %in% names(FR)))

for (i in seq_len(nrow(TARGET))) {
  fr <- FR[[TARGET$biotype[i]]]
  stopifnot(is.data.table(fr$base), TARGET$feature[i] %in% colnames(fr$scores))
  d  <- copy(fr$base)
  d[, ME := as.numeric(fr$scores[, TARGET$feature[i]])]
  for (sp in names(KIRC_SPECS))
    collect(sex_models(d, KIRC_SPECS[[sp]], "TCGA-KIRC", sp, "prognostic_modules",
                       TARGET$biotype[i], TARGET$module[i]))
}
# The additive HR under the principal specification must reproduce stage 03.
chk <- rbindlist(INT)[cohort == "TCGA-KIRC" & spec == "principal"]
chk <- merge(chk[, .(biotype, module, HR_additive)], TARGET[, .(biotype, module, HR_full)],
             by = c("biotype", "module"))
msg("Check against 03 (principal, additive): max |HR - HR_full| = ",
    signif(max(abs(chk$HR_additive - chk$HR_full)), 2))
stopifnot(max(abs(chk$HR_additive - chk$HR_full)) < 1e-6)

# ---- STAR metrics and the lncRNA axis as exposures (clinical set) ----
# As in stage 07, each exposure is standardised over all samples where it is
# available (528 for STAR metrics, 511 for the axis) before the complete-case
# restriction. log10_pct_noFeature is a transform of pct_noFeature, so it is
# flagged and excluded from the BH family.
kirc_exposure <- function(samples, x, biotype, name, is_transform = FALSE) {
  d <- kirc_base(samples)
  d[, ME := per_sd(x)]
  collect(sex_models(d, KIRC_SPECS$clinical, "TCGA-KIRC", "clinical", "quality_axis",
                     biotype, name, is_transform = is_transform))
}
kirc_exposure(cohort$sample_barcode, cohort$pct_noFeature, "quality_metric", "pct_noFeature")
kirc_exposure(cohort$sample_barcode, log10(pmax(cohort$pct_noFeature, 1e-3)),
              "quality_metric", "log10_pct_noFeature", is_transform = TRUE)
kirc_exposure(cohort$sample_barcode, cohort$pct_multimapping, "quality_metric",
              "pct_multimapping")
kirc_exposure(cohort$sample_barcode, log10(cohort$assigned_reads), "quality_metric",
              "log10_assigned_reads")
kirc_exposure(AXIS$sample_barcode, AXIS$lnc_axis, "axis", "lnc_axis")
# The additive HRs of the exposures modelled in stage 07 must match it to 3 dp.
chk_q <- rbindlist(INT)[family == "quality_axis"]
nf07  <- fread(file.path(RESULTS_DIR, "07_noFeature_as_exposure.tsv"))
ax07  <- fread(file.path(RESULTS_DIR, "07_axis_nested_adjustment.tsv"))
ref07 <- c(pct_noFeature       = nf07[exposure_scale == "linear, per SD" &
                                      grepl("^\\(b\\)", model), HR][1],
           log10_pct_noFeature = nf07[exposure_scale == "log10, per SD" &
                                      grepl("^\\(b\\)", model), HR][1],
           lnc_axis            = ax07[grepl("age, sex", model), HR][1])
for (v in names(ref07))
  msg("Check against 07 (", v, ", clinical set): 33 HR = ",
      round(chk_q[module == v, HR_additive], 4), " vs 07 HR = ", ref07[[v]],
      " (|diff| after rounding = ",
      signif(abs(round(chk_q[module == v, HR_additive], 3) - ref07[[v]]), 2), ")")
stopifnot(all(abs(chk_q[match(names(ref07), module), HR_additive] -
                  as.numeric(ref07)) < 0.005))

# ---- 2. CPTAC-3 (locked-model scores) and TCGA-KIRP (subtype cache) ----
banner("2 | CPTAC-3 and TCGA-KIRP")
vcl <- as.data.table(L$vcl); Ev <- L$Ev; Xv <- L$X$clinical$Xv
stopifnot(identical(rownames(Ev), vcl$sample_barcode),
          identical(rownames(Xv), vcl$sample_barcode))
d_cptac <- data.table(sample_barcode = vcl$sample_barcode,
                      os_time = vcl$os_time, os_event = vcl$os_event,
                      age = Xv[, "age"], male = Xv[, "male"], T_stage = Xv[, "T_stage"],
                      M1 = Xv[, "M1"], grade = Xv[, "grade"])
CPTAC_COVS <- c("age", "male", "T_stage", "M1", "grade")
msg("CPTAC-3 evaluated set: ", nrow(d_cptac), " patients, ", sum(d_cptac$os_event),
    " deaths; women ", sum(d_cptac$male == 0), " (", sum(d_cptac$os_event[d_cptac$male == 0]),
    " deaths), men ", sum(d_cptac$male == 1), " (", sum(d_cptac$os_event[d_cptac$male == 1]),
    " deaths)")
for (i in seq_len(nrow(TARGET))) {
  if (!TARGET$feature[i] %in% colnames(Ev)) next
  d <- copy(d_cptac); d[, ME := per_sd(Ev[, TARGET$feature[i]])]
  collect(sex_models(d, CPTAC_COVS, "CPTAC-3", "parsimonious_09", "prognostic_modules",
                     TARGET$biotype[i], TARGET$module[i]))
}
chk <- rbindlist(INT)[cohort == "CPTAC-3"]
rep_par <- fread(file.path(RESULTS_DIR, "09_module_replication_parsimonious.tsv"))
chk <- merge(chk[, .(biotype, module, HR_additive)], rep_par[, .(biotype, module, HR_cptac)],
             by = c("biotype", "module"))
msg("Check against 09 (parsimonious, additive): max |HR - HR_cptac| = ",
    signif(max(abs(chk$HR_additive - chk$HR_cptac)), 2))

# ---- subtype cohorts, processed as in stage 11 ----
# Sub-stages are collapsed to the roman numeral, as in stages 01 and 11.
parse_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "NOT REPORTED", "[NOT AVAILABLE]", "[UNKNOWN]")] <- NA
  rn <- sub("^STAGE\\s*", "", x); rn <- sub("[A-C]$", "", rn)
  out <- rep(NA_integer_, length(x))
  out[rn == "I"] <- 1L; out[rn == "II"] <- 2L; out[rn == "III"] <- 3L; out[rn == "IV"] <- 4L
  out
}
# Cohort with stage (cached by 11) and the M1 reconciliation of stages 01/11.
# log2(FPKM + 1) on network transcripts is residualised within cohort and
# projected onto the discovery loadings.
prep_subtype <- function(tag) {
  obj <- readRDS(file.path(CACHE_DIR, paste0("subtype_", tag, ".rds")))
  co  <- as.data.table(copy(obj$cohort))
  f_st <- file.path(CACHE_DIR, paste0("subtype_stage_", tag, ".rds"))
  if (file.exists(f_st)) {
    st <- readRDS(f_st)
    co[, stage_raw := st$stage_raw[match(patient, st$patient)]]
  } else {
    msg(tag, ": no cached stage table (run 11 first); no M reconciliation possible")
    co[, stage_raw := NA_character_]
  }
  co[, stage_num := parse_stage(stage_raw)]
  co[, M1_as_coded := as.integer(M1)]
  co[, M1_imputed  := !grepl("^M[01]", toupper(trimws(pM)))]
  co[, M1_stage_reconciled := as.integer(M1_imputed & stage_num %in% 4L &
                                         T_stage %in% 1:3 & M1_as_coded == 0L)]
  co[M1_stage_reconciled == 1L, M1 := 1L]
  g_net <- unique(c(colnames(nets$mrna$expr), colnames(nets$lnc$expr)))
  g     <- intersect(rownames(obj$fpkm), g_net)
  E_obs <- t(log2(obj$fpkm[g, , drop = FALSE] + 1))
  q     <- as.data.table(obj$qc)[match(rownames(E_obs), sample_barcode)]
  E     <- remove_technical(E_obs, tech_covariates(q))
  M     <- score_modules(E, LOAD)
  co    <- co[match(rownames(M), sample_barcode)]
  stopifnot(identical(co$sample_barcode, rownames(M)), !anyNA(co$sex))
  msg(obj$project, ": ", nrow(M), " scored (", sum(co$os_event), " deaths; women ",
      sum(co$sex == "female"), " with ", sum(co$os_event[co$sex == "female"]),
      " deaths, men ", sum(co$sex == "male"), " with ", sum(co$os_event[co$sex == "male"]),
      " deaths); M1 reconciled from stage IV n = ", sum(co$M1_stage_reconciled),
      "; T missing n = ", sum(is.na(co$T_stage)))
  list(project = obj$project, cohort = co, M = M, q = q, E_obs = E_obs)
}
SUB <- list(KIRP = prep_subtype("KIRP"), KICH = prep_subtype("KICH"))

# KIRP: the parsimonious covariate set of stage 11 on complete cases. Age,
# non-feature fraction and score are per SD over the analysed set.
kp <- SUB$KIRP
X0 <- data.table(sample_barcode = kp$cohort$sample_barcode,
                 os_time = kp$cohort$os_time, os_event = kp$cohort$os_event,
                 age = as.numeric(kp$cohort$age), male = as.numeric(kp$cohort$sex == "male"),
                 T_stage = as.numeric(kp$cohort$T_stage), M1 = as.numeric(kp$cohort$M1),
                 noFeat = as.numeric(kp$q$pct_noFeature))
cc <- complete.cases(X0) & is.finite(X0$os_time)
d_kirp <- X0[cc]
d_kirp[, `:=`(age = per_sd(age), noFeat = per_sd(noFeat))]
KIRP_COVS <- c("age", "male", "T_stage", "M1", "noFeat")
msg("TCGA-KIRP analysed set: ", nrow(d_kirp), " complete cases of ", nrow(X0), ", ",
    sum(d_kirp$os_event), " deaths; women ", sum(d_kirp$male == 0), " (",
    sum(d_kirp$os_event[d_kirp$male == 0]), " deaths), men ", sum(d_kirp$male == 1), " (",
    sum(d_kirp$os_event[d_kirp$male == 1]), " deaths)")
for (i in seq_len(nrow(TARGET))) {
  if (!TARGET$feature[i] %in% colnames(kp$M)) next
  d <- copy(d_kirp); d[, ME := per_sd(kp$M[cc, TARGET$feature[i]])]
  collect(sex_models(d, KIRP_COVS, "TCGA-KIRP", "parsimonious_11", "prognostic_modules",
                     TARGET$biotype[i], TARGET$module[i]))
}
chk <- rbindlist(INT)[cohort == "TCGA-KIRP"]
sub_eff <- fread(file.path(RESULTS_DIR, "11_subtype_module_effects.tsv"))[cohort == "TCGA-KIRP"]
chk <- merge(chk[, .(biotype, module, HR_additive)], sub_eff[, .(biotype, module, HR)],
             by = c("biotype", "module"))
msg("Check against 11 (KIRP, additive): max |HR - HR| = ",
    signif(max(abs(chk$HR_additive - chk$HR)), 2))

# ---- assemble, with BH within cohort x specification x family ----
# Rows flagged is_transform are left out of p.adjust() and carry NA.
int_tbl <- rbindlist(INT); str_tbl <- rbindlist(STR); con_tbl <- rbindlist(CON)
int_tbl[is_transform == FALSE,
        `:=`(fdr_interaction     = p.adjust(p_interaction,     "BH"),
             fdr_interaction_LRT = p.adjust(p_interaction_LRT, "BH"),
             n_tests_fdr         = .N), by = .(cohort, spec, family)]
# The sex-specific fits get their own BH over the same family. The pair-level
# contrast is repeated on both sex rows, so correcting within sex covers it.
str_tbl[is_transform == FALSE,
        `:=`(fdr_within_sex         = p.adjust(p, "BH"),
             fdr_sex_difference     = p.adjust(p_sex_difference,     "BH"),
             fdr_sex_difference_LRT = p.adjust(p_sex_difference_LRT, "BH"),
             n_tests_fdr            = .N),
        by = .(cohort, spec, family, sex)]
con_tbl[is_transform == FALSE,
        `:=`(strat_fdr_sex_difference     = p.adjust(strat_p_sex_difference, "BH"),
             strat_fdr_sex_difference_LRT = p.adjust(strat_p_sex_difference_LRT, "BH"),
             int_fdr_interaction          = p.adjust(int_p_interaction, "BH"),
             int_fdr_interaction_LRT      = p.adjust(int_p_interaction_LRT, "BH"),
             n_tests_fdr                  = .N), by = .(cohort, spec, family)]
# Transform rows carry the size of their family in n_tests_fdr.
int_tbl[, n_tests_fdr := max(n_tests_fdr, na.rm = TRUE), by = .(cohort, spec, family)]
str_tbl[, n_tests_fdr := max(n_tests_fdr, na.rm = TRUE), by = .(cohort, spec, family, sex)]
con_tbl[, n_tests_fdr := max(n_tests_fdr, na.rm = TRUE), by = .(cohort, spec, family)]
num_round <- function(tb, hr_cols, p_cols) {
  for (v in intersect(hr_cols, names(tb))) tb[, (v) := round(get(v), 4)]
  for (v in intersect(p_cols,  names(tb))) tb[, (v) := signif(get(v), 3)]
  if ("epv" %in% names(tb)) tb[, epv := round(epv, 2)]
  tb
}
int_tbl <- num_round(int_tbl,
  c("HR_additive", "lo_additive", "hi_additive",
    "HR_female_interaction_model", "lo_female_interaction_model", "hi_female_interaction_model",
    "HR_male_interaction_model", "lo_male_interaction_model", "hi_male_interaction_model",
    "HR_interaction", "lo_interaction", "hi_interaction"),
  c("p_additive", "p_female_interaction_model", "p_male_interaction_model",
    "p_interaction", "p_interaction_LRT", "fdr_interaction", "fdr_interaction_LRT"))
str_tbl <- num_round(str_tbl,
  c("HR", "lo", "hi", "HR_per_stratum_SD", "lo_per_stratum_SD", "hi_per_stratum_SD",
    "HR_ratio_male_over_female", "ratio_lo", "ratio_hi"),
  c("p", "fdr_within_sex", "p_sex_difference", "p_sex_difference_LRT",
    "fdr_sex_difference", "fdr_sex_difference_LRT"))
str_tbl[, sd_score_in_stratum := round(sd_score_in_stratum, 4)]
setcolorder(int_tbl, c("cohort", "spec", "family", "biotype", "module", "model",
                       "is_transform", "covariates", "n", "events", "n_female",
                       "events_female", "n_male", "events_male"))
setcolorder(str_tbl, c("cohort", "spec", "family", "biotype", "module", "model",
                       "is_transform", "sex", "covariates", "n", "events",
                       "HR", "lo", "hi", "p"))
save_tsv(int_tbl, "33_sex_interaction_models.tsv")
save_tsv(str_tbl, "33_sex_stratified_hr.tsv")
print(int_tbl[, .(cohort, spec, biotype, module, events_female, events_male,
                  HR_additive, HR_female_interaction_model, HR_male_interaction_model,
                  HR_interaction, p_interaction, fdr_interaction, epv, epv_below_min)],
      row.names = FALSE)
print(str_tbl[, .(cohort, spec, biotype, module, sex, n, events, HR, lo, hi, p,
                  epv, epv_below_min, implausible, covariates_dropped)], row.names = FALSE)

# ---- reconciliation table: both models side by side ----
# tests_agree records whether the two sex-difference tests agree in
# significance at FDR_ALPHA on raw p. They test different models, so they can
# disagree.
con_tbl[, tests_agree := (strat_p_sex_difference < FDR_ALPHA) ==
                         (int_p_interaction < FDR_ALPHA)]
con_tbl[, tests_agree_LRT := (strat_p_sex_difference_LRT < FDR_ALPHA) ==
                             (int_p_interaction_LRT < FDR_ALPHA)]
con_tbl[, disagreement := fifelse(tests_agree, "",
  fifelse(strat_p_sex_difference < FDR_ALPHA,
          "sex-specific fits show a difference, the common-covariate interaction does not",
          "common-covariate interaction shows a difference, the sex-specific fits do not"))]
stopifnot(max(con_tbl$strat_lrt_loglik_gap) < 1e-4)
con_tbl <- num_round(con_tbl,
  c("strat_HR_female", "strat_lo_female", "strat_hi_female",
    "strat_HR_male", "strat_lo_male", "strat_hi_male",
    "strat_HR_ratio_male_over_female", "strat_ratio_lo", "strat_ratio_hi",
    "int_HR_female", "int_lo_female", "int_hi_female",
    "int_HR_male", "int_lo_male", "int_hi_male",
    "int_HR_interaction", "int_lo_interaction", "int_hi_interaction"),
  c("strat_p_female", "strat_p_male", "strat_p_sex_difference",
    "strat_p_sex_difference_LRT", "strat_fdr_sex_difference",
    "strat_fdr_sex_difference_LRT", "int_p_female", "int_p_male",
    "int_p_interaction", "int_p_interaction_LRT", "int_fdr_interaction",
    "int_fdr_interaction_LRT"))
con_tbl[, strat_lrt_loglik_gap := signif(strat_lrt_loglik_gap, 3)]
setcolorder(con_tbl, c("cohort", "spec", "family", "biotype", "module", "is_transform",
                       "n", "events", "n_female", "events_female", "n_male", "events_male",
                       "tests_agree", "disagreement"))
save_tsv(con_tbl, "33_sex_contrast_reconciliation.tsv")
print(con_tbl[, .(cohort, spec, biotype, module,
                  strat_HR_female, strat_p_female, strat_HR_male, strat_p_male,
                  strat_p_sex_difference, int_HR_interaction, int_p_interaction,
                  tests_agree)], row.names = FALSE)
dis <- con_tbl[tests_agree == FALSE]
msg("Sex-difference tests disagreeing in significance at ", FDR_ALPHA, ": ", nrow(dis),
    " of ", nrow(con_tbl), " model rows")
if (nrow(dis)) print(dis[, .(cohort, spec, biotype, module, strat_p_sex_difference,
                             int_p_interaction, disagreement)], row.names = FALSE)

# ---- 3. association with sex: module scores, STAR metrics and the axis ----
banner("3 | Correlation with sex in every cohort")
# Point-biserial r is Pearson with the male indicator (positive = higher in
# men). The Wilcoxon rank-sum test does not assume normality.
sex_assoc <- function(x, male, cohort_name, family, biotype, variable) {
  ok <- is.finite(x) & !is.na(male); x <- x[ok]; male <- male[ok]
  if (sum(male == 1) < 3 || sum(male == 0) < 3) return(NULL)
  ct <- suppressWarnings(cor.test(x, male))
  wt <- suppressWarnings(wilcox.test(x[male == 1], x[male == 0], exact = FALSE))
  data.table(cohort = cohort_name, family = family, biotype = biotype, variable = variable,
             n = length(x), n_female = sum(male == 0), n_male = sum(male == 1),
             median_female = median(x[male == 0]), median_male = median(x[male == 1]),
             mean_female = mean(x[male == 0]), mean_male = mean(x[male == 1]),
             r_pb = unname(ct$estimate), r_lo = ct$conf.int[1], r_hi = ct$conf.int[2],
             p_pearson = ct$p.value, W = unname(wt$statistic), p_wilcoxon = wt$p.value)
}
# Leading axis as in stages 07/09/11: column-centred, unscaled PC1,
# sign-aligned to mean expression.
pc1_of <- function(E) {
  x <- scale(E, center = TRUE, scale = FALSE)
  s <- svd(x, nu = 1, nv = 1)
  pc <- s$u[, 1] * s$d[1]
  if (stats::cor(pc, rowMeans(E)) < 0) pc <- -pc
  setNames(pc, rownames(E))
}
cohort_assoc <- function(cohort_name, M_list, qc, axis, male_lookup) {
  rows <- list()
  for (M in M_list) for (cn in colnames(M)) {
    bt  <- if (startsWith(cn, "mRNA_")) "mRNA" else "lncRNA"
    rows[[length(rows) + 1]] <- sex_assoc(M[, cn], male_lookup[rownames(M)], cohort_name,
                                          "module", bt, sub("^(mRNA|lnc)_ME", "", cn))
  }
  # The STAR metrics and the axis form one BH family, as in the model table.
  qm <- list(pct_noFeature = qc$pct_noFeature, pct_multimapping = qc$pct_multimapping,
             log10_assigned_reads = log10(qc$assigned_reads))
  for (v in names(qm))
    rows[[length(rows) + 1]] <- sex_assoc(qm[[v]], male_lookup[qc$sample_barcode],
                                          cohort_name, "quality_axis", "quality_metric", v)
  rows[[length(rows) + 1]] <- sex_assoc(axis, male_lookup[names(axis)], cohort_name,
                                        "quality_axis", "lncRNA", "lnc_axis")
  rbindlist(rows)
}
male_of <- function(cl) setNames(as.numeric(cl$sex == "male"), cl$sample_barcode)

# TCGA-KIRC: discovery scores, metrics for all 528, and the axis from stage 07.
cor_kirc <- cohort_assoc("TCGA-KIRC", list(SC$mrna, SC$lnc),
                         cohort[, .(sample_barcode, pct_noFeature, pct_multimapping, assigned_reads)],
                         setNames(AXIS$lnc_axis, AXIS$sample_barcode), male_of(cohort))
# CPTAC-3: Ev on the evaluated set, metrics and the axis on all clear-cell
# patients (observed log2(FPKM + 1) on network lncRNAs, as in stage 09).
vd  <- readRDS(file.path(CACHE_DIR, "validation_dataset.rds"))
vco <- as.data.table(vd$cohort)
vqc <- as.data.table(vd$qc)[match(vco$sample_barcode, sample_barcode)]
gl  <- intersect(colnames(nets$lnc$expr), rownames(vd$fpkm))
Vl_obs <- t(log2(vd$fpkm[gl, vco$sample_barcode, drop = FALSE] + 1))
cor_cptac <- cohort_assoc("CPTAC-3", list(Ev), vqc, pc1_of(Vl_obs), male_of(vco))
msg("CPTAC-3 clear-cell cohort: ", nrow(vco), " patients (women ", sum(vco$sex == "female"),
    ", men ", sum(vco$sex == "male"), "); lncRNA axis on ", length(gl), " transcripts")
# KIRP and KICH
cor_sub <- rbindlist(lapply(SUB, function(s) {
  gl_s <- intersect(colnames(s$E_obs), colnames(nets$lnc$expr))
  cohort_assoc(s$project, list(s$M), s$q, pc1_of(s$E_obs[, gl_s, drop = FALSE]),
               male_of(s$cohort))
}))
cor_tbl <- rbindlist(list(cor_kirc, cor_cptac, cor_sub))
cor_tbl[, prognostic := family == "module" &
                        paste(biotype, variable) %in% TARGET[, paste(biotype, module)]]
cor_tbl[, direction := fifelse(r_pb > 0, "higher in men", "higher in women")]
cor_tbl[, fdr_wilcoxon := p.adjust(p_wilcoxon, "BH"), by = .(cohort, family)]
cor_tbl[, n_tests_fdr := .N, by = .(cohort, family)]
for (v in c("median_female", "median_male", "mean_female", "mean_male", "r_pb", "r_lo", "r_hi"))
  cor_tbl[, (v) := round(get(v), 4)]
for (v in c("p_pearson", "p_wilcoxon", "fdr_wilcoxon")) cor_tbl[, (v) := signif(get(v), 3)]
cor_tbl[, cohort := factor(cohort, levels = c("TCGA-KIRC", "CPTAC-3", "TCGA-KIRP", "TCGA-KICH"))]
setorderv(cor_tbl, c("cohort", "family", "prognostic", "biotype", "variable"),
          c(1L, 1L, -1L, 1L, 1L), na.last = TRUE)
cor_tbl[, cohort := as.character(cohort)]
setcolorder(cor_tbl, c("cohort", "family", "biotype", "variable", "prognostic"))
save_tsv(cor_tbl, "33_sex_correlations.tsv")
print(cor_tbl[prognostic == TRUE | family != "module",
              .(cohort, family, biotype, variable, n_female, n_male, r_pb, r_lo, r_hi,
                p_wilcoxon, fdr_wilcoxon, direction)], row.names = FALSE)
sig_sex <- cor_tbl[family == "module" & fdr_wilcoxon < FDR_ALPHA,
                   .(mods = paste(paste(biotype, variable), collapse = ", ")), by = cohort]
msg("Module scores associated with sex at FDR < ", FDR_ALPHA, " (per cohort): ",
    if (nrow(sig_sex)) paste(sig_sex[, paste0(cohort, ": ", mods)], collapse = " | ")
    else "none")
# Summary sentence per cohort, generated from cor_tbl.
sex_sentence <- function(coh) {
  m <- cor_tbl[cohort == coh & family == "module"]
  if (!nrow(m)) return(paste0(coh, ": no module scores"))
  s <- m[fdr_wilcoxon < FDR_ALPHA][order(fdr_wilcoxon, -abs(r_pb))]
  paste0(coh, ": ", nrow(s), " of the ", nrow(m),
         " module scores are associated with sex at BH < ", FDR_ALPHA, ", including ",
         sum(s$prognostic), " of the ", sum(m$prognostic), " prognostic modules",
         if (nrow(s)) paste0(" [", paste0(s$biotype, " ", s$variable, " r_pb ", s$r_pb,
                                          " BH ", s$fdr_wilcoxon, collapse = "; "), "]")
         else "")
}
SEX_SENTENCE <- vapply(levels(factor(cor_tbl$cohort,
                                     levels = c("TCGA-KIRC", "CPTAC-3", "TCGA-KIRP",
                                                "TCGA-KICH"))),
                       sex_sentence, character(1))
for (s in SEX_SENTENCE) msg(s)

# ---- 4. X- and Y-encoded members of the lncRNA modules ----
banner("4 | chrX / chrY members of the lncRNA modules")
genes <- fread(file.path(RESULTS_DIR, "02_lncRNA_module_genes.tsv"))
msg("lncRNA network: ", nrow(genes), " transcripts in ", uniqueN(genes$module),
    " modules (grey included)")

# Gene records from the GENCODE v36 GTF, parsed once and cached. GDC STAR
# counts use this annotation, so versioned gene IDs match.
GTF_GZ   <- file.path(CACHE_DIR, "28_gencode_v36.gtf.gz")
GTF_RDS  <- file.path(CACHE_DIR, "33_gencode_v36_genes.rds")
parse_gencode_genes <- function(gz) {
  con <- gzfile(gz, open = "rt"); on.exit(close(con))
  out <- list(); n_lines <- 0
  repeat {
    ln <- readLines(con, n = 400000L)
    if (!length(ln)) break
    n_lines <- n_lines + length(ln)
    g <- ln[grepl("\tgene\t", ln, fixed = TRUE)]
    if (length(g)) {
      f <- tstrsplit(g, "\t", fixed = TRUE)
      a <- f[[9]]
      out[[length(out) + 1]] <- data.table(
        gene_id   = sub('.*gene_id "([^"]+)".*',   "\\1", a),
        gene_name = sub('.*gene_name "([^"]+)".*', "\\1", a),
        gene_type = sub('.*gene_type "([^"]+)".*', "\\1", a),
        chr = f[[1]], start = as.integer(f[[4]]), end = as.integer(f[[5]]), strand = f[[7]])
    }
  }
  msg("  GTF: ", n_lines, " lines read")
  rbindlist(out)
}
ann <- NULL; ann_source <- NA_character_
if (file.exists(GTF_RDS)) {
  ann <- readRDS(GTF_RDS); ann_source <- "GENCODE v36 GTF (cache/28_gencode_v36.gtf.gz)"
  msg("Gene coordinates from cached parse of the GENCODE v36 GTF: ", nrow(ann), " genes")
} else if (file.exists(GTF_GZ)) {
  ann <- tryCatch(parse_gencode_genes(GTF_GZ), error = function(e) {
    msg("  GTF parse failed: ", conditionMessage(e)); NULL })
  if (!is.null(ann) && nrow(ann) > 50000) {
    saveRDS(ann, GTF_RDS)
    ann_source <- "GENCODE v36 GTF (cache/28_gencode_v36.gtf.gz)"
    msg("Gene coordinates parsed from the GENCODE v36 GTF: ", nrow(ann), " genes")
  } else ann <- NULL
}
if (is.null(ann)) {
  # Fallback: Ensembl REST lookup/id in batches of 1000 unversioned IDs, with
  # up to five attempts per batch (the service returns 5xx intermittently).
  msg("GENCODE GTF not available; querying the Ensembl REST API")
  ENS_RDS <- file.path(CACHE_DIR, "33_ensembl_chr_lookup.rds")
  ensembl_lookup <- function(ids) {
    ch <- split(ids, ceiling(seq_along(ids) / 1000))
    res <- list()
    for (k in seq_along(ch)) {
      j <- NULL
      for (a in 1:5) {
        r <- tryCatch(httr::POST("https://rest.ensembl.org/lookup/id",
                                 body = jsonlite::toJSON(list(ids = ch[[k]])),
                                 httr::content_type_json(), httr::accept_json(),
                                 httr::timeout(180)), error = function(e) NULL)
        if (!is.null(r) && httr::status_code(r) == 200) {
          j <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"),
                                  simplifyVector = FALSE); break
        }
        msg("  batch ", k, " attempt ", a, " failed (",
            if (is.null(r)) "no response" else httr::status_code(r), ")")
        Sys.sleep(5 * a)
      }
      if (is.null(j)) return(NULL)
      res[[k]] <- rbindlist(lapply(names(j), function(id) {
        x <- j[[id]]
        if (is.null(x) || is.null(x$seq_region_name))
          return(data.table(gene_id_unv = id, gene_name = NA_character_,
                            gene_type = NA_character_, chr = NA_character_,
                            start = NA_integer_, end = NA_integer_, strand = NA_character_))
        data.table(gene_id_unv = id,
                   gene_name = if (is.null(x$display_name)) NA_character_ else x$display_name,
                   gene_type = if (is.null(x$biotype)) NA_character_ else x$biotype,
                   chr = paste0("chr", x$seq_region_name),
                   start = as.integer(x$start), end = as.integer(x$end),
                   strand = if (identical(x$strand, 1L) || identical(x$strand, 1)) "+" else "-")
      }))
      msg("  batch ", k, "/", length(ch), " done")
      Sys.sleep(1)
    }
    rbindlist(res)
  }
  ens <- if (file.exists(ENS_RDS)) readRDS(ENS_RDS) else NULL
  if (is.null(ens)) {
    ens <- ensembl_lookup(unique(sub("\\..*$", "", genes$gene_id)))
    if (!is.null(ens)) saveRDS(ens, ENS_RDS)
  }
  if (!is.null(ens)) {
    ann <- copy(ens); ann[, gene_id := genes$gene_id[match(gene_id_unv, sub("\\..*$", "", genes$gene_id))]]
    ann[, gene_id_unv := NULL]
    ann_source <- "Ensembl REST lookup/id (GRCh38, current release)"
  }
}

if (is.null(ann)) {
  writeLines(c("33_sex_stratified.R: neither cache/28_gencode_v36.gtf.gz nor the Ensembl REST API",
               "was available; chromosome columns of 33_xy_lncRNA_members.tsv are empty.",
               format(Sys.time())), file.path(RESULTS_DIR, "33_NOTE_annotation_unreachable.txt"))
  xy <- genes[0][, `:=`(chr = character(), start = integer(), end = integer(),
                        strand = character(), gene_type_annotation = character())]
  counts <- data.table(module = unique(genes$module), n_genes = NA_integer_)
  named  <- data.table()
} else {
  gm <- merge(genes, ann[, .(gene_id, chr, start, end, strand, gene_type_annotation = gene_type,
                             gene_name_annotation = gene_name)],
              by = "gene_id", all.x = TRUE)
  msg("Network transcripts with a chromosome: ", sum(!is.na(gm$chr)), "/", nrow(gm),
      " (source: ", ann_source, ")")
  gm[, prognostic_module := module %in% TARGET[biotype == "lncRNA", module]]
  gm[, named_group := fifelse(gene_name %in% c("XIST", "TSIX"), gene_name,
                       fifelse(chr == "chrY", "Y-linked", ""))]
  xy <- gm[chr %in% c("chrX", "chrY")]
  xy[, source := ann_source]
  setorder(xy, chr, module, -kME)
  xy <- xy[, .(gene_id, gene_name, chr, start, end, strand, gene_type_annotation,
               module, kME = round(kME, 4), prognostic_module, named_group, source)]

  # counts per module, with the annotation-wide lncRNA denominators
  counts <- gm[, .(n_genes = .N, n_with_chr = sum(!is.na(chr)),
                   n_chrX = sum(chr == "chrX", na.rm = TRUE),
                   n_chrY = sum(chr == "chrY", na.rm = TRUE)), by = module]
  counts[, n_XY := n_chrX + n_chrY]
  counts[, pct_XY := round(100 * n_XY / n_genes, 2)]
  counts[, prognostic_module := module %in% TARGET[biotype == "lncRNA", module]]
  counts[, row_type := "module"]
  counts[, xy_genes := vapply(module, function(m)
    paste(xy[module == m][order(chr, -kME), gene_name], collapse = ";"), character(1))]
  setorderv(counts, c("prognostic_module", "n_XY"), c(-1L, -1L), na.last = TRUE)
  tot <- gm[, .(module = "all network transcripts", n_genes = .N, n_with_chr = sum(!is.na(chr)),
                n_chrX = sum(chr == "chrX", na.rm = TRUE), n_chrY = sum(chr == "chrY", na.rm = TRUE))]
  # Summary rows: prognostic_module FALSE rather than NA, row_type "summary".
  tot[, `:=`(n_XY = n_chrX + n_chrY, pct_XY = round(100 * (n_chrX + n_chrY) / n_genes, 2),
             prognostic_module = FALSE, row_type = "summary", xy_genes = "")]
  if (identical(ann_source, "GENCODE v36 GTF (cache/28_gencode_v36.gtf.gz)")) {
    al <- ann[gene_type == "lncRNA" & !grepl("_PAR_Y$", gene_id)]
    tot <- rbind(tot, data.table(
      module = "GENCODE v36 lncRNA genes (annotation)", n_genes = nrow(al), n_with_chr = nrow(al),
      n_chrX = sum(al$chr == "chrX"), n_chrY = sum(al$chr == "chrY"),
      n_XY = sum(al$chr %in% c("chrX", "chrY")),
      pct_XY = round(100 * mean(al$chr %in% c("chrX", "chrY")), 2),
      prognostic_module = FALSE, row_type = "summary", xy_genes = ""))
  }
  counts <- rbind(counts, tot)
  counts[, source := ann_source]

  # ---- X/Y enrichment per module ----
  # Fisher's exact test (small counts) against GENCODE v36 lncRNAs minus the
  # tested genes, so the two arms of the 2x2 table do not overlap.
  if (identical(ann_source, "GENCODE v36 GTF (cache/28_gencode_v36.gtf.gz)")) {
    bg_all   <- ann[gene_type == "lncRNA" & !grepl("_PAR_Y$", gene_id)]
    bg_n     <- nrow(bg_all); bg_xy <- sum(bg_all$chr %in% c("chrX", "chrY"))
    set_ids  <- function(m) if (m == "all network transcripts") gm$gene_id
                            else if (m == "GENCODE v36 lncRNA genes (annotation)") character(0)
                            else gm[module == m, gene_id]
    fisher_row <- function(m, n_set, xy_set) {
      ids <- set_ids(m)
      if (!length(ids)) return(data.table(n_genes_ref = NA_integer_, n_XY_ref = NA_integer_,
                                          pct_XY_ref = NA_real_, OR_fisher = NA_real_,
                                          or_lo = NA_real_, or_hi = NA_real_,
                                          p_fisher = NA_real_))
      inb <- bg_all[gene_id %in% ids]
      n_ref  <- bg_n - nrow(inb); xy_ref <- bg_xy - sum(inb$chr %in% c("chrX", "chrY"))
      ft <- fisher.test(matrix(c(xy_set, n_set - xy_set, xy_ref, n_ref - xy_ref), nrow = 2))
      data.table(n_genes_ref = n_ref, n_XY_ref = xy_ref,
                 pct_XY_ref = round(100 * xy_ref / n_ref, 2),
                 OR_fisher = round(unname(ft$estimate), 3),
                 or_lo = round(ft$conf.int[1], 3), or_hi = round(ft$conf.int[2], 3),
                 p_fisher = signif(ft$p.value, 3))
    }
    counts <- cbind(counts, rbindlist(lapply(seq_len(nrow(counts)), function(i)
      fisher_row(counts$module[i], counts$n_genes[i], counts$n_XY[i]))))
    counts[, ref_set := fifelse(is.na(p_fisher), "not applicable (this row is the reference)",
                                "GENCODE v36 lncRNA genes not in the tested set")]
    counts[, fdr_fisher := NA_real_]
    counts[row_type == "module",
           fdr_fisher := signif(p.adjust(p_fisher, "BH"), 3)]
    msg("X/Y enrichment (Fisher, vs the GENCODE lncRNA background): network as a whole ",
        counts[module == "all network transcripts",
               paste0(n_XY, "/", n_genes, " = ", pct_XY, "% vs ", pct_XY_ref,
                      "%, OR ", OR_fisher, " (", or_lo, "-", or_hi, "), p ", p_fisher)])
    for (m in TARGET[biotype == "lncRNA", module])
      msg("  ", m, ": ", counts[module == m,
          paste0(n_XY, "/", n_genes, " = ", pct_XY, "%, OR ", OR_fisher,
                 " (", or_lo, "-", or_hi, "), p ", p_fisher, ", BH ", fdr_fisher)])
  }

  # the named genes: XIST, TSIX, every chrY lncRNA in the annotation that
  # entered the network, and the X/Y members of the prognostic modules
  named_ids <- unique(c(ann[gene_name %in% c("XIST", "TSIX"), gene_id],
                        xy[chr == "chrY", gene_id], xy[prognostic_module == TRUE, gene_id]))
  named <- ann[gene_id %in% named_ids, .(gene_id, gene_name, chr, gene_type_annotation = gene_type)]
  named <- merge(named, genes[, .(gene_id, module, kME)], by = "gene_id", all.x = TRUE)
  named[, in_network := !is.na(module)]
  named[, status := fifelse(is.na(module), "not in network (below the expression or MAD filter of 01/02)",
                     fifelse(module == "grey", "in network, unassigned (grey)",
                             paste0("module ", module)))]
  named[, prognostic_module := module %in% TARGET[biotype == "lncRNA", module]]
  named[, group := fifelse(gene_name %in% c("XIST", "TSIX"), gene_name,
                    fifelse(chr == "chrY", "Y-linked", "X-linked, prognostic module"))]
  named[, kME := round(kME, 4)]
  setorderv(named, c("group", "chr", "in_network", "kME"), c(1L, 1L, -1L, -1L),
            na.last = TRUE)
}
save_tsv(xy,     "33_xy_lncRNA_members.tsv")
save_tsv(counts, "33_xy_lncRNA_counts_by_module.tsv")
save_tsv(named,  "33_xy_named_genes.tsv")
print(counts[, .SD, .SDcols = intersect(c("module", "row_type", "n_genes", "n_chrX", "n_chrY",
                                          "n_XY", "pct_XY", "pct_XY_ref", "OR_fisher",
                                          "p_fisher", "fdr_fisher", "prognostic_module"),
                                        names(counts))], row.names = FALSE)
if (nrow(named)) print(named[, .(group, gene_name, chr, status, kME)], row.names = FALSE)

# ---- 5. STAR metrics and the lncRNA axis by sex in TCGA-KIRC ----
# 33_noFeature_by_sex.tsv holds all five exposures modelled in section 1.
banner("5 | STAR metrics and lncRNA axis by sex (TCGA-KIRC)")
by_sex <- function(x, male, variable) {
  ok <- is.finite(x) & !is.na(male); x <- x[ok]; male <- male[ok]
  qf <- quantile(x[male == 0], c(0.25, 0.5, 0.75)); qm <- quantile(x[male == 1], c(0.25, 0.5, 0.75))
  wt <- suppressWarnings(wilcox.test(x[male == 1], x[male == 0], exact = FALSE, conf.int = TRUE))
  ct <- suppressWarnings(cor.test(x, male))
  data.table(cohort = "TCGA-KIRC", variable = variable, n = length(x),
             n_female = sum(male == 0), n_male = sum(male == 1),
             median_female = qf[[2]], q25_female = qf[[1]], q75_female = qf[[3]],
             median_male = qm[[2]], q25_male = qm[[1]], q75_male = qm[[3]],
             mean_female = mean(x[male == 0]), sd_female = sd(x[male == 0]),
             mean_male = mean(x[male == 1]), sd_male = sd(x[male == 1]),
             diff_median_male_minus_female = qm[[2]] - qf[[2]],
             hodges_lehmann_shift = unname(wt$estimate),
             hl_lo = wt$conf.int[1], hl_hi = wt$conf.int[2],
             W = unname(wt$statistic), p_wilcoxon = wt$p.value,
             r_pb = unname(ct$estimate), r_lo = ct$conf.int[1], r_hi = ct$conf.int[2],
             p_pearson = ct$p.value)
}
mk <- male_of(cohort)
nf_tbl <- rbind(
  by_sex(cohort$pct_noFeature, mk[cohort$sample_barcode], "pct_noFeature"),
  by_sex(log10(pmax(cohort$pct_noFeature, 1e-3)), mk[cohort$sample_barcode], "log10_pct_noFeature"),
  by_sex(cohort$pct_multimapping, mk[cohort$sample_barcode], "pct_multimapping"),
  by_sex(log10(cohort$assigned_reads), mk[cohort$sample_barcode], "log10_assigned_reads"),
  by_sex(AXIS$lnc_axis, mk[AXIS$sample_barcode], "lnc_axis"))
num_cols <- setdiff(names(nf_tbl), c("cohort", "variable", "n", "n_female", "n_male", "W",
                                     "p_wilcoxon", "p_pearson"))
for (v in num_cols) nf_tbl[, (v) := round(get(v), 4)]
for (v in c("p_wilcoxon", "p_pearson")) nf_tbl[, (v) := signif(get(v), 3)]
# log10_pct_noFeature re-expresses pct_noFeature (identical Wilcoxon p).
nf_tbl[, is_transform := variable == "log10_pct_noFeature"]
setcolorder(nf_tbl, c("cohort", "variable", "is_transform"))
save_tsv(nf_tbl, "33_noFeature_by_sex.tsv")
print(nf_tbl[, .(variable, is_transform, n_female, n_male, median_female, median_male,
                 hodges_lehmann_shift, hl_lo, hl_hi, p_wilcoxon, r_pb)], row.names = FALSE)

# ---- 6. sensitivity: X/Y members deleted, and recorded versus marker-based sex ----
banner("6 | Sensitivity: X/Y members deleted; recorded vs marker-based sex")
SENS <- list()
# One TCGA-KIRC principal-specification model as a single row, built from the
# reconciliation block (strat_* and int_* columns). These rows are not part
# of any BH family.
sens_model <- function(d, biotype, module, analysis, variant, extra = list()) {
  r <- sex_models(d, KIRC_SPECS$principal, "TCGA-KIRC", "principal",
                  paste0("sensitivity: ", analysis), biotype, module)
  x <- copy(r$contrast)
  x[, `:=`(analysis = analysis, variant = variant)]
  a <- sex_assoc(d$ME, d$male, "TCGA-KIRC", "sensitivity", biotype, module)
  x[, `:=`(r_pb = round(a$r_pb, 4), p_wilcoxon = signif(a$p_wilcoxon, 3))]
  ii <- r$interaction
  x[, `:=`(covariates = ii$covariates,
           HR_additive = ii$HR_additive, lo_additive = ii$lo_additive,
           hi_additive = ii$hi_additive, p_additive = ii$p_additive,
           n_parameters = ii$n_parameters, epv = ii$epv,
           epv_below_min = ii$epv_below_min,
           covariates_dropped = ii$covariates_dropped)]
  for (nm in names(extra)) x[, (nm) := extra[[nm]]]
  SENS[[length(SENS) + 1]] <<- x
  invisible(x)
}

# ---- 6a. leave-out of chrX/chrY members ----
# Delete every chrX/chrY member, refit the module loadings on the remaining
# genes (fit_module_loadings) and re-score.
if (!is.null(ann)) {
  xy_ids <- ann[chr %in% c("chrX", "chrY"), gene_id]
  for (i in seq_len(nrow(TARGET))) {
    bt  <- TARGET$biotype[i]; mo <- TARGET$module[i]
    net <- nets[[NET_TAG[[bt]]]]
    gt  <- as.data.table(net$gene_tbl)[module == mo, .(gene_id, module)]
    kp2 <- gt[!gene_id %in% xy_ids]
    n_drop <- nrow(gt) - nrow(kp2)
    d_full <- copy(FR[[bt]]$base)
    d_full[, ME := as.numeric(FR[[bt]]$scores[, TARGET$feature[i]])]
    sens_model(d_full, bt, mo, "drop_XY_members", "full module (baseline)",
               list(n_genes_used = nrow(gt), n_genes_dropped = 0L,
                    cor_with_full_score = 1, n_discordant_sex_labels = 0L))
    if (n_drop == 0L) {
      msg(bt, " ", mo, ": no chrX/chrY members, leave-out not applicable")
      next
    }
    L2 <- fit_module_loadings(net$expr, kp2, prefix = "drop_")
    M2 <- score_modules(net$expr, L2)
    stopifnot(identical(rownames(M2), rownames(FR[[bt]]$scores)))
    d_drop <- copy(FR[[bt]]$base); d_drop[, ME := as.numeric(M2[, 1])]
    rho <- cor(d_drop$ME, d_full$ME)
    sens_model(d_drop, bt, mo, "drop_XY_members", "chrX/chrY members deleted",
               list(n_genes_used = nrow(kp2), n_genes_dropped = n_drop,
                    cor_with_full_score = round(rho, 4), n_discordant_sex_labels = 0L))
    msg(bt, " ", mo, ": ", n_drop, " X/Y member(s) deleted, r(full, reduced) = ",
        round(rho, 4))
  }
}

# ---- 6b. recorded sex against XIST and the chrY lncRNAs ----
# Fixed cuts on log2(FPKM + 1) with an indeterminate band: XIST on above 1 and off
# below 0.5, mean chrY lncRNA on above 0.2 and off below 0.05. A sample is called only
# when the two markers agree. Otherwise it is "ambiguous" and not counted as discordant.
lab_tbl <- data.table(); disc <- character(0)
if (!is.null(ann)) {
  E_lnc   <- obs_expr(nets$lnc)
  xist_id <- intersect(ann[gene_name == "XIST", gene_id], colnames(E_lnc))
  y_ids   <- intersect(ann[chr == "chrY", gene_id], colnames(E_lnc))
  if (length(xist_id) == 1L && length(y_ids) >= 1L) {
    XIST_ON <- 1.0; XIST_OFF <- 0.5; CHRY_ON <- 0.2; CHRY_OFF <- 0.05
    lab_tbl <- data.table(
      sample_barcode = rownames(E_lnc),
      recorded_sex   = ifelse(mk[rownames(E_lnc)] == 1, "male", "female"),
      XIST_log2fpkm  = round(as.numeric(E_lnc[, xist_id]), 4),
      chrY_lnc_mean_log2fpkm = round(rowMeans(E_lnc[, y_ids, drop = FALSE]), 4),
      n_chrY_lncRNA_used = length(y_ids),
      XIST_on_cut = XIST_ON, XIST_off_cut = XIST_OFF,
      chrY_on_cut = CHRY_ON, chrY_off_cut = CHRY_OFF)
    lab_tbl[, XIST_call := fifelse(XIST_log2fpkm > XIST_ON, "on",
                            fifelse(XIST_log2fpkm < XIST_OFF, "off", "indeterminate"))]
    lab_tbl[, chrY_call := fifelse(chrY_lnc_mean_log2fpkm > CHRY_ON, "on",
                            fifelse(chrY_lnc_mean_log2fpkm < CHRY_OFF, "off", "indeterminate"))]
    lab_tbl[, marker_sex := fifelse(XIST_call == "on"  & chrY_call == "off", "female",
                             fifelse(XIST_call == "off" & chrY_call == "on", "male",
                                     "ambiguous"))]
    lab_tbl[, concordance := fifelse(marker_sex == "ambiguous", "ambiguous",
                             fifelse(marker_sex == recorded_sex, "concordant", "discordant"))]
    setorderv(lab_tbl, c("concordance", "recorded_sex", "XIST_log2fpkm"))
    save_tsv(lab_tbl, "33_sex_label_check.tsv")
    disc <- lab_tbl[concordance == "discordant", sample_barcode]
    msg("Recorded sex vs XIST and ", length(y_ids), " chrY lncRNAs on ", nrow(lab_tbl),
        " discovery samples: ", sum(lab_tbl$concordance == "concordant"), " concordant, ",
        length(disc), " discordant, ", sum(lab_tbl$concordance == "ambiguous"), " ambiguous")
    if (length(disc))
      print(lab_tbl[concordance == "discordant",
                    .(sample_barcode, recorded_sex, marker_sex, XIST_log2fpkm,
                      chrY_lnc_mean_log2fpkm)], row.names = FALSE)
    # Refit the prognostic-module models with the discordant labels corrected,
    # and again with those samples removed.
    if (length(disc)) {
      mk2 <- mk
      mk2[disc] <- as.numeric(lab_tbl[concordance == "discordant"][
        match(disc, sample_barcode), marker_sex] == "male")
      for (i in seq_len(nrow(TARGET))) {
        bt <- TARGET$biotype[i]; mo <- TARGET$module[i]
        ng <- nrow(as.data.table(nets[[NET_TAG[[bt]]]]$gene_tbl)[module == mo])
        d0 <- copy(FR[[bt]]$base)
        d0[, ME := as.numeric(FR[[bt]]$scores[, TARGET$feature[i]])]
        nd <- sum(d0$sample_barcode %in% disc)
        base_extra <- list(n_genes_used = ng, n_genes_dropped = 0L, cor_with_full_score = 1)
        sens_model(d0, bt, mo, "sex_label_check", "recorded sex (baseline)",
                   c(base_extra, list(n_discordant_sex_labels = 0L)))
        d1 <- copy(d0); d1[, male := as.numeric(mk2[sample_barcode])]
        sens_model(d1, bt, mo, "sex_label_check", "marker-based sex where discordant",
                   c(base_extra, list(n_discordant_sex_labels = nd)))
        d2 <- d0[!sample_barcode %in% disc]
        sens_model(d2, bt, mo, "sex_label_check", "discordant samples dropped",
                   c(base_extra, list(n_discordant_sex_labels = nd)))
      }
    }
  } else msg("XIST or the chrY lncRNAs are not in the network matrix; label check skipped")
}
if (length(SENS)) {
  sens_tbl <- rbindlist(SENS, fill = TRUE)
  sens_tbl <- num_round(sens_tbl,
    c("HR_additive", "lo_additive", "hi_additive",
      "strat_HR_female", "strat_lo_female", "strat_hi_female",
      "strat_HR_male", "strat_lo_male", "strat_hi_male",
      "strat_HR_ratio_male_over_female", "strat_ratio_lo", "strat_ratio_hi",
      "int_HR_female", "int_lo_female", "int_hi_female",
      "int_HR_male", "int_lo_male", "int_hi_male",
      "int_HR_interaction", "int_lo_interaction", "int_hi_interaction"),
    c("p_additive", "strat_p_female", "strat_p_male", "strat_p_sex_difference",
      "strat_p_sex_difference_LRT", "int_p_female", "int_p_male",
      "int_p_interaction", "int_p_interaction_LRT"))
  sens_tbl[, strat_lrt_loglik_gap := signif(strat_lrt_loglik_gap, 3)]
  drop_cols <- intersect(c("family", "is_transform"), names(sens_tbl))
  if (length(drop_cols)) sens_tbl[, (drop_cols) := NULL]
  # Columns grouped by model: sex-specific fits, then the interaction model.
  setcolorder(sens_tbl, c("analysis", "variant", "cohort", "spec", "biotype", "module",
                          "n", "events", "n_female", "events_female", "n_male", "events_male",
                          "r_pb", "p_wilcoxon", "strat_model",
                          "strat_HR_female", "strat_lo_female", "strat_hi_female",
                          "strat_p_female",
                          "strat_HR_male", "strat_lo_male", "strat_hi_male", "strat_p_male",
                          "strat_HR_ratio_male_over_female", "strat_ratio_lo",
                          "strat_ratio_hi", "strat_p_sex_difference",
                          "strat_p_sex_difference_LRT", "int_model",
                          "int_HR_female", "int_lo_female", "int_hi_female", "int_p_female",
                          "int_HR_male", "int_lo_male", "int_hi_male", "int_p_male",
                          "int_HR_interaction", "int_lo_interaction", "int_hi_interaction",
                          "int_p_interaction", "int_p_interaction_LRT"))
  save_tsv(sens_tbl, "33_sensitivity_analyses.tsv")
  print(sens_tbl[, .(analysis, variant, biotype, module, n, events, r_pb, p_wilcoxon,
                     strat_HR_female, strat_p_female, strat_HR_male, strat_p_male,
                     strat_p_sex_difference, int_HR_interaction, int_p_interaction)],
        row.names = FALSE)
}

# ---- summary ----
banner("Summary")
msg("Common-covariate interaction, BH over the five prognostic modules, below ", FDR_ALPHA,
    ": ", sum(int_tbl[family == "prognostic_modules", fdr_interaction < FDR_ALPHA]), " of ",
    nrow(int_tbl[family == "prognostic_modules"]), " cohort x specification x module tests")
msg("Sex-specific fits, their OWN difference test (Wald contrast), BH over the same family, below ",
    FDR_ALPHA, ": ",
    sum(con_tbl[family == "prognostic_modules", strat_fdr_sex_difference < FDR_ALPHA]),
    " of ", nrow(con_tbl[family == "prognostic_modules"]),
    " (smallest raw p ", min(con_tbl[family == "prognostic_modules", strat_p_sex_difference]),
    ", smallest BH ",
    min(con_tbl[family == "prognostic_modules", strat_fdr_sex_difference]), ")")
msg("Models with fewer than ", MIN_EPV, " events per parameter: ",
    sum(int_tbl$epv_below_min), " interaction, ", sum(str_tbl$epv_below_min),
    " stratified (of ", nrow(int_tbl), " and ", nrow(str_tbl), ")")
for (s in SEX_SENTENCE) msg(s)
gy <- int_tbl[biotype == "lncRNA" & module == "greenyellow"]
msg("lncRNA greenyellow: interaction HR below 1 in ", sum(gy$HR_interaction < 1), " of ",
    nrow(gy), " rows and in all ", uniqueN(gy[HR_interaction < 1, cohort]),
    " cohorts (", paste0(gy$cohort, " ", gy$spec, " HR_int ", gy$HR_interaction,
                         ", p ", gy$p_interaction, collapse = "; "),
    "); smallest interaction p in the whole table is ",
    min(int_tbl$p_interaction), " and no test survives BH")
gyk <- int_tbl[cohort == "TCGA-KIRC" & spec == "principal" & module == "greenyellow"]
gys <- str_tbl[cohort == "TCGA-KIRC" & spec == "principal" & module == "greenyellow"]
gyc <- con_tbl[cohort == "TCGA-KIRC" & spec == "principal" & module == "greenyellow"]
# lncRNA greenyellow in TCGA-KIRC: each HR with the p-value of its own model.
msg("Headline (lncRNA greenyellow, TCGA-KIRC principal). SEX-SPECIFIC FITS, per SD of the ",
    "score over the ", nrow(FR$lncRNA$base), "-sample discovery set: ",
    paste0(gys$sex, " HR ", gys$HR, " (", gys$lo, "-", gys$hi, "), within-stratum p ", gys$p,
           collapse = "; "),
    "; difference between them, from those same fits: HR ratio male/female ",
    gyc$strat_HR_ratio_male_over_female, " (", gyc$strat_ratio_lo, "-", gyc$strat_ratio_hi,
    "), Wald p ", gyc$strat_p_sex_difference, ", likelihood-ratio p ",
    gyc$strat_p_sex_difference_LRT, ", BH ", gyc$strat_fdr_sex_difference)
msg("  per within-stratum SD (SD ", paste(gys$sd_score_in_stratum, collapse = " vs "), "): ",
    paste0(gys$sex, " ", gys$HR_per_stratum_SD, " (", gys$lo_per_stratum_SD, "-",
           gys$hi_per_stratum_SD, ")", collapse = ", "))
msg("  COMMON-COVARIATE INTERACTION MODEL, a different model and a different pair: female ",
    gyk$HR_female_interaction_model, " (p ", gyk$p_female_interaction_model, "), male ",
    gyk$HR_male_interaction_model, " (p ", gyk$p_male_interaction_model,
    "), interaction HR ", gyk$HR_interaction, " (", gyk$lo_interaction, "-",
    gyk$hi_interaction, "), p ", gyk$p_interaction, ", BH ", gyk$fdr_interaction)
msg("  The two tests of the same hypothesis disagree in significance here (",
    gyc$strat_p_sex_difference, " against ", gyc$int_p_interaction,
    "): the sex-specific fits let every covariate differ by sex, the interaction model ",
    "holds them common. Neither survives BH over the five prognostic modules.")
msg("Elapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
write_session_info("33_sex_stratified")
banner("33 | done")
