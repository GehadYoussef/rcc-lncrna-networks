# 15_sensitivity_tissue_source_site.R: tissue source site as a confounder
#
# Tissue source site is TCGA barcode field 2. On the discovery cohort:
#   1. variance in each quality metric, the lncRNA axis, the ESTIMATE scores and
#      the clinical variables explained by site (ANOVA R2, Kruskal-Wallis test)
#   2. the axis hazard ratio with and without site strata, under the reduced and
#      full covariate sets, with and without the STAR metrics
#   3. whether site is prognostic after the clinical terms (likelihood ratio)
# Inputs: 01_patient_level_data.tsv, cache/lnc_global_axis.rds and
# cache/estimate_scores.rds (07). Outputs: 15_tss_variance.tsv,
# 15_site_models.tsv (also saved as 20_site_models.tsv), 15_site_lrt.tsv,
# 15_tss_merged.tsv.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(survival) })
banner("15 | Sensitivity: tissue source site")

MIN_SITE_N <- 10

# Site factor for one analysis set: sites with fewer than MIN_SITE_N patients
# in that set are pooled into "other".
pool_sites <- function(tss) {
  sz <- table(tss)
  factor(ifelse(tss %in% names(sz)[sz >= MIN_SITE_N], tss, "other"))
}

pt   <- fread(file.path(RESULTS_DIR, "01_patient_level_data.tsv"))
axis <- as.data.table(readRDS(file.path(CACHE_DIR, "lnc_global_axis.rds")))
est  <- as.data.table(readRDS(file.path(CACHE_DIR, "estimate_scores.rds")))

d <- merge(merge(pt, axis[, .(sample_barcode, lnc_axis)], by = "sample_barcode"),
           est[, .(sample_barcode, StromalScore, ImmuneScore)], by = "sample_barcode")
if (!"tss" %in% names(d)) d[, tss := tstrsplit(sample_barcode, "-", keep = 2)[[1]]]
msg("Patients with axis, ESTIMATE and STAR metrics: ", nrow(d),
    "; distinct sites: ", uniqueN(d$tss))

sz <- d[, .N, by = tss][order(-N)]
msg("Site sizes: max ", max(sz$N), ", median ", median(sz$N),
    ", sites with n >= ", MIN_SITE_N, ": ", sum(sz$N >= MIN_SITE_N))
# Pooling on the merged set, for 15_tss_merged.tsv only.
d[, site := as.character(pool_sites(tss))]

# ---- 1. Variance explained by site ----
# Each variable on its own complete cases, with site pooling recomputed there.
site_var <- function(v) {
  ok <- !is.na(d[[v]]) & !is.na(d$tss)
  y  <- d[[v]][ok]; g <- pool_sites(d$tss[ok])
  a  <- anova(lm(y ~ g))
  data.table(variable = v, n = sum(ok), n_sites = nlevels(g),
             R2_anova = a[1, "Sum Sq"] / sum(a[, "Sum Sq"]),
             p_kruskal = kruskal.test(y ~ g)$p.value)
}
VARS <- c("pct_noFeature", "lnc_axis", "pct_multimapping", "StromalScore",
          "ImmuneScore", "age", "T_stage", "grade_num")
tss_var <- rbindlist(lapply(VARS, site_var))
tss_var[, `:=`(R2_anova = round(R2_anova, 4), p_kruskal = signif(p_kruskal, 3))]
print(tss_var, row.names = FALSE)
save_tsv(tss_var, "15_tss_variance.tsv")

# ---- 2. Axis models with and without site strata ----
# One complete-case set for all six models. Covariates are scaled per SD over
# the merged cohort before the complete-case restriction, as in 07.
d[, `:=`(axis_z  = as.numeric(scale(lnc_axis)),
         male    = as.numeric(sex == "male"),
         stromal = as.numeric(scale(StromalScore)),
         immune  = as.numeric(scale(ImmuneScore)),
         nf      = as.numeric(scale(pct_noFeature)),
         mm      = as.numeric(scale(pct_multimapping)),
         dep     = as.numeric(scale(log10(assigned_reads))))]
need <- c("os_time", "os_event", "axis_z", "age", "male", "T_stage", "N_pos", "M1",
          "grade_num", "stromal", "immune", "nf", "mm", "dep", "tss")
dc <- d[complete.cases(d[, ..need])]
dc[, site := as.character(pool_sites(tss))]
msg("Complete cases for the site models: ", nrow(dc), " patients, ",
    sum(dc$os_event), " events; site levels after pooling: ", uniqueN(dc$site))
y <- Surv(dc$os_time, dc$os_event)

mods <- list(
  `parsimonious: age sex T M1 grade`            = "axis_z + age + male + T_stage + M1 + grade_num",
  `parsimonious + site stratum`                 = "axis_z + age + male + T_stage + M1 + grade_num + strata(site)",
  `full principal set (no quality metrics)`     = "axis_z + age + male + T_stage + N_pos + M1 + grade_num + stromal + immune",
  `full principal set + site stratum`           = "axis_z + age + male + T_stage + N_pos + M1 + grade_num + stromal + immune + strata(site)",
  `full principal set + quality metrics`        = "axis_z + age + male + T_stage + N_pos + M1 + grade_num + stromal + immune + nf + mm + dep",
  `full principal set + quality + site stratum` = "axis_z + age + male + T_stage + N_pos + M1 + grade_num + stromal + immune + nf + mm + dep + strata(site)")
site_tbl <- rbindlist(lapply(names(mods), function(nm) {
  f <- coxph(as.formula(paste("y ~", mods[[nm]])), data = dc)
  s <- summary(f)
  zp <- tryCatch(cox.zph(f)$table["axis_z", "p"], error = function(e) NA_real_)
  data.table(model = nm,
             HR = round(s$conf.int["axis_z", 1], 3), lo = round(s$conf.int["axis_z", 3], 3),
             hi = round(s$conf.int["axis_z", 4], 3), p = signif(s$coefficients["axis_z", 5], 3),
             n = s$n, events = s$nevent, C = round(s$concordance[1], 3),
             ph_p = signif(zp, 3), site_stratified = grepl("strata", mods[[nm]]))
}))
print(site_tbl, row.names = FALSE)
save_tsv(site_tbl, "15_site_models.tsv")
save_tsv(site_tbl, "20_site_models.tsv")

# ---- 3. Is site itself prognostic? Likelihood-ratio test of site as a factor ----
lrt_base <- list(
  `age sex T M1 grade`                     = "age + male + T_stage + M1 + grade_num",
  `age sex T N M1 grade + ESTIMATE`        = "age + male + T_stage + N_pos + M1 + grade_num + stromal + immune",
  `age sex T N M1 grade + ESTIMATE + STAR` = "age + male + T_stage + N_pos + M1 + grade_num + stromal + immune + nf + mm + dep")
site_lrt <- rbindlist(lapply(names(lrt_base), function(nm) {
  f0 <- coxph(as.formula(paste("y ~", lrt_base[[nm]])), data = dc)
  f1 <- coxph(as.formula(paste("y ~", lrt_base[[nm]], "+ factor(site)")), data = dc)
  a  <- anova(f0, f1)
  data.table(base_model = nm, n = f1$n, events = f1$nevent,
             df = a$Df[2], chisq = round(a$Chisq[2], 3), p_lrt = signif(a$`Pr(>|Chi|)`[2], 3),
             n_site_levels = uniqueN(dc$site))
}))
print(site_lrt, row.names = FALSE)
save_tsv(site_lrt, "15_site_lrt.tsv")

save_tsv(d[, .(sample_barcode, tss, site, pct_noFeature, pct_multimapping,
               assigned_reads, lnc_axis, StromalScore, ImmuneScore)],
         "15_tss_merged.tsv")

write_session_info("15_sensitivity_tissue_source_site")
banner("15 | done")
