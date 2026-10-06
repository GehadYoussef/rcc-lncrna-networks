# 45_cptac_handling_time.R: does specimen handling time track the non-feature
# fraction in CPTAC-3?
#
# The handling-time reading of the metric (Discussion) holds that a longer time
# from devascularisation to freezing raises the non-feature fraction and, in
# larger and more complex nephrectomies, also marks higher-risk tumours.
# TCGA-KIRC records no handling time. CPTAC-3 records, per sample, the time from
# excision to freezing (all validation patients) and from clamping to freezing
# (a minority) in the GDC sample records. Its protocol caps excision to freezing
# at 30 minutes, so the range is narrow and a null result is weak evidence.
#
# Expected directions, recorded before the first run:
#   handling time vs non-feature fraction: positive (handling-time reading)
#   handling time vs projected lncRNA axis: positive
#   handling time vs T category and stage: positive (size and complexity)
#   handling time vs overall survival: hazard ratio above 1
# Tests: Spearman with 2,000-resample bootstrap intervals; Benjamini-Hochberg
# across the correlation tests of each sample set. Survival: Cox per SD of
# handling time, alone and with age, sex, T category, M1 and ordinal grade.
# Post hoc analysis, specified after the main analyses were complete.
# Inputs: results/08_validation_patient_level_data.tsv,
#   results/25_paired_scores_per_sample.tsv, results/25_normal_library_metrics.tsv,
#   GDC API sample records (cached in data/derived/cache).
# Outputs: results/45_cptac_handling_time_per_sample.tsv,
#   45_cptac_handling_time_correlations.tsv, 45_cptac_handling_time_survival.tsv.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(httr); library(survival)
})
banner("45 | CPTAC-3 specimen handling time")

gdc <- function(endpoint, body) {
  r <- httr::POST(paste0("https://api.gdc.cancer.gov/", endpoint),
                  body = jsonlite::toJSON(body, auto_unbox = TRUE),
                  httr::content_type_json(), httr::accept_json(), httr::timeout(300))
  httr::stop_for_status(r)
  jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"),
                     simplifyDataFrame = TRUE)$data$hits
}

# ---- 1. handling times per GDC sample ----------------------------------------
ht_rds <- file.path(CACHE_DIR, "45_cptac_handling_time.rds")
if (file.exists(ht_rds)) {
  ht <- readRDS(ht_rds); msg("Handling times from cache")
} else {
  filters <- list(op = "and", content = list(
    list(op = "in", content = list(field = "project.project_id", value = list("CPTAC-3"))),
    list(op = "in", content = list(field = "primary_site", value = list("Kidney")))))
  hits <- gdc("cases", list(filters = filters, format = "JSON", size = "5000",
                            fields = paste("submitter_id", "samples.submitter_id",
                                           "samples.sample_type",
                                           "samples.time_between_excision_and_freezing",
                                           "samples.time_between_clamping_and_freezing",
                                           sep = ",")))
  num <- function(x) if (is.null(x)) NA_real_ else suppressWarnings(as.numeric(x))
  ht <- rbindlist(lapply(seq_len(nrow(hits)), function(i) {
    s <- hits$samples[[i]]
    if (is.null(s) || !nrow(s)) return(NULL)
    data.table(patient = hits$submitter_id[i], sample_barcode = s$submitter_id,
               sample_type = s$sample_type,
               excision_to_freezing_min = num(s$time_between_excision_and_freezing),
               clamping_to_freezing_min = num(s$time_between_clamping_and_freezing))
  }), fill = TRUE)
  saveRDS(ht, ht_rds)
}
msg("GDC kidney samples: ", nrow(ht))

# ---- 2. join to the sequenced libraries ---------------------------------------
val <- fread(file.path(RESULTS_DIR, "08_validation_patient_level_data.tsv"))
ps  <- fread(file.path(RESULTS_DIR, "25_paired_scores_per_sample.tsv"))[cohort == "CPTAC-3"]
nm  <- fread(file.path(RESULTS_DIR, "25_normal_library_metrics.tsv"))[cohort == "CPTAC-3"]

tum <- merge(val, ps[tissue == "tumour", .(sample_barcode, lnc_axis_proj)],
             by = "sample_barcode", all.x = TRUE)
tum <- merge(tum, ht[, .(sample_barcode, excision_to_freezing_min, clamping_to_freezing_min)],
             by = "sample_barcode", all.x = TRUE)
tum[, `:=`(tissue = "tumour", log_noFeature = log10(pct_noFeature),
           log_depth = log10(assigned_reads))]
stopifnot(nrow(tum) == nrow(val), !anyDuplicated(tum$sample_barcode))

nor <- nm[!(failed_reads %in% TRUE), .(patient, sample_barcode, pct_noFeature,
                                       pct_multimapping, assigned_reads)]
nor <- merge(nor, ps[tissue == "normal", .(sample_barcode, lnc_axis_proj)],
             by = "sample_barcode", all.x = TRUE)
nor <- merge(nor, ht[, .(sample_barcode, excision_to_freezing_min, clamping_to_freezing_min)],
             by = "sample_barcode", all.x = TRUE)
nor[, `:=`(tissue = "normal", log_noFeature = log10(pct_noFeature),
           log_depth = log10(assigned_reads))]

per <- rbind(tum[, .(tissue, patient, sample_barcode, excision_to_freezing_min,
                     clamping_to_freezing_min, pct_noFeature, pct_multimapping,
                     assigned_reads, lnc_axis_proj, T_stage, stage_num, grade_num,
                     os_time, os_event)],
             nor[, .(tissue, patient, sample_barcode, excision_to_freezing_min,
                     clamping_to_freezing_min, pct_noFeature, pct_multimapping,
                     assigned_reads, lnc_axis_proj)], fill = TRUE)
save_tsv(per, "45_cptac_handling_time_per_sample.tsv")
msg("Tumours with excision time: ", tum[is.finite(excision_to_freezing_min), .N], " of ", nrow(tum),
    "; with clamping time: ", tum[is.finite(clamping_to_freezing_min), .N])
msg("Normals with excision time: ", nor[is.finite(excision_to_freezing_min), .N], " of ", nrow(nor))
print(summary(tum$excision_to_freezing_min))

# ---- 3. correlations -----------------------------------------------------------
set.seed(SEED)
sp_boot <- function(x, y, B = 2000L) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 10) return(list(n = n, rho = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  bs <- replicate(B, { i <- sample.int(n, n, TRUE); cor(x[i], y[i], method = "spearman") })
  list(n = n, rho = unname(ct$estimate), lo = quantile(bs, 0.025, na.rm = TRUE),
       hi = quantile(bs, 0.975, na.rm = TRUE), p = ct$p.value)
}
EXPECT <- c(pct_noFeature = "positive", lnc_axis_proj = "positive",
            pct_multimapping = "none stated", log_depth = "none stated",
            T_stage = "positive", stage_num = "positive", grade_num = "none stated")
run_set <- function(d, set, time_var, targets) {
  out <- rbindlist(lapply(targets, function(v) {
    r <- sp_boot(d[[time_var]], d[[v]])
    data.table(sample_set = set, handling_time = time_var, variable = v,
               expected = EXPECT[[v]], n = r$n, spearman_rho = r$rho,
               boot_lo = r$lo, boot_hi = r$hi, p = r$p)
  }))
  out[, fdr_within_set := p.adjust(p, "BH")]
  out
}
tvars <- c("pct_noFeature", "lnc_axis_proj", "pct_multimapping", "log_depth",
           "T_stage", "stage_num", "grade_num")
cors <- rbind(
  run_set(tum, "CPTAC-3 tumours", "excision_to_freezing_min", tvars),
  run_set(tum, "CPTAC-3 tumours", "clamping_to_freezing_min", tvars),
  run_set(nor, "CPTAC-3 normal tissue", "excision_to_freezing_min",
          c("pct_noFeature", "lnc_axis_proj", "pct_multimapping", "log_depth")))
cors[, direction_as_expected := fifelse(expected == "positive",
                                        fifelse(spearman_rho > 0, "yes", "no"), "--")]
save_tsv(cors, "45_cptac_handling_time_correlations.tsv"); print(cors)

# ---- 4. handling time and overall survival -------------------------------------
tum[, ht_z := (excision_to_freezing_min - mean(excision_to_freezing_min, na.rm = TRUE)) /
              sd(excision_to_freezing_min, na.rm = TRUE)]
tum[, male := as.integer(sex == "male")]
fit_one <- function(fml, label) {
  d <- tum[complete.cases(tum[, all.vars(fml), with = FALSE])]
  f <- coxph(fml, data = d)
  s <- summary(f)$coefficients["ht_z", , drop = FALSE]; ci <- summary(f)$conf.int["ht_z", , drop = FALSE]
  data.table(model = label, n = nrow(d), events = sum(d$os_event),
             HR = ci[, "exp(coef)"], lo = ci[, "lower .95"], hi = ci[, "upper .95"],
             p = s[, "Pr(>|z|)"], epv = sum(d$os_event) / length(coef(f)))
}
surv <- rbind(
  fit_one(Surv(os_time, os_event) ~ ht_z, "excision-to-freezing time alone"),
  fit_one(Surv(os_time, os_event) ~ ht_z + age + male + T_stage + M1 + grade_num,
          "+ age, sex, T category, M1, ordinal grade"))
save_tsv(surv, "45_cptac_handling_time_survival.tsv"); print(surv)

write_session_info("45_cptac_handling_time")
msg("45 done")
