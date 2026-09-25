# 30_endpoint_sensitivity.R: disease-specific and progression-free endpoints from the TCGA Clinical Data Resource
# Repeats the main survival analyses on the TCGA-CDR endpoints (Liu et al. 2018, Cell):
# OS, DSS (approximate for KIRC), PFI and DFI (few events, reported for coverage only),
# censored at OS_CENSOR_DAYS, with the pipeline OS as the reference row. Discovery
# cohort only (CPTAC-3 has no DSS or PFI).
# Sections: 1 CDR fetch and OS concordance, 2 prognostic modules, 3 non-feature fraction,
#   4 lncRNA axis nested models, 5 locked models re-scored (apparent), 6 repeated
#   10 x 10-fold CV increment, 7 fold-wise network rebuild for DSS and OS.
# Inputs: caches (dataset, networks, ESTIMATE, STAR QC, lncRNA axis, locked model).
# Outputs: 30_*.tsv, caches 30_TCGA-CDR.xlsx and 30_foldwise_rebuild.rds.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet)
  library(httr); library(readxl); library(parallel); library(WGCNA)
})
banner("30 | Endpoint sensitivity: TCGA-CDR overall, disease-specific and progression-free")
set.seed(SEED)
t_start <- Sys.time()

# Bootstrap seeds swept in section 6, and fold assignments for the section 7 rebuild
# (each assignment rebuilds FOLDWISE_K x 2 networks).
BOOT_SEED_SWEEP    <- 200L
FOLDWISE30_N_ASSIGN <- 5L

ds     <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets   <- readRDS(file.path(CACHE_DIR, "networks.rds"))
est    <- as.data.table(readRDS(file.path(CACHE_DIR, "estimate_scores.rds")))
QC     <- as.data.table(readRDS(file.path(CACHE_DIR, "star_qc.rds")))
axis   <- as.data.table(readRDS(file.path(CACHE_DIR, "lnc_global_axis.rds")))
L      <- readRDS(file.path(CACHE_DIR, "locked_model.rds"))
stopifnot(identical(L$version, "v9"))
cohort <- as.data.table(ds$cohort_full)
stopifnot(all(c("patient", "sample_barcode", "os_time", "os_event", "age", "sex", "T_stage",
                "N_pos", "M1", "grade_num", "pct_noFeature", "pct_multimapping", "libsize")
              %in% names(cohort)))
msg("Discovery cohort: ", nrow(cohort), " patients, ", sum(cohort$os_event), " deaths (pipeline OS)")

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
ph_term <- function(fit, term) tryCatch(cox.zph(fit)$table[term, "p"], error = function(e) NA_real_)
# One Cox fit on the complete cases of its own covariates, summarising the `term`
# coefficient (per SD of the term as supplied).
cox_row <- function(d, term, terms, time = "time", event = "event") {
  vars <- c(time, event, term, terms)
  dc <- d[complete.cases(d[, ..vars])]
  f  <- coxph(as.formula(paste0("Surv(", time, ", ", event, ") ~ ",
                                paste(c(term, terms), collapse = " + "))), data = dc)
  s  <- summary(f)
  list(fit = f, data = dc,
       row = data.table(n = s$n, events = s$nevent,
                        HR = round(s$conf.int[term, "exp(coef)"], 3),
                        lo = round(s$conf.int[term, "lower .95"], 3),
                        hi = round(s$conf.int[term, "upper .95"], 3),
                        p  = signif(s$coefficients[term, "Pr(>|z|)"], 3),
                        C  = round(unname(s$concordance["C"]), 3),
                        ph_p = signif(ph_term(f, term), 3)))
}

# ---- 1. TCGA Clinical Data Resource: fetch, restrict to KIRC, merge, compare ----
banner("1 | TCGA-CDR: fetch, coverage and overall-survival concordance")
CDR_URL  <- "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
CDR_FILE <- file.path(CACHE_DIR, "30_TCGA-CDR.xlsx")
if (!file.exists(CDR_FILE) || file.size(CDR_FILE) < 1e6) {
  msg("Downloading TCGA-CDR-SupplementalTableS1.xlsx from ", CDR_URL)
  got <- tryCatch({
    r <- httr::GET(CDR_URL, httr::write_disk(CDR_FILE, overwrite = TRUE), httr::timeout(600))
    httr::stop_for_status(r)
    msg("  HTTP ", httr::status_code(r), "; ",
        httr::headers(r)[["content-disposition"]], "; ", file.size(CDR_FILE), " bytes")
    TRUE
  }, error = function(e) e)
  if (inherits(got, "error")) {
    unlink(CDR_FILE)
    stop("TCGA-CDR download failed: ", conditionMessage(got),
         ". Nothing in this script can run without the workbook; retry when ",
         "https://gdc.cancer.gov/about-data/publications/pancanatlas is reachable.")
  }
} else msg("TCGA-CDR workbook from cache: ", CDR_FILE, " (", file.size(CDR_FILE), " bytes)")

sheets <- readxl::excel_sheets(CDR_FILE)
stopifnot("TCGA-CDR" %in% sheets)
cdr_all <- as.data.table(readxl::read_excel(CDR_FILE, sheet = "TCGA-CDR", guess_max = 20000))
cdr <- cdr_all[type == "KIRC"]
msg("TCGA-CDR: ", nrow(cdr_all), " patients across ", uniqueN(cdr_all$type),
    " cancer types; KIRC rows: ", nrow(cdr))
for (v in c("OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time",
            "last_contact_days_to", "death_days_to", "new_tumor_event_dx_days_to"))
  set(cdr, j = v, value = num(cdr[[v]]))

# The CDR's own event tallies over all KIRC patients in the resource.
cdr_tally <- function(e) c(events = sum(cdr[[e]] == 1, na.rm = TRUE),
                           censored = sum(cdr[[e]] == 0, na.rm = TRUE),
                           not_available = sum(is.na(cdr[[e]])))
# Per-endpoint recommendations for KIRC, transcribed from Liu et al. 2018,
# Table 3 (check = recommended, check* = with caution, app. = approximate).
CDR_RECOMMENDATION <- c(
  OS_cdr = "recommended (check); Table 3 KIRC: 177 events, 360 censored",
  PFI    = "recommended (check); Table 3 KIRC: 162 events, 375 censored; the CDR notes sheet prefers PFI over OS given TCGA's short follow-up",
  DFI    = "recommended with caution (check*), explanation 'number of events is small'; Table 3 KIRC: 15 events, 102 censored",
  DSS    = "recommended, approximate (check, app.); Table 3 KIRC: 110 events, 415 censored; the notes sheet states DSS is approximated for all types except CESC, PAAD and UVM")
CDR_GENERAL_NOTE <- "TCGA-CDR_Notes sheet: 'For clinical outcome endpoints, we recommend the use of PFI for progression-free interval, and OS for overall survival. Both endpoints are relatively accurate. Given the relatively short follow-up time, PFI is preferred over OS.'"

# ---- merge with the cohort ---------------------------------------------------
m <- merge(cohort[, .(patient, sample_barcode, os_time, os_event, vital_pipeline = vital)],
           cdr[, .(patient = bcr_patient_barcode, in_cdr = TRUE, vital_status, tumor_status,
                   cause_of_death, last_contact_days_to, death_days_to,
                   new_tumor_event_type, new_tumor_event_dx_days_to,
                   OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time, Redaction)],
           by = "patient", all.x = TRUE)
m[is.na(in_cdr), in_cdr := FALSE]
stopifnot(nrow(m) == nrow(cohort), !anyDuplicated(m$patient))
msg("Cohort patients found in the CDR: ", sum(m$in_cdr), " of ", nrow(m))

# ---- redaction guard ---------------------------------------------------------
# Patients whose records the CDR marks as redacted are excluded and counted.
n_redacted_resource <- sum(!is.na(cdr_all$Redaction))
n_redacted_kirc     <- sum(!is.na(cdr$Redaction))
redacted_patients   <- m$patient[!is.na(m$Redaction)]
n_cohort_before_redaction <- nrow(m)
m <- m[is.na(Redaction)]
msg("CDR redaction flag: ", n_redacted_resource, " patients across the resource, ",
    n_redacted_kirc, " in TCGA-KIRC, ", length(redacted_patients),
    " in this cohort; excluded here: ", n_cohort_before_redaction - nrow(m))
stopifnot(nrow(m) == n_cohort_before_redaction - length(redacted_patients),
          !any(!is.na(m$Redaction)))
n_analysed_cohort <- nrow(m)

# ---- endpoint table, one row per patient x endpoint --------------------------
# Censored at OS_CENSOR_DAYS like the pipeline OS. Non-positive times cannot
# enter a Cox model and are counted as unusable.
ENDPOINT_LABELS <- c(
  OS_pipeline = "overall survival, pipeline (01)",
  OS_cdr      = "overall survival, TCGA-CDR",
  DSS         = "disease-specific survival, TCGA-CDR (approximate)",
  PFI         = "progression-free interval, TCGA-CDR",
  DFI         = "disease-free interval, TCGA-CDR")
ANALYSED <- c("OS_pipeline", "OS_cdr", "DSS", "PFI")
mk_ep <- function(name, time, event) {
  time <- num(time); event <- num(event)
  ok <- is.finite(time) & is.finite(event)
  data.table(patient = m$patient, sample_barcode = m$sample_barcode, endpoint = name,
             time_raw = time, event_raw = event,
             time  = ifelse(ok, pmin(time, OS_CENSOR_DAYS), NA_real_),
             event = ifelse(ok, ifelse(time > OS_CENSOR_DAYS, 0, event), NA_real_),
             available = ok, usable = ok & time > 0,
             # pipeline OS status, to show whether missingness sits in decedents
             os_death = as.numeric(m$os_event))
}
EP <- rbindlist(list(mk_ep("OS_pipeline", m$os_time, m$os_event),
                     mk_ep("OS_cdr", m$OS.time, m$OS),
                     mk_ep("DSS", m$DSS.time, m$DSS),
                     mk_ep("PFI", m$PFI.time, m$PFI),
                     mk_ep("DFI", m$DFI.time, m$DFI)))
EP[, endpoint := factor(endpoint, levels = names(ENDPOINT_LABELS))]
ep_of <- function(name) EP[endpoint == name & usable == TRUE, .(sample_barcode, time, event)]

rev_km_median <- function(time, event) {
  f <- survfit(Surv(time, 1 - event) ~ 1)
  unname(summary(f)$table["median"])
}
coverage <- EP[, {
  u <- usable
  list(n_cohort = .N, n_in_cdr = sum(m$in_cdr),
    n_available = sum(available), n_time_nonpositive = sum(available & time_raw <= 0),
    n_analysable = sum(u), events_raw = sum(event_raw[u] == 1),
    n_censored_at_10y = sum(time_raw[u] > OS_CENSOR_DAYS),
    events_10y = sum(event[u] == 1),
    median_followup_days_10y = round(rev_km_median(time[u], event[u])),
    median_time_to_event_days = round(median(time_raw[u & event_raw == 1])),
    # ---- missingness structure (which patients this endpoint cannot score) ---
    n_missing_status          = sum(is.na(event_raw)),
    n_missing_time_with_status= sum(!is.na(event_raw) & is.na(time_raw)),
    n_missing_any             = sum(!u),
    n_missing_among_os_deaths = sum(!u & os_death == 1),
    n_missing_among_os_alive  = sum(!u & os_death == 0),
    n_events_lost             = sum(event_raw == 1 & !u, na.rm = TRUE))
}, by = endpoint]
coverage[, cdr_kirc_events_all537 := sapply(as.character(endpoint), function(e)
  if (e == "OS_pipeline") NA_integer_ else cdr_tally(sub("_cdr$", "", e))[["events"]])]
coverage[, cdr_kirc_censored_all537 := sapply(as.character(endpoint), function(e)
  if (e == "OS_pipeline") NA_integer_ else cdr_tally(sub("_cdr$", "", e))[["censored"]])]
coverage[, cdr_kirc_not_available_all537 := sapply(as.character(endpoint), function(e)
  if (e == "OS_pipeline") NA_integer_ else cdr_tally(sub("_cdr$", "", e))[["not_available"]])]
coverage[, description := ENDPOINT_LABELS[as.character(endpoint)]]
coverage[, cdr_recommendation_KIRC := ifelse(as.character(endpoint) == "OS_pipeline",
                                             "not applicable (pipeline endpoint)",
                                             CDR_RECOMMENDATION[as.character(endpoint)])]
coverage[, analysed_here := as.character(endpoint) %in% ANALYSED]
coverage[, epv_principal_12_parameters := round(events_10y / 12, 1)]

# ---- missingness structure, in words ----------------------------------------
# Which patients an endpoint cannot score, and whether they are decedents.
tumour_status_breakdown <- function(patients) {
  ts <- m$tumor_status[match(patients, m$patient)]
  ts[is.na(ts)] <- "not recorded"
  tb <- sort(table(ts), decreasing = TRUE)
  paste(sprintf("%s n=%d", names(tb), as.integer(tb)), collapse = "; ")
}
coverage[, missingness_note := vapply(as.character(endpoint), function(e) {
  mm <- EP[as.character(endpoint) == e & usable == FALSE]
  if (!nrow(mm)) return("every cohort patient is scoreable for this endpoint")
  sprintf(paste0("%d of %d patients unscoreable (%d with no status, %d with a status but ",
                 "no time, %d with a non-positive time); pipeline overall-survival status of ",
                 "the unscoreable: dead %d, alive %d; CDR tumour status of the ",
                 "unscoreable: %s; recorded events lost with them: %d"),
          nrow(mm), n_analysed_cohort, sum(is.na(mm$event_raw)),
          sum(!is.na(mm$event_raw) & is.na(mm$time_raw)),
          sum(mm$available & mm$time_raw <= 0),
          sum(mm$os_death == 1), sum(mm$os_death == 0),
          tumour_status_breakdown(mm$patient),
          sum(mm$event_raw == 1, na.rm = TRUE))
}, character(1))]
setcolorder(coverage, c("endpoint", "description"))
print(coverage[, .(endpoint, n_available, n_analysable, events_10y, n_censored_at_10y,
                   median_followup_days_10y, epv_principal_12_parameters, analysed_here)])

# ---- concordance of the CDR overall survival with the pipeline overall survival
cc <- m[in_cdr == TRUE & is.finite(OS.time) & is.finite(OS)]
cc[, `:=`(cdr_time_10y = pmin(OS.time, OS_CENSOR_DAYS),
          cdr_event_10y = ifelse(OS.time > OS_CENSOR_DAYS, 0, OS))]
cc[, `:=`(diff_days = os_time - cdr_time_10y,
          event_differs = os_event != cdr_event_10y)]
cc[, time_differs := abs(diff_days) > 0]
sp_raw <- cor.test(cc$os_time, cc$OS.time, method = "spearman", exact = FALSE)
sp_10y <- cor.test(cc$os_time, cc$cdr_time_10y, method = "spearman", exact = FALSE)
disc <- cc[event_differs | time_differs]
conc <- data.table(quantity = c(
  "n_compared", "n_events_pipeline", "n_events_cdr_raw", "n_events_cdr_10y",
  "event_agree_raw_n", "event_agree_raw_pct", "event_agree_10y_n", "event_agree_10y_pct",
  "event_differs_n", "spearman_rho_time_raw", "spearman_p_time_raw",
  "spearman_rho_time_10y", "spearman_p_time_10y",
  "time_identical_n", "time_differs_n", "time_differs_gt30d_n", "time_differs_gt365d_n",
  "time_differs_median_abs_days", "time_differs_max_abs_days",
  "time_differs_cdr_longer_n", "time_differs_pipeline_longer_n",
  "pipeline_censored_at_10y_n", "event_or_time_differs_n"),
  value = c(
  nrow(cc), sum(cc$os_event), sum(cc$OS), sum(cc$cdr_event_10y),
  sum(cc$os_event == cc$OS), round(100 * mean(cc$os_event == cc$OS), 2),
  sum(!cc$event_differs), round(100 * mean(!cc$event_differs), 2),
  sum(cc$event_differs), round(unname(sp_raw$estimate), 4), signif(sp_raw$p.value, 3),
  round(unname(sp_10y$estimate), 4), signif(sp_10y$p.value, 3),
  sum(!cc$time_differs), sum(cc$time_differs), sum(abs(cc$diff_days) > 30),
  sum(abs(cc$diff_days) > 365),
  if (any(cc$time_differs)) median(abs(cc$diff_days[cc$time_differs])) else 0,
  max(abs(cc$diff_days)),
  sum(cc$diff_days < 0), sum(cc$diff_days > 0),
  sum(cc$os_time == OS_CENSOR_DAYS), nrow(disc)))
conc[, note := ""]
conc[quantity == "n_compared", note := "cohort patients with a CDR OS record; both endpoints censored at OS_CENSOR_DAYS for the comparison"]
conc[quantity == "spearman_rho_time_raw", note := "pipeline OS time (censored at 10 y) vs CDR OS.time as supplied"]
conc[quantity == "pipeline_censored_at_10y_n", note := "patients whose pipeline time is the 10-year cap; they differ from the raw CDR time by construction"]
conc[quantity == "time_differs_cdr_longer_n", note := "the CDR follow-up files record a later last contact or death than the GDC API / BCR XML used in 01"]
print(conc)

# ---- how far the endpoints are from being independent tests -------------------
# The BH adjustments in sections 5 and 6 treat correlated endpoints as a family.
# Phi coefficients between event indicators, on patients both endpoints can
# score, are quoted with the adjusted p values.
ep_ind <- dcast(EP[usable == TRUE], sample_barcode ~ endpoint, value.var = "event")
phi_pair <- function(a, b) {
  k <- is.finite(ep_ind[[a]]) & is.finite(ep_ind[[b]])
  # stats::cor explicitly: WGCNA masks cor() and returns a matrix
  c(phi = as.numeric(stats::cor(ep_ind[[a]][k], ep_ind[[b]][k])), n = sum(k))
}
EP_PAIRS <- combn(ANALYSED, 2, simplify = FALSE)
ep_dep <- rbindlist(lapply(EP_PAIRS, function(p) {
  z <- phi_pair(p[1], p[2])
  data.table(pair = paste(p, collapse = " vs "), phi = round(z[["phi"]], 4),
             n_both_scoreable = as.integer(z[["n"]]))
}))
phi_of <- function(a, b) ep_dep[pair == paste(a, b, sep = " vs "), phi]
ENDPOINT_DEPENDENCE <- sprintf(
  paste0("the endpoints are not independent: the pipeline and CDR overall survival are the ",
         "same endpoint (Spearman %s on times, %s%% event agreement), and the event ",
         "indicators correlate at phi %.2f between overall and disease-specific survival, ",
         "%.2f between disease-specific survival and the progression-free interval and ",
         "%.2f between overall survival and the progression-free interval, so the ",
         "adjustment is conservative-crude rather than exact; the CDR endpoints were ",
         "prespecified as sensitivity analyses, not as independent primary tests"),
  conc[quantity == "spearman_rho_time_raw", value],
  conc[quantity == "event_agree_10y_pct", value],
  phi_of("OS_pipeline", "DSS"), phi_of("DSS", "PFI"), phi_of("OS_pipeline", "PFI"))
print(ep_dep)

cov_chr <- coverage[, lapply(.SD, as.character)]
cov_long <- rbindlist(list(
  melt(cov_chr, id.vars = "endpoint", variable.name = "quantity", value.name = "value",
       variable.factor = FALSE)[, .(section = "coverage", endpoint = endpoint,
                                    quantity, value, note = "")],
  data.table(section = "consort", endpoint = "all",
             quantity = c("cohort_patients_entering_30",
                          "cdr_redaction_flagged_whole_resource",
                          "cdr_redaction_flagged_TCGA_KIRC",
                          "cdr_redaction_flagged_in_this_cohort",
                          "excluded_cdr_redaction",
                          "analysed_after_redaction_exclusion"),
             value = as.character(c(n_cohort_before_redaction, n_redacted_resource,
                                    n_redacted_kirc, length(redacted_patients),
                                    n_cohort_before_redaction - n_analysed_cohort,
                                    n_analysed_cohort)),
             note = c("cohort_full of 01, after the 01 exclusions",
                      "patients the CDR marks Redacted anywhere in the resource",
                      "of those, patients in TCGA-KIRC",
                      "of those, patients in this analysis cohort",
                      "the redaction flag is acted on, not merely reported; no patient in this cohort carries it, so the exclusion removes nobody",
                      "every count in every 30_ file is on this set")),
  data.table(section = "missingness", endpoint = as.character(coverage$endpoint),
             quantity = "missingness_structure", value = coverage$missingness_note,
             note = "patients this endpoint cannot score, and whether they are decedents; see also n_missing_* and n_events_lost in the coverage section"),
  data.table(section = "endpoint_dependence", endpoint = ep_dep$pair,
             quantity = "phi_event_indicators", value = as.character(ep_dep$phi),
             note = paste0("phi (Pearson on 0/1 event indicators) among the ",
                           ep_dep$n_both_scoreable, " patients both endpoints can score")),
  data.table(section = "endpoint_dependence", endpoint = "all",
             quantity = "multiplicity_caveat", value = ENDPOINT_DEPENDENCE,
             note = "applies to the BH-adjusted p columns of 30_endpoint_cv_increment.tsv and 30_endpoint_locked_model.tsv"),
  data.table(section = "CDR_recommendation_KIRC", endpoint = c("OS_cdr", "PFI", "DFI", "DSS"),
             quantity = "Liu_2018_Table3", value = unname(CDR_RECOMMENDATION[c("OS_cdr", "PFI", "DFI", "DSS")]),
             note = "transcribed from Table 3 of Liu et al. 2018 Cell 173:400; counts are the CDR's own for all 537 KIRC patients"),
  data.table(section = "CDR_recommendation_general", endpoint = "all", quantity = "TCGA-CDR_Notes",
             value = CDR_GENERAL_NOTE, note = "sheet TCGA-CDR_Notes of the workbook"),
  data.table(section = "OS_concordance", endpoint = "OS_cdr vs OS_pipeline",
             quantity = conc$quantity, value = as.character(conc$value), note = conc$note)))
save_tsv(cov_long, "30_cdr_coverage_and_os_concordance.tsv")
msg("Missingness: ", coverage[endpoint == "DSS", missingness_note])
msg("Missingness: ", coverage[endpoint == "PFI", missingness_note])

disc_out <- disc[, .(patient, os_time_pipeline = os_time, os_event_pipeline = os_event,
                     vital_pipeline, OS_time_cdr_raw = OS.time, OS_cdr_raw = OS,
                     OS_time_cdr_10y = cdr_time_10y, OS_cdr_10y = cdr_event_10y,
                     diff_days_pipeline_minus_cdr = diff_days, event_differs, time_differs,
                     cdr_vital_status = vital_status, cdr_last_contact_days_to = last_contact_days_to,
                     cdr_death_days_to = death_days_to, cdr_tumor_status = tumor_status,
                     DSS, DSS.time, PFI, PFI.time)][order(-abs(diff_days_pipeline_minus_cdr))]
# A later data freeze can extend follow-up but cannot move a recorded death by a
# year. Patients dead in both sources with times more than a year apart are
# flagged as not explained by the freeze. They are recorded, not corrected.
disc_out[, dead_in_both := os_event_pipeline == 1 & OS_cdr_10y == 1]
disc_out[, explained_by_later_freeze := !(dead_in_both &
             abs(diff_days_pipeline_minus_cdr) > 365)]
disc_out[, flag := fifelse(
  !explained_by_later_freeze,
  paste0("NOT explained by a later freeze: death recorded in both sources with the times ",
         "more than a year apart; check whether 01 selected an earlier follow-up record"),
  fifelse(diff_days_pipeline_minus_cdr < 0,
          "benign: the CDR records later follow-up than the GDC API / BCR XML used in 01 (later freeze)",
          "benign: the pipeline records later follow-up than the CDR for this patient"))]
save_tsv(disc_out, "30_cdr_os_discrepant_patients.tsv")
msg("Patients whose OS event or (10-year-censored) time differs between the pipeline and the CDR: ",
    nrow(disc_out), "; not explained by a later freeze: ",
    sum(!disc_out$explained_by_later_freeze),
    if (any(!disc_out$explained_by_later_freeze))
      paste0(" (", paste(disc_out$patient[!disc_out$explained_by_later_freeze], collapse = ", "), ")")
    else "")

# Cross-tabulation of the CDR endpoints for the log (OS deaths that are DSS events,
# PFI events that are deaths).
ep_w <- dcast(EP[usable == TRUE], sample_barcode ~ endpoint, value.var = "event")
msg("Cohort: OS_cdr deaths ", sum(ep_w$OS_cdr == 1, na.rm = TRUE),
    "; of these DSS events ", sum(ep_w$OS_cdr == 1 & ep_w$DSS == 1, na.rm = TRUE),
    ", DSS censored ", sum(ep_w$OS_cdr == 1 & ep_w$DSS == 0, na.rm = TRUE),
    ", DSS not available ", sum(ep_w$OS_cdr == 1 & is.na(ep_w$DSS), na.rm = TRUE),
    "; PFI events ", sum(ep_w$PFI == 1, na.rm = TRUE),
    " of which without death ", sum(ep_w$PFI == 1 & ep_w$OS_cdr == 0, na.rm = TRUE))

# ---- 2. The five prognostic modules per endpoint ----
banner("2 | Prognostic modules under each endpoint")
SC <- discovery_scores(nets, discovery_loadings(nets))
MODULES <- data.table(biotype = c("protein_coding", "protein_coding", "lncRNA", "lncRNA", "lncRNA"),
                      network = c("mrna", "mrna", "lnc", "lnc", "lnc"),
                      module  = c("green", "purple", "blue", "greenyellow", "turquoise"))
MODULES[, column := paste0(ifelse(network == "mrna", "mRNA_ME", "lnc_ME"), module)]
stopifnot(all(MODULES[network == "mrna", column] %in% colnames(SC$mrna)),
          all(MODULES[network == "lnc",  column] %in% colnames(SC$lnc)))

# Model frames as in stage 03: covariates standardised over all network samples
# before complete-case restriction. Module score is the unit-variance eigengene.
frame_for <- function(S) {
  cl <- cohort[match(rownames(S), sample_barcode)]
  stopifnot(identical(cl$sample_barcode, rownames(S)))
  e  <- est[match(cl$sample_barcode, sample_barcode)]
  stopifnot(!anyNA(e$StromalScore))
  data.table(sample_barcode = cl$sample_barcode,
             age = cl$age, sex = factor(cl$sex, levels = c("female", "male")),
             T_stage = as.numeric(cl$T_stage), N_pos = as.numeric(cl$N_pos),
             M1 = as.numeric(cl$M1), grade_o = as.numeric(cl$grade_num),
             stromal  = as.numeric(scale(e$StromalScore)),
             immune   = as.numeric(scale(e$ImmuneScore)),
             noFeat   = as.numeric(scale(cl$pct_noFeature)),
             multimap = as.numeric(scale(cl$pct_multimapping)),
             depth    = as.numeric(scale(log10(cl$libsize))))
}
FR <- list(mrna = list(scores = SC$mrna, base = frame_for(SC$mrna)),
           lnc  = list(scores = SC$lnc,  base = frame_for(SC$lnc)))
SPECS <- list(
  unadjusted = character(0),
  clinical   = c("age", "sex", "T_stage", "N_pos", "M1", "grade_o"),
  principal  = c("age", "sex", "T_stage", "N_pos", "M1", "grade_o",
                 "stromal", "immune", "noFeat", "multimap", "depth"))
SPEC_LABEL <- c(unadjusted = "module alone",
                clinical   = "module + age, sex, T, N, M1, ordinal grade",
                principal  = "module + clinical + ESTIMATE + non-feature, multimapping, log10 depth")

mod_tbl <- rbindlist(lapply(ANALYSED, function(ep) {
  y <- ep_of(ep)
  rbindlist(lapply(seq_len(nrow(MODULES)), function(i) {
    fr <- FR[[MODULES$network[i]]]
    d  <- merge(copy(fr$base)[, ME := as.numeric(fr$scores[, MODULES$column[i]])],
                y, by = "sample_barcode")
    rbindlist(lapply(names(SPECS), function(sp) {
      r <- cox_row(d, "ME", SPECS[[sp]])$row
      n_par <- length(SPECS[[sp]]) + 1L
      cbind(data.table(endpoint = ep, biotype = MODULES$biotype[i], module = MODULES$module[i],
                       specification = sp, specification_terms = SPEC_LABEL[[sp]]),
            r, data.table(n_parameters = n_par, epv = round(r$events / n_par, 1),
                          adequate_epv = r$events / n_par >= MIN_EPV))
    }))
  }))
}))
mod_tbl[, fdr_5modules := signif(p.adjust(p, "BH"), 3), by = .(endpoint, specification)]
setcolorder(mod_tbl, c("endpoint", "biotype", "module", "specification", "n", "events",
                       "HR", "lo", "hi", "p", "fdr_5modules", "C", "ph_p"))
save_tsv(mod_tbl, "30_endpoint_modules.tsv")
print(dcast(mod_tbl[specification == "principal"], biotype + module ~ endpoint,
            value.var = "HR"))
print(mod_tbl[specification == "principal",
              .(endpoint, biotype, module, n, events, HR, lo, hi, p, fdr_5modules, epv)])

# ---- 3. The non-feature fraction as the exposure ----
banner("3 | Non-feature fraction as exposure, per endpoint")
# Standardised over the whole cohort before complete-case restriction, as in
# stage 07. Both models use one complete-case set per endpoint.
cov <- merge(cohort[, .(sample_barcode, age, sex, T_stage, N_pos, M1, grade_num,
                        pct_noFeature, pct_multimapping, libsize)],
             est[, .(sample_barcode, StromalScore, ImmuneScore)], by = "sample_barcode")
cov[, `:=`(male    = as.numeric(sex == "male"),
           stromal = as.numeric(scale(StromalScore)),
           immune  = as.numeric(scale(ImmuneScore)),
           nf      = as.numeric(scale(pct_noFeature)),
           nf_log  = as.numeric(scale(log10(pct_noFeature))),
           mm      = as.numeric(scale(pct_multimapping)),
           dep     = as.numeric(scale(log10(libsize))))]
CLIN_RHS <- c("age", "male", "T_stage", "N_pos", "M1", "grade_num")
NF_MODELS <- list(`(a) pct_noFeature alone` = character(0),
                  `(b) + age, sex, T, N, M1, ordinal grade` = CLIN_RHS)
exp_tbl <- rbindlist(lapply(ANALYSED, function(ep) {
  d  <- merge(cov, ep_of(ep), by = "sample_barcode")
  need <- c("time", "event", "nf", "nf_log", CLIN_RHS)
  dc <- d[complete.cases(d[, ..need])]
  rbindlist(lapply(names(NF_MODELS), function(nm) rbindlist(lapply(
    c("nf", "nf_log"), function(term) {
      r <- cox_row(dc, term, NF_MODELS[[nm]])$row
      cbind(data.table(endpoint = ep,
                       exposure_scale = if (term == "nf") "linear, per SD" else "log10, per SD",
                       model = nm), r)
    }))))
}))
exp_tbl[, exposure_scale := factor(exposure_scale, levels = c("linear, per SD", "log10, per SD"))]
setorder(exp_tbl, exposure_scale, endpoint, model)
exp_tbl[, exposure_scale := as.character(exposure_scale)]
save_tsv(exp_tbl, "30_endpoint_exposure.tsv")
print(exp_tbl)

# ---- 4. The lncRNA axis: nested sequence per endpoint ----
banner("4 | lncRNA axis nested sequence, per endpoint")
NESTED <- list(
  `Axis alone`                               = character(0),
  `+ age, sex, T, N, M1, ordinal grade`      = c(CLIN_RHS),
  `+ ESTIMATE stromal and immune scores`     = c(CLIN_RHS, "stromal", "immune"),
  `+ non-feature, multimapping, log10 depth` = c(CLIN_RHS, "stromal", "immune", "nf", "mm", "dep"))
NEED_AX <- c("time", "event", "axis_z", CLIN_RHS, "stromal", "immune", "nf", "mm", "dep")
vif_of <- function(d, term, others) {
  if (!length(others)) return(1)
  r2 <- summary(lm(as.formula(paste(term, "~", paste(others, collapse = " + "))), data = d))$r.squared
  1 / (1 - r2)
}
ax <- merge(axis, cov, by = "sample_barcode")
ax[, axis_z := as.numeric(scale(lnc_axis))]           # per SD over the axis samples, as in 07
axis_tbl <- rbindlist(lapply(ANALYSED, function(ep) {
  d  <- merge(ax, ep_of(ep), by = "sample_barcode")
  dc <- d[complete.cases(d[, ..NEED_AX])]
  rbindlist(lapply(names(NESTED), function(nm) {
    r <- cox_row(dc, "axis_z", NESTED[[nm]])$row
    cbind(data.table(endpoint = ep, model = nm), r,
          data.table(vif_axis = round(vif_of(dc, "axis_z", NESTED[[nm]]), 3)))
  }))
}))
save_tsv(axis_tbl, "30_endpoint_axis_nested.tsv")
print(axis_tbl)

# ---- 5. Locked-model increment in discovery (apparent), per endpoint ----
banner("5 | Locked models of 09 re-scored against each endpoint (discovery, apparent)")
# Apparent, in-sample: the stage 09 coefficients were fitted on the pipeline OS
# of these same patients and are re-scored against each endpoint.
EVALUATION_LABEL <- paste0(
  "apparent (in-sample): coefficients estimated on the pipeline overall survival of ",
  "these same patients")
INTERVAL_LABEL <- "paired patient bootstrap 95%, predictors fixed, in-sample"
lp_of <- function(X, b) as.numeric(X[, names(b), drop = FALSE] %*% b)
locked_tbl <- rbindlist(lapply(names(L$models), function(cmp) {
  M  <- L$models[[cmp]]; Xt <- L$X[[cmp]]$Xt
  stopifnot(identical(rownames(Xt), L$tsamp), identical(rownames(L$Et), L$tsamp))
  lp_ref <- lp_of(Xt, M$b_clin)
  lp_new <- lp_of(cbind(Xt, L$Et), M$b_full)
  rbindlist(lapply(ANALYSED, function(ep) {
    y  <- ep_of(ep)
    i  <- match(L$tsamp, y$sample_barcode); ok <- !is.na(i)
    yy <- Surv(y$time[i[ok]], y$event[i[ok]])
    c_ref <- cindex_ci(yy, lp_ref[ok]); c_new <- cindex_ci(yy, lp_new[ok])
    b  <- paired_boot_delta_c(yy, lp_new[ok], lp_ref[ok])
    data.table(endpoint = ep, comparator = cmp,
               n_locked_set = length(L$tsamp), n = sum(ok), events = sum(yy[, 2]),
               n_modules_retained = length(M$sel_me),
               C_comparator = round(c_ref[["C"]], 4), C_comparator_lo = round(c_ref[["lo"]], 4),
               C_comparator_hi = round(c_ref[["hi"]], 4),
               C_combined = round(c_new[["C"]], 4), C_combined_lo = round(c_new[["lo"]], 4),
               C_combined_hi = round(c_new[["hi"]], 4),
               delta_C = round(b[["delta"]], 4), lo = round(b[["lo"]], 4), hi = round(b[["hi"]], 4),
               p_boot = signif(b[["p_boot"]], 3), prop_positive = round(b[["prop_positive"]], 3),
               n_boot = b[["n_boot"]],
               evaluation = EVALUATION_LABEL, interval_type = INTERVAL_LABEL,
               coefficients = if (ep == "OS_pipeline")
                 "fitted on this endpoint in these patients (apparent, as 09)"
               else paste0("fitted on the pipeline overall survival of these same patients ",
                           "and re-scored against a different endpoint in the same patients; ",
                           "not external validation"))
  }))
}))
# BH across all bootstrap p values in this table.
locked_tbl[, p_boot_BH_8rows := signif(p.adjust(p_boot, "BH"), 3)]
locked_tbl[, multiplicity_note := paste0(
  "Benjamini-Hochberg across the ", nrow(locked_tbl),
  " bootstrap p-values in this file (", length(ANALYSED), " endpoints x ",
  length(L$models), " comparators); ", ENDPOINT_DEPENDENCE,
  "; the two comparators are nested and share the same patients, which makes the ",
  "adjustment cruder still")]
save_tsv(locked_tbl, "30_endpoint_locked_model.tsv")
print(locked_tbl[, .(endpoint, comparator, n, events, C_comparator, C_combined, delta_C, lo, hi,
                     p_boot, p_boot_BH_8rows)])

# ---- 6. Repeated cross-validation: clinical vs clinical + all module eigengenes ----
banner("6 | 10 x 10 cross-validation of the eigengene increment, per endpoint")
# Analysis set as in stage 05: samples in both networks with complete clinical
# inputs, then per endpoint the rows with a usable endpoint.
CLIN_INPUTS <- c("age", "sex", "T_stage", "N_pos", "M1", "grade_num")
in_nets <- intersect(rownames(nets$lnc$expr), rownames(nets$mrna$expr))
cand    <- cohort[sample_barcode %in% in_nets]
cand    <- cand[complete.cases(cand[, ..CLIN_INPUTS])]
Xe_all  <- cbind(SC$mrna[match(cand$sample_barcode, rownames(SC$mrna)), , drop = FALSE],
                 SC$lnc [match(cand$sample_barcode, rownames(SC$lnc)),  , drop = FALSE])
rownames(Xe_all) <- cand$sample_barcode
stopifnot(!anyNA(Xe_all))
msg("Candidate set (both networks, complete clinical inputs): ", nrow(cand),
    "; module scores as predictors: ", ncol(Xe_all))

# Penalised fit as in stage 05: elastic net, truncated lambda path, comparator
# columns unpenalised, lambda by inner 5-fold cv.glmnet. inner_foldid() draws
# the inner split from the repeat and fold index and restores the RNG state,
# so lambda.min does not depend on how many models were fitted before.
INNER_NFOLDS <- 5L
inner_foldid <- function(n_train, seed, nfolds = INNER_NFOLDS) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv))
           get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  fid <- sample(rep(seq_len(nfolds), length.out = n_train))
  if (is.null(old)) rm(".Random.seed", envir = .GlobalEnv)
  else assign(".Random.seed", old, envir = .GlobalEnv)
  fid
}
fit_cv <- function(x, yy, pfac, fid) {
  cv.glmnet(x, yy, family = "cox", alpha = ML_ALPHA, foldid = fid,
            nlambda = ML_NLAMBDA, lambda.min.ratio = ML_LAMBDA_MIN_RATIO,
            thresh = ML_THRESH, maxit = ML_MAXIT, penalty.factor = pfac)
}
lp_cox <- function(fit, Xte) {
  b <- coef(fit); b[is.na(b)] <- 0        # NA only for a column constant within the fold
  as.numeric(Xte[, names(b), drop = FALSE] %*% b)
}
lp_net <- function(fit, Xte) as.numeric(predict(fit, newx = Xte, s = "lambda.min"))
cv_repeat <- function(r, y, Xc, Xe, ev) {
  set.seed(SEED + r)
  n <- nrow(Xc)
  out <- matrix(NA_real_, n, 2, dimnames = list(rownames(Xc), c("clinical", "clinical_eig")))
  folds <- sample(rep(seq_len(ML_N_FOLDS), length.out = n))
  Xce <- cbind(Xc, Xe); pf <- c(rep(0, ncol(Xc)), rep(1, ncol(Xe)))
  for (k in seq_len(ML_N_FOLDS)) {
    tr <- folds != k; te <- !tr
    if (sum(ev[tr]) < 5) next
    out[te, "clinical"] <- lp_cox(coxph(y[tr] ~ ., data = as.data.frame(Xc[tr, , drop = FALSE])),
                                  Xc[te, , drop = FALSE])
    fid <- inner_foldid(sum(tr), SEED + 1000L * r + k)
    out[te, "clinical_eig"] <- lp_net(fit_cv(Xce[tr, ], y[tr], pf, fid), Xce[te, , drop = FALSE])
  }
  out
}
par_cl <- makeCluster(min(10L, N_THREADS))
invisible(clusterEvalQ(par_cl, { library(survival); library(glmnet) }))
clusterExport(par_cl, c("fit_cv", "lp_cox", "lp_net", "inner_foldid", "INNER_NFOLDS",
                        "SEED", "ML_N_FOLDS", "ML_ALPHA",
                        "ML_NLAMBDA", "ML_LAMBDA_MIN_RATIO", "ML_THRESH", "ML_MAXIT",
                        "cindex", "paired_boot_delta_c", "BOOT_B"))
cv_res <- lapply(ANALYSED, function(ep) {
  y0 <- ep_of(ep)
  cl <- cand[sample_barcode %in% y0$sample_barcode]
  y0 <- y0[match(cl$sample_barcode, sample_barcode)]
  Xc <- clinical_design(cl, "clinical"); stopifnot(!anyNA(Xc))
  Xe <- Xe_all[cl$sample_barcode, , drop = FALSE]
  y  <- Surv(y0$time, y0$event); ev <- y0$event; n <- nrow(Xc)
  msg("[", ep, "] analysis set n = ", n, ", events = ", sum(ev), "; running ",
      ML_N_REPEATS, " x ", ML_N_FOLDS, "-fold CV ...")
  t0 <- Sys.time()
  reps <- parLapply(par_cl, seq_len(ML_N_REPEATS), cv_repeat, y = y, Xc = Xc, Xe = Xe, ev = ev)
  msg("[", ep, "] done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  per_rep <- rbindlist(lapply(seq_along(reps), function(r) data.table(
    endpoint = ep, repeat_id = r, n = n, events = sum(ev),
    C_clinical = cindex(y, reps[[r]][, "clinical"]),
    C_clinical_eig = cindex(y, reps[[r]][, "clinical_eig"]))))
  per_rep[, delta_clinical_eig_vs_clinical := C_clinical_eig - C_clinical]
  lp_avg <- Reduce(`+`, reps) / length(reps)
  c_ref <- cindex_ci(y, lp_avg[, "clinical"]); c_new <- cindex_ci(y, lp_avg[, "clinical_eig"])
  bt <- paired_boot_delta_c(y, lp_avg[, "clinical_eig"], lp_avg[, "clinical"])
  d  <- per_rep$delta_clinical_eig_vs_clinical

  # ---- how much of the interval is the bootstrap seed? ----------------------
  # The paired bootstrap is repeated over BOOT_SEED_SWEEP seeds on the same out-of-fold
  # predictors to show the spread of the lower bound and p. SEED is the first seed.
  msg("[", ep, "] bootstrap-seed sweep, ", BOOT_SEED_SWEEP, " seeds x ", BOOT_B,
      " resamples ...")
  t1 <- Sys.time()
  sweep_seeds <- SEED + seq_len(BOOT_SEED_SWEEP) - 1L
  sw <- as.data.table(do.call(rbind, parLapply(
    par_cl, sweep_seeds,
    function(s, yy, a, b) {
      z <- paired_boot_delta_c(yy, a, b, seed = s)
      c(seed = s, lo = z[["lo"]], hi = z[["hi"]], p = z[["p_boot"]])
    },
    yy = y, a = lp_avg[, "clinical_eig"], b = lp_avg[, "clinical"])))
  msg("[", ep, "] seed sweep done in ",
      round(difftime(Sys.time(), t1, units = "mins"), 1), " min; lower bound ",
      signif(min(sw$lo), 2), " to ", signif(max(sw$lo), 2), ", p ",
      signif(min(sw$p), 2), " to ", signif(max(sw$p), 2), ", ",
      round(100 * mean(sw$lo > 0)), "% of seeds exclude zero")
  sw_out <- cbind(data.table(endpoint = ep), sw)

  summ <- data.table(endpoint = ep, comparison = "clinical_eig - clinical", n = n, events = sum(ev),
                     n_module_scores = ncol(Xe),
                     C_clinical_cv_mean = round(mean(per_rep$C_clinical), 4),
                     C_clinical_eig_cv_mean = round(mean(per_rep$C_clinical_eig), 4),
                     delta_mean = round(mean(d), 4), range_lo = round(min(d), 4), range_hi = round(max(d), 4),
                     pct2.5 = round(unname(quantile(d, 0.025)), 4),
                     pct97.5 = round(unname(quantile(d, 0.975)), 4),
                     n_repeats_positive = sum(d > 0), n_repeats = length(d),
                     C_avgscore_clinical = round(c_ref[["C"]], 4),
                     C_avgscore_clinical_lo = round(c_ref[["lo"]], 4),
                     C_avgscore_clinical_hi = round(c_ref[["hi"]], 4),
                     C_avgscore_clinical_eig = round(c_new[["C"]], 4),
                     C_avgscore_clinical_eig_lo = round(c_new[["lo"]], 4),
                     C_avgscore_clinical_eig_hi = round(c_new[["hi"]], 4),
                     boot_delta = round(bt[["delta"]], 4), boot_lo = round(bt[["lo"]], 4),
                     boot_hi = round(bt[["hi"]], 4), boot_p = signif(bt[["p_boot"]], 3),
                     boot_prop_positive = round(bt[["prop_positive"]], 3), n_boot = bt[["n_boot"]],
                     boot_seed_pipeline = SEED, boot_seed_n = nrow(sw),
                     boot_lo_median   = round(median(sw$lo), 4),
                     boot_lo_pct2.5   = round(unname(quantile(sw$lo, 0.025)), 4),
                     boot_lo_pct97.5  = round(unname(quantile(sw$lo, 0.975)), 4),
                     boot_lo_min      = round(min(sw$lo), 4),
                     boot_lo_max      = round(max(sw$lo), 4),
                     boot_hi_median   = round(median(sw$hi), 4),
                     boot_p_median    = signif(median(sw$p), 3),
                     boot_p_min       = signif(min(sw$p), 3),
                     boot_p_max       = signif(max(sw$p), 3),
                     boot_prop_seeds_excluding_zero = round(mean(sw$lo > 0), 3),
                     boot_seed_note = sprintf(paste0(
                       "boot_lo/boot_hi/boot_p are the single pipeline seed %d; over %d seeds ",
                       "on the same fixed predictors the lower bound runs %.4f to %.4f ",
                       "(median %.4f) and p runs %.3f to %.3f, and %.0f%% of seeds gave an ",
                       "interval excluding zero"),
                       SEED, nrow(sw), min(sw$lo), max(sw$lo), median(sw$lo),
                       min(sw$p), max(sw$p), 100 * mean(sw$lo > 0)),
                     increment_note = sprintf(paste0(
                       "the increment that carries the interval (boot_delta = %+.4f) is that of ",
                       "the predictor averaged over the %d cross-validation repeats; the mean of ",
                       "the %d single-run increments is delta_mean = %+.4f (range %+.4f to ",
                       "%+.4f) and carries no interval"),
                       bt[["delta"]], length(d), length(d), mean(d), min(d), max(d)))
  list(per_rep = per_rep, summary = summ, seed_sweep = sw_out)
})
stopCluster(par_cl)
cv_per_rep <- rbindlist(lapply(cv_res, `[[`, "per_rep"))
cv_tbl     <- rbindlist(lapply(cv_res, `[[`, "summary"))
cv_sweep   <- rbindlist(lapply(cv_res, `[[`, "seed_sweep"))
# BH across the endpoint increments in this table.
cv_tbl[, boot_p_BH_4endpoints := signif(p.adjust(boot_p, "BH"), 3)]
cv_tbl[, multiplicity_note := paste0(
  "Benjamini-Hochberg across the ", nrow(cv_tbl),
  " endpoint bootstrap p-values in this file; ", ENDPOINT_DEPENDENCE)]
save_tsv(cv_per_rep, "30_endpoint_cv_per_repeat.tsv")
save_tsv(cv_sweep, "30_endpoint_cv_boot_seed_sweep.tsv")
save_tsv(cv_tbl, "30_endpoint_cv_increment.tsv")
print(cv_tbl[, .(endpoint, n, events, C_clinical_cv_mean, C_clinical_eig_cv_mean, delta_mean,
                 n_repeats_positive, boot_delta, boot_lo, boot_hi, boot_p,
                 boot_p_BH_4endpoints)])
print(cv_tbl[, .(endpoint, boot_lo, boot_lo_median, boot_lo_pct2.5, boot_lo_pct97.5,
                 boot_p, boot_p_median, boot_p_min, boot_p_max,
                 boot_prop_seeds_excluding_zero)])

# ---- 7. Fold-wise network rebuild for the disease-specific endpoint ----
banner("7 | Modules recomputed inside every training fold (DSS, with OS as reference)")
# Inside each training fold: technical residualisation, both networks (stage 02
# parameters) and module loadings are fitted on the training fold and applied to the
# held-out fold, then an elastic-net Cox (comparator unpenalised) is scored out of fold.
# The same folds score the full-cohort ("fixed") modules, so rebuilt and fixed increments
# differ only in whether modules saw the held-out patients. OS is scored on the same
# patients and modules as the endpoint-matched reference.
enableWGCNAThreads(N_THREADS)

FW_ENDPOINTS <- c("DSS", "OS_pipeline")
fw_usable <- Reduce(intersect, lapply(FW_ENDPOINTS, function(e) ep_of(e)$sample_barcode))
f_samp <- L$tsamp[L$tsamp %in% fw_usable]
f_cl   <- cohort[match(f_samp, sample_barcode)]
f_qc   <- QC[match(f_samp, sample_barcode)]
stopifnot(identical(f_cl$sample_barcode, f_samp), identical(f_qc$sample_barcode, f_samp),
          !anyNA(f_qc$pct_noFeature), !anyNA(f_qc$assigned_reads))
Xf <- clinical_design(f_cl, "clinical")
stopifnot(!anyNA(Xf), identical(rownames(Xf), f_samp))
y_fw <- lapply(setNames(FW_ENDPOINTS, FW_ENDPOINTS), function(e) {
  yy <- ep_of(e)[match(f_samp, sample_barcode)]
  stopifnot(identical(yy$sample_barcode, f_samp), !anyNA(yy$time), !anyNA(yy$event))
  Surv(yy$time, yy$event)
})
cov_f   <- tech_covariates(f_qc)
E_obs_f <- list(mrna = obs_expr(nets$mrna)[f_samp, , drop = FALSE],
                lnc  = obs_expr(nets$lnc)[f_samp, , drop = FALSE])
S_fixed <- cbind(SC$mrna[f_samp, , drop = FALSE], SC$lnc[f_samp, , drop = FALSE])
stopifnot(!anyNA(S_fixed), identical(rownames(S_fixed), f_samp))
NET_PAR <- list(
  mrna = list(prefix = "mRNA_ME", power = nets$mrna$power, deep_split = DEEP_SPLIT,
              merge_cut = MERGE_CUT_HEIGHT, min_kme = MIN_KME_TO_STAY),
  lnc  = list(prefix = "lnc_ME", power = LNC_POWER,
              deep_split = if (is.na(LNC_DEEPSPLIT)) DEEP_SPLIT else LNC_DEEPSPLIT,
              merge_cut  = if (is.na(LNC_MERGE)) MERGE_CUT_HEIGHT else LNC_MERGE,
              min_kme    = if (is.na(LNC_MINKME)) MIN_KME_TO_STAY else LNC_MINKME))
for (e in FW_ENDPOINTS)
  msg("Fold-wise set: n = ", length(f_samp), ", ", e, " events = ", sum(y_fw[[e]][, 2]))
msg("K = ", FOLDWISE_K, " folds x ", FOLDWISE30_N_ASSIGN, " assignments; ",
    ncol(S_fixed), " full-cohort module scores as the fixed comparison")

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

# Each assignment (FOLDWISE_K x 2 network rebuilds) is cached when it finishes
# and reused on a re-run.
FW_CACHE <- file.path(CACHE_DIR, "30_foldwise_rebuild.rds")
fw_key <- list(samples = f_samp, K = FOLDWISE_K, endpoints = FW_ENDPOINTS,
               net_par = NET_PAR, n_fixed_scores = ncol(S_fixed))
fw <- if (file.exists(FW_CACHE)) readRDS(FW_CACHE) else NULL
if (!is.null(fw) && !identical(fw$key, fw_key)) {
  msg("Fold-wise cache does not match the current analysis set or parameters; recomputing")
  fw <- NULL
}
if (is.null(fw)) fw <- list(key = fw_key, per_fold = list(), per_assign = list())

for (a in seq_len(FOLDWISE30_N_ASSIGN)) {
  ka <- as.character(a)
  if (!is.null(fw$per_assign[[ka]])) {
    msg("  assignment ", a, "/", FOLDWISE30_N_ASSIGN, " from cache"); next
  }
  t_assign <- Sys.time()
  set.seed(SEED + a)
  folds <- sample(rep(seq_len(FOLDWISE_K), length.out = length(f_samp)))
  oof <- lapply(setNames(FW_ENDPOINTS, FW_ENDPOINTS), function(e)
    list(comparator = rep(NA_real_, length(f_samp)),
         rebuilt    = rep(NA_real_, length(f_samp)),
         fixed      = rep(NA_real_, length(f_samp))))
  fold_rows <- list()
  for (k in seq_len(FOLDWISE_K)) {
    t_fold <- Sys.time()
    tr <- folds != k; te <- !tr
    S_tr <- list(); S_te <- list(); n_mod <- setNames(integer(2), names(E_obs_f))
    for (nm in names(E_obs_f)) {
      tf  <- fit_technical(E_obs_f[[nm]][tr, , drop = FALSE], cov_f[tr, , drop = FALSE])
      Etr <- apply_technical(E_obs_f[[nm]][tr, , drop = FALSE], tf, cov_f[tr, , drop = FALSE])
      Ete <- apply_technical(E_obs_f[[nm]][te, , drop = FALSE], tf, cov_f[te, , drop = FALSE])
      gt  <- rebuild_modules(Etr, NET_PAR[[nm]])
      n_mod[nm] <- length(setdiff(unique(gt$module), "grey"))
      Ltr <- fit_module_loadings(Etr, gt, NET_PAR[[nm]]$prefix)
      S_tr[[nm]] <- score_modules(Etr, Ltr, gene_standardise = "discovery",
                                  score_standardise = "discovery")
      S_te[[nm]] <- score_modules(Ete, Ltr, gene_standardise = "discovery",
                                  score_standardise = "discovery")
    }
    S_tr <- do.call(cbind, S_tr); S_te <- do.call(cbind, S_te)
    stopifnot(identical(colnames(S_tr), colnames(S_te)))
    MODSETS <- list(rebuilt = list(tr = S_tr, te = S_te),
                    fixed   = list(tr = S_fixed[tr, , drop = FALSE],
                                   te = S_fixed[te, , drop = FALSE]))
    # One inner fold split per outer fold, shared by both module sets and both
    # endpoints, so the module definition is the only thing that differs.
    fid <- inner_foldid(sum(tr), SEED + 500000L + 1000L * a + k, nfolds = 10L)
    for (e in FW_ENDPOINTS) {
      y_tr <- y_fw[[e]][tr]
      Xtr <- Xf[tr, , drop = FALSE]; Xte <- Xf[te, , drop = FALSE]; dropped <- ""
      fc  <- coxph(y_tr ~ ., data = as.data.frame(Xtr))
      # A comparator term constant within a training fold is dropped from every
      # model of that fold, keeping the comparison paired, and is recorded.
      bad <- is.na(coef(fc))
      if (any(bad)) {
        dropped <- paste(names(coef(fc))[bad], collapse = ";")
        warning("Fold-wise ", e, ": dropping inestimable term(s) ", dropped,
                " in assignment ", a, ", fold ", k, call. = FALSE, immediate. = TRUE)
        Xtr <- Xtr[, !bad, drop = FALSE]; Xte <- Xte[, !bad, drop = FALSE]
        fc  <- coxph(y_tr ~ ., data = as.data.frame(Xtr))
        stopifnot(!anyNA(coef(fc)))
      }
      oof[[e]]$comparator[te] <- as.numeric(Xte[, names(coef(fc)), drop = FALSE] %*% coef(fc))
      n_ret <- setNames(integer(length(MODSETS)), names(MODSETS))
      for (ms in names(MODSETS)) {
        Str <- MODSETS[[ms]]$tr; Ste <- MODSETS[[ms]]$te
        if (!ncol(Str)) {
          oof[[e]][[ms]][te] <- oof[[e]]$comparator[te]; n_ret[ms] <- 0L; next
        }
        fit <- cv.glmnet(cbind(Xtr, Str), y_tr, family = "cox", alpha = ML_ALPHA,
                         foldid = fid, nlambda = ML_NLAMBDA,
                         lambda.min.ratio = ML_LAMBDA_MIN_RATIO, thresh = ML_THRESH,
                         maxit = ML_MAXIT,
                         penalty.factor = c(rep(0, ncol(Xtr)), rep(1, ncol(Str))))
        newx <- cbind(Xte, Ste)
        stopifnot(identical(colnames(newx), colnames(cbind(Xtr, Str))))
        oof[[e]][[ms]][te] <- as.numeric(predict(fit, newx = newx, s = "lambda.min"))
        bf <- as.matrix(coef(fit, s = "lambda.min"))[, 1]
        n_ret[ms] <- sum(bf[colnames(Str)] != 0)
      }
      fold_rows[[length(fold_rows) + 1]] <- data.table(
        assignment = a, fold = k, endpoint = e, comparator = "clinical",
        n_train = sum(tr), n_test = sum(te), events_train = sum(y_tr[, 2]),
        n_modules_mRNA_rebuilt = n_mod[["mrna"]], n_modules_lncRNA_rebuilt = n_mod[["lnc"]],
        n_scores_rebuilt = ncol(S_tr), n_scores_fixed = ncol(S_fixed),
        n_retained_rebuilt = n_ret[["rebuilt"]], n_retained_fixed = n_ret[["fixed"]],
        dropped_terms = dropped,
        elapsed_min = round(as.numeric(difftime(Sys.time(), t_fold, units = "mins")), 2))
    }
    msg(sprintf("  assignment %d/%d, fold %d/%d: %d mRNA + %d lncRNA modules rebuilt (%d scores), %.1f min",
                a, FOLDWISE30_N_ASSIGN, k, FOLDWISE_K, n_mod[["mrna"]], n_mod[["lnc"]],
                ncol(S_tr), as.numeric(difftime(Sys.time(), t_fold, units = "mins"))))
  }
  el_a <- as.numeric(difftime(Sys.time(), t_assign, units = "mins"))
  fw$per_fold[[ka]] <- rbindlist(fold_rows)
  fw$per_assign[[ka]] <- rbindlist(lapply(FW_ENDPOINTS, function(e) {
    yy <- y_fw[[e]]
    stopifnot(!anyNA(oof[[e]]$comparator), !anyNA(oof[[e]]$rebuilt), !anyNA(oof[[e]]$fixed))
    Cc <- cindex(yy, oof[[e]]$comparator)
    Cr <- cindex(yy, oof[[e]]$rebuilt)
    Cx <- cindex(yy, oof[[e]]$fixed)
    data.table(assignment = a, seed = SEED + a, endpoint = e, comparator = "clinical",
               n = length(f_samp), events = sum(yy[, 2]),
               C_comparator_oof = round(Cc, 4),
               C_combined_rebuilt_oof = round(Cr, 4),
               C_combined_fixed_oof = round(Cx, 4),
               delta_rebuilt = round(Cr - Cc, 4), delta_fixed = round(Cx - Cc, 4),
               elapsed_min = round(el_a, 1))
  }))
  saveRDS(fw, FW_CACHE)
  msg(sprintf("  assignment %d done in %.1f min", a, el_a))
  print(fw$per_assign[[ka]])
}

keys          <- as.character(seq_len(FOLDWISE30_N_ASSIGN))
per_fold_fw   <- rbindlist(fw$per_fold[keys])
per_assign_fw <- rbindlist(fw$per_assign[keys])
save_tsv(per_fold_fw, "30_endpoint_foldwise_per_fold.tsv")
save_tsv(per_assign_fw, "30_endpoint_foldwise_per_assignment.tsv")
print(per_assign_fw)

# Apparent locked-model increment on the same patients, so the three levels of
# optimism control (none, fixed modules cross-validated, modules rebuilt in
# fold) share one set.
f_idx    <- match(f_samp, L$tsamp)
Xt_f     <- L$X$clinical$Xt[f_idx, , drop = FALSE]
Et_f     <- L$Et[f_idx, , drop = FALSE]
stopifnot(identical(rownames(Xt_f), f_samp), identical(rownames(Et_f), f_samp))
lp_ref_f <- lp_of(Xt_f, L$models$clinical$b_clin)
lp_new_f <- lp_of(cbind(Xt_f, Et_f), L$models$clinical$b_full)

fw_tbl <- rbindlist(lapply(FW_ENDPOINTS, function(e) {
  pa  <- per_assign_fw[endpoint == e]
  yy  <- y_fw[[e]]
  bt  <- paired_boot_delta_c(yy, lp_new_f, lp_ref_f)
  cvr <- cv_tbl[endpoint == e]
  pr  <- cv_per_rep[endpoint == e]
  optimism <- mean(pa$delta_fixed) - mean(pa$delta_rebuilt)
  rbind(
    data.table(
      analysis = "modules fixed on full cohort (apparent, locked coefficients of 09)",
      endpoint = e, comparator = "clinical", n = length(f_samp), events = sum(yy[, 2]),
      C_clinical = round(cindex(yy, lp_ref_f), 4),
      C_combined = round(cindex(yy, lp_new_f), 4),
      delta_C = round(bt[["delta"]], 4), delta_sd = NA_real_,
      delta_lo = round(bt[["lo"]], 4), delta_hi = round(bt[["hi"]], 4),
      n_assignments = NA_integer_,
      interval_type = "paired patient bootstrap 95% (in-sample)",
      note = paste0("coefficients estimated on the pipeline overall survival of these same ",
                    "patients; no cross-validation of any kind")),
    data.table(
      analysis = "modules fixed on full cohort (10 x 10 CV, section 6)",
      endpoint = e, comparator = "clinical", n = cvr$n, events = cvr$events,
      C_clinical = cvr$C_avgscore_clinical, C_combined = cvr$C_avgscore_clinical_eig,
      delta_C = cvr$boot_delta,
      delta_sd = round(sd(pr$delta_clinical_eig_vs_clinical), 4),
      delta_lo = cvr$boot_lo, delta_hi = cvr$boot_hi,
      n_assignments = cvr$n_repeats,
      interval_type = "paired patient bootstrap 95% on the repeat-averaged out-of-fold predictor",
      note = paste0("section 6 of this script; n = ", cvr$n,
                    " on its own analysis set, which is ",
                    ifelse(cvr$n == length(f_samp), "the same as",
                           paste0("wider than")), " the ", length(f_samp),
                    " patients of the fold-wise rows")),
    data.table(
      analysis = sprintf("modules fixed on full cohort (%d x %d assignments, identical folds)",
                         FOLDWISE_K, FOLDWISE30_N_ASSIGN),
      endpoint = e, comparator = "clinical", n = length(f_samp), events = sum(yy[, 2]),
      C_clinical = round(mean(pa$C_comparator_oof), 4),
      C_combined = round(mean(pa$C_combined_fixed_oof), 4),
      delta_C = round(mean(pa$delta_fixed), 4), delta_sd = round(sd(pa$delta_fixed), 4),
      delta_lo = round(min(pa$delta_fixed), 4), delta_hi = round(max(pa$delta_fixed), 4),
      n_assignments = nrow(pa), interval_type = "range across fold assignments",
      note = paste0("the matched reference for the row below: same patients, same folds, ",
                    "same inner fold split, same model class; only the module definition differs")),
    data.table(
      analysis = sprintf("modules recomputed inside each fold (%d x %d assignments)",
                         FOLDWISE_K, FOLDWISE30_N_ASSIGN),
      endpoint = e, comparator = "clinical", n = length(f_samp), events = sum(yy[, 2]),
      C_clinical = round(mean(pa$C_comparator_oof), 4),
      C_combined = round(mean(pa$C_combined_rebuilt_oof), 4),
      delta_C = round(mean(pa$delta_rebuilt), 4), delta_sd = round(sd(pa$delta_rebuilt), 4),
      delta_lo = round(min(pa$delta_rebuilt), 4), delta_hi = round(max(pa$delta_rebuilt), 4),
      n_assignments = nrow(pa), interval_type = "range across fold assignments",
      note = sprintf(paste0("both networks, both residualisations and all loadings rebuilt ",
                            "inside each training fold; held-out patients projected with the ",
                            "training fold's centring, scaling and rotation. Module-definition ",
                            "optimism on this endpoint and these patients = fixed minus rebuilt ",
                            "= %+.4f"), optimism)))
}))
save_tsv(fw_tbl, "30_endpoint_foldwise.tsv")
print(fw_tbl[, .(analysis, endpoint, n, events, C_clinical, C_combined, delta_C, delta_sd,
                 delta_lo, delta_hi, n_assignments)])
for (e in FW_ENDPOINTS) {
  pa <- per_assign_fw[endpoint == e]
  msg(sprintf(paste0("[%s] fixed modules %+.4f (sd %.4f, %.4f to %.4f) vs rebuilt %+.4f ",
                     "(sd %.4f, %.4f to %.4f) over %d assignments; optimism %+.4f; ",
                     "assignments with a positive rebuilt increment: %d of %d"),
              e, mean(pa$delta_fixed), sd(pa$delta_fixed), min(pa$delta_fixed),
              max(pa$delta_fixed), mean(pa$delta_rebuilt), sd(pa$delta_rebuilt),
              min(pa$delta_rebuilt), max(pa$delta_rebuilt), nrow(pa),
              mean(pa$delta_fixed) - mean(pa$delta_rebuilt),
              sum(pa$delta_rebuilt > 0), nrow(pa)))
}

banner("Summary")
cat("Endpoints analysed:", paste(ANALYSED, collapse = ", "), "\n")
cat("Coverage (cohort of", nrow(cohort), "):\n")
print(coverage[, .(endpoint, n_analysable, events_10y, epv_principal_12_parameters, analysed_here)])
cat("\nPrincipal-specification hazard ratios per SD, five modules x endpoints:\n")
print(dcast(mod_tbl[specification == "principal"], biotype + module ~ endpoint, value.var = "HR"))
cat("\nNon-feature fraction (linear, per SD) + clinical terms:\n")
print(exp_tbl[exposure_scale == "linear, per SD" & grepl("^\\(b\\)", model),
              .(endpoint, n, events, HR, lo, hi, p)])
cat("\nLocked-model increment (clinical comparator; APPARENT, in-sample):\n")
print(locked_tbl[comparator == "clinical",
                 .(endpoint, n, events, delta_C, lo, hi, p_boot, p_boot_BH_8rows)])
cat("\nCross-validated increment (boot_delta is the repeat-averaged predictor):\n")
print(cv_tbl[, .(endpoint, n, events, delta_mean, boot_delta, boot_lo, boot_hi, boot_p,
                 boot_p_BH_4endpoints)])
cat("\nBootstrap-seed stability of the cross-validated interval:\n")
print(cv_tbl[, .(endpoint, boot_lo, boot_lo_median, boot_lo_min, boot_lo_max,
                 boot_p_median, boot_p_min, boot_p_max, boot_prop_seeds_excluding_zero)])
cat("\nFold-wise network rebuild (clinical comparator, ", length(f_samp), " patients):\n", sep = "")
print(fw_tbl[, .(endpoint, analysis, n, events, delta_C, delta_sd, delta_lo, delta_hi)])
msg("Total runtime: ", round(difftime(Sys.time(), t_start, units = "mins"), 1), " min")

write_session_info("30")
banner("30 | done")
