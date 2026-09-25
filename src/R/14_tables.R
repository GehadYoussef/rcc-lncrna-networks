# 14_tables.R: main tables, built from the result files and written as markdown.
#
# Reads result tables from stages 01, 03, 05, 08, 09, 11 and 12 and the cached
# subtype cohorts. Writes results/14_manuscript_tables.md. Run after stage 12.
# The primary comparator is the clinical model (age, sex, T, N, M1, grade). The
# augmented comparator adds composition and quality metrics. Module scores
# standardised within the scored cohort are primary. Filters on optional
# columns are guarded, so a partial results directory also works.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(survival) })
banner("14 | Manuscript tables")

R_  <- function(f) fread(file.path(RESULTS_DIR, f))
Rx  <- function(f) if (file.exists(file.path(RESULTS_DIR, f))) R_(f) else NULL
has <- function(x, col) col %in% names(x)
# Keep rows where `col` is in `vals`. Returns `x` unchanged if the column is
# missing, or with a note if no row matches.
keep_if <- function(x, col, vals) {
  if (!has(x, col)) return(x)
  y <- x[x[[col]] %in% vals]
  if (nrow(y)) return(y)
  msg("  note: no rows with ", col, " in {", paste(vals, collapse = ", "),
      "}; using all rows")
  x
}
col_or <- function(x, col, default = NA) if (has(x, col)) x[[col]] else rep(default, nrow(x))

# ---- formatting -------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = d))
ci_txt <- function(lo, hi, d = 3, sep = "-")
  ifelse(is.na(lo) | is.na(hi), NA_character_, paste0(fmt(lo, d), sep, fmt(hi, d)))
est_ci <- function(est, lo, hi, d = 3, sep = "-")
  ifelse(is.na(est), NA_character_, paste0(fmt(est, d), " (", ci_txt(lo, hi, d, sep), ")"))
p_txt <- function(p) ifelse(is.na(p), NA_character_,
                     ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                            as.character(signif(p, 2))))

out <- c()
add <- function(...) out <<- c(out, ...)
md <- function(dt) {
  dt <- as.data.table(dt)[, lapply(.SD, function(z) {
    z <- as.character(z); z[is.na(z) | z == "NA"] <- "--"; z })]
  hdr <- paste0("| ", paste(names(dt), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|")
  rows <- apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(hdr, sep, rows)
}

# ---- cohort characteristics ----
# Discovery and CPTAC-3 rows come from stages 01 and 08. Subtype rows use the
# stage 11 summaries and the cached cohorts for the remaining columns.
t1d <- R_("01_table1_cohort.tsv")
t1v <- R_("08_validation_table1.tsv")
qm  <- R_("11_cohort_quality_metrics.tsv")
ss  <- R_("11_subtype_cohort_summary.tsv")

t1_row <- function(name, n, deaths, fu, age, male, t34, m1_all, m1_ass, g34,
                   nf_med, nf_lo, nf_hi)
  data.table(
    Cohort = name, Patients = n, Deaths = deaths,
    `Median follow-up, reverse KM (d)` = round(fu),
    `Median age (y)` = fmt(age, 1),
    `Male (%)` = fmt(male, 1),
    `T3-T4 (%)` = fmt(t34, 1),
    `M1 (%), all patients (assessed)` =
      ifelse(is.na(m1_all), NA_character_,
             paste0(fmt(m1_all, 1), " (", ifelse(is.na(m1_ass), "--", fmt(m1_ass, 1)), ")")),
    `Grade 3-4 (%)` = fmt(g34, 1),
    `Non-feature reads (%), median (IQR)` =
      ifelse(is.na(nf_med), NA_character_,
             paste0(fmt(nf_med, 1), " (", fmt(nf_lo, 1), "-", fmt(nf_hi, 1), ")")))

from_table1 <- function(x, name)
  t1_row(name, x$n, x$deaths, x$median_fu_reverse_km, x$age_median, x$male_pct,
         x$T3_T4_pct, x$M1_pct_all, x$M1_pct_assessed, x$grade_G3_4_pct,
         x$noFeature_median, x$noFeature_IQR_lo, x$noFeature_IQR_hi)

subtype_row <- function(tag, name) {
  cn <- paste0("TCGA-", tag)
  q <- qm[cohort == cn]; s <- ss[cohort == cn]
  n      <- if (nrow(q)) q$n else if (nrow(s)) s$n_scored else NA
  deaths <- if (nrow(s)) s$events else NA
  m1_all <- if (nrow(s) && all(has(s, c("n_M1_as_coded", "n_M1_reconciled"))))
    100 * (s$n_M1_as_coded + s$n_M1_reconciled) / s$n_scored else NA
  nf <- if (nrow(q)) c(q$noFeature_median, q$noFeature_q25, q$noFeature_q75) else rep(NA, 3)
  fu <- age <- male <- t34 <- m1_ass <- g34 <- NA
  f <- file.path(CACHE_DIR, paste0("subtype_", tag, ".rds"))
  if (file.exists(f)) {
    co <- as.data.table(readRDS(f)$cohort)
    fu   <- unname(summary(survfit(Surv(os_time, 1 - os_event) ~ 1, data = co))$table["median"])
    age  <- median(co$age, na.rm = TRUE)
    male <- 100 * mean(co$sex == "male", na.rm = TRUE)
    t34  <- 100 * mean(co$T_stage >= 3, na.rm = TRUE)
    if (any(!is.na(co$grade_num))) g34 <- 100 * mean(co$grade_num >= 3, na.rm = TRUE)
    assessed <- !is.na(co$pM) & co$pM %in% c("M0", "M1")
    if (any(assessed)) m1_ass <- 100 * mean(co$pM[assessed] == "M1")
    if (is.na(m1_all)) m1_all <- 100 * mean(co$M1 == 1, na.rm = TRUE)
    if (is.na(n)) n <- nrow(co)
    if (is.na(deaths)) deaths <- sum(co$os_event)
  }
  t1_row(name, n, deaths, fu, age, male, t34, m1_all, m1_ass, g34, nf[1], nf[2], nf[3])
}

t1 <- rbind(from_table1(t1d, "TCGA-KIRC (discovery, ccRCC)"),
            from_table1(t1v, "CPTAC-3 (validation, ccRCC)"),
            subtype_row("KIRP", "TCGA-KIRP (papillary)"),
            subtype_row("KICH", "TCGA-KICH (chromophobe)"))
add("**Table 1. Characteristics of the four cohorts.** All were quantified",
    "through the identical GDC STAR-Counts pipeline; TCGA libraries are",
    "poly(A)-selected and CPTAC-3 libraries ribo-depleted total RNA, which is",
    "why the non-feature read fraction differs by an order of magnitude and is",
    "summarised per cohort. Follow-up is the reverse Kaplan-Meier median. M1",
    "is given as a percentage of all patients (stage IV without T4 counted as",
    "M1) with the percentage of patients with an assessed M category in",
    "parentheses. Grade is not standardly applicable to papillary or",
    "chromophobe tumours in these data. Nodal status was assessed in",
    sprintf("%d of %d discovery and %d of %d CPTAC-3 patients and is reported",
            t1d$N_assessed, t1d$n, t1v$N_assessed, t1v$n),
    "in the text rather than tabulated. For TCGA-KIRP and TCGA-KICH, deaths",
    "are counted among the patients entering the adjusted models",
    sprintf("(complete T category: %s).",
            paste(ss[cohort != "TCGA-KIRC", paste0(cohort, " ", n_complete, " of ", n_scored)],
                  collapse = "; ")),
    "", md(t1), "")

# ---- module associations, discovery and external validation ----
sv <- rbind(R_("03_mRNA_module_survival.tsv")[, biotype := "mRNA"],
            R_("03_lncRNA_module_survival.tsv")[, biotype := "lncRNA"], fill = TRUE)
for (cc in c("n_full", "events_full", "fdr_full_joint",
             "HR_clin", "HR_clin_lo", "HR_clin_hi", "fdr_clin"))
  if (!has(sv, cc)) sv[, (cc) := NA_real_]
sig <- sv[fdr_full < FDR_ALPHA]
disc <- sig[, .(biotype, module, Genes = n_genes, n = n_full, Events = events_full,
                `HR, principal (95% CI)` = est_ci(HR_full, HR_full_lo, HR_full_hi, 2),
                `FDR, principal (within network)` = p_txt(fdr_full),
                `FDR, principal (joint, all modules)` = p_txt(fdr_full_joint),
                `HR, clinical comparator (95% CI)` = est_ci(HR_clin, HR_clin_lo, HR_clin_hi, 2),
                `FDR, clinical comparator` = p_txt(fdr_clin),
                .fdr = fdr_full)]

rep12 <- R_("09_module_replication.tsv")
repP  <- Rx("09_module_replication_parsimonious.tsv")
rep_key <- intersect(c("biotype", "module"), names(rep12))
r12 <- rep12[, c(rep_key, "HR_cptac", "lo", "hi", "p_cptac", "same_direction"), with = FALSE]
r12[, `HR CPTAC-3, 12-parameter (95% CI)` := est_ci(HR_cptac, lo, hi, 2)]
r12[, `p (12-parameter)` := p_txt(p_cptac)]
r12[, `Direction replicated` := ifelse(same_direction, "yes", "NO")]
r12 <- r12[, c(rep_key, "HR CPTAC-3, 12-parameter (95% CI)", "p (12-parameter)",
               "Direction replicated"), with = FALSE]
t2 <- merge(disc, r12, by = rep_key, all.x = TRUE)
if (!is.null(repP)) {
  rp <- repP[, c(rep_key, "HR_cptac", "lo", "hi", "p_cptac"), with = FALSE]
  rp[, `HR CPTAC-3, 6-parameter (95% CI)` := est_ci(HR_cptac, lo, hi, 2)]
  rp[, `p (6-parameter)` := p_txt(p_cptac)]
  t2 <- merge(t2, rp[, c(rep_key, "HR CPTAC-3, 6-parameter (95% CI)", "p (6-parameter)"),
                     with = FALSE], by = rep_key, all.x = TRUE)
}
setorder(t2, biotype, .fdr)
t2[, .fdr := NULL]
t2[biotype == "mRNA", biotype := "protein-coding"]
setnames(t2, c("biotype", "module"), c("Biotype", "Module"))
setcolorder(t2, c("Biotype", "Module", "Genes", "n", "Events"))
add("**Table 2. Modules associated with overall survival in discovery and",
    "their external replication.** Discovery rows are the modules with",
    sprintf("within-network FDR < %s under the principal specification (module", FDR_ALPHA),
    "eigengene + age, sex, T category, nodal status, M1, ordinal grade, ESTIMATE",
    "stromal and immune scores and the three STAR library-quality metrics), with",
    "the joint FDR over all modules of both networks and the estimate under the",
    "clinical comparator alone (age, sex, T, N, M1, grade). Hazard ratios are",
    "per 1 SD of the module eigengene. In CPTAC-3, module membership and",
    "eigengene loadings were fixed in discovery and scores standardised within",
    "CPTAC-3; the 12-parameter model uses the principal covariate set and the",
    "6-parameter model the module eigengene + age, sex, T, M1 and grade",
    sprintf("(%s).", if (has(rep12, "n")) paste0(unique(rep12$n), " patients, ",
                                                 unique(rep12$events), " events") else "CPTAC-3"),
    "", md(t2), "")

# ---- model performance ----
# (a) discrimination of the locked models, every comparator and standardisation
ci <- R_("09_validation_cindex.tsv")
t3a <- data.table(Cohort = ci$cohort,
                  Comparator = col_or(ci, "comparator", "clinical"),
                  Model = ci$model,
                  Standardisation = col_or(ci, "standardisation", "cohort"),
                  n = ci$n, Events = ci$events,
                  `C-index (95% CI)` = est_ci(ci$C, ci$lo, ci$hi, 3))

# (b) incremental value of the module eigengenes over each comparator
rows <- list()
b_row <- function(estimate, comparator, delta, lo = NA, hi = NA, p = NA,
                  spread = NA_character_, n = NA, events = NA)
  data.table(Estimate = estimate, Comparator = comparator,
             `Delta C` = fmt(delta, 4), `95% CI` = ci_txt(lo, hi, 4, " to "),
             p = p_txt(p), `Spread across repeats` = spread, n = n, Events = events)

# (i) apparent increment in discovery, with the stage 09 bootstrap interval
app <- ci[grepl("discovery|KIRC", cohort)]
app <- keep_if(app, "standardisation", "cohort")
if (!has(app, "comparator")) app[, comparator := "clinical"]
dv <- Rx("09_delta_cindex_validation.tsv")
if (!is.null(dv) && !has(dv, "comparator")) dv[, comparator := "clinical"]
if (!is.null(dv)) dv <- keep_if(dv, "standardisation", "cohort")
for (cmp in unique(app$comparator)) {
  a  <- app[comparator == cmp]
  c0 <- a[model == "comparator only", C]
  c1 <- a[grepl("module", model), C]
  if (length(c0) != 1L || length(c1) != 1L) next
  d <- if (!is.null(dv)) dv[grepl("discovery|KIRC", cohort) & comparator == cmp] else dv[0]
  rows[[length(rows) + 1L]] <- b_row(
    "Apparent, discovery (locked model scored in-sample; paired bootstrap)", cmp,
    c1 - c0,
    lo = if (nrow(d) == 1L) d$lo else NA, hi = if (nrow(d) == 1L) d$hi else NA,
    p  = if (nrow(d) == 1L) d$p_boot else NA, n = a$n[1], events = a$events[1])
}

# (ii) eigengenes and (iii) hub lncRNAs in repeated 10 x 10 cross-validation
dl <- R_("05_delta_cindex.tsv")
cv_row <- function(cmpn, label, cmp) {
  r <- dl[comparison == cmpn]
  if (!nrow(r)) return(NULL)
  spread <- if (all(has(r, c("range_lo", "range_hi"))))
    paste0("min-max ", fmt(r$range_lo, 4), " to ", fmt(r$range_hi, 4),
           if (all(has(r, c("n_repeats_positive", "n_repeats"))))
             paste0(" (", r$n_repeats_positive, "/", r$n_repeats, " repeats > 0)") else "")
    else NA_character_
  b_row(label, cmp, r$delta_mean,
        lo = col_or(r, "boot_lo"), hi = col_or(r, "boot_hi"), p = col_or(r, "boot_p"),
        spread = spread, n = col_or(r, "n"), events = col_or(r, "events"))
}
rows <- c(rows, list(
  cv_row("clinical_eig - clinical",
         sprintf("Repeated %d x %d CV, module eigengenes (modules fixed on full cohort; paired bootstrap of the repeat-averaged predictor)",
                 ML_N_REPEATS, ML_N_FOLDS), "clinical"),
  cv_row("augmented_eig - augmented",
         sprintf("Repeated %d x %d CV, module eigengenes (modules fixed on full cohort; paired bootstrap of the repeat-averaged predictor)",
                 ML_N_REPEATS, ML_N_FOLDS), "augmented"),
  cv_row("clinical_hub - clinical",
         sprintf("Repeated %d x %d CV, hub lncRNAs selected inside each fold", ML_N_REPEATS, ML_N_FOLDS),
         "clinical"),
  cv_row("augmented_hub - augmented",
         sprintf("Repeated %d x %d CV, hub lncRNAs selected inside each fold", ML_N_REPEATS, ML_N_FOLDS),
         "augmented")))

# (iv) networks and loadings rebuilt inside each fold
fw <- R_("12_foldwise_module_sensitivity.tsv")
fwr <- fw[grepl("recomputed", analysis)]
if (!has(fwr, "comparator")) fwr[, comparator := "clinical"]
for (i in seq_len(nrow(fwr))) {
  r <- fwr[i]
  rows[[length(rows) + 1L]] <- b_row(
    paste0("Internal CV, networks and loadings rebuilt inside each fold",
           if (has(r, "n_assignments")) sprintf(" (%d folds x %d assignments)", FOLDWISE_K, r$n_assignments) else ""),
    r$comparator, r$delta_C,
    lo = col_or(r, "delta_lo"), hi = col_or(r, "delta_hi"),
    spread = if (has(r, "delta_sd")) paste0("SD across assignments ", fmt(r$delta_sd, 4)) else NA_character_)
}

# (v) external validation, locked model, both standardisations
if (!is.null(dv)) {
  dvv <- R_("09_delta_cindex_validation.tsv")
  if (!has(dvv, "comparator")) dvv[, comparator := "clinical"]
  if (!has(dvv, "standardisation")) dvv[, standardisation := "cohort"]
  ext <- dvv[grepl("CPTAC|validation", cohort)]
  for (i in seq_len(nrow(ext))) {
    r <- ext[i]
    rows[[length(rows) + 1L]] <- b_row(
      paste0("External validation (CPTAC-3), locked model, ", r$standardisation,
             "-standardised scores (paired bootstrap)"),
      r$comparator, r$delta_C, lo = r$lo, hi = r$hi, p = r$p_boot,
      n = col_or(r, "n"), events = col_or(r, "events"))
  }
}
t3b <- rbindlist(rows)

# (c) Brier score, clinical comparator, cohort-standardised scores
br <- R_("12_brier_scores.tsv")
br <- keep_if(br, "comparator", "clinical")
br <- keep_if(br, "standardisation", "cohort")
t3c <- br[, .(Cohort = cohort, Years = years, `Clinical model` = brier_clinical,
              `Clinical + modules` = brier_combined, Null = brier_null)]

# (d) calibration of the clinical + modules model in CPTAC-3 at both horizons
cs <- Rx("12_calibration_summary.tsv")
t3d <- NULL
if (!is.null(cs)) {
  cs <- keep_if(cs, "comparator", "clinical")
  cs <- keep_if(cs, "standardisation", "cohort")
  if (has(cs, "model")) cs <- cs[grepl("module", model)]
  cs <- cs[grepl("CPTAC|validation", cohort) & years %in% c(PRIMARY_HORIZON_YR, SECONDARY_HORIZON_YR)]
  if (nrow(cs))
    t3d <- cs[, .(Cohort = cohort, Years = years,
                  `Calibration slope (95% CI)` = est_ci(slope, slope_lo, slope_hi, 2),
                  `Observed risk` = fmt(observed, 3), `Expected risk` = fmt(expected, 3),
                  `O/E` = fmt(OE, 2), `n at risk` = n_at_risk)]
}

add("**Table 3. Model performance.** (a) Discrimination of the locked models",
    "in discovery and CPTAC-3 for each comparator, with module scores",
    "standardised within the scored cohort (primary) or with the discovery",
    "constants. (b) Increment in concordance from adding the module eigengenes",
    "(or fold-selected hub lncRNAs) to each comparator: the apparent in-sample",
    "increment, repeated cross-validation with the modules fixed on the full",
    "cohort, cross-validation with both networks and all loadings rebuilt",
    "inside each training fold, and external validation with the locked model.",
    "Confidence intervals are paired patient bootstraps",
    sprintf("(%d resamples) unless stated; the spread column gives the min-max", BOOT_B),
    "across cross-validation repeats or the SD across fold assignments.",
    "(c) Brier score of the clinical comparator, lower being better; the null",
    "model assigns every patient the Kaplan-Meier risk.",
    if (!is.null(t3d)) sprintf("(d) Calibration of the clinical + modules model in CPTAC-3 at %d and %d years, absolute risk from the discovery baseline.",
                               PRIMARY_HORIZON_YR, SECONDARY_HORIZON_YR) else NULL,
    "", "*(a) Discrimination*", "", md(t3a), "",
    "*(b) Incremental value*", "", md(t3b), "",
    "*(c) Brier score*", "", md(t3c), "")
if (!is.null(t3d)) add("*(d) Calibration in CPTAC-3*", "", md(t3d), "") else
  add("*(d) Calibration in CPTAC-3: 12_calibration_summary.tsv not yet written; re-run 14 after stage 12.*", "")

writeLines(out, file.path(RESULTS_DIR, "14_manuscript_tables.md"))
msg("Wrote results/14_manuscript_tables.md")
write_session_info("14_tables")
banner("14 | done")
