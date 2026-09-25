# 31_supplementary_data_and_index.R: case list with file identifiers and results index
# Records which expression files entered each cohort and which script wrote each
# results file. Feeds no other analysis. Run last.
#   1. Case list (Supplementary Data 1): one row per STAR-Counts file of
#      TCGA-KIRC, CPTAC-3, TCGA-KIRP and TCGA-KICH, with GDC metadata, a local
#      md5 check, retention and the exclusion reason reconstructed from 01, 08
#      and 11, plus a column dictionary.
#   2. Results index: every writer call in src/R, found in the parse tree and
#      matched to files in results/ and figures/, with a currency test.
#   3. write_session_info() coverage and installed package versions.
# Writes: SupplementaryData1_*.tsv, 31_*.tsv, and 31_NOTE_gdc_unreachable.txt
#   when a GDC request failed. Caches 31_*.rds (delete them to force a refetch).

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(httr)
})
banner("31 | Supplementary Data 1 case list and results index")
set.seed(SEED)
t_start <- Sys.time()

GDC_API    <- "https://api.gdc.cancer.gov/"
GDC_BATCH  <- 200L                       # file_ids per /files request
# Fields requested per file. /files has no per-file data_release field, so the
# release is read once from /status. Changing this string invalidates the
# metadata cache.
GDC_FILE_FIELDS <- paste("file_id", "file_name", "md5sum", "file_size", "state",
                         "analysis.workflow_version", "created_datetime",
                         "updated_datetime", sep = ",")
VALID_ROOT <- file.path(DOWNLOAD_DIR, "validation_data")
SUBTYPES   <- c("TCGA-KIRP", "TCGA-KICH")
COHORT_ORDER <- c("TCGA-KIRC", "CPTAC-3", "TCGA-KIRP", "TCGA-KICH")
REASON_NOT_SEPARABLE <- "excluded by the cohort build (reason not separable)"
REASON_ORDER <- c("non-primary sample", "file not present locally", "failed library",
                  "duplicate aliquot", "no clinical record", "non-clear-cell histology",
                  "non-positive or missing survival time", "missing age or sex",
                  REASON_NOT_SEPARABLE)
# Written where a column does not apply to a cohort, to keep it distinct from FALSE.
NOT_APPLICABLE <- "not applicable"
# Two writes of one script more than this many minutes apart belong to different runs.
# It must exceed the longest silence inside one run (30_endpoint_sensitivity.R writes
# nothing for about 227 minutes during its fold-wise rebuild) and fall short of the
# smallest gap between two runs of a stage. Both margins are measured in section 5e
# (31_run_gap_diagnostics.tsv). Runs less than RUN_GAP_MIN apart are reported as one.
RUN_GAP_MIN <- 330
# A stage whose most recent run predates the pipeline's most recent run by more
# than this was not re-run with the pipeline. Coarse by design: it separates
# pipeline runs, not stages within one run.
PIPELINE_GAP_MIN <- 7L * 24L * 60L
# Status of a claiming script, in words. Retired scripts (in a subdirectory of
# the script directory) are exempt from the pipeline test, so they get their
# own status and an empty logical.
RUN_STATUS_IN         <- "in the pipeline's most recent run"
RUN_STATUS_NOT_RERUN  <- "live script not re-run with the pipeline's most recent run"
RUN_STATUS_SUPERSEDED <- paste("superseded: moved to analysis/R/_superseded and part of no",
                               "pipeline run, so the pipeline-currency test does not apply")
RUN_STATUS_NO_OUTPUT  <- "no output of this script is on disk, so neither currency test was applied"
# Whether the merge-side margin exists for a script at all (section 5e).
MERGE_MEASURED       <- "measured"
MERGE_NOT_MEASURABLE <- paste("not measurable: every matched output of this script was written",
                              "in one run, so there is no earlier run to be separated from")
gdc_errors <- character(0)

# ---- guarded GDC access -----------------------------------------------------
gdc_post <- function(endpoint, body, timeout = 300) {
  r <- httr::POST(paste0(GDC_API, endpoint),
                  body = jsonlite::toJSON(body, auto_unbox = TRUE),
                  httr::content_type_json(), httr::timeout(timeout))
  httr::stop_for_status(r)
  jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"),
                     simplifyDataFrame = TRUE)
}
# Evaluate a GDC request. On any error, record the message and return NULL so the
# stage completes with empty columns instead of aborting.
gdc_try <- function(what, expr) {
  out <- tryCatch(expr, error = function(e) e)
  if (inherits(out, "error")) {
    gdc_errors <<- c(gdc_errors, paste0(what, ": ", conditionMessage(out)))
    msg("  GDC request FAILED (", what, "): ", conditionMessage(out))
    return(NULL)
  }
  out
}
# /files manifest hits -> the file map used by 01/08/11
hits_to_file_map <- function(hits) {
  rbindlist(lapply(seq_len(nrow(hits)), function(i) {
    cs <- hits$cases[[i]]; sm <- cs$samples[[1]]
    data.table(file_id = hits$file_id[i], file_name = hits$file_name[i],
               patient = cs$submitter_id[1], sample_barcode = sm$submitter_id[1],
               sample_type = sm$sample_type[1])
  }))
}
# /cases hits -> survival, age and sex, with the rules of 11.
parse_cases <- function(h) {
  n <- length(h$submitter_id)
  dx1 <- function(l, c_) {
    if (is.null(l)) return(rep(NA_character_, n))
    vapply(l, function(d) {
      if (is.null(d) || !is.data.frame(d) || !c_ %in% names(d) || nrow(d) == 0)
        return(NA_character_)
      v <- as.character(d[[c_]])
      v <- v[!is.na(v) & nzchar(v) & !tolower(v) %in% c("not reported", "unknown")]
      if (!length(v)) NA_character_ else v[1] }, character(1))
  }
  dxm <- function(l, c_) {
    if (is.null(l)) return(rep(NA_real_, n))
    vapply(l, function(d) {
      if (is.null(d) || !is.data.frame(d) || !c_ %in% names(d) || nrow(d) == 0)
        return(NA_real_)
      v <- suppressWarnings(as.numeric(d[[c_]]))
      if (!any(is.finite(v))) NA_real_ else max(v, na.rm = TRUE) }, numeric(1))
  }
  gt <- function(df, c_) if (!is.null(df) && c_ %in% names(df)) df[[c_]] else rep(NA, n)
  g1 <- as.character(gt(h$demographic, "gender"))
  g2 <- as.character(gt(h$demographic, "sex_at_birth"))
  cl <- data.table(
    patient  = h$submitter_id,
    sex_chr  = ifelse(!is.na(g1) & nzchar(g1), g1, g2),
    vital    = as.character(gt(h$demographic, "vital_status")),
    dtd      = suppressWarnings(as.numeric(gt(h$demographic, "days_to_death"))),
    dtlf1    = dxm(h$diagnoses, "days_to_last_follow_up"),
    dtlf2    = dxm(h$follow_ups, "days_to_follow_up"),
    age_days = dxm(h$diagnoses, "age_at_diagnosis"))
  cl[, `:=`(age = age_days / 365.25,
            sex = factor(tolower(sex_chr), levels = c("female", "male")),
            os_event = as.integer(vital == "Dead"))]
  cl[, os_time := ifelse(!is.na(dtd), dtd, pmax(dtlf1, dtlf2, na.rm = TRUE))]
  cl[!is.na(os_time) & os_time > OS_CENSOR_DAYS, os_event := 0L]
  cl[!is.na(os_time), os_time := pmin(os_time, OS_CENSOR_DAYS)]
  cl
}
# Clear-cell histology rule of 08.
is_clear_cell <- function(h) {
  h0   <- tolower(trimws(as.character(h)))
  excl <- grepl("papillary|chromophobe|collecting duct|leiomyomatosis|hlrcc", h0) |
          (grepl("sarcomatoid", h0) & !grepl("clear cell", h0))
  keep <- !is.na(h) & (h %in% c("Renal cell carcinoma, NOS",
                                "Clear cell adenocarcinoma, NOS") |
                       grepl("clear cell", h0))
  keep & !excl
}
# Every file of a project must be either retained or carry exactly one reason.
check_consistency <- function(D, label) {
  stopifnot(all(is.na(D$reason) == D$retained))
  msg("  ", label, ": ", nrow(D), " files, ", sum(D$retained), " retained; ",
      paste(sprintf("%s = %d", names(table(D$reason)), as.integer(table(D$reason))),
            collapse = "; "))
}
compare_consort <- function(D, consort, map, label) {
  for (r in names(map)) {
    got <- sum(D$reason %in% r); exp <- consort[step == map[[r]], n]
    if (length(exp) && got != exp)
      warning(label, ": reconstructed '", r, "' = ", got, " but ", map[[r]],
              " = ", exp, " in the consort table")
  }
}

# ---- 1. Files per cohort and the consort logic that excluded them -----------
banner("1 | Files per cohort and the consort logic")
L_locked <- readRDS(file.path(CACHE_DIR, "locked_model.rds"))

# ---- TCGA-KIRC (01) ----------------------------------------------------------
ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
kco  <- as.data.table(ds$cohort_full); kadj <- as.data.table(ds$cohort_adj)
ksi  <- as.data.table(readRDS(file.path(CACHE_DIR, "expr_raw.rds"))$sample_info)
kfm  <- as.data.table(readRDS(file.path(CACHE_DIR, "gdc_file_map.rds")))
kcl  <- as.data.table(readRDS(file.path(CACHE_DIR, "gdc_clinical.rds")))
stopifnot(setequal(kfm$file_id, ksi$file_id), !anyDuplicated(kfm$file_id))
K <- merge(kfm[, .(file_id, file_name, patient, sample_barcode, sample_type)],
           ksi[, .(file_id, libsize, pct_noFeature)], by = "file_id")
K[, local_path := file.path(EXPR_DIR, file_id, file_name)]
K[, reason := NA_character_]
# 01 section 3, in order: sample type, library failure (per file, before the
# one-per-patient rule), largest library per patient, then section 6: clinical
# record present, positive survival time.
K[sample_type != "Primary Tumor", reason := "non-primary sample"]
K[is.na(reason) & (libsize < MIN_ASSIGNED_READS | pct_noFeature > MAX_NOFEATURE_PCT),
  reason := "failed library"]
cand <- K[is.na(reason)]; setorder(cand, patient, -libsize)
K[file_id %in% cand[duplicated(patient), file_id], reason := "duplicate aliquot"]
K[is.na(reason) & !(patient %in% kcl[!is.na(vital_status), patient]),
  reason := "no clinical record"]
K[is.na(reason) & !(file_id %in% kco$file_id),
  reason := "non-positive or missing survival time"]
K[, retained := file_id %in% kco$file_id]
K[, in_adjusted_models := retained & sample_barcode %in% kadj$sample_barcode]
K[, in_prediction_set  := retained & sample_barcode %in% L_locked$tsamp]
K[, `:=`(cohort = "TCGA-KIRC", role = "discovery")]
check_consistency(K, "TCGA-KIRC")
kcons <- fread(file.path(RESULTS_DIR, "01_consort_flow.tsv"))
stopifnot(sum(K$retained) == kcons[step == "cohort_survival_analysis", n])
compare_consort(K, kcons, list(
  "non-primary sample" = "excluded_not_primary_tumour",
  "failed library" = "excluded_failed_library",
  "duplicate aliquot" = "excluded_duplicate_patient_samples",
  "no clinical record" = "excluded_no_clinical_record",
  "non-positive or missing survival time" = "excluded_nonpositive_os_time"), "TCGA-KIRC")

# ---- CPTAC-3 (08) ------------------------------------------------------------
vd   <- readRDS(file.path(CACHE_DIR, "validation_dataset.rds"))
vco  <- as.data.table(vd$cohort)
vsi  <- as.data.table(readRDS(file.path(CACHE_DIR, "valid_expr_raw.rds"))$sample_info)
vfm  <- as.data.table(readRDS(file.path(CACHE_DIR, "valid_file_map.rds")))
vcl  <- as.data.table(readRDS(file.path(CACHE_DIR, "valid_clinical.rds")))
stopifnot(setequal(vfm$file_id, vsi$file_id), !anyDuplicated(vfm$file_id),
          !anyDuplicated(vcl$patient))
V <- merge(vfm[, .(file_id, file_name, patient, sample_barcode, sample_type)],
           vsi[, .(file_id, libsize)], by = "file_id")
V[, local_path := file.path(VALID_ROOT, "CPTAC-3", file_id, file_name)]
V[, reason := NA_character_]
# 08 section 5, in order: primary tumour, largest library per patient, clinical
# record (inner merge), clear-cell histology, usable survival.
V[sample_type != "Primary Tumor", reason := "non-primary sample"]
cand <- V[is.na(reason)]; setorder(cand, patient, -libsize)
V[file_id %in% cand[duplicated(patient), file_id], reason := "duplicate aliquot"]
V[is.na(reason) & !(patient %in% vcl$patient), reason := "no clinical record"]
V[, histology := vcl$histology[match(patient, vcl$patient)]]
V[is.na(reason) & !is_clear_cell(histology), reason := "non-clear-cell histology"]
V[is.na(reason) & !(file_id %in% vco$file_id),
  reason := "non-positive or missing survival time"]
V[, retained := file_id %in% vco$file_id]
vcc <- vco[complete.cases(vco[, .(age, sex, T_stage, N_pos, M1, grade_num)]), sample_barcode]
V[, in_adjusted_models := retained & sample_barcode %in% vcc]
V[, in_prediction_set  := retained & sample_barcode %in% L_locked$vsamp]
V[, `:=`(cohort = "CPTAC-3", role = "validation")]
check_consistency(V, "CPTAC-3")
vcons <- fread(file.path(RESULTS_DIR, "08_validation_consort.tsv"))
stopifnot(sum(V$retained) == vcons[step == "usable_survival", n])
if (sum(V$reason %in% "non-clear-cell histology") != vcons[step == "excluded_non_clear_cell", n] ||
    sum(V$reason %in% "duplicate aliquot") != vcons[step == "dropped_duplicate_aliquots", n] ||
    sum(V$reason %in% "non-primary sample") != vcons[step == "files_total", n] - vcons[step == "primary_tumour_samples", n])
  warning("CPTAC-3: reconstructed exclusion counts differ from 08_validation_consort.tsv")

# ---- TCGA-KIRP and TCGA-KICH (11) --------------------------------------------
build_subtype_files <- function(proj) {
  tag <- sub("TCGA-", "", proj); dir <- file.path(VALID_ROOT, proj)
  sobj <- readRDS(file.path(CACHE_DIR, paste0("subtype_", tag, ".rds")))
  co   <- as.data.table(sobj$cohort)
  stopifnot(!anyDuplicated(co$libsize))     # library size identifies the analysed file

  # (a) file manifest, all sample types, fetched as 11 does (cached under 31_)
  fm_f <- file.path(CACHE_DIR, paste0("31_subtype_file_map_", tag, ".rds"))
  fm <- if (file.exists(fm_f)) readRDS(fm_f) else NULL
  if (is.null(fm)) {
    msg(proj, ": querying GDC /files for the STAR-Counts manifest ...")
    hits <- gdc_try(paste0(proj, " files manifest"), {
      filt <- list(op = "and", content = list(
        list(op = "in", content = list(field = "cases.project.project_id", value = list(proj))),
        list(op = "in", content = list(field = "data_category", value = list("Transcriptome Profiling"))),
        list(op = "in", content = list(field = "data_type", value = list("Gene Expression Quantification"))),
        list(op = "in", content = list(field = "analysis.workflow_type", value = list("STAR - Counts"))),
        list(op = "in", content = list(field = "access", value = list("open")))))
      j <- gdc_post("files", list(filters = filt, format = "JSON", size = "5000",
                                  fields = paste("file_id", "file_name", "cases.submitter_id",
                                                 "cases.samples.submitter_id",
                                                 "cases.samples.sample_type", sep = ",")))
      if (is.null(j$data$hits) || !nrow(j$data$hits)) stop("no files returned")
      j$data$hits })
    if (!is.null(hits)) { fm <- hits_to_file_map(hits); saveRDS(fm, fm_f) }
  }
  fm_known <- !is.null(fm)

  # local files: 11 downloaded primary tumours only, one directory per file
  loc <- data.table(file_id = list.dirs(dir, full.names = FALSE, recursive = FALSE))
  loc[, file_name := vapply(file_id, function(d) {
    f <- list.files(file.path(dir, d), pattern = "\\.tsv$"); if (length(f)) f[1] else NA_character_
  }, character(1))]
  loc <- loc[!is.na(file_name)]
  if (!fm_known) {
    msg(proj, ": file map unavailable; falling back to the local listing ",
        "(primary tumours only were downloaded by 11)")
    fm <- loc[, .(file_id, file_name, patient = NA_character_,
                  sample_barcode = NA_character_, sample_type = "Primary Tumor")]
  }
  S <- merge(as.data.table(fm)[, .(file_id, file_name, patient, sample_barcode, sample_type)],
             loc[, .(file_id, on_disk = TRUE)], by = "file_id", all.x = TRUE)
  S[is.na(on_disk), on_disk := FALSE]
  S[, local_path := file.path(dir, file_id, file_name)]

  # (b) STAR summary of every local file (cached): assigned reads is the
  #     libsize 11 stored, so it identifies the analysed aliquot
  ss_f <- file.path(CACHE_DIR, paste0("31_subtype_star_summary_", tag, ".rds"))
  ss <- if (file.exists(ss_f)) as.data.table(readRDS(ss_f)) else data.table(file_id = character(0))
  need <- S[on_disk == TRUE & !(file_id %in% ss$file_id)]
  if (nrow(need)) {
    msg(proj, ": reading STAR summaries for ", nrow(need), " local files ...")
    add <- rbindlist(lapply(seq_len(nrow(need)), function(i)
      cbind(data.table(file_id = need$file_id[i]),
            read_star_summary(winlong(need$local_path[i])))))
    ss <- rbindlist(list(ss, add), fill = TRUE); saveRDS(ss, ss_f)
  }
  stopifnot("assigned_reads" %in% names(ss))
  S[, libsize := ss$assigned_reads[match(file_id, ss$file_id)]]

  # (c) retained = the local file whose assigned reads equal the cached libsize
  m <- match(S$libsize, co$libsize)
  S[, retained := on_disk & !is.na(m)]
  stopifnot(sum(S$retained) == nrow(co), !anyDuplicated(m[!is.na(m)]))
  if (fm_known) {
    stopifnot(identical(S$sample_barcode[S$retained], co$sample_barcode[m[S$retained]]))
  } else {
    # `m` is indexed by row of S, so select by full row numbers. An i-subset
    # of the retained column would recycle m to nrow(S).
    ret <- which(S$retained)
    S[ret, `:=`(sample_barcode = co$sample_barcode[m[ret]],
                patient        = co$patient[m[ret]])]
  }

  # (d) clinical records as 11 read them (cached under 31_)
  cs_f <- file.path(CACHE_DIR, paste0("31_subtype_cases_", tag, ".rds"))
  cl <- if (file.exists(cs_f)) readRDS(cs_f) else NULL
  if (is.null(cl)) {
    msg(proj, ": querying GDC /cases ...")
    h <- gdc_try(paste0(proj, " cases"), {
      j <- gdc_post("cases", list(
        filters = list(op = "in", content = list(field = "project.project_id", value = list(proj))),
        fields = "submitter_id", expand = "demographic,diagnoses,follow_ups",
        format = "JSON", size = "3000"))
      if (is.null(j$data$hits) || !length(j$data$hits$submitter_id)) stop("no cases returned")
      j$data$hits })
    if (!is.null(h)) { cl <- parse_cases(h); saveRDS(cl, cs_f) }
  }

  # (e) exclusion reasons in the order of 11's build_cohort()
  S[, reason := NA_character_]
  S[sample_type != "Primary Tumor", reason := "non-primary sample"]
  S[is.na(reason) & !on_disk, reason := "file not present locally"]
  if (fm_known) {
    cand <- S[is.na(reason)]; setorder(cand, patient, -libsize)
    S[file_id %in% cand[duplicated(patient), file_id], reason := "duplicate aliquot"]
    if (!is.null(cl)) {
      S[is.na(reason) & !(patient %in% cl$patient), reason := "no clinical record"]
      S[, c("os_time", "os_event", "age", "sex") :=
          cl[match(S$patient, cl$patient), .(os_time, os_event, age, as.character(sex))]]
      S[is.na(reason) & (is.na(os_time) | os_time <= 0 | is.na(os_event)),
        reason := "non-positive or missing survival time"]
      S[is.na(reason) & (is.na(age) | is.na(sex)), reason := "missing age or sex"]
    }
  }
  S[is.na(reason) & !retained, reason := REASON_NOT_SEPARABLE]
  S[, in_adjusted_models := FALSE]
  cc <- co[complete.cases(co[, .(age, sex, T_stage, M1)]), sample_barcode]  # 11 covariate set
  S[retained == TRUE, in_adjusted_models := sample_barcode %in% cc]
  S[, in_prediction_set := NA]
  S[, `:=`(cohort = proj, role = "subtype")]
  check_consistency(S, proj)
  S
}
sub_tables <- lapply(SUBTYPES, build_subtype_files)
sub_summary <- fread(file.path(RESULTS_DIR, "11_subtype_cohort_summary.tsv"))
for (S in sub_tables) {
  ex <- sub_summary[cohort == S$cohort[1]]
  if (nrow(ex) && (sum(S$retained) != ex$n_scored || sum(S$in_adjusted_models) != ex$n_complete))
    warning(S$cohort[1], ": retained / complete counts differ from 11_subtype_cohort_summary.tsv")
}

COLS <- c("cohort", "role", "patient", "sample_barcode", "sample_type", "file_id",
          "file_name", "local_path", "libsize", "retained", "reason",
          "in_adjusted_models", "in_prediction_set")
ALL <- rbindlist(c(list(K[, ..COLS], V[, ..COLS]), lapply(sub_tables, function(S) S[, ..COLS])))
stopifnot(!anyDuplicated(ALL$file_id))
msg("Files listed: ", nrow(ALL), " (", sum(ALL$retained), " retained)")

# ---- 2. GDC md5sum, size and data release for every file --------------------
banner("2 | GDC file metadata")
meta_f <- file.path(CACHE_DIR, "31_gdc_file_metadata.rds")
meta0 <- list(meta = data.table(file_id = character(0)), query_release = NA_character_,
              fetched = NA, fields = GDC_FILE_FIELDS)
meta <- if (file.exists(meta_f)) readRDS(meta_f) else meta0
# A cache written with a different field set is discarded rather than merged.
if (!identical(meta$fields, GDC_FILE_FIELDS)) {
  msg("Cached GDC metadata was fetched with a different field set; refetching")
  meta <- meta0
}
todo <- setdiff(ALL$file_id, meta$meta$file_id)
if (length(todo)) {
  msg("Fetching GDC metadata for ", length(todo), " files in batches of ", GDC_BATCH, " ...")
  rel <- gdc_try("status", {
    r <- httr::GET(paste0(GDC_API, "status"), httr::timeout(60)); httr::stop_for_status(r)
    jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"))$data_release })
  if (!is.null(rel)) msg("  GDC data release at query time: ", rel)
  chunks <- split(todo, ceiling(seq_along(todo) / GDC_BATCH)); got <- list(); warned <- FALSE
  for (i in seq_along(chunks)) {
    ids <- chunks[[i]]
    j <- gdc_try(paste0("files batch ", i, "/", length(chunks)), {
      j <- gdc_post("files", list(
        filters = list(op = "in", content = list(field = "file_id", value = as.list(ids))),
        fields  = GDC_FILE_FIELDS,
        format  = "JSON", size = as.character(length(ids))))
      if (is.null(j$data$hits) || !nrow(j$data$hits)) stop("no hits returned")
      j })
    if (is.null(j)) break
    if (!warned && length(j$warnings)) {
      msg("  GDC warning (once): ", paste(unlist(j$warnings), collapse = "; ")); warned <- TRUE
    }
    h <- j$data$hits
    # analysis.workflow_version comes back as a nested data.frame column.
    wv <- if (is.data.frame(h$analysis) && "workflow_version" %in% names(h$analysis))
            as.character(h$analysis$workflow_version) else rep(NA_character_, nrow(h))
    got[[i]] <- cbind(as.data.table(h[, setdiff(names(h), "analysis"), drop = FALSE]),
                      workflow_version = wv)
    msg("  batch ", i, "/", length(chunks), ": ", nrow(got[[i]]), " of ", length(ids), " files")
  }
  if (length(got)) {
    new <- rbindlist(got, fill = TRUE)
    for (cc in c("file_name", "md5sum", "file_size", "state", "workflow_version",
                 "created_datetime", "updated_datetime"))
      if (!cc %in% names(new)) new[[cc]] <- NA
    new <- new[, .(file_id, gdc_file_name = as.character(file_name), gdc_md5sum = as.character(md5sum),
                   gdc_file_size = as.numeric(file_size), gdc_state = as.character(state),
                   gdc_workflow_version = as.character(workflow_version),
                   gdc_created_datetime = as.character(created_datetime),
                   gdc_updated_datetime = as.character(updated_datetime))]
    meta$meta <- rbindlist(list(meta$meta, new), fill = TRUE)[!duplicated(file_id)]
    if (!is.null(rel)) meta$query_release <- rel
    meta$fetched <- Sys.time(); meta$fields <- GDC_FILE_FIELDS; saveRDS(meta, meta_f)
  }
} else msg("GDC metadata for all ", nrow(ALL), " files from cache (fetched ", format(meta$fetched), ")")
ALL <- merge(ALL, meta$meta, by = "file_id", all.x = TRUE, sort = FALSE)
# Guarantee the GDC columns exist even if every request failed and nothing was
# ever cached, so the table below has a fixed shape.
for (cc in c("gdc_file_name", "gdc_md5sum", "gdc_state", "gdc_workflow_version",
             "gdc_created_datetime", "gdc_updated_datetime"))
  if (!cc %in% names(ALL)) ALL[[cc]] <- NA_character_
if (!"gdc_file_size" %in% names(ALL)) ALL[["gdc_file_size"]] <- NA_real_
# One release for the whole fetch (see GDC_FILE_FIELDS above).
ALL[, gdc_data_release := as.character(meta$query_release)]
ALL[, gdc_metadata_fetched := format(meta$fetched, "%Y-%m-%d")]
n_meta <- sum(!is.na(ALL$gdc_md5sum))
msg("GDC md5sum available for ", n_meta, " / ", nrow(ALL), " files")
if (any(!is.na(ALL$gdc_file_name) & ALL$gdc_file_name != ALL$file_name))
  warning("GDC file_name differs from the cached file map for ",
          sum(!is.na(ALL$gdc_file_name) & ALL$gdc_file_name != ALL$file_name), " files")

# ---- 3. Local copies: presence, size, date and md5sum -----------------------
banner("3 | Local copies")
lp <- winlong(ALL$local_path)
ALL[, local_file_present := file.exists(lp)]
ALL[, local_file_size := ifelse(local_file_present, file.size(lp), NA_real_)]
ALL[, local_file_date := ifelse(local_file_present, format(file.mtime(lp), "%Y-%m-%d"), NA_character_)]
md5_f <- file.path(CACHE_DIR, "31_local_md5.rds")
md5c <- if (file.exists(md5_f)) readRDS(md5_f) else
        data.table(local_path = character(0), size = numeric(0), date = character(0), md5 = character(0))
key_all <- paste(ALL$local_path, ALL$local_file_size, ALL$local_file_date)
key_c   <- paste(md5c$local_path, md5c$size, md5c$date)
need <- ALL[local_file_present & !(key_all %in% key_c)]
if (nrow(need)) {
  msg("Computing md5sum of ", nrow(need), " local files ...")
  h <- unname(tools::md5sum(winlong(need$local_path)))
  md5c <- rbindlist(list(md5c, data.table(local_path = need$local_path, size = need$local_file_size,
                                          date = need$local_file_date, md5 = h)))
  md5c <- md5c[!duplicated(paste(local_path, size, date), fromLast = TRUE)]
  saveRDS(md5c, md5_f)
}
ALL[, local_md5sum := md5c$md5[match(key_all, paste(md5c$local_path, md5c$size, md5c$date))]]
ALL[, local_md5_matches_gdc := ifelse(is.na(gdc_md5sum) | is.na(local_md5sum), NA,
                                      local_md5sum == gdc_md5sum)]
msg("Local files present: ", sum(ALL$local_file_present), "; md5 verified against GDC: ",
    sum(ALL$local_md5_matches_gdc %in% TRUE), "; mismatches: ",
    sum(ALL$local_md5_matches_gdc %in% FALSE))
if (any(ALL$local_md5_matches_gdc %in% FALSE))
  warning("md5 mismatch for: ", paste(ALL[local_md5_matches_gdc %in% FALSE, file_id], collapse = ", "))

# ---- GDC note --------------------------------------------------------------
note_file <- file.path(RESULTS_DIR, "31_NOTE_gdc_unreachable.txt")
if (length(gdc_errors)) {
  writeLines(c("31_supplementary_data_and_index.R: one or more GDC API requests failed.",
               "Columns that depend on them are empty for the affected files; re-run to retry",
               "(successful batches are cached in cache/31_*.rds and are not refetched).",
               paste0("Files without GDC md5sum / size: ", sum(is.na(ALL$gdc_md5sum)), " of ", nrow(ALL)),
               "Errors:", paste0("  ", gdc_errors), paste0("Time: ", format(Sys.time()))), note_file)
  msg("Wrote ", basename(note_file))
} else unlink(note_file)

# ---- 4. Supplementary Data 1 and its summary --------------------------------
banner("4 | Supplementary Data 1")
ALL[, cohort := factor(cohort, levels = COHORT_ORDER)]
setorder(ALL, cohort, -retained, patient, sample_barcode, file_id)
# in_prediction_set: membership of the two sample sets in locked_model.rds
# (09_validate.R). tsamp holds the discovery patients every locked model is
# fitted on and vsamp the CPTAC-3 patients they are scored on. Subtype files
# get NOT_APPLICABLE.
SD1 <- ALL[, .(cohort = as.character(cohort), role, patient, sample_barcode, sample_type,
               file_id, file_name,
               gdc_md5sum, gdc_file_size, gdc_data_release, gdc_metadata_fetched,
               gdc_state, gdc_workflow_version, gdc_created_datetime,
               gdc_updated_datetime,
               local_file_present, local_file_size, local_file_date, local_md5sum,
               local_md5_matches_gdc,
               assigned_reads = libsize,
               retained_in_analysed_cohort = retained,
               exclusion_reason = reason,
               in_adjusted_models,
               in_prediction_set = fifelse(role == "subtype", NOT_APPLICABLE,
                                           as.character(in_prediction_set)))]
save_tsv(SD1, "SupplementaryData1_case_list.tsv")

# ---- data dictionary for the case list --------------------------------------
# Every column is defined with the cohorts it applies to. The locked-model definition
# is read from locked_model.rds. Per comparator set, b_clin is the unpenalised Cox model
# of the comparator terms, b_full an elastic net penalising only the module scores in Et,
# and sel_me the retained modules. The assertions check that all fits use tsamp and
# vsamp and contain only comparator terms and module scores.
stopifnot(identical(L_locked$version, "v9"),
          identical(colnames(L_locked$Et), colnames(L_locked$Ev)),
          setequal(colnames(L_locked$Et), names(L_locked$loadings)),
          identical(names(L_locked$models), names(L_locked$X)))
n_me     <- ncol(L_locked$Et)
n_me_pc  <- sum(grepl("^mRNA_", colnames(L_locked$Et)))
n_me_lnc <- sum(grepl("^lnc_",  colnames(L_locked$Et)))
stopifnot(n_me_pc + n_me_lnc == n_me)
CMP <- names(L_locked$models)
cmp_terms <- lapply(CMP, function(cm) colnames(L_locked$X[[cm]]$Xt))
names(cmp_terms) <- CMP
for (cm in CMP) {
  M <- L_locked$models[[cm]]
  stopifnot(nrow(L_locked$X[[cm]]$Xt) == length(L_locked$tsamp),
            nrow(L_locked$X[[cm]]$Xv) == length(L_locked$vsamp),
            setequal(names(M$b_clin), cmp_terms[[cm]]),
            all(M$sel_me %in% colnames(L_locked$Et)),
            # no term outside the comparator set and the module scores
            !length(setdiff(names(M$b_full), c(cmp_terms[[cm]], colnames(L_locked$Et)))))
}
# One phrase per comparator set: its terms and the number of retained module scores.
cmp_phrase <- paste(vapply(CMP, function(cm) sprintf(
  "%s (%s), of which %d module scores are retained", cm,
  paste(cmp_terms[[cm]], collapse = ", "), length(L_locked$models[[cm]]$sel_me)),
  character(1)), collapse = "; ")
LOCKED_DEF <- sprintf(paste0(
  "TRUE if the sample is in the locked prediction set written by 09_validate.R ",
  "(cache/locked_model.rds): the %d TCGA-KIRC patients every model in that file ",
  "is fitted on, or the %d CPTAC-3 patients they are then scored on. The file ",
  "holds %d comparator sets and two Cox fits of each on those same patients -- ",
  "the comparator terms alone, unpenalised, and an elastic net that leaves ",
  "those terms unpenalised and offers the same %d module scores (%d ",
  "protein-coding, %d lncRNA) to the penalty. The comparator sets are: %s. No ",
  "fit contains an individual hub lncRNA gene. Written as '%s' for TCGA-KIRP ",
  "and TCGA-KICH, for which no locked model exists."),
  length(L_locked$tsamp), length(L_locked$vsamp), length(CMP),
  n_me, n_me_pc, n_me_lnc, cmp_phrase, NOT_APPLICABLE)
msg("Locked models of 09_validate.R: ", length(CMP), " comparator sets (",
    paste(CMP, collapse = ", "), "), each fitted twice -- comparator terms alone ",
    "and comparator terms plus an elastic net over ", n_me, " module scores (",
    n_me_pc, " protein-coding, ", n_me_lnc, " lncRNA); modules retained -- ",
    paste(vapply(CMP, function(cm)
            sprintf("%s %d", cm, length(L_locked$models[[cm]]$sel_me)),
          character(1)), collapse = ", "),
    "; all fitted on ", length(L_locked$tsamp), " discovery and scored on ",
    length(L_locked$vsamp), " CPTAC-3 patients")
DEFS <- rbindlist(list(
  list("cohort", "GDC project the expression file belongs to (TCGA-KIRC, CPTAC-3, TCGA-KIRP or TCGA-KICH).", "all"),
  list("role", "Role of that cohort in the study: discovery, validation or subtype.", "all"),
  list("patient", "Submitter (case) identifier of the patient the aliquot came from.", "all"),
  list("sample_barcode", "Submitter identifier of the sample (aliquot) the library was made from.", "all"),
  list("sample_type", "GDC sample type of the aliquot; only Primary Tumor samples can enter a cohort.", "all"),
  list("file_id", "GDC file UUID of the STAR-Counts file.", "all"),
  list("file_name", "GDC file name of the STAR-Counts file.", "all"),
  list("gdc_md5sum", "md5 checksum of the file as recorded by the GDC /files endpoint.", "all"),
  list("gdc_file_size", "Size in bytes of the file as recorded by the GDC /files endpoint.", "all"),
  list("gdc_data_release", "GDC data release current when this metadata was read (/status endpoint). One value for the whole fetch: /files carries no per-file release field.", "all"),
  list("gdc_metadata_fetched", "Date on which this script read the GDC metadata.", "all"),
  list("gdc_state", "GDC file state (released for every file used here).", "all"),
  list("gdc_workflow_version", "Version of the GDC mRNA analysis workflow that produced the file.", "all"),
  list("gdc_created_datetime", "Datetime the file was created in the GDC.", "all"),
  list("gdc_updated_datetime", "Datetime the file record was last updated in the GDC.", "all"),
  list("local_file_present", "TRUE if the file is present in this project's local download tree. Only primary tumours were downloaded for the subtype cohorts, so the non-primary subtype files are absent by design.", "all"),
  list("local_file_size", "Size in bytes of the local copy; empty when the file was never downloaded.", "all"),
  list("local_file_date", "Modification date of the local copy; empty when the file was never downloaded.", "all"),
  list("local_md5sum", "md5 checksum recomputed from the local copy; empty when the file was never downloaded.", "all"),
  list("local_md5_matches_gdc", "TRUE if the local md5 equals the GDC md5; empty when either checksum is unavailable.", "all"),
  list("assigned_reads", "STAR assigned reads: the sum of unstranded counts over ENSG genes, which is the library size used throughout the pipeline. Empty when the file was never downloaded.", "all"),
  list("retained_in_analysed_cohort", "TRUE if this file supplied the expression profile of one patient in the analysed cohort.", "all"),
  list("exclusion_reason", "Why the file was excluded, reconstructed from the cohort build of 01 (TCGA-KIRC), 08 (CPTAC-3) or 11 (subtypes). Empty for retained files.", "all"),
  list("in_adjusted_models", "TRUE if the patient contributed to the covariate-adjusted survival models: complete age, sex, T stage, N, M1 and grade in TCGA-KIRC and CPTAC-3; complete age, sex, T stage and M1 in the subtype cohorts.", "all"),
  list("in_prediction_set", LOCKED_DEF, "TCGA-KIRC and CPTAC-3")))
setnames(DEFS, c("column", "definition", "applicable_cohorts"))
DEFS[, type := vapply(SD1, function(x) class(x)[1], character(1))[column]]
setcolorder(DEFS, c("column", "type", "definition", "applicable_cohorts"))
stopifnot(identical(DEFS$column, names(SD1)))
save_tsv(DEFS, "SupplementaryData1_column_definitions.tsv")

ALL[, reason_lab := fifelse(retained, "(retained)", reason)]
ALL[, reason_lab := factor(reason_lab, levels = c("(retained)", REASON_ORDER))]
# n_patients counts every patient with a file in the group, n_patients_retained the
# analysed count.
SUMM_COUNTS <- quote(.(n_files = .N, n_patients = uniqueN(patient),
                       n_patients_retained = uniqueN(patient[retained]),
                       n_local_files = sum(local_file_present),
                       n_gdc_metadata = sum(!is.na(gdc_md5sum)),
                       n_md5_verified = sum(local_md5_matches_gdc %in% TRUE),
                       n_in_adjusted_models = sum(in_adjusted_models %in% TRUE),
                       n_in_prediction_set = sum(in_prediction_set %in% TRUE)))
summ <- ALL[, eval(SUMM_COUNTS),
            by = .(cohort, role, retained_in_analysed_cohort = retained, exclusion_reason = reason_lab)]
setorder(summ, cohort, -retained_in_analysed_cohort, exclusion_reason)
summ[, `:=`(cohort = as.character(cohort), exclusion_reason = as.character(exclusion_reason))]
summ[exclusion_reason == "(retained)", exclusion_reason := NA_character_]
# Total row across cohorts, which also gives the size of the locked prediction set.
tot <- ALL[, eval(SUMM_COUNTS)]
summ <- rbind(summ, cbind(data.table(cohort = "ALL COHORTS", role = "total",
                                     retained_in_analysed_cohort = NA,
                                     exclusion_reason = NA_character_), tot))
# Blank, not zero, where the locked model does not apply.
summ[role == "subtype", n_in_prediction_set := NA_integer_]
save_tsv(summ, "31_supplementary_data1_summary.tsv")
print(summ)
msg("Locked prediction set: ", sum(ALL$in_prediction_set %in% TRUE), " files (",
    paste(vapply(c("TCGA-KIRC", "CPTAC-3"), function(cc)
      sprintf("%s %d of %d retained", cc, sum(ALL$in_prediction_set %in% TRUE & ALL$cohort == cc),
              sum(ALL$retained & ALL$cohort == cc)), character(1)), collapse = "; "), ")")
stopifnot(sum(ALL$in_prediction_set %in% TRUE) ==
            length(L_locked$tsamp) + length(L_locked$vsamp))

# ---- 5. Results index, session-info coverage and package versions -----------
banner("5 | Results index")
# Writers and the argument that carries the file name (name, else position).
WRITERS <- list(save_tsv = list("file", 2L), save_fig = list("name", 2L),
                save_fig_base = list("name", 1L), fwrite = list("file", 2L),
                write.table = list("file", 2L), write.csv = list("file", 2L),
                writeLines = list("con", 2L), write.xlsx = list("file", 2L),
                saveWorkbook = list("file", 2L), ggsave = list("filename", 1L))
FIG_WRITERS <- c("save_fig", "save_fig_base")
EXT_RE <- "\\.(tsv|txt|csv|xlsx|svg|png|pdf|docx|rds|gct|json|md)$"

# script_files is the live pipeline (top level of the script directory) and is
# used for the session-info and package tables. index_files adds any
# subdirectories, so a file written only by a retired script is attributed.
script_files <- list.files(R_DIR, pattern = "\\.R$", full.names = TRUE)
script_files <- script_files[basename(script_files) != "00_config.R"]
script_files <- script_files[order(!grepl("^[0-9]", basename(script_files)), basename(script_files))]
all_r <- list.files(R_DIR, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
all_r <- all_r[basename(all_r) != "00_config.R"]   # helpers, not a writer of results
index_files <- c(script_files, setdiff(all_r, script_files))
script_dir_of <- function(f) {
  d <- sub(paste0("^", regex_escape(R_DIR), "/?"), "", dirname(normalizePath(f, winslash = "/")))
  ifelse(nzchar(d), d, ".")
}
UNTRACKED_LIVE <- setNames(character(0), character(0))

parse_pd <- function(f) {
  ex <- tryCatch(parse(f, keep.source = TRUE), error = function(e) NULL)
  if (is.null(ex)) return(NULL)
  as.data.table(utils::getParseData(ex, includeText = TRUE))
}
strip_q <- function(s) sub("^(['\"])(.*)\\1$", "\\2", s)
regex_escape <- function(s) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", s)
descendants <- function(pd, id) {
  ids <- id
  repeat { nw <- pd[parent %in% ids & !(id %in% ids), id]; if (!length(nw)) break; ids <- c(ids, nw) }
  ids
}
# Arguments of a call node: positional and named, in source order.
call_args <- function(pd, call_id, fn_expr_id) {
  ch <- pd[parent == call_id][order(line1, col1)]
  args <- list(); pending <- ""
  for (k in seq_len(nrow(ch))) {
    if (ch$token[k] == "SYMBOL_SUB") pending <- ch$text[k]
    else if (ch$token[k] == "expr" && ch$id[k] != fn_expr_id) {
      args[[length(args) + 1]] <- list(name = pending, id = ch$id[k], text = ch$text[k]); pending <- ""
    }
  }
  args
}
pick_arg <- function(args, name, pos) {
  nm <- vapply(args, `[[`, "", "name")
  if (any(nm == name)) return(args[[which(nm == name)[1]]])
  un <- args[nm == ""]
  if (length(un) >= pos) un[[pos]] else NULL
}
# String literals a bare symbol can stand for in the same script: either the
# head of a for() loop it is the variable of (for (p in c("GO.db", ...))) or
# the right-hand side of an assignment (note_file <- file.path(RESULTS_DIR, "...")).
str_const <- function(pd, root) strip_q(
  pd[id %in% descendants(pd, root) & token == "STR_CONST"][order(line1, col1), text])
symbol_strings <- function(pd, sym) {
  for (fid in pd[token == "forcond", id]) {
    ch <- pd[parent == fid]
    if (any(ch$token == "SYMBOL" & ch$text == sym)) {
      s <- str_const(pd, fid)
      if (length(s)) return(s)
    }
  }
  for (p in pd[token == "SYMBOL" & text == sym, parent]) {
    pr <- pd[id == p, parent]; if (!length(pr)) next
    asg <- pd[parent == pr][order(line1, col1)]
    if (nrow(asg) >= 3 && asg$id[1] == p && asg$token[2] == "LEFT_ASSIGN")
      return(str_const(pd, asg$id[3]))
  }
  character(0)
}
# String literals inside an argument (resolving a bare symbol via str_const()), and
# whether the built name is anchored at each end. paste0("03_", lab, "_KM_", m)
# continues past its last literal, so its pattern must end in a wildcard.
arg_literals <- function(pd, arg) {
  kids <- pd[parent == arg$id]
  if (nrow(kids) == 1 && kids$token[1] == "SYMBOL")   # resolved name, complete
    return(list(lits = symbol_strings(pd, kids$text[1]), start = TRUE, end = TRUE))
  tk <- pd[id %in% descendants(pd, arg$id) &
           token %in% c("STR_CONST", "SYMBOL")][order(line1, col1)]
  list(lits  = strip_q(tk[token == "STR_CONST", text]),
       start = nrow(tk) > 0 && tk$token[1] == "STR_CONST",
       end   = nrow(tk) > 0 && tk$token[nrow(tk)] == "STR_CONST")
}
# Strip leading hashes. Banner rules and numbered section headers become empty.
# The hyphen must be last in a bracket expression ("[=\\- ]" is an invalid TRE
# range). Double quotes become single quotes because save_tsv() writes with
# quote = FALSE.
clean_comment <- function(x) {
  x <- trimws(sub("^#+\\s?", "", trimws(x)))
  x <- gsub("^[-= ]+|[-= ]+$", "", x)   # "---- hub gene tables ----" -> the text
  x <- gsub("[\"“”]", "'", x)
  x[grepl("^[= -]*$", x)] <- ""
  x[grepl("^[0-9]+[a-z]?[.)]\\s", x)] <- ""            # "8. Persist"
  x[grepl("^(step|section|part)\\s+[0-9A-Za-z]+[.:)]?\\s*$", x, ignore.case = TRUE)] <- ""
  x
}
# Trim function words and connector punctuation left by cutting a sentence.
DANGLING <- paste0("\\s+(a|an|and|are|as|at|be|but|by|for|from|in|is|it|it is|its|",
                   "of|on|or|that|the|their|then|this|to|was|were|which|with)$")
tidy_description <- function(s, max_chars = 160) {
  if (is.na(s) || !nzchar(s)) return(NA_character_)
  s <- gsub("\\s+", " ", trimws(s))
  # keep the first sentence when the block runs on
  cut <- regexpr("[.;] ", s)
  if (cut > 20) s <- substr(s, 1, cut - 1)
  if (nchar(s) > max_chars)
    s <- sub("\\s+\\S*$", "", substr(s, 1, max_chars))
  # A cut inside a parenthesis ends the text before the parenthesis opened.
  op <- gregexpr("(", s, fixed = TRUE)[[1]]; cl <- gregexpr(")", s, fixed = TRUE)[[1]]
  n_op <- if (op[1] == -1) 0L else length(op); n_cl <- if (cl[1] == -1) 0L else length(cl)
  if (n_op > n_cl) s <- trimws(substr(s, 1, op[n_cl + 1] - 1))
  # "... i.e" / "... e.g" left by cutting at the abbreviation's own full stop.
  s <- sub("[,;:]?\\s*(i\\.e|e\\.g|cf|viz|etc)\\.?$", "", s, ignore.case = TRUE)
  repeat {
    s2 <- sub("[,;:]$", "", sub(DANGLING, "", s, ignore.case = TRUE))
    if (identical(s2, s)) break
    s <- s2
  }
  s <- trimws(sub("[.,;:]$", "", s))
  if (!nzchar(s)) NA_character_ else s
}
# Trailing comment on the call line, else the nearest contiguous comment block
# above it that starts in the same column as the call. Section headers and rules
# are blanked. The source line of the description is returned with it.
nearest_comment <- function(pd, line, col = 1L, window = 15) {
  none <- list(text = NA_character_, line = NA_integer_)
  cm <- pd[token == "COMMENT"]
  same <- cm[line1 == line]
  if (nrow(same)) return(list(text = tidy_description(clean_comment(same$text[1])),
                              line = line))
  cm <- cm[col1 == col]
  prev <- cm[line1 < line & line1 >= line - window][order(-line1)]
  if (!nrow(prev)) return(none)
  ln <- prev$line1[1]; block <- ln
  while ((ln - 1) %in% cm$line1) { ln <- ln - 1; block <- c(ln, block) }
  txt <- clean_comment(cm[line1 %in% block][order(line1), text])
  # Drop leading blanks (a rule line), then read the block as one sentence.
  txt <- txt[cumsum(nzchar(txt)) > 0]
  if (!length(txt)) return(none)
  d <- tidy_description(paste(txt[nzchar(txt)], collapse = " "))
  if (is.na(d)) none else list(text = d, line = block[1])
}
# Where a writer call can land. Only results/ and figures/ are covered by the
# unclaimed listing. The top level of the cache and the manuscript directory
# are also searched, so a call writing an intermediate is not reported missing.
OUT_DIRS <- list(results = list(d = RESULTS_DIR, recursive = TRUE,  unclaimed = TRUE),
                 figures = list(d = FIG_DIR,     recursive = TRUE,  unclaimed = TRUE),
                 cache   = list(d = CACHE_DIR,   recursive = FALSE, unclaimed = FALSE),
                 manuscript = list(d = file.path(ANALYSIS_DIR, "manuscript"),
                                   recursive = FALSE, unclaimed = FALSE))
list_outputs <- function() {
  rbindlist(lapply(names(OUT_DIRS), function(nm) {
    o <- OUT_DIRS[[nm]]
    if (!dir.exists(o$d)) return(NULL)
    f <- list.files(o$d, recursive = o$recursive, full.names = TRUE)
    f <- f[!dir.exists(f)]
    data.table(location = nm, file = basename(f),
               rel = file.path(nm, sub(paste0("^", regex_escape(o$d), "/"), "", f)),
               size = file.size(f), mtime = format(file.mtime(f), "%Y-%m-%d %H:%M"),
               mtime_num = as.numeric(file.mtime(f)),
               unclaimed_scope = o$unclaimed)
  }))
}
# Start of the newest cluster of timestamps: walk down and stop at the first gap
# longer than gap_min minutes. Applied to one script's output mtimes
# (RUN_GAP_MIN) to find its most recent run, and to the run starts of all live
# scripts (PIPELINE_GAP_MIN) to find the pipeline's most recent run.
last_run_start <- function(times, gap_min = RUN_GAP_MIN) {
  t <- sort(unique(times[is.finite(times)]), decreasing = TRUE)
  if (!length(t)) return(NA_real_)
  lo <- t[1]
  for (k in seq_along(t)[-1]) {
    if (lo - t[k] <= gap_min * 60) lo <- t[k] else break
  }
  lo
}
# The two margins of the timing rules for one set of timestamps: the largest
# gap inside the newest cluster (must stay below gap_min) and the separation
# from the write before it (must stay above gap_min). NA where undefined.
run_margins <- function(times, gap_min = RUN_GAP_MIN) {
  t <- sort(unique(times[is.finite(times)]))
  if (!length(t)) return(list(run_start = NA_real_, run_end = NA_real_, n_times = 0L,
                              span_min = NA_real_, max_internal_gap_min = NA_real_,
                              gap_to_previous_min = NA_real_))
  rs <- last_run_start(t, gap_min)
  tr <- t[t >= rs]
  list(run_start = rs, run_end = max(tr), n_times = length(tr),
       span_min = (max(tr) - min(tr)) / 60,
       max_internal_gap_min = if (length(tr) > 1) max(diff(tr)) / 60 else 0,
       gap_to_previous_min = if (any(t < rs)) (rs - max(t[t < rs])) / 60 else NA_real_)
}
fmt_t <- function(x) if (is.na(x)) NA_character_ else
  format(as.POSIXct(x, origin = "1970-01-01"), "%Y-%m-%d %H:%M")
# `sess` maps script name to the mtime of its session-info file, the end of its
# most recent run. Output newer than that record means the script was still
# running when the index was built, or records provenance before its last write
# (as 31 and 30_endpoint_sensitivity.R do). Reported, not treated as superseded output.
build_index <- function(sess = numeric(0)) {
  listing <- list_outputs(); rows <- list(); hits <- list()
  # Time of this listing, written into all three index tables.
  snapshot <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  for (f in index_files) {
    sc <- basename(f); sdir <- script_dir_of(f); pd <- parse_pd(f)
    if (is.null(pd)) {
      # A script that fails to parse is recorded, not skipped.
      rows[[length(rows) + 1]] <- data.table(
        file = NA_character_, script = sc, script_dir = sdir, exists = NA,
        description = "PARSE ERROR (script did not parse when the index was built)",
        writer = NA_character_)
      next
    }
    fn <- pd[token == "SYMBOL_FUNCTION_CALL" & text %in% names(WRITERS)][order(line1, col1)]
    for (k in seq_len(nrow(fn))) {
      w <- fn$text[k]; fn_expr <- fn$parent[k]; call_id <- pd[id == fn_expr, parent]
      line <- pd[id == call_id, line1]; col <- pd[id == call_id, col1]
      nc <- nearest_comment(pd, line, col)
      a <- pick_arg(call_args(pd, call_id, fn_expr), WRITERS[[w]][[1]], WRITERS[[w]][[2]])
      al <- if (is.null(a)) list(lits = character(0), start = FALSE, end = FALSE)
            else arg_literals(pd, a)
      lits <- al$lits
      is_fig <- w %in% FIG_WRITERS
      ntype <- if (!length(lits)) "unresolved" else
               if (length(lits) == 1 && al$start && al$end) "literal" else "dynamic"
      if (ntype == "unresolved") {
        shown <- if (is.null(a)) NA_character_ else a$text; pat <- NA_character_; open <- FALSE
      } else {
        core <- paste(regex_escape(lits), collapse = ".*"); glob <- paste(lits, collapse = "*")
        # A wildcard replaces whatever the expression contributes before the
        # first and after the last literal. A name that already ends in a known
        # extension needs no tail.
        head_re <- if (al$start) "" else ".*"; head_gl <- if (al$start) "" else "*"
        tail_open <- !al$end || (!is_fig && !grepl(EXT_RE, lits[length(lits)]))
        tail_re <- if (tail_open) ".*" else ""; tail_gl <- if (tail_open) "*" else ""
        if (is_fig) {
          pat <- paste0("^", head_re, core, tail_re,
                        "\\.(", paste(FIG_FORMATS, collapse = "|"), ")$")
          shown <- paste0(head_gl, glob, tail_gl, ".", FIG_FORMATS[1])
        } else {
          pat <- paste0("^", head_re, core, tail_re, "$")
          shown <- paste0(head_gl, glob, tail_gl)
        }
        # A pattern is open when a wildcard stands for something computed at run
        # time. FIG_FORMATS is an enumeration, not a wildcard.
        open <- nzchar(head_re) || nzchar(tail_re) || length(lits) > 1
      }
      # Match the basename and the path relative to the output root, so a call
      # that names a subdirectory also matches.
      m <- if (is.na(pat)) listing[0] else listing[grepl(pat, file) | grepl(pat, rel)]
      row_id <- length(rows) + 1L
      if (nrow(m))
        hits[[length(hits) + 1L]] <- data.table(row_id = row_id, script = sc, script_dir = sdir,
                                                line = line, rel = m$rel, file = m$file,
                                                mtime_num = m$mtime_num, open = open)
      rows[[row_id]] <- data.table(
        file = shown, script = sc, script_dir = sdir, exists = nrow(m) > 0,
        description = nc$text, description_line = nc$line,
        writer = w, line = line, name_type = ntype,
        matched_by = if (ntype == "unresolved") NA_character_
                     else if (open) "wildcard" else "literal",
        location = paste(sort(unique(m$location)), collapse = ";"), n_matching_files = nrow(m),
        matched_files = paste(head(sort(m$file), 8), collapse = "; "),
        argument = if (is.null(a)) NA_character_ else substr(gsub("\\s+", " ", a$text), 1, 120))
    }
  }
  idx <- rbindlist(rows, fill = TRUE)
  H <- if (length(hits)) rbindlist(hits) else
       data.table(row_id = integer(0), script = character(0), script_dir = character(0),
                  line = integer(0), rel = character(0), file = character(0),
                  mtime_num = numeric(0), open = logical(0))
  # A script's most recent run, from the mtimes of all files its writer calls
  # match, with the two margins of that rule.
  runs <- H[, { p <- run_margins(mtime_num)
                .(run_start = p$run_start, run_end = p$run_end,
                  n_run_times = p$n_times, run_span_min = p$span_min,
                  max_internal_gap_min = p$max_internal_gap_min,
                  gap_to_previous_min = p$gap_to_previous_min) },
            by = .(script, script_dir)]
  # The pipeline's most recent run: the newest cluster of live script runs.
  # Output of a script outside it comes from an earlier run.
  pipe <- run_margins(runs[script_dir == ".", run_start], PIPELINE_GAP_MIN)
  pipe_start <- if (is.na(pipe$run_start)) -Inf else pipe$run_start
  runs[, run_in_pipeline_run := run_start >= pipe_start]
  # Status in words, and an empty logical where the pipeline test does not
  # apply. The internal logical still exempts retired scripts.
  runs[, superseded_script := script_dir != "."]
  runs[, run_status := fifelse(superseded_script, RUN_STATUS_SUPERSEDED,
                       fifelse(run_in_pipeline_run, RUN_STATUS_IN, RUN_STATUS_NOT_RERUN))]
  runs[, reported_in_pipeline_run := fifelse(superseded_script, NA, run_in_pipeline_run)]
  # Output newer than the script's session-info record. The 10-minute tolerance
  # covers the delay between a stage's last write and write_session_info().
  runs[, session_info_time := unname(sess[script])]
  runs[, output_after_session_info := is.finite(session_info_time) & run_end > session_info_time + 600]
  # Two runs closer than RUN_GAP_MIN appear as one and leave no
  # gap_to_previous_min. The signal is a session-info record inside the
  # apparent run: the script finished a run and wrote again. Flagged, not
  # corrected.
  runs[, apparent_run_merges_two_runs :=
         output_after_session_info & is.finite(session_info_time) &
         session_info_time >= run_start]
  H[runs, on = .(script, script_dir),
    `:=`(run_start = i.run_start, run_end = i.run_end,
         run_in_pipeline_run = i.run_in_pipeline_run,
         output_after_session_info = i.output_after_session_info)]
  # Currency test for every match: the file was written in the claiming
  # script's most recent run and, for a live script, that run belongs to the
  # pipeline's most recent run. Retired scripts are exempt from the second test.
  H[, in_pipeline := script_dir != "." | run_in_pipeline_run]
  H[, current := mtime_num >= run_start & in_pipeline]
  H[, pre_pipeline := !in_pipeline]
  live <- H[script_dir == "."]
  claimed <- unique(live[current == TRUE, rel])
  stale <- live[current == FALSE & !(rel %in% claimed)][order(rel, row_id)]
  stale <- stale[!duplicated(rel), .(rel, stale_script = script, stale_line = line,
                                     stale_run = run_start, stale_pre = pre_pipeline,
                                     stale_open = open, stale_running = output_after_session_info)]
  stale[, stale_reason := fifelse(stale_pre,
    paste("written before the pipeline's most recent run: the script that claims it",
          "was not re-run with the pipeline, so the file is output of an earlier",
          "version of the analysis"),
    fifelse(stale_running %in% TRUE,
      paste("the script that claims it has output newer than its own session-info",
            "record, so it was running as this index was built or records provenance",
            "before its last write, and this file is from the earlier of the two",
            "apparent runs; re-run the index when the pipeline is quiet"),
      fifelse(stale_open,
        paste("superseded run: matched only by a dynamic (wildcard) file name,",
              "and older than that script's most recent run"),
        paste("superseded run: older than the most recent run of the script whose",
              "writer call names it"))))]
  sup <- unique(H[script_dir != "." & current == TRUE, .(rel, sup_script = script)])
  # Per-call counts, so a call matching stale files is visible in the index.
  cur_n <- H[, .(n_claimed_files = sum(current), n_stale_files = sum(!current),
                 stale_files = paste(head(sort(file[!current]), 8), collapse = "; ")),
             by = row_id]
  idx[, row_id := .I]
  idx[cur_n, on = "row_id", `:=`(n_claimed_files = i.n_claimed_files,
                                 n_stale_files = i.n_stale_files, stale_files = i.stale_files)]
  idx[is.na(n_claimed_files), `:=`(n_claimed_files = 0L, n_stale_files = 0L, stale_files = "")]
  idx[, row_id := NULL]
  # exists: something on disk matches the name. exists_current: some of it
  # passed the currency test.
  idx[, exists_current := n_claimed_files > 0]
  idx[runs, on = .(script, script_dir),
      `:=`(claiming_run_in_pipeline_run = i.reported_in_pipeline_run,
           claiming_run_status = i.run_status)]
  # A script with no matching output has no run to test and gets
  # RUN_STATUS_NO_OUTPUT.
  idx[is.na(claiming_run_status), claiming_run_status := RUN_STATUS_NO_OUTPUT]
  idx[, index_snapshot_time := snapshot]

  unc <- listing[unclaimed_scope == TRUE & !(rel %in% claimed)][order(location, file)]
  unc[, `:=`(reason = "no writer call in analysis/R (scanned recursively) produces this name",
             written_by = NA_character_, writing_script_last_run = NA_character_)]
  unc[stale, on = "rel", `:=`(
    reason = i.stale_reason,
    written_by = paste0(i.stale_script, ":", i.stale_line),
    writing_script_last_run = vapply(i.stale_run, fmt_t, character(1)))]
  unc[sup, on = "rel", `:=`(
    reason = "written only by a superseded script in analysis/R/_superseded",
    written_by = i.sup_script)]
  unc[file %in% names(UNTRACKED_LIVE), reason := unname(UNTRACKED_LIVE[file])]
  unc[, `:=`(unclaimed_scope = NULL, mtime_num = NULL)]
  unc[, index_snapshot_time := snapshot]
  # Per-script timing and the headroom of each rule in this tree.
  cnt <- H[, .(n_matched_files = uniqueN(rel), n_current_files = uniqueN(rel[current])),
           by = .(script, script_dir)]
  gaps <- merge(copy(runs), cnt, by = c("script", "script_dir"), all.x = TRUE)
  gaps[, run_in_pipeline_run := reported_in_pipeline_run]
  gaps[, `:=`(superseded_script = NULL, reported_in_pipeline_run = NULL)]
  gaps[, `:=`(run_start = vapply(run_start, fmt_t, character(1)),
              run_end = vapply(run_end, fmt_t, character(1)),
              session_info_time = vapply(session_info_time, fmt_t, character(1)),
              run_span_min = round(run_span_min, 2),
              max_internal_gap_min = round(max_internal_gap_min, 2),
              gap_to_previous_min = round(gap_to_previous_min, 2),
              run_gap_min = RUN_GAP_MIN)]
  gaps[, `:=`(split_headroom_min = round(RUN_GAP_MIN - max_internal_gap_min, 2),
              merge_margin_min = round(gap_to_previous_min - RUN_GAP_MIN, 2))]
  # The split margin exists for every row, the merge margin only where an
  # earlier run of the script is still on disk.
  gaps[, merge_margin_status := fifelse(is.na(merge_margin_min),
                                        MERGE_NOT_MEASURABLE, MERGE_MEASURED)]
  gaps[, index_snapshot_time := snapshot]
  gcol <- c("script", "script_dir", "run_start", "run_end", "n_run_times",
            "run_span_min", "max_internal_gap_min", "gap_to_previous_min",
            "run_in_pipeline_run", "run_status", "session_info_time",
            "output_after_session_info", "apparent_run_merges_two_runs",
            "n_matched_files", "n_current_files",
            "run_gap_min", "split_headroom_min", "merge_margin_min",
            "merge_margin_status", "index_snapshot_time")
  setcolorder(gaps, c(intersect(gcol, names(gaps)), setdiff(names(gaps), gcol)))
  setorder(gaps, -max_internal_gap_min)
  list(index = idx, unclaimed = unc, gaps = gaps, pipeline = pipe, snapshot = snapshot)
}

# Written now so the coverage table includes 31's own session-info file. The call at
# the end overwrites it with identical content.
write_session_info("31_supplementary_data_and_index")

# ---- 5b. write_session_info coverage ----------------------------------------
cov <- rbindlist(lapply(script_files, function(f) {
  sc <- basename(f); pd <- parse_pd(f); tag <- NA_character_; calls <- FALSE
  if (!is.null(pd)) {
    fn <- pd[token == "SYMBOL_FUNCTION_CALL" & text == "write_session_info"]
    calls <- nrow(fn) > 0
    if (calls) {
      call_id <- pd[id == fn$parent[1], parent]
      a <- pick_arg(call_args(pd, call_id, fn$parent[1]), "tag", 1L)
      lits <- if (is.null(a)) character(0) else arg_literals(pd, a)$lits
      if (length(lits)) tag <- lits[1]
    }
  }
  sif <- if (is.na(tag)) NA_character_ else paste0("sessionInfo_", tag, ".txt")
  data.table(script = sc, numbered = grepl("^[0-9]", sc), parsed = !is.null(pd),
             calls_write_session_info = calls, tag = tag, sessionInfo_file = sif,
             sessionInfo_file_exists = !is.na(sif) && file.exists(file.path(LOG_DIR, sif)),
             v9_log_exists = file.exists(file.path(LOG_DIR, paste0("v9_", sub("\\.R$", "", sc), ".log"))))
}))
save_tsv(cov, "31_session_info_coverage.tsv")
# Denominator: every live script except 00_config.R.
msg("Scripts never calling write_session_info(): ",
    sum(cov$calls_write_session_info == FALSE), " of ", nrow(cov), " live scripts (",
    paste(cov[calls_write_session_info == FALSE, script], collapse = ", "), ")")

# ---- 5c. package versions ----------------------------------------------------
pkg_rows <- rbindlist(lapply(c(file.path(R_DIR, "00_config.R"), script_files), function(f) {
  sc <- basename(f); pd <- parse_pd(f); if (is.null(pd)) return(NULL)
  out <- list()
  ns <- unique(pd[token == "SYMBOL_PACKAGE", text])
  if (length(ns)) out[[1]] <- data.table(package = ns, how = "namespace (::)")
  fn <- pd[token == "SYMBOL_FUNCTION_CALL" & text %in% c("library", "require", "requireNamespace")]
  for (k in seq_len(nrow(fn))) {
    call_id <- pd[id == fn$parent[k], parent]
    a <- pick_arg(call_args(pd, call_id, fn$parent[k]), "package", 1L)
    if (is.null(a)) next
    kids <- pd[parent == a$id]
    if (nrow(kids) != 1 || !kids$token[1] %in% c("SYMBOL", "STR_CONST")) next
    # library(data.table) names the package directly. requireNamespace(p) in a
    # for (p in c(...)) loop names the loop variable, resolved to its strings.
    nm <- if (kids$token[1] == "STR_CONST") strip_q(kids$text[1]) else {
      s <- symbol_strings(pd, kids$text[1]); if (length(s)) s else kids$text[1]
    }
    out[[length(out) + 1]] <- data.table(package = nm, how = fn$text[k])
  }
  if (!length(out)) return(NULL)
  cbind(rbindlist(out), script = sc)
}))
pkgs <- pkg_rows[, .(how = paste(sort(unique(how)), collapse = "; "),
                     n_scripts = uniqueN(script),
                     scripts = paste(unique(script), collapse = "; ")), by = package][order(package)]
pkgs[, version := vapply(package, function(p)
  tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), character(1))]
pkgs[, installed := !is.na(version)]
pkgs <- rbind(data.table(package = "R", how = "interpreter", n_scripts = nrow(cov) + 1L,
                         scripts = "all", version = paste(R.version$major, R.version$minor, sep = "."),
                         installed = TRUE), pkgs)
pkgs[, platform := R.version$platform]
save_tsv(pkgs, "31_package_versions.tsv")
msg("Packages referenced: ", nrow(pkgs) - 1L, "; not installed: ",
    paste(pkgs[installed == FALSE, package], collapse = ", "))

# ---- 5d. results index (two passes so that the index reports its own files) --
# Session-info mtimes: the recorded end of each live script's most recent run.
sess_time <- setNames(vapply(cov$sessionInfo_file, function(f) {
  if (is.na(f)) return(NA_real_)
  p <- file.path(LOG_DIR, f)
  if (file.exists(p)) as.numeric(file.mtime(p)) else NA_real_
}, numeric(1)), cov$script)
for (pass in 1:2) {
  ix <- build_index(sess_time)
  save_tsv(ix$index, "31_results_index.tsv")
  save_tsv(ix$unclaimed, "31_results_index_unclaimed.tsv")
  save_tsv(ix$gaps, "31_run_gap_diagnostics.tsv")
}
msg("Results index: ", nrow(ix$index), " writer calls in ", uniqueN(ix$index$script),
    " scripts; ", sum(ix$index$exists %in% TRUE), " resolve to an existing file, ",
    sum(ix$index$exists %in% FALSE), " do not; ", sum(ix$index$name_type %in% "unresolved"),
    " unresolved; ", sum(ix$index$matched_by %in% "wildcard"), " matched by a wildcard")
msg("Writer calls whose file is absent: ",
    paste(ix$index[exists %in% FALSE, paste0(script, ":", line, " ", file)], collapse = "; "))
# Unclaimed files by reason: no writer, a superseded run of a live script, or a
# stage not re-run with the pipeline.
unc_by <- ix$unclaimed[, .(n_files = .N,
                           files = paste(head(sort(file), 4), collapse = "; ")), by = reason]
msg("Files in results/ and figures/ not claimed by the live pipeline: ", nrow(ix$unclaimed),
    " (", sum(ix$unclaimed$location == "figures"), " figure files, ",
    sum(ix$unclaimed$location == "results"), " results files)")
for (r in seq_len(nrow(unc_by)))
  msg("  ", unc_by$n_files[r], ": ", unc_by$reason[r], " [", unc_by$files[r], " ...]")
stale_calls <- ix$index[n_stale_files > 0]
if (nrow(stale_calls))
  msg("Writer calls also matching output the current pipeline did not write: ",
      paste(stale_calls[, paste0(script, ":", line, " (", n_stale_files, " not current of ",
                                 n_matching_files, ")")], collapse = "; "))
# Live means the top level of the script directory. Retired scripts are
# reported on their own line.
not_rerun <- unique(ix$index[script_dir == "." & claiming_run_in_pipeline_run %in% FALSE, script])
msg("Live scripts outside the pipeline's most recent run: ",
    if (length(not_rerun)) paste(not_rerun, collapse = ", ") else "none",
    " (of ", uniqueN(ix$index[script_dir == ".", script]), " live scripts with writer calls)")
sup_scripts <- sort(unique(ix$index[script_dir != ".", script]))
msg("Superseded scripts in analysis/R/_superseded, exempt from that test and reported as such: ",
    if (length(sup_scripts)) paste(sup_scripts, collapse = ", ") else "none")

# ---- 5e. the margins the claiming rules stand on ----------------------------
# RUN_GAP_MIN must exceed every gap inside a stage run and fall short of every
# separation between two runs of the same stage. Each margin is the minimum
# over the rows of 31_run_gap_diagnostics.tsv on which it is defined:
#   split = min(split_headroom_min) over all scripts
#   merge = min(merge_margin_min, na.rm = TRUE) over scripts with an earlier
#           run on disk (merge_margin_status "measured")
# Stages re-run while the index is built are kept, because they show the
# tightest separations.
gp <- ix$gaps
running   <- gp[output_after_session_info %in% TRUE, script]
measured  <- gp[merge_margin_status == MERGE_MEASURED]
worst_split <- gp[which.min(split_headroom_min)]
worst_merge <- measured[which.min(merge_margin_min)]
msg("Run-gap rule: RUN_GAP_MIN = ", RUN_GAP_MIN,
    " min, tested on both sides against 31_run_gap_diagnostics.tsv (", nrow(gp), " scripts).")
if (nrow(worst_split))
  msg("  Split side (the rule must exceed every gap inside one run): largest gap inside ",
      "one stage's most recent run ", worst_split$max_internal_gap_min, " min (",
      worst_split$script, ", ", worst_split$n_matched_files, " files spanning ",
      worst_split$run_span_min, " min), leaving ", worst_split$split_headroom_min,
      " min -- that is min(split_headroom_min) over all ", nrow(gp), " rows.")
msg("  Merge side (the rule must fall short of every separation between two runs of one ",
    "stage): measurable for ", nrow(measured), " of ", nrow(gp),
    " scripts; the other ", nrow(gp) - nrow(measured),
    " overwrote every matched output in one run and have no separation to measure.")
if (nrow(worst_merge)) {
  msg("    Smallest separation between two runs of one stage: ",
      worst_merge$gap_to_previous_min, " min (", worst_merge$script,
      "), leaving ", worst_merge$merge_margin_min,
      " min -- that is min(merge_margin_min, na.rm = TRUE), over ",
      paste(sort(measured$script), collapse = ", "), ".")
} else {
  msg("    No script has an earlier run on disk, so the merge margin is not measurable ",
      "in this tree at all and RUN_GAP_MIN is untested on that side.")
}
msg("  Stages with output newer than their own session-info record, so either still ",
    "running as this index was built or recording provenance before their last write: ",
    if (length(running)) paste(running, collapse = ", ") else "none",
    ". Their timing is reported like any other; the binding row above ",
    if (worst_split$output_after_session_info %in% TRUE ||
        (nrow(worst_merge) && worst_merge$output_after_session_info %in% TRUE))
      "IS one of them, which is what a tree written to while the index is built looks like."
    else "is not one of them.")
# A separation below RUN_GAP_MIN leaves no gap_to_previous_min, so rows flagged
# here are merged runs whose outputs are all reported as current.
merged <- gp[apparent_run_merges_two_runs %in% TRUE]
msg("  Stages whose apparent most recent run MERGES two real runs (their session-info ",
    "record lies inside it, so they finished a run and wrote again less than RUN_GAP_MIN ",
    "later): ",
    if (nrow(merged)) paste(sprintf("%s (%s to %s, largest internal gap %s min)",
                                    merged$script, merged$run_start, merged$run_end,
                                    merged$max_internal_gap_min), collapse = "; ") else "none",
    ". A separation below RUN_GAP_MIN leaves no gap_to_previous_min to measure, so the ",
    "merge margin above is a minimum over the separations the rule separated correctly; ",
    "any row flagged here is one it did not, and every output of both runs is reported ",
    "as current.")
msg("  Pipeline run: started ", fmt_t(ix$pipeline$run_start), ", ", ix$pipeline$n_times,
    " stage runs, largest gap between consecutive stage runs ",
    round(ix$pipeline$max_internal_gap_min, 1), " min against PIPELINE_GAP_MIN = ",
    PIPELINE_GAP_MIN, " min; previous run ",
    if (is.na(ix$pipeline$gap_to_previous_min)) "not present"
    else paste0(round(ix$pipeline$gap_to_previous_min / 60 / 24, 1), " days earlier"), ".")
tight_split <- nrow(worst_split) && worst_split$split_headroom_min < 0.25 * RUN_GAP_MIN
tight_merge <- nrow(worst_merge) && worst_merge$merge_margin_min < 0.25 * RUN_GAP_MIN
if (isTRUE(tight_split) || isTRUE(tight_merge))
  warning("RUN_GAP_MIN = ", RUN_GAP_MIN, " min has under a quarter of its value in hand ",
          "on one side (split headroom ",
          if (nrow(worst_split)) worst_split$split_headroom_min else NA, " min from ",
          if (nrow(worst_split)) worst_split$script else NA, ", merge margin ",
          if (nrow(worst_merge)) worst_merge$merge_margin_min else NA, " min from ",
          if (nrow(worst_merge)) worst_merge$script else NA,
          "): check 31_run_gap_diagnostics.tsv before quoting the index as current")

# The index is a snapshot of results/ and figures/ at index_snapshot_time.
msg("Index snapshot of results/ and figures/ taken at ", ix$snapshot,
    "; written into index_snapshot_time in 31_results_index.tsv, ",
    "31_results_index_unclaimed.tsv and 31_run_gap_diagnostics.tsv. Re-run once every ",
    "other stage has settled: any stage re-run after this moment is not in it.")
msg("31 finished in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
write_session_info("31_supplementary_data_and_index")
banner("31 | done")
