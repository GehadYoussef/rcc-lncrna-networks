# 08_validation_data.R: build the CPTAC-3 clear-cell RCC external validation cohort.
#
# CPTAC-3 kidney RNA-seq uses the same GDC STAR-Counts pipeline as TCGA-KIRC,
# but a ribo-depleted total RNA protocol (TCGA used polyA selection). Its
# non-feature fraction is therefore much higher, so quality metrics are
# summarised and residualised within cohort (see 00_config.R).
# Inputs: GDC API (open STAR-Counts files and clinical records), cached in data/derived/cache.
# Outputs: validation_dataset.rds (FPKM, one primary-tumour aliquot per clear-cell
#   patient, harmonised clinical data as in 01, STAR metrics from the same file)
#   and results 08_validation_{consort,histology,table1,patient_level_data}.tsv.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(httr); library(survival)
})
banner("08 | Building CPTAC-3 ccRCC validation cohort")

VALID_PROJECT <- "CPTAC-3"
VALID_DIR     <- file.path(DOWNLOAD_DIR, "validation_data", VALID_PROJECT)
dir.create(VALID_DIR, recursive = TRUE, showWarnings = FALSE)

gdc <- function(endpoint, body) {
  r <- httr::POST(paste0("https://api.gdc.cancer.gov/", endpoint),
                  body = jsonlite::toJSON(body, auto_unbox = TRUE),
                  httr::content_type_json(), httr::timeout(300))
  httr::stop_for_status(r)
  jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"),
                     simplifyDataFrame = TRUE)$data$hits
}

# -----------------------------------------------------------------------------
# 1. File manifest
# -----------------------------------------------------------------------------
fm_rds <- file.path(CACHE_DIR, "valid_file_map.rds")
if (file.exists(fm_rds)) {
  file_map <- readRDS(fm_rds); msg("Validation file map from cache")
} else {
  filters <- list(op = "and", content = list(
    list(op = "in", content = list(field = "cases.project.project_id",
                                   value = list(VALID_PROJECT))),
    list(op = "in", content = list(field = "cases.primary_site",
                                   value = list("Kidney"))),
    list(op = "in", content = list(field = "data_category",
                                   value = list("Transcriptome Profiling"))),
    list(op = "in", content = list(field = "data_type",
                                   value = list("Gene Expression Quantification"))),
    list(op = "in", content = list(field = "analysis.workflow_type",
                                   value = list("STAR - Counts"))),
    list(op = "in", content = list(field = "access", value = list("open")))))
  hits <- gdc("files", list(filters = filters, format = "JSON", size = "5000",
                            fields = paste("file_id", "file_name",
                                           "cases.submitter_id",
                                           "cases.samples.submitter_id",
                                           "cases.samples.sample_type", sep = ",")))
  file_map <- rbindlist(lapply(seq_len(nrow(hits)), function(i) {
    cs <- hits$cases[[i]]; smp <- cs$samples[[1]]
    data.table(file_id = hits$file_id[i], file_name = hits$file_name[i],
               patient = cs$submitter_id[1],
               sample_barcode = smp$submitter_id[1],
               sample_type = smp$sample_type[1])
  }))
  saveRDS(file_map, fm_rds)
}
msg(VALID_PROJECT, ": ", nrow(file_map), " files")
print(table(file_map$sample_type))

# -----------------------------------------------------------------------------
# 2. Download (chunked, resumable)
# -----------------------------------------------------------------------------
file_map[, path := winlong(file.path(VALID_DIR, file_id, file_name))]
todo <- file_map[!file.exists(path)]
msg("Already on disk: ", nrow(file_map) - nrow(todo), " / ", nrow(file_map))

if (nrow(todo)) {
  # GDC /data returns a tarball for several ids but the raw file for one id,
  # so a trailing one-file chunk is folded into the previous chunk.
  chunks <- split(todo$file_id, ceiling(seq_len(nrow(todo)) / 40))
  if (length(chunks) > 1 && length(chunks[[length(chunks)]]) == 1) {
    chunks[[length(chunks) - 1]] <- c(chunks[[length(chunks) - 1]],
                                      chunks[[length(chunks)]])
    chunks[[length(chunks)]] <- NULL
  }
  for (i in seq_along(chunks)) {
    ids <- chunks[[i]]
    msg("Downloading chunk ", i, " / ", length(chunks),
        " (", length(ids), " files) ...")
    tf <- tempfile(fileext = ".tar.gz")
    ok <- FALSE
    for (attempt in 1:3) {
      r <- try(httr::POST("https://api.gdc.cancer.gov/data",
                          body = jsonlite::toJSON(list(ids = ids),
                                                  auto_unbox = TRUE),
                          httr::content_type_json(),
                          httr::write_disk(tf, overwrite = TRUE),
                          httr::timeout(1800)), silent = TRUE)
      if (!inherits(r, "try-error") && httr::status_code(r) == 200) { ok <- TRUE; break }
      msg("  attempt ", attempt, " failed; retrying")
    }
    if (!ok) stop("Download failed for chunk ", i)
    if (length(ids) == 1) {
      row <- file_map[file_id == ids]
      dir.create(file.path(VALID_DIR, row$file_id), showWarnings = FALSE)
      file.copy(tf, winlong(file.path(VALID_DIR, row$file_id, row$file_name)),
                overwrite = TRUE)
    } else {
      untar(tf, exdir = VALID_DIR)
    }
    unlink(tf)
  }
}
file_map[, on_disk := file.exists(path)]
msg("Present locally: ", sum(file_map$on_disk), " / ", nrow(file_map))
stopifnot(all(file_map$on_disk))

# -----------------------------------------------------------------------------
# 3. Expression matrix
# -----------------------------------------------------------------------------
vexpr_rds <- file.path(CACHE_DIR, "valid_expr_raw.rds")
if (file.exists(vexpr_rds)) {
  vraw <- readRDS(vexpr_rds); msg("Validation expression from cache")
} else {
  msg("Reading ", nrow(file_map), " STAR count files ...")
  first <- fread(file_map$path[1], skip = 1, showProgress = FALSE)
  kr <- grepl("^ENSG", first$gene_id)
  ann <- data.table(gene_id = first$gene_id[kr], gene_name = first$gene_name[kr],
                    gene_type = first$gene_type[kr])
  fp <- matrix(NA_real_, nrow(ann), nrow(file_map),
               dimnames = list(ann$gene_id, file_map$sample_barcode))
  lib <- numeric(nrow(file_map))
  for (i in seq_len(nrow(file_map))) {
    d <- fread(file_map$path[i], skip = 1, showProgress = FALSE,
               select = c("gene_id", "unstranded", "fpkm_unstranded"))
    d <- d[grepl("^ENSG", gene_id)]
    if (!identical(d$gene_id, ann$gene_id)) d <- d[match(ann$gene_id, d$gene_id)]
    fp[, i] <- d$fpkm_unstranded; lib[i] <- sum(d$unstranded, na.rm = TRUE)
    if (i %% 50 == 0) msg("  ", i, " / ", nrow(file_map))
  }
  vraw <- list(fpkm = fp, gene_ann = ann,
               sample_info = cbind(file_map, libsize = lib))
  saveRDS(vraw, vexpr_rds)
}
msg("Validation matrix: ", nrow(vraw$fpkm), " genes x ", ncol(vraw$fpkm), " samples")

# -----------------------------------------------------------------------------
# 4. Clinical
# -----------------------------------------------------------------------------
vclin_rds <- file.path(CACHE_DIR, "valid_clinical.rds")
if (file.exists(vclin_rds)) {
  vclin <- readRDS(vclin_rds)
} else {
  f <- list(op = "and", content = list(
    list(op = "in", content = list(field = "project.project_id",
                                   value = list(VALID_PROJECT))),
    list(op = "in", content = list(field = "primary_site",
                                   value = list("Kidney")))))
  h <- gdc("cases", list(filters = f, fields = "submitter_id,disease_type",
                         expand = "demographic,diagnoses,follow_ups",
                         format = "JSON", size = "2000"))
  dx1 <- function(l, col) vapply(l, function(d) {
    if (is.null(d) || !is.data.frame(d) || !col %in% names(d) || nrow(d) == 0)
      return(NA_character_)
    v <- as.character(d[[col]])
    v <- v[!is.na(v) & nzchar(v) &
             !tolower(v) %in% c("not reported", "unknown", "not applicable")]
    if (!length(v)) NA_character_ else v[1] }, character(1))
  dxm <- function(l, col) vapply(l, function(d) {
    if (is.null(d) || !is.data.frame(d) || !col %in% names(d) || nrow(d) == 0)
      return(NA_real_)
    v <- suppressWarnings(as.numeric(d[[col]]))
    if (!any(is.finite(v))) NA_real_ else max(v, na.rm = TRUE) }, numeric(1))
  gt <- function(df, col) if (!is.null(df) && col %in% names(df)) df[[col]]
                          else rep(NA, nrow(h))
  # CPTAC-3 records sex under `sex_at_birth`, TCGA under `gender`.
  g1 <- as.character(gt(h$demographic, "gender"))
  g2 <- as.character(gt(h$demographic, "sex_at_birth"))
  sex_chr <- ifelse(!is.na(g1) & nzchar(g1), g1, g2)
  vclin <- data.table(
    patient   = h$submitter_id,
    histology = dx1(h$diagnoses, "primary_diagnosis"),
    gender    = sex_chr,
    vital     = as.character(gt(h$demographic, "vital_status")),
    dtd       = suppressWarnings(as.numeric(gt(h$demographic, "days_to_death"))),
    dtlf_dx   = dxm(h$diagnoses, "days_to_last_follow_up"),
    dtlf_fu   = dxm(h$follow_ups, "days_to_follow_up"),
    stage_raw = dx1(h$diagnoses, "ajcc_pathologic_stage"),
    grade_raw = dx1(h$diagnoses, "tumor_grade"),
    pT        = dx1(h$diagnoses, "ajcc_pathologic_t"),
    pN        = dx1(h$diagnoses, "ajcc_pathologic_n"),
    pM        = dx1(h$diagnoses, "ajcc_pathologic_m"),
    age_days  = dxm(h$diagnoses, "age_at_diagnosis"))
  saveRDS(vclin, vclin_rds)
}
msg("Clinical records: ", nrow(vclin))
print(head(sort(table(vclin$histology), decreasing = TRUE), 8))

# -----------------------------------------------------------------------------
# 5. Harmonise with the same rules as TCGA-KIRC in 01
# -----------------------------------------------------------------------------
# Sub-stages (IIIA, IVB, ...) are collapsed to the numeral.
parse_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "NOT REPORTED")] <- NA
  rn <- sub("[A-C]$", "", sub("^STAGE\\s*", "", x))
  out <- rep(NA_integer_, length(x))
  out[rn == "I"] <- 1L; out[rn == "II"] <- 2L
  out[rn == "III"] <- 3L; out[rn == "IV"] <- 4L
  out
}
parse_grade <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[x == "G1"] <- 1L; out[x == "G2"] <- 2L
  out[x == "G3"] <- 3L; out[x == "G4"] <- 4L
  out                                            # GX stays missing
}
parse_T <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[grepl("^T1", x)] <- 1L; out[grepl("^T2", x)] <- 2L
  out[grepl("^T3", x)] <- 3L; out[grepl("^T4", x)] <- 4L
  out
}
vclin[, `:=`(stage_num = parse_stage(stage_raw), grade_num = parse_grade(grade_raw),
             age = age_days / 365.25,
             sex = factor(tolower(gender), levels = c("female", "male")),
             os_event = as.integer(vital == "Dead"))]
vclin[, os_time := ifelse(!is.na(dtd), dtd, pmax(dtlf_dx, dtlf_fu, na.rm = TRUE))]
vclin[!is.na(os_time) & os_time > OS_CENSOR_DAYS, os_event := 0L]
vclin[!is.na(os_time), os_time := pmin(os_time, OS_CENSOR_DAYS)]

# T comes from pT only and is not imputed from stage. Missing T drops out of
# the adjusted models as an incomplete case.
vclin[, T_stage := parse_T(pT)]
vclin[, N_pos := ifelse(grepl("^N[1-9]", toupper(trimws(pN))), 1L,
                 ifelse(grepl("^N0", toupper(trimws(pN))), 0L, NA_integer_))]
vclin[, M1_as_coded := ifelse(grepl("^M1", toupper(trimws(pM))), 1L,
                       ifelse(grepl("^M0", toupper(trimws(pM))), 0L, NA_integer_))]
vclin[, `:=`(N_pos_imputed = is.na(N_pos), M1_imputed = is.na(M1_as_coded))]
vclin[is.na(N_pos),       N_pos       := 0L]        # NX -> N0
vclin[is.na(M1_as_coded), M1_as_coded := 0L]        # MX -> M0, as coded
# Stage IV without a T4 primary implies distant metastasis, so MX patients at
# stage IV with T < 4 are set to M1. M1_as_coded keeps the unreconciled value
# for the sensitivity analyses in 09.
vclin[, M1 := M1_as_coded]
vclin[, M1_stage_reconciled := as.integer(M1_imputed & !is.na(stage_num) &
                                          stage_num == 4L & !is.na(T_stage) &
                                          T_stage < 4L)]
vclin[M1_stage_reconciled == 1L, M1 := 1L]

# Descriptive groupings only. The models use the ordinal and binary TNM terms.
vclin[, stage_group := factor(fifelse(is.na(stage_num), NA_character_,
                              fifelse(stage_num <= 2, "I-II", "III-IV")),
                              levels = c("I-II", "III-IV"))]
vclin[, grade_group := factor(fifelse(is.na(grade_num), NA_character_,
                              fifelse(grade_num <= 2, "G1-2", "G3-4")),
                              levels = c("G1-2", "G3-4"))]

# One primary-tumour aliquot per patient (deepest library). file_id is kept so
# expression and STAR metrics come from the same file.
samp_all <- as.data.table(vraw$sample_info)
stopifnot(nrow(samp_all) == ncol(vraw$fpkm))
flow <- list(files_total = nrow(samp_all))
samp <- samp_all[sample_type == "Primary Tumor"]
flow$primary_tumour_samples <- nrow(samp)
setorder(samp, patient, -libsize)
flow$dropped_duplicate_aliquots <- sum(duplicated(samp$patient))
samp <- samp[!duplicated(patient)]

v <- merge(samp[, .(patient, sample_barcode, file_id, file_name, libsize)],
           vclin, by = "patient")
flow$with_clinical <- nrow(v)

# Clear-cell histology only. "Renal cell carcinoma, NOS" is kept because GDC
# codes the CPTAC ccRCC cohort that way. Other subtypes and missing histology
# are excluded.
is_clear_cell <- function(h) {
  h0   <- tolower(trimws(as.character(h)))
  excl <- grepl("papillary|chromophobe|collecting duct|leiomyomatosis|hlrcc", h0) |
          (grepl("sarcomatoid", h0) & !grepl("clear cell", h0))
  keep <- !is.na(h) & (h %in% c("Renal cell carcinoma, NOS",
                                "Clear cell adenocarcinoma, NOS") |
                       grepl("clear cell", h0))
  keep & !excl
}
v[, clear_cell := is_clear_cell(histology)]
hist_tbl <- v[, .(n = .N, retained = clear_cell[1]),
              by = .(histology = fifelse(is.na(histology), "(not reported)", histology))]
setorder(hist_tbl, -retained, -n)
save_tsv(hist_tbl, "08_validation_histology.tsv")
print(hist_tbl)
flow$excluded_non_clear_cell <- sum(!v$clear_cell)
v <- v[clear_cell == TRUE][, clear_cell := NULL]
flow$clear_cell_or_rcc_nos <- nrow(v)

v <- v[!is.na(os_time) & os_time > 0 & !is.na(os_event)]
flow$usable_survival <- nrow(v)
flow$deaths <- sum(v$os_event)
flow$T_missing <- sum(is.na(v$T_stage))
flow$grade_missing <- sum(is.na(v$grade_num))
flow$N_assessed <- sum(!v$N_pos_imputed)
flow$M_assessed <- sum(!v$M1_imputed)
flow$m1_reconciled_from_stage_iv <- sum(v$M1_stage_reconciled)
cc_clin <- complete.cases(v[, .(age, sex, T_stage, N_pos, M1, grade_num)])
flow$complete_clinical_set <- sum(cc_clin)
flow$deaths_complete_clinical_set <- sum(v$os_event[cc_clin])

flow_dt <- data.frame(step = names(flow), n = unlist(flow), row.names = NULL)
save_tsv(flow_dt, "08_validation_consort.tsv")
print(flow_dt)
stopifnot(!anyDuplicated(v$sample_barcode), !anyDuplicated(v$file_id))

# -----------------------------------------------------------------------------
# 6. STAR alignment metrics for every retained sample
# -----------------------------------------------------------------------------
# Keyed on file_id because one sample can have several aliquots. The cache is
# recomputed if its (sample_barcode, file_id) pairs differ from the cohort.
vqc_f <- file.path(CACHE_DIR, "valid_star_qc.rds")
qc_ok <- FALSE
if (file.exists(vqc_f)) {
  vqc <- as.data.table(readRDS(vqc_f))
  qc_ok <- all(c("sample_barcode", "file_id", "pct_unmapped", "pct_multimapping",
                 "pct_noFeature", "pct_ambiguous", "assigned_reads") %in% names(vqc)) &&
           nrow(vqc) == nrow(v) &&
           setequal(paste(vqc$sample_barcode, vqc$file_id),
                    paste(v$sample_barcode, v$file_id))
  msg("STAR QC cache ", if (qc_ok) "matches the cohort" else "stale; recomputing")
}
if (!qc_ok) {
  paths <- winlong(file.path(DOWNLOAD_DIR, "validation_data", VALID_PROJECT,
                             v$file_id, v$file_name))
  stopifnot(all(file.exists(paths)))
  msg("Reading STAR alignment summaries for ", nrow(v), " CPTAC-3 files ...")
  vqc <- rbindlist(lapply(seq_len(nrow(v)), function(i) {
    cbind(data.table(sample_barcode = v$sample_barcode[i], file_id = v$file_id[i]),
          read_star_summary(paths[i]))
  }))
  saveRDS(vqc, vqc_f)
}
vqc <- vqc[match(v$sample_barcode, sample_barcode)]
stopifnot(identical(vqc$file_id, v$file_id), !anyNA(vqc$pct_noFeature))
# assigned_reads and libsize are both the sum of unstranded ENSG counts, so a
# mismatch means they were read from different files.
if (any(abs(vqc$assigned_reads - v$libsize) > 1))
  warning("assigned_reads differs from libsize for ",
          sum(abs(vqc$assigned_reads - v$libsize) > 1), " samples")
qc <- vqc[, .(sample_barcode, pct_unmapped, pct_multimapping, pct_noFeature,
              pct_ambiguous, assigned_reads)]
msg("CPTAC-3 non-feature fraction: median ", round(median(qc$pct_noFeature), 1),
    "% (IQR ", round(quantile(qc$pct_noFeature, 0.25), 1), "-",
    round(quantile(qc$pct_noFeature, 0.75), 1), "%)")

# -----------------------------------------------------------------------------
# 7. Save the validation dataset
# -----------------------------------------------------------------------------
# Columns are selected by file_id so they match the QC source file.
idx <- match(v$file_id, samp_all$file_id)
stopifnot(!anyNA(idx))
vfpkm <- vraw$fpkm[, idx, drop = FALSE]
colnames(vfpkm) <- v$sample_barcode
saveRDS(list(fpkm = vfpkm, gene_ann = vraw$gene_ann, cohort = v, qc = qc),
        file.path(CACHE_DIR, "validation_dataset.rds"))
msg("Validation dataset: ", nrow(vfpkm), " genes x ", ncol(vfpkm), " samples, ",
    sum(v$os_event), " deaths")

# -----------------------------------------------------------------------------
# 8. Cohort characteristics table and patient-level export
# -----------------------------------------------------------------------------
rkm <- survfit(Surv(os_time, 1 - os_event) ~ 1, data = v)
median_fu_rkm <- unname(summary(rkm)$table[["median"]])
vq <- merge(v, qc, by = "sample_barcode")
tbl <- vq[, .(n = .N, deaths = sum(os_event),
              median_fu = round(median(os_time)),
              median_fu_reverse_km = round(median_fu_rkm),
              age_median = round(median(age, na.rm = TRUE), 1),
              male_pct = round(100 * mean(sex == "male", na.rm = TRUE), 1),
              stage_III_IV_pct = round(100 * mean(stage_group == "III-IV", na.rm = TRUE), 1),
              grade_G3_4_pct = round(100 * mean(grade_group == "G3-4", na.rm = TRUE), 1),
              T_available = sum(!is.na(T_stage)),
              T3_T4_pct = round(100 * mean(T_stage >= 3, na.rm = TRUE), 1),
              M1_pct_all = round(100 * mean(M1), 1),
              M_assessed = sum(!M1_imputed),
              M1_pct_assessed = round(100 * mean(M1_as_coded[!M1_imputed]), 1),
              n_M1_reconciled = sum(M1_stage_reconciled),
              N_assessed = sum(!N_pos_imputed),
              N_pos_pct_assessed = round(100 * mean(N_pos[!N_pos_imputed]), 1),
              grade_available = sum(!is.na(grade_num)),
              noFeature_median = round(median(pct_noFeature), 2),
              noFeature_IQR_lo = round(quantile(pct_noFeature, 0.25), 2),
              noFeature_IQR_hi = round(quantile(pct_noFeature, 0.75), 2),
              multimap_median = round(median(pct_multimapping), 2),
              assigned_reads_median = round(median(assigned_reads)))]
save_tsv(tbl, "08_validation_table1.tsv")
print(tbl)

pld <- vq[, .(patient, sample_barcode, age = round(age, 2), sex = as.character(sex),
              T_stage, N_pos, M1, M1_as_coded, M1_stage_reconciled,
              N_pos_imputed = as.integer(N_pos_imputed),
              M1_imputed = as.integer(M1_imputed),
              grade_num, stage_num, os_time, os_event,
              pct_noFeature, pct_multimapping, pct_unmapped, pct_ambiguous,
              assigned_reads, histology)]
setorder(pld, patient)
save_tsv(pld, "08_validation_patient_level_data.tsv")

write_session_info("08_validation_data")
banner("08 | done")
