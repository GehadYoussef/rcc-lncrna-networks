# 01_build_data.R: build the TCGA-KIRC expression matrices and clinical cohort.
#
# Inputs: the STAR gene-count files and BCR clinical XML from 00, plus GDC API
# file and clinical metadata. Primary tumours only, one aliquot per patient
# (largest library). Caches expr_raw.rds and dataset.rds in data/derived/cache
# and writes the CONSORT flow, cohort summary, TNM completeness and
# patient-level tables to results.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(httr); library(xml2)
  library(survival)
})
banner("01 | Building dataset")

# ---- 1. file UUID to barcode and sample-type map (GDC API, cached) ----
gdc_post <- function(endpoint, ..., size = 5000L) {
  extra <- list(...)
  if (endpoint == "files") {
    filters <- list(op = "and", content = list(
      list(op = "in", content = list(field = "cases.project.project_id",
                                     value = list(PROJECT_ID))),
      list(op = "in", content = list(field = "data_category",
                                     value = list("Transcriptome Profiling"))),
      list(op = "in", content = list(field = "data_type",
                                     value = list("Gene Expression Quantification"))),
      list(op = "in", content = list(field = "analysis.workflow_type",
                                     value = list("STAR - Counts")))
    ))
  } else {
    filters <- list(op = "in", content = list(field = "project.project_id",
                                              value = list(PROJECT_ID)))
  }
  body <- c(list(filters = filters, format = "JSON", size = as.character(size)),
            extra)
  r <- httr::POST(paste0("https://api.gdc.cancer.gov/", endpoint),
                  body = jsonlite::toJSON(body, auto_unbox = TRUE),
                  httr::content_type_json(), httr::timeout(180))
  httr::stop_for_status(r)
  jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"),
                     simplifyDataFrame = TRUE)$data$hits
}

file_map_rds <- file.path(CACHE_DIR, "gdc_file_map.rds")
if (file.exists(file_map_rds)) {
  file_map <- readRDS(file_map_rds); msg("File map loaded from cache")
} else {
  msg("Querying GDC /files for barcode + sample-type mapping ...")
  hits <- gdc_post("files", fields = paste(
    "file_id", "file_name", "cases.submitter_id",
    "cases.samples.submitter_id", "cases.samples.sample_type", sep = ","))
  file_map <- rbindlist(lapply(seq_len(nrow(hits)), function(i) {
    cs <- hits$cases[[i]]
    smp <- cs$samples[[1]]
    data.table(file_id      = hits$file_id[i],
               file_name    = hits$file_name[i],
               patient      = cs$submitter_id[1],
               sample_barcode = smp$submitter_id[1],
               sample_type  = smp$sample_type[1])
  }))
  saveRDS(file_map, file_map_rds)
}
msg("GDC reports ", nrow(file_map), " STAR-Counts files for ", PROJECT_ID)
print(table(file_map$sample_type))

# Restrict to files on disk. winlong() adds the Windows extended-length prefix.
file_map[, path := winlong(file.path(EXPR_DIR, file_id, file_name))]
file_map[, on_disk := file.exists(path)]
msg("Present locally: ", sum(file_map$on_disk), " / ", nrow(file_map))

consort <- list(
  gdc_files_total            = nrow(file_map),
  files_present_locally      = sum(file_map$on_disk)
)

# ---- 2. read expression matrices (cached) ----
expr_rds <- file.path(CACHE_DIR, "expr_raw.rds")
if (file.exists(expr_rds)) {
  expr_raw <- readRDS(expr_rds); msg("Expression matrices loaded from cache")
} else {
  fm <- file_map[on_disk == TRUE]
  msg("Reading ", nrow(fm), " STAR count files ...")
  first <- fread(fm$path[1], skip = 1, showProgress = FALSE)
  keep_rows <- grepl("^ENSG", first$gene_id)
  gene_ann  <- data.table(gene_id   = first$gene_id[keep_rows],
                          gene_name = first$gene_name[keep_rows],
                          gene_type = first$gene_type[keep_rows])
  n_g <- nrow(gene_ann)

  fpkm_mat  <- matrix(NA_real_, nrow = n_g, ncol = nrow(fm),
                      dimnames = list(gene_ann$gene_id, fm$sample_barcode))
  libsize   <- numeric(nrow(fm))
  # STAR alignment summary rows. Stage 02 regresses pct_noFeature out before
  # network construction.
  qc_cols   <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
  qc_mat    <- matrix(NA_real_, nrow(fm), length(qc_cols),
                      dimnames = list(fm$sample_barcode, qc_cols))

  for (i in seq_len(nrow(fm))) {
    d <- fread(fm$path[i], skip = 1, showProgress = FALSE,
               select = c("gene_id", "unstranded", "fpkm_unstranded"))
    nn <- setNames(as.numeric(d$unstranded), d$gene_id)
    qc_mat[i, ] <- nn[qc_cols]
    d <- d[grepl("^ENSG", gene_id)]
    if (!identical(d$gene_id, gene_ann$gene_id)) {
      d <- d[match(gene_ann$gene_id, d$gene_id)]
    }
    fpkm_mat[, i] <- d$fpkm_unstranded
    libsize[i]    <- sum(d$unstranded, na.rm = TRUE)
    if (i %% 50 == 0) msg("  ", i, " / ", nrow(fm))
  }
  tot <- libsize + rowSums(qc_mat, na.rm = TRUE)
  qc_pct <- as.data.table(100 * qc_mat / tot)
  setnames(qc_pct, c("pct_unmapped", "pct_multimapping",
                     "pct_noFeature", "pct_ambiguous"))
  expr_raw <- list(fpkm = fpkm_mat, gene_ann = gene_ann,
                   sample_info = cbind(fm, libsize = libsize, qc_pct))
  saveRDS(expr_raw, expr_rds)
  msg("Expression matrices cached")
}

fpkm     <- expr_raw$fpkm
gene_ann <- expr_raw$gene_ann
samp     <- as.data.table(expr_raw$sample_info)
msg("Expression matrix: ", nrow(fpkm), " genes x ", ncol(fpkm), " samples")

# ---- 3. sample-level filtering ----
samp[, keep := sample_type %in% KEEP_SAMPLE_TYPES]
consort$excluded_not_primary_tumour <- sum(!samp$keep)
msg("Excluded ", sum(!samp$keep), " non-primary-tumour samples")

samp_t <- samp[keep == TRUE]

# Library failure: too few assigned reads or an extreme non-feature fraction.
# Applied per file before the one-sample-per-patient rule, so a patient with a
# second usable aliquot keeps it.
samp_t[, failed_library := libsize < MIN_ASSIGNED_READS | pct_noFeature > MAX_NOFEATURE_PCT]
consort$excluded_failed_library <- sum(samp_t$failed_library)
if (any(samp_t$failed_library))
  msg("Excluded ", sum(samp_t$failed_library), " failed librar",
      if (sum(samp_t$failed_library) == 1) "y: " else "ies: ",
      paste(sprintf("%s (%.1f M assigned reads, %.1f%% non-feature)",
                    samp_t[failed_library == TRUE, sample_barcode],
                    samp_t[failed_library == TRUE, libsize] / 1e6,
                    samp_t[failed_library == TRUE, pct_noFeature]), collapse = "; "))
samp_t <- samp_t[failed_library == FALSE]

# One sample per patient: keep the aliquot with the largest library size.
setorder(samp_t, patient, -libsize)
dup_n <- sum(duplicated(samp_t$patient))
samp_t <- samp_t[!duplicated(patient)]
consort$excluded_duplicate_patient_samples <- dup_n
msg("Excluded ", dup_n, " duplicate samples from the same patient")
consort$tumour_samples_one_per_patient <- nrow(samp_t)

# ---- 4. clinical data: GDC harmonised API and local BCR XML ----
clin_rds <- file.path(CACHE_DIR, "gdc_clinical.rds")
if (file.exists(clin_rds)) {
  clin_api <- readRDS(clin_rds); msg("Clinical data loaded from cache")
} else {
  msg("Querying GDC /cases for harmonised clinical data ...")
  # `expand` is used because requesting demographic.* by name returns only a
  # subset of the sub-fields.
  h <- gdc_post("cases", fields = "submitter_id",
                expand = "demographic,diagnoses,follow_ups")

  # A case can carry several diagnosis records. Take the first record that
  # reports a value, and the latest follow-up across all records.
  dx1 <- function(lst, col) vapply(lst, function(d) {
    if (is.null(d) || !is.data.frame(d) || !col %in% names(d) || nrow(d) == 0)
      return(NA_character_)
    v <- as.character(d[[col]])
    v <- v[!is.na(v) & nzchar(v) & !tolower(v) %in%
             c("not reported", "unknown", "not applicable")]
    if (!length(v)) NA_character_ else v[1]
  }, character(1))
  dxmax <- function(lst, col) vapply(lst, function(d) {
    if (is.null(d) || !is.data.frame(d) || !col %in% names(d) || nrow(d) == 0)
      return(NA_real_)
    v <- suppressWarnings(as.numeric(d[[col]]))
    if (!any(is.finite(v))) NA_real_ else max(v, na.rm = TRUE)
  }, numeric(1))
  dm  <- h$demographic
  get <- function(df, col) if (!is.null(df) && col %in% names(df)) df[[col]]
                           else rep(NA, nrow(h))

  clin_api <- data.table(
    patient       = h$submitter_id,
    gender        = as.character(get(dm, "gender")),
    vital_status  = as.character(get(dm, "vital_status")),
    days_to_death = suppressWarnings(as.numeric(get(dm, "days_to_death"))),
    race          = as.character(get(dm, "race")),
    age_days      = suppressWarnings(as.numeric(dx1(h$diagnoses, "age_at_diagnosis"))),
    stage_raw     = dx1(h$diagnoses, "ajcc_pathologic_stage"),
    grade_api     = dx1(h$diagnoses, "tumor_grade"),
    dtlf          = dxmax(h$diagnoses, "days_to_last_follow_up"),
    dtlf_fu       = dxmax(h$follow_ups, "days_to_follow_up"),
    pT            = dx1(h$diagnoses, "ajcc_pathologic_t"),
    pN            = dx1(h$diagnoses, "ajcc_pathologic_n"),
    pM            = dx1(h$diagnoses, "ajcc_pathologic_m"),
    prior_malig   = dx1(h$diagnoses, "prior_malignancy")
  )
  saveRDS(clin_api, clin_rds)
}
msg("Clinical records from API: ", nrow(clin_api))

# ---- grade from the BCR XML ----
# The harmonised `tumor_grade` field is often "not reported" for KIRC, so
# neoplasm_histologic_grade in the BCR XML is the primary source.
grade_rds <- file.path(CACHE_DIR, "grade_xml.rds")
if (file.exists(grade_rds)) {
  grade_xml <- readRDS(grade_rds)
} else {
  xmls <- list.files(CLIN_XML_DIR, pattern = "\\.xml$", full.names = TRUE,
                     recursive = TRUE)
  msg("Parsing ", length(xmls), " BCR clinical XML files for tumour grade ...")
  # The BCR files use many namespaces, so nodes are selected on local-name().
  ln <- function(x, nm) {
    v <- xml_text(xml_find_first(x, sprintf(".//*[local-name()='%s']", nm)))
    if (length(v) == 0) NA_character_ else trimws(v)
  }
  # Follow-up fields recur once per follow-up record, so take the latest contact.
  ln_max <- function(x, nm) {
    v <- suppressWarnings(as.numeric(xml_text(
      xml_find_all(x, sprintf(".//*[local-name()='%s']", nm)))))
    if (!any(is.finite(v))) NA_real_ else max(v, na.rm = TRUE)
  }
  ln_all <- function(x, nm) {
    v <- trimws(xml_text(xml_find_all(x, sprintf(".//*[local-name()='%s']", nm))))
    v[nzchar(v)]
  }
  grade_xml <- rbindlist(lapply(xmls, function(f) {
    x  <- read_xml(f)
    vs <- ln_all(x, "vital_status")
    data.table(patient    = ln(x, "bcr_patient_barcode"),
               grade_xml  = ln(x, "neoplasm_histologic_grade"),
               stage_xml  = ln(x, "pathologic_stage"),
               sex_xml    = ln(x, "gender"),
               age_xml    = ln(x, "age_at_initial_pathologic_diagnosis"),
               vital_xml  = if (any(toupper(vs) == "DEAD")) "Dead"
                            else if (length(vs)) "Alive" else NA_character_,
               dtd_xml    = ln_max(x, "days_to_death"),
               dtlf_xml   = ln_max(x, "days_to_last_followup"))
  }), fill = TRUE)
  grade_xml <- grade_xml[!is.na(patient) & nzchar(patient)]
  grade_xml <- grade_xml[!duplicated(patient)]
  saveRDS(grade_xml, grade_rds)
}
msg("Grade recovered from XML for ", sum(!is.na(grade_xml$grade_xml)), " patients")

clin <- merge(clin_api, grade_xml, by = "patient", all.x = TRUE)

# ---- 5. derive analysis variables ----
# Stage: sub-stages (IIIA, IIC ...) are collapsed to the Roman numeral.
parse_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "NOT REPORTED", "[NOT AVAILABLE]", "[UNKNOWN]")] <- NA
  rn <- sub("^STAGE\\s*", "", x)
  rn <- sub("[A-C]$", "", rn)          # strip sub-stage letter only at the end
  out <- rep(NA_integer_, length(x))
  out[rn == "I"]   <- 1L
  out[rn == "II"]  <- 2L
  out[rn == "III"] <- 3L
  out[rn == "IV"]  <- 4L
  out
}
clin[, stage_num := parse_stage(stage_raw)]
clin[is.na(stage_num), stage_num := parse_stage(stage_xml)]

parse_grade <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[x %in% c("G1")] <- 1L
  out[x %in% c("G2")] <- 2L
  out[x %in% c("G3")] <- 3L
  out[x %in% c("G4")] <- 4L
  out                                    # GX / not reported stay NA
}
clin[, grade_num := parse_grade(grade_xml)]
clin[is.na(grade_num), grade_num := parse_grade(grade_api)]

# ---- TNM ----
# T category, nodal status and distant metastasis are modelled as separate
# terms in the clinical baseline.
parse_T <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[grepl("^T1", x)] <- 1L; out[grepl("^T2", x)] <- 2L
  out[grepl("^T3", x)] <- 3L; out[grepl("^T4", x)] <- 4L
  out
}
clin[, T_stage := parse_T(pT)]
clin[, N_pos   := ifelse(grepl("^N[1-9]", toupper(trimws(pN))), 1L,
                  ifelse(grepl("^N0", toupper(trimws(pN))), 0L, NA_integer_))]
clin[, M1      := ifelse(grepl("^M1", toupper(trimws(pM))), 1L,
                  ifelse(grepl("^M0", toupper(trimws(pM))), 0L, NA_integer_))]
# NX and MX ("not assessed") are coded 0, as in TCGA-KIRC they mostly mean
# clinically node-negative or non-metastatic. The imputed counts are reported.
clin[, `:=`(N_pos_imputed = is.na(N_pos), M1_imputed = is.na(M1))]
clin[is.na(N_pos), N_pos := 0L]
clin[is.na(M1),    M1    := 0L]
# M1_as_coded keeps pM as coded. Stage IV without T4 implies M1, so MX at
# stage IV with T < 4 is recoded M1. Explicit M0 is never overridden.
clin[, M1_as_coded := M1]
clin[, M1_stage_reconciled := as.integer(!is.na(stage_num) & stage_num == 4L &
                                         !is.na(T_stage) & T_stage < 4L &
                                         M1_imputed)]
clin[M1_stage_reconciled == 1L, M1 := 1L]

clin[, age := age_days / 365.25]
clin[is.na(age) & !is.na(age_xml), age := suppressWarnings(as.numeric(age_xml))]
clin[, sex_chr := tolower(gender)]
clin[is.na(sex_chr) | !nzchar(sex_chr), sex_chr := tolower(sex_xml)]
clin[, sex := factor(sex_chr, levels = c("female", "male"))]
# Vital status from the API, with the BCR XML as fallback. A death in either
# source counts as an event.
clin[, vital := vital_status]
clin[is.na(vital) | !nzchar(vital), vital := vital_xml]
clin[!is.na(vital_xml) & vital_xml == "Dead", vital := "Dead"]
clin[, os_event := as.integer(vital == "Dead")]

# Follow-up time. days_to_last_follow_up is missing from the GDC diagnoses
# records for most living patients, so the latest contact is the maximum over
# GDC diagnoses, GDC follow_ups and the BCR XML follow-up records.
clin[, dtd_any  := pmax(days_to_death, dtd_xml, na.rm = TRUE)]
clin[, dtlf_any := pmax(dtlf, dtlf_fu, dtlf_xml, na.rm = TRUE)]
clin[, os_time  := ifelse(os_event == 1,
                          pmax(dtd_any, dtlf_any, na.rm = TRUE),
                          dtlf_any)]

# Administrative censoring at 10 years
clin[!is.na(os_time) & os_time > OS_CENSOR_DAYS, os_event := 0L]
clin[!is.na(os_time), os_time := pmin(os_time, OS_CENSOR_DAYS)]

clin[, stage_group := factor(ifelse(is.na(stage_num), NA,
                             ifelse(stage_num <= 2, "I-II", "III-IV")),
                             levels = c("I-II", "III-IV"))]
clin[, grade_group := factor(ifelse(is.na(grade_num), NA,
                             ifelse(grade_num <= 2, "G1-2", "G3-4")),
                             levels = c("G1-2", "G3-4"))]

# ---- 6. assemble the analysis cohort ----
dat <- merge(samp_t[, .(patient, sample_barcode, file_id, libsize,
                        pct_unmapped, pct_multimapping, pct_noFeature,
                        pct_ambiguous)],
             clin, by = "patient", all.x = TRUE)
# assigned_reads = sum of unstranded counts over ENSG genes. tss = tissue source
# site (barcode field 2). The plate is in the aliquot barcode, which the sample
# barcode lacks. Stage 23 fills it from the GDC biospecimen records.
dat[, assigned_reads := libsize]
dat[, tss := tstrsplit(sample_barcode, "-", keep = 2)[[1]]]
dat[, plate := NA_character_]
# Stage 02 residualises on the STAR summary, so it must be complete.
stopifnot(!anyNA(dat$pct_noFeature), !anyNA(dat$pct_multimapping),
          !anyNA(dat$pct_unmapped), !anyNA(dat$pct_ambiguous))

n0 <- nrow(dat)
drop_no_clin <- sum(is.na(dat$vital_status))
dat <- dat[!is.na(vital_status)]

drop_bad_time <- sum(is.na(dat$os_time) | dat$os_time <= 0)
dat <- dat[!is.na(os_time) & os_time > 0]

drop_no_stage <- sum(is.na(dat$stage_num))
drop_no_T     <- sum(is.na(dat$T_stage))
drop_no_grade <- sum(is.na(dat$grade_num))
drop_no_age   <- sum(is.na(dat$age))
drop_no_sex   <- sum(is.na(dat$sex))

consort$excluded_no_clinical_record  <- drop_no_clin
consort$excluded_nonpositive_os_time <- drop_bad_time
consort$missing_stage                <- drop_no_stage
consort$missing_T                    <- drop_no_T
consort$missing_grade                <- drop_no_grade
consort$missing_age                  <- drop_no_age
consort$missing_sex                  <- drop_no_sex
consort$m1_reconciled_from_stage_iv  <- sum(dat$M1_stage_reconciled == 1L)

# cohort_full: complete survival data. cohort_adj: also complete on the
# clinical-model inputs, used for every adjusted model. N_pos and M1 are never
# missing, so the restriction acts through T_stage, grade_num, age and sex.
cohort_full <- copy(dat)
cohort_adj  <- dat[complete.cases(dat[, .(os_time, os_event, age, sex,
                                         T_stage, N_pos, M1, grade_num)])]

consort$cohort_survival_analysis <- nrow(cohort_full)
consort$cohort_adjusted_models   <- nrow(cohort_adj)
msg("M1 reconciled from stage IV with T < 4: ",
    consort$m1_reconciled_from_stage_iv, " patients")

msg("Cohort (survival): n = ", nrow(cohort_full),
    ", events = ", sum(cohort_full$os_event))
msg("Cohort (adjusted): n = ", nrow(cohort_adj),
    ", events = ", sum(cohort_adj$os_event))

# ---- 7. split expression by biotype, filter, transform ----
# Raw columns are one per file, and some sample barcodes carry two files.
# Columns are selected by file_id so each matches the aliquot kept above.
.col_idx <- match(cohort_full$file_id, samp$file_id)
stopifnot(!anyNA(.col_idx), !anyDuplicated(.col_idx))
fpkm <- fpkm[, .col_idx, drop = FALSE]
colnames(fpkm) <- cohort_full$sample_barcode

split_and_filter <- function(biotype, min_fpkm, min_frac, top_n) {
  idx <- which(gene_ann$gene_type == biotype)
  m   <- fpkm[idx, , drop = FALSE]
  ann <- copy(gene_ann[idx])
  n_start <- nrow(m)

  expressed <- rowMeans(m >= min_fpkm, na.rm = TRUE) >= min_frac
  m <- m[expressed, , drop = FALSE]; ann <- ann[expressed]

  lg <- log2(m + 1)                                  # variance stabilisation
  if (is.finite(top_n) && nrow(lg) > top_n) {
    mad_v <- matrixStats::rowMads(lg)
    sel   <- order(mad_v, decreasing = TRUE)[seq_len(top_n)]
    lg <- lg[sel, , drop = FALSE]; ann <- ann[sel]
  }
  # Resolve duplicate symbols by keeping the most variable Ensembl entry
  ann[, mad_v := matrixStats::rowMads(lg)]
  setorder(ann, gene_name, -mad_v)
  keep <- !duplicated(ann$gene_name) & !is.na(ann$gene_name) & nzchar(ann$gene_name)
  lg <- lg[ann$gene_id[keep], , drop = FALSE]
  ann <- ann[keep]
  setkey(ann, gene_id); ann <- ann[rownames(lg)]

  msg(biotype, ": ", n_start, " -> ", nrow(lg), " genes after filtering")
  list(expr = lg, ann = ann)
}

mrna <- split_and_filter("protein_coding", MRNA_MIN_FPKM, MRNA_MIN_FRAC, MRNA_TOP_N_MAD)
lnc  <- split_and_filter("lncRNA",         LNC_MIN_FPKM,  LNC_MIN_FRAC,  LNC_TOP_N_MAD)

consort$genes_mrna_analysed   <- nrow(mrna$expr)
consort$genes_lncrna_analysed <- nrow(lnc$expr)

# ---- 8. save ----
saveRDS(list(mrna = mrna, lnc = lnc,
             cohort_full = cohort_full, cohort_adj = cohort_adj),
        file.path(CACHE_DIR, "dataset.rds"))

consort_df <- data.frame(step = names(consort),
                         n    = unlist(consort), row.names = NULL)
save_tsv(consort_df, "01_consort_flow.tsv")
print(consort_df)

# Median follow-up by reverse Kaplan-Meier (event indicator inverted).
fu_rkm <- survfit(Surv(os_time, 1 - os_event) ~ 1, data = cohort_full)
median_fu_rkm <- unname(summary(fu_rkm)$table["median"])

tbl1 <- cohort_full[, .(
  n              = .N,
  deaths         = sum(os_event),
  median_fu_days = round(median(os_time)),
  median_fu_reverse_km = round(median_fu_rkm),
  age_median     = round(median(age, na.rm = TRUE), 1),
  male_pct       = round(100 * mean(sex == "male", na.rm = TRUE), 1),
  stage_III_IV_pct = round(100 * mean(stage_group == "III-IV", na.rm = TRUE), 1),
  T3_T4_pct      = round(100 * mean(T_stage >= 3, na.rm = TRUE), 1),
  M1_pct_all     = round(100 * mean(M1 == 1), 1),
  M1_pct_assessed = round(100 * mean(M1_as_coded[!M1_imputed] == 1), 1),
  N_assessed     = sum(!N_pos_imputed),
  N_pos_pct_assessed = round(100 * mean(N_pos[!N_pos_imputed] == 1), 1),
  grade_G3_4_pct   = round(100 * mean(grade_group == "G3-4", na.rm = TRUE), 1),
  grade_available  = sum(!is.na(grade_num)),
  noFeature_median = round(median(pct_noFeature, na.rm = TRUE), 2),
  noFeature_IQR_lo = round(quantile(pct_noFeature, 0.25, na.rm = TRUE), 2),
  noFeature_IQR_hi = round(quantile(pct_noFeature, 0.75, na.rm = TRUE), 2),
  multimap_median  = round(median(pct_multimapping, na.rm = TRUE), 2),
  assigned_reads_median = round(median(assigned_reads, na.rm = TRUE))
)]
save_tsv(tbl1, "01_table1_cohort.tsv")
print(tbl1)

tnm <- cohort_full[, .(T_available = sum(!is.na(T_stage)),
                       N_assessed  = sum(!N_pos_imputed),
                       N_positive  = sum(N_pos == 1),
                       M_assessed  = sum(!M1_imputed),
                       M1_as_coded_count = sum(M1_as_coded == 1),
                       M1_count    = sum(M1 == 1),
                       m1_reconciled_from_stage_iv = sum(M1_stage_reconciled == 1L),
                       grade_available = sum(!is.na(grade_num)))]
save_tsv(tnm, "01_tnm_completeness.tsv")
print(tnm)

save_tsv(cohort_full[, .(patient, sample_barcode, age, sex, stage_group,
                         grade_group, T_stage, N_pos, M1, grade_num,
                         os_time, os_event,
                         assigned_reads, pct_noFeature, pct_multimapping,
                         pct_unmapped, pct_ambiguous, tss, M1_as_coded,
                         M1_stage_reconciled, N_pos_imputed, M1_imputed,
                         stage_num)],
         "01_patient_level_data.tsv")

write_session_info("01_build_data")
banner("01 | done")
