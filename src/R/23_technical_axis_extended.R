# 23_technical_axis_extended.R: biospecimen measurements behind the lncRNA axis.
#
# Relates the leading lncRNA axis and the non-feature read fraction to RIN,
# A260/A280, slide necrosis and composition, plate, batch and site. Sections:
#   1  GDC biospecimen fetch and BCR batch numbers
#   2  correlations and one-way variance by plate, batch and site
#   3  nested Cox models with RIN, necrosis and a plate stratum
#   4  the stage 12 with/without adjustment contrast, with RIN and necrosis
#   5  the lncRNA network rebuilt on the observed matrix
#   6  PC1 against per-sample mean expression
# Inputs (cache): dataset.rds, networks.rds, estimate_scores.rds,
# lnc_global_axis.rds. Caches biospecimen_kirc.rds and writes plate into
# dataset.rds. If the GDC API is unreachable, 23_NOTE_gdc_unreachable.txt
# records why.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(httr); library(xml2)
  library(survival); library(WGCNA)
})
banner("23 | Technical axis: biospecimen measurements and network demonstration")
set.seed(SEED)
enableWGCNAThreads(N_THREADS)
t_start <- Sys.time()

ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
est  <- as.data.table(readRDS(file.path(CACHE_DIR, "estimate_scores.rds")))
axis <- as.data.table(readRDS(file.path(CACHE_DIR, "lnc_global_axis.rds")))
cohort <- as.data.table(ds$cohort_full)
if (!"tss" %in% names(cohort))
  cohort[, tss := tstrsplit(sample_barcode, "-", keep = 2)[[1]]]
if (!"assigned_reads" %in% names(cohort)) cohort[, assigned_reads := libsize]
msg("Discovery cohort: ", nrow(cohort), " patients, ", sum(cohort$os_event), " deaths")

# Pool levels with fewer than min_n samples into "other".
pool_levels <- function(x, min_n = 10) {
  x <- as.character(x); tb <- table(x)
  x[x %in% names(tb)[tb < min_n]] <- "other"
  factor(x)
}
spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(c(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}
# PC1 of a samples x genes matrix: centred, unscaled, sign-aligned to mean
# expression (as in stage 07).
pc1_of <- function(E) {
  x <- scale(E, center = TRUE, scale = FALSE)
  s <- svd(x, nu = 1, nv = 1)
  pc <- s$u[, 1] * s$d[1]
  if (stats::cor(pc, rowMeans(E)) < 0) pc <- -pc
  list(pc = pc, var_share = s$d[1]^2 / sum(x^2))
}

# ---- 1. biospecimen measurements ----
banner("1 | GDC biospecimen fetch and BCR batch numbers")
bio_rds   <- file.path(CACHE_DIR, "biospecimen_kirc.rds")
bio_cache <- if (file.exists(bio_rds)) readRDS(bio_rds) else list()
gdc_error <- NULL

BIO_FIELDS <- paste(
  "submitter_id", "samples.submitter_id", "samples.sample_type",
  "samples.portions.submitter_id",
  "samples.portions.analytes.analyte_type",
  "samples.portions.analytes.submitter_id",
  "samples.portions.analytes.rna_integrity_number",
  "samples.portions.analytes.a260_a280_ratio",
  "samples.portions.analytes.aliquots.submitter_id",
  "samples.portions.slides.section_location",
  "samples.portions.slides.percent_necrosis",
  "samples.portions.slides.percent_tumor_nuclei",
  "samples.portions.slides.percent_stromal_cells",
  "samples.portions.slides.percent_normal_cells", sep = ",")

# The raw JSON is walked explicitly because the nesting (case > sample >
# portion > analyte > aliquot) defeats data-frame simplification.
gdc_raw <- function(endpoint, body) {
  r <- httr::POST(paste0("https://api.gdc.cancer.gov/", endpoint),
                  body = jsonlite::toJSON(body, auto_unbox = TRUE),
                  httr::content_type_json(), httr::timeout(300))
  httr::stop_for_status(r)
  j <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"),
                          simplifyVector = FALSE)
  if (!is.null(j$data$pagination$total) &&
      j$data$pagination$total > length(j$data$hits))
    warning("GDC returned ", length(j$data$hits), " of ",
            j$data$pagination$total, " records from /", endpoint, "; increase size")
  j$data$hits
}
fetch_biospecimen <- function() {
  filt <- list(op = "in", content = list(field = "project.project_id",
                                         value = list(PROJECT_ID)))
  body <- list(filters = filt, fields = BIO_FIELDS, format = "JSON", size = "600")
  hits <- tryCatch(gdc_raw("cases", body), error = function(e) e)
  if (inherits(hits, "error")) {
    # An unknown field fails the whole request, so retry with expand (all fields).
    msg("  field-list request failed (", conditionMessage(hits),
        "); retrying with expand")
    body <- list(filters = filt, fields = "submitter_id,samples.submitter_id,samples.sample_type",
                 expand = paste("samples", "samples.portions", "samples.portions.analytes",
                                "samples.portions.analytes.aliquots",
                                "samples.portions.slides", sep = ","),
                 format = "JSON", size = "600")
    hits <- gdc_raw("cases", body)
  }
  hits
}
# Sequenced aliquot of each STAR-Counts file (associated_entities), so plate
# and RIN describe that aliquot.
fetch_aliquot_map <- function(file_ids) {
  body <- list(filters = list(op = "in", content = list(field = "file_id",
                                                        value = as.list(unique(file_ids)))),
               fields = "file_id,associated_entities.entity_submitter_id,associated_entities.entity_type",
               format = "JSON", size = "600")
  hits <- gdc_raw("files", body)
  rbindlist(lapply(hits, function(h) {
    ae <- h$associated_entities
    if (!length(ae)) return(NULL)
    data.table(file_id = chr1(h$file_id),
               aliquot_id  = vapply(ae, function(e) chr1(e$entity_submitter_id), character(1)),
               entity_type = vapply(ae, function(e) chr1(e$entity_type), character(1)))
  }))
}

num1 <- function(x) if (is.null(x)) NA_real_ else suppressWarnings(as.numeric(x[[1]]))
chr1 <- function(x) if (is.null(x)) NA_character_ else as.character(x[[1]])
mean_na <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
max_na  <- function(x) if (any(is.finite(x))) max(x[is.finite(x)])  else NA_real_
modal   <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA_character_
                         else names(sort(table(x), decreasing = TRUE))[1] }
# Field k of a hyphenated TCGA barcode (aliquot: field 5 = portion + analyte
# letter, 6 = plate, 7 = centre).
bc_field <- function(x, k) vapply(strsplit(as.character(x), "-", fixed = TRUE), function(z)
  if (length(z) >= k) z[k] else NA_character_, character(1))
# Match 16-character sample barcodes, falling back to the first 15 characters.
match_sample <- function(x, y) {
  i <- match(x, y)
  if (all(is.na(i))) i <- match(substr(x, 1, 15), substr(y, 1, 15))
  i
}

# Returns one row per aliquot (analyte type, RIN, A260/A280) and one row per
# sample (slide percentages pooled over portions).
parse_biospecimen <- function(hits) {
  an <- list(); sl <- list()
  for (h in hits) for (s in h$samples) {
    sb <- chr1(s$submitter_id); pt <- chr1(h$submitter_id); st <- chr1(s$sample_type)
    slides <- list()
    for (p in s$portions) {
      for (a in p$analytes) {
        al <- vapply(a$aliquots, function(q) chr1(q$submitter_id), character(1))
        if (!length(al)) al <- NA_character_
        an[[length(an) + 1]] <- data.table(
          sample_barcode = sb, patient = pt, sample_type = st,
          analyte_id = chr1(a$submitter_id), analyte_type = chr1(a$analyte_type),
          rin = num1(a$rna_integrity_number), a260_a280 = num1(a$a260_a280_ratio),
          aliquot_id = al)
      }
      for (d in p$slides) slides[[length(slides) + 1]] <- data.table(
        section = chr1(d$section_location),
        pct_necrosis = num1(d$percent_necrosis),
        pct_tumor_nuclei = num1(d$percent_tumor_nuclei),
        pct_stromal = num1(d$percent_stromal_cells),
        pct_normal = num1(d$percent_normal_cells))
    }
    slides <- rbindlist(slides)
    sl[[length(sl) + 1]] <- data.table(
      sample_barcode = sb, patient = pt, sample_type = st, n_slides = nrow(slides),
      pct_necrosis_mean = if (nrow(slides)) mean_na(slides$pct_necrosis) else NA_real_,
      pct_necrosis_max  = if (nrow(slides)) max_na(slides$pct_necrosis)  else NA_real_,
      pct_tumor_nuclei_mean = if (nrow(slides)) mean_na(slides$pct_tumor_nuclei) else NA_real_,
      pct_stromal_mean = if (nrow(slides)) mean_na(slides$pct_stromal) else NA_real_,
      pct_normal_mean  = if (nrow(slides)) mean_na(slides$pct_normal)  else NA_real_)
  }
  list(analytes = rbindlist(an), slides = rbindlist(sl))
}

if (is.null(bio_cache$analytes)) {
  msg("Querying GDC /cases for biospecimen fields ...")
  got <- tryCatch({
    hits <- fetch_biospecimen()
    b <- parse_biospecimen(hits)
    if (!length(hits) || !nrow(b$analytes) || !nrow(b$slides))
      stop("GDC returned no biospecimen records")
    msg("  ", length(hits), " cases, ", nrow(b$slides), " samples, ",
        nrow(b$analytes), " aliquots returned")
    b
  }, error = function(e) e)
  if (inherits(got, "error")) {
    gdc_error <- conditionMessage(got)
    msg("GDC biospecimen fetch FAILED: ", gdc_error)
  } else {
    bio_cache$analytes <- got$analytes; bio_cache$slides <- got$slides
    bio_cache$gdc_fetched <- Sys.time()
    saveRDS(bio_cache, bio_rds)
  }
} else msg("Biospecimen tables from cache (fetched ", format(bio_cache$gdc_fetched), ")")

# On failure, plate and RIN fall back to the sample-level rule.
if (!is.null(bio_cache$analytes) && is.null(bio_cache$aliquot_map)) {
  msg("Querying GDC /files for the aliquot of each expression file ...")
  got <- tryCatch({
    am <- fetch_aliquot_map(cohort$file_id)
    if (!nrow(am)) stop("GDC returned no associated entities")
    msg("  aliquot recovered for ", uniqueN(am$file_id), " of ", uniqueN(cohort$file_id), " files")
    am
  }, error = function(e) e)
  if (inherits(got, "error")) {
    msg("GDC aliquot lookup FAILED (", conditionMessage(got),
        "): plate and RIN fall back to the sample-level rule")
  } else {
    bio_cache$aliquot_map <- got
    saveRDS(bio_cache, bio_rds)
  }
}

# Processing batch from the BCR clinical XML supplements (admin:batch_number).
if (is.null(bio_cache$batch)) {
  xmls <- list.files(CLIN_XML_DIR, pattern = "\\.xml$", full.names = TRUE,
                     recursive = TRUE)
  msg("Parsing ", length(xmls), " BCR clinical XML files for batch_number ...")
  ln <- function(x, nm) {
    v <- xml_text(xml_find_first(x, sprintf(".//*[local-name()='%s']", nm)))
    if (length(v) == 0 || is.na(v)) NA_character_ else trimws(v)
  }
  batch_xml <- rbindlist(lapply(xmls, function(f) {
    x <- tryCatch(read_xml(f), error = function(e) NULL)
    if (is.null(x)) return(NULL)
    data.table(patient = ln(x, "bcr_patient_barcode"), batch = ln(x, "batch_number"),
               is_clinical = grepl("clinical\\.", basename(f)))
  }))
  if (!nrow(batch_xml))
    batch_xml <- data.table(patient = character(0), batch = character(0),
                            is_clinical = logical(0))
  batch_xml <- batch_xml[!is.na(patient) & nzchar(patient) & !is.na(batch) & nzchar(batch)]
  setorder(batch_xml, patient, -is_clinical)
  bio_cache$batch <- batch_xml[!duplicated(patient), .(patient, batch)]
  saveRDS(bio_cache, bio_rds)
}
msg("Batch number recovered for ", nrow(bio_cache$batch), " patients")

# One row per cohort sample. Plate and RIN come from the sequenced aliquot
# where known. Otherwise the fallback is the mean RIN and A260/A280 over RNA
# analytes and the modal plate of aliquots with analyte letter R and centre 07.
bio_chr  <- c("aliquot_id", "plate", "plate_source", "rin_source")
bio_cols <- c(bio_chr, "n_plates", "n_rna_analytes", "n_rna_aliquots", "rin", "a260_a280",
              "n_slides", "pct_necrosis_mean", "pct_necrosis_max",
              "pct_tumor_nuclei_mean", "pct_stromal_mean", "pct_normal_mean")
slide_cols <- c("n_slides", "pct_necrosis_mean", "pct_necrosis_max",
                "pct_tumor_nuclei_mean", "pct_stromal_mean", "pct_normal_mean")
bio <- cohort[, .(sample_barcode, patient, tss, file_id)]
if (!is.null(bio_cache$analytes)) {
  an <- copy(bio_cache$analytes)
  an[, `:=`(analyte_letter = sub("^[0-9]+", "", bc_field(aliquot_id, 5)),
            plate_id = bc_field(aliquot_id, 6), centre = bc_field(aliquot_id, 7))]
  rna <- an[!is.na(analyte_type) & startsWith(analyte_type, "RNA")]
  rseq <- rna$analyte_letter %in% "R" & rna$centre %in% "07"
  fb <- rna[, .(n_rna_analytes = uniqueN(analyte_id),
                n_rna_aliquots = sum(!is.na(aliquot_id)),
                rin = mean_na(rin[!duplicated(analyte_id)]),
                a260_a280 = mean_na(a260_a280[!duplicated(analyte_id)]),
                plate = modal(plate_id[rseq[.I]]),
                n_plates = uniqueN(plate_id[rseq[.I]], na.rm = TRUE)),
            by = sample_barcode]
  if (!nrow(fb))
    fb <- data.table(sample_barcode = character(0), n_rna_analytes = integer(0),
                     n_rna_aliquots = integer(0), rin = numeric(0), a260_a280 = numeric(0),
                     plate = character(0), n_plates = integer(0))
  seq_al <- rep(NA_character_, nrow(bio))
  if (!is.null(bio_cache$aliquot_map) && nrow(bio_cache$aliquot_map)) {
    am <- bio_cache$aliquot_map[is.na(entity_type) | entity_type == "aliquot"]
    am <- am[!is.na(aliquot_id)][!duplicated(file_id)]
    seq_al <- am$aliquot_id[match(bio$file_id, am$file_id)]
  }
  ex <- an[match(seq_al, an$aliquot_id)]
  exact <- !is.na(seq_al) & !is.na(ex$aliquot_id)
  i_fb <- match_sample(bio$sample_barcode, fb$sample_barcode)
  i_sl <- match_sample(bio$sample_barcode, bio_cache$slides$sample_barcode)
  rin_ex <- ex$rin; a260_ex <- ex$a260_a280
  bio[, aliquot_id := seq_al]
  bio[, plate := fifelse(exact, bc_field(seq_al, 6), fb$plate[i_fb])]
  bio[, plate_source := fifelse(exact, "file_aliquot",
                        fifelse(!is.na(fb$plate[i_fb]), "fallback", NA_character_))]
  bio[, n_plates := fb$n_plates[i_fb]]
  bio[, n_rna_analytes := fb$n_rna_analytes[i_fb]]
  bio[, n_rna_aliquots := fb$n_rna_aliquots[i_fb]]
  bio[, rin := fifelse(exact & is.finite(rin_ex), rin_ex, fb$rin[i_fb])]
  bio[, rin_source := fifelse(exact & is.finite(rin_ex), "file_analyte",
                      fifelse(is.finite(fb$rin[i_fb]), "rna_analyte_mean", NA_character_))]
  bio[, a260_a280 := fifelse(exact & is.finite(a260_ex), a260_ex, fb$a260_a280[i_fb])]
  bio[, (slide_cols) := bio_cache$slides[i_sl, ..slide_cols]]
  msg("Plate from the sequenced aliquot for ", sum(exact), " samples, from the sample-level rule for ",
      sum(!exact & !is.na(bio$plate)), "; slides matched for ", sum(!is.na(i_sl)), " of ", nrow(bio))
} else {
  for (cc in bio_cols) set(bio, j = cc, value = if (cc %in% bio_chr) NA_character_ else NA_real_)
}
bio[, batch := bio_cache$batch$batch[match(patient, bio_cache$batch$patient)]]
have_gdc <- !is.null(bio_cache$analytes) && sum(is.finite(bio$rin)) >= 20
if (!have_gdc) msg("RIN unavailable for the cohort: GDC-dependent analyses are skipped")
save_tsv(bio, "23_biospecimen_kirc.tsv")

# Stage 01 leaves cohort_full$plate NA. Fill only that column.
if (!is.null(bio_cache$analytes) && any(!is.na(bio$plate))) {
  cf <- ds$cohort_full
  pl <- bio$plate[match(cf$sample_barcode, bio$sample_barcode)]
  if (is.data.table(cf)) { cf <- copy(cf); set(cf, j = "plate", value = pl) } else cf$plate <- pl
  ds$cohort_full <- cf
  saveRDS(ds, file.path(CACHE_DIR, "dataset.rds"))
  msg("plate written to cohort_full of cache/dataset.rds for ", sum(!is.na(pl)), " samples")
}

completeness <- rbindlist(lapply(
  c("aliquot_id", "rin", "a260_a280", "plate", "batch", "pct_necrosis_mean", "pct_necrosis_max",
    "pct_tumor_nuclei_mean", "pct_stromal_mean", "pct_normal_mean"), function(v) {
  x <- bio[[v]]; ok <- !is.na(x) & (!is.character(x) | nzchar(x))
  data.table(measure = v, n_cohort = nrow(bio), n_available = sum(ok),
             pct_available = round(100 * mean(ok), 1),
             n_levels = if (is.character(x)) length(unique(x[ok])) else NA_integer_,
             median = if (is.numeric(x)) round(median(x[ok]), 2) else NA_real_,
             iqr_lo = if (is.numeric(x)) round(quantile(x[ok], 0.25), 2) else NA_real_,
             iqr_hi = if (is.numeric(x)) round(quantile(x[ok], 0.75), 2) else NA_real_)
}))
save_tsv(completeness, "23_biospecimen_completeness.tsv"); print(completeness)

note_file <- file.path(RESULTS_DIR, "23_NOTE_gdc_unreachable.txt")
if (!is.null(gdc_error)) {
  writeLines(c("23_technical_axis_extended.R: the GDC cases endpoint could not be reached,",
               "so RIN, A260/A280, slide percentages and plate are missing and the",
               "analyses that depend on them were skipped.",
               paste0("Error: ", gdc_error), paste0("Time: ", format(Sys.time()))),
             note_file)
} else unlink(note_file)

# ---- 2. correlations, plate, batch and site ----
banner("2 | Correlations with biospecimen measurements; plate, batch and site")
# The placeholder plate column is dropped. The measured plate comes from bio.
d <- merge(cohort[, setdiff(names(cohort), "plate"), with = FALSE],
           est[, .(sample_barcode, StromalScore, ImmuneScore)],
           by = "sample_barcode", all.x = TRUE)
d <- merge(d, axis[, .(sample_barcode, lnc_axis)], by = "sample_barcode", all.x = TRUE)
d <- merge(d, bio[, c("sample_barcode", "batch", bio_cols), with = FALSE],
           by = "sample_barcode", all.x = TRUE)
d[, log_depth := log10(assigned_reads)]

measures <- c("rin", "a260_a280", "pct_necrosis_mean", "pct_tumor_nuclei_mean",
              "pct_stromal_mean", "StromalScore", "ImmuneScore")
variables <- c("pct_noFeature", "pct_multimapping", "log_depth", "lnc_axis")
cor_tbl <- rbindlist(lapply(variables, function(v) rbindlist(lapply(measures, function(m) {
  s <- spearman(d[[v]], d[[m]])
  data.table(variable = v, measure = m, n = s[["n"]],
             spearman_rho = round(s[["rho"]], 3), p = signif(s[["p"]], 3))
}))))
save_tsv(cor_tbl, "23_axis_vs_biospecimen.tsv")
print(dcast(cor_tbl, variable ~ measure, value.var = "spearman_rho"))

batch_var <- rbindlist(lapply(c("pct_noFeature", "lnc_axis"), function(v)
  rbindlist(lapply(c("plate", "batch", "tss"), function(f) {
    ok <- is.finite(d[[v]]) & !is.na(d[[f]]) & nzchar(d[[f]])
    if (sum(ok) < 20) return(data.table(variable = v, factor = f, n = sum(ok)))
    g <- pool_levels(d[[f]][ok]); y <- d[[v]][ok]
    if (nlevels(g) < 2) return(data.table(variable = v, factor = f, n = sum(ok),
                                          n_levels = nlevels(g)))
    data.table(variable = v, factor = f, n = sum(ok), n_levels = nlevels(g),
               n_pooled_other = sum(g == "other"),
               R2 = round(summary(lm(y ~ g))$r.squared, 4),
               anova_p = signif(anova(lm(y ~ g))[["Pr(>F)"]][1], 3),
               kruskal_p = signif(kruskal.test(y ~ g)$p.value, 3))
  }), fill = TRUE)), fill = TRUE)
save_tsv(batch_var, "23_batch_variance.tsv"); print(batch_var)

# ---- 3. nested Cox models and determinants of the non-feature fraction ----
banner("3 | Nested Cox models with RIN, necrosis and plate")
d[, `:=`(axis_z = as.numeric(scale(lnc_axis)), male = as.numeric(sex == "male"),
         stromal = as.numeric(scale(StromalScore)), immune = as.numeric(scale(ImmuneScore)),
         nf = as.numeric(scale(pct_noFeature)), nf_log = as.numeric(scale(log10(pct_noFeature))),
         mm = as.numeric(scale(pct_multimapping)),
         dep = as.numeric(scale(log_depth)),
         rin_z = as.numeric(scale(rin)), necrosis_z = as.numeric(scale(pct_necrosis_mean)))]
d[, plate_f := pool_levels(plate)]
d[, tss_f   := pool_levels(tss)]

# One complete-case set for all models, so added terms do not change the patients.
need <- c("os_time", "os_event", "axis_z", "age", "male", "T_stage", "N_pos", "M1",
          "grade_num", "stromal", "immune", "nf", "mm", "dep")
if (have_gdc) need <- c(need, "rin_z", "necrosis_z", "plate_f")
dc <- d[complete.cases(d[, ..need])]
msg("Complete cases: ", nrow(dc), " patients, ", sum(dc$os_event), " events",
    if (have_gdc) " (RIN, necrosis and plate available)" else " (STAR metrics only)")

fit_row <- function(label, rhs, term, dat) {
  fit <- coxph(as.formula(paste("Surv(os_time, os_event) ~", rhs)), data = dat)
  s  <- summary(fit)
  zp <- tryCatch(cox.zph(fit)$table[term, "p"], error = function(e) NA_real_)
  data.table(model = label, term = term, n = s$n, events = s$nevent,
             HR = round(s$conf.int[term, 1], 3), lo = round(s$conf.int[term, 3], 3),
             hi = round(s$conf.int[term, 4], 3), p = signif(s$coefficients[term, 5], 3),
             C = round(s$concordance[1], 3), ph_p = signif(zp, 3))
}
clin <- "age + male + T_stage + N_pos + M1 + grade_num"
axis_models <- list(
  `Axis alone`                                = "axis_z",
  `+ age, sex, T, N, M1, ordinal grade`       = paste("axis_z +", clin),
  `+ ESTIMATE stromal and immune scores`      = paste("axis_z +", clin, "+ stromal + immune"),
  `+ non-feature, multimapping, log10 depth`  = paste("axis_z +", clin, "+ stromal + immune + nf + mm + dep"))
if (have_gdc) axis_models <- c(axis_models, list(
  `+ RIN`                                     = paste("axis_z +", clin, "+ stromal + immune + rin_z"),
  `+ RIN + necrosis`                          = paste("axis_z +", clin, "+ stromal + immune + rin_z + necrosis_z"),
  `+ RIN + necrosis + STAR metrics`           = paste("axis_z +", clin, "+ stromal + immune + rin_z + necrosis_z + nf + mm + dep"),
  `+ RIN + necrosis + STAR metrics, stratified by plate` =
    paste("axis_z +", clin, "+ stromal + immune + rin_z + necrosis_z + nf + mm + dep + strata(plate_f)")))
axis_ext <- rbindlist(lapply(names(axis_models), function(n)
  fit_row(n, axis_models[[n]], "axis_z", dc)))
save_tsv(axis_ext, "23_axis_nested_extended.tsv"); print(axis_ext, row.names = FALSE)

nf_models <- list(
  `Non-feature fraction alone`                = "nf",
  `+ age, sex, T, N, M1, ordinal grade`       = paste("nf +", clin),
  `+ ESTIMATE stromal and immune scores`      = paste("nf +", clin, "+ stromal + immune"),
  `+ multimapping, log10 depth`               = paste("nf +", clin, "+ stromal + immune + mm + dep"),
  `+ ESTIMATE, stratified by tissue source site` =
    paste("nf +", clin, "+ stromal + immune + strata(tss_f)"))
if (have_gdc) nf_models <- c(nf_models, list(
  `+ RIN`                                     = paste("nf +", clin, "+ stromal + immune + rin_z"),
  `+ RIN + necrosis`                          = paste("nf +", clin, "+ stromal + immune + rin_z + necrosis_z"),
  `+ RIN + necrosis + multimapping, log10 depth` =
    paste("nf +", clin, "+ stromal + immune + rin_z + necrosis_z + mm + dep"),
  `+ RIN + necrosis + multimapping, log10 depth, stratified by plate` =
    paste("nf +", clin, "+ stromal + immune + rin_z + necrosis_z + mm + dep + strata(plate_f)")))
# Linear and log10 scale of the exposure (the fraction is right-skewed).
nf_models_log <- lapply(nf_models, function(r) trimws(sub("^nf ", "nf_log ", paste0(r, " "))))
nf_ext <- rbind(
  cbind(exposure_scale = "linear, per SD",
        rbindlist(lapply(names(nf_models), function(n) fit_row(n, nf_models[[n]], "nf", dc)))),
  cbind(exposure_scale = "log10, per SD",
        rbindlist(lapply(names(nf_models_log), function(n) fit_row(n, nf_models_log[[n]], "nf_log", dc)))))
save_tsv(nf_ext, "23_noFeature_exposure_extended.tsv"); print(nf_ext, row.names = FALSE)

# Non-feature fraction on RIN, plate and site: one-way and drop-one partial R2.
det_terms <- c(rin = "rin", plate = "plate_f", tss = "tss_f")
dd <- d[complete.cases(d[, .(pct_noFeature, rin, plate_f, tss_f)])]
dd[, `:=`(plate_f = droplevels(plate_f), tss_f = droplevels(tss_f))]
if (have_gdc && nrow(dd) >= 30) {
  full <- lm(pct_noFeature ~ rin + plate_f + tss_f, data = dd)
  r2_full <- summary(full)$r.squared
  determinants <- rbindlist(c(lapply(names(det_terms), function(tm) {
    alone <- lm(as.formula(paste("pct_noFeature ~", det_terms[[tm]])), data = dd)
    red   <- lm(as.formula(paste("pct_noFeature ~",
                                 paste(setdiff(det_terms, det_terms[[tm]]), collapse = " + "))),
                data = dd)
    a <- anova(red, full)
    data.table(term = tm, n = nrow(dd),
               df = if (tm == "rin") 1L else nlevels(dd[[det_terms[[tm]]]]) - 1L,
               R2_alone = round(summary(alone)$r.squared, 4),
               partial_R2 = round((deviance(red) - deviance(full)) / deviance(red), 4),
               p_alone = signif(anova(alone)[["Pr(>F)"]][1], 3),
               p_partial = signif(a[["Pr(>F)"]][2], 3))
  }), list(data.table(term = "rin + plate + tss", n = nrow(dd), df = full$rank - 1L,
                      R2_alone = round(r2_full, 4), partial_R2 = round(r2_full, 4),
                      p_alone = signif(pf(summary(full)$fstatistic[1], summary(full)$fstatistic[2],
                                          summary(full)$fstatistic[3], lower.tail = FALSE), 3),
                      p_partial = NA_real_))), fill = TRUE)
  save_tsv(determinants, "23_noFeature_determinants.tsv"); print(determinants)
} else {
  msg("Determinants of the non-feature fraction skipped (RIN/plate unavailable)")
}

# ---- 4. with/without technical adjustment, with RIN and necrosis ----
banner("4 | With/without adjustment contrast with RIN and necrosis as covariates")
# As in stage 12 section B: same module membership, eigengenes scored on the
# observed or the residualised matrix. Covariates are the clinical set, then
# the clinical set plus RIN and mean necrosis, on one complete-case set.
LOAD <- discovery_loadings(nets)
module_cox_arm <- function(net, tag, prefix, arm, cl_all) {
  E <- if (arm == "residualised") net$expr else obs_expr(net)
  Lm <- if (arm == "residualised") LOAD[grep(paste0("^", prefix), names(LOAD))]
        else fit_module_loadings(E, net$gene_tbl, prefix)
  M  <- score_modules(E, Lm)
  cl <- cl_all[match(rownames(M), sample_barcode)]
  Xc <- clinical_design(cl, set = "clinical")
  Xr <- cbind(Xc, rin = as.numeric(scale(cl$rin)),
              necrosis = as.numeric(scale(cl$pct_necrosis_mean)))
  ok <- complete.cases(Xr) & is.finite(cl$os_time) & is.finite(cl$os_event)
  y  <- Surv(cl$os_time[ok], cl$os_event[ok])
  rbindlist(lapply(c("clinical", "clinical_rin_necrosis"), function(cs) {
    X <- if (cs == "clinical") Xc[ok, , drop = FALSE] else Xr[ok, , drop = FALSE]
    rbindlist(lapply(colnames(M), function(nm) {
      dd <- data.frame(ME = as.numeric(scale(M[ok, nm])), X)
      s  <- summary(coxph(y ~ ., data = dd))
      data.table(biotype = tag, module = sub(paste0("^", prefix), "", nm),
                 covariate_set = cs, arm = arm, n = s$n, events = s$nevent,
                 HR = s$conf.int["ME", 1], lo = s$conf.int["ME", 3],
                 hi = s$conf.int["ME", 4], p = s$coefficients["ME", 5])
    }))
  }))
}
if (have_gdc) {
  cl_bio <- merge(cohort, bio[, .(sample_barcode, rin, pct_necrosis_mean)],
                  by = "sample_barcode", all.x = TRUE)
  sens <- rbindlist(lapply(c("residualised", "observed"), function(a) rbind(
    module_cox_arm(nets$mrna, "mRNA",   "mRNA_ME", a, cl_bio),
    module_cox_arm(nets$lnc,  "lncRNA", "lnc_ME",  a, cl_bio))))
  sens[, fdr := p.adjust(p, "BH"), by = .(biotype, covariate_set, arm)]
  sens[, arm := fifelse(arm == "residualised", "adjusted", "unadjusted")]
  wide <- dcast(sens, biotype + module + covariate_set ~ arm,
                value.var = c("HR", "lo", "hi", "p", "fdr"))
  wide <- merge(wide, sens[arm == "adjusted", .(biotype, module, covariate_set, n, events)],
                by = c("biotype", "module", "covariate_set"))
  wide[, status := fifelse(fdr_adjusted < FDR_ALPHA & fdr_unadjusted < FDR_ALPHA, "robust",
                   fifelse(fdr_adjusted >= FDR_ALPHA & fdr_unadjusted < FDR_ALPHA, "LOST on adjustment",
                   fifelse(fdr_adjusted < FDR_ALPHA & fdr_unadjusted >= FDR_ALPHA, "gained on adjustment",
                           "not significant")))]
  setorder(wide, covariate_set, biotype, fdr_adjusted)
  out4 <- wide[, .(biotype, module, covariate_set, n, events,
                   HR_unadjusted = round(HR_unadjusted, 3), lo_unadjusted = round(lo_unadjusted, 3),
                   hi_unadjusted = round(hi_unadjusted, 3), p_unadjusted = signif(p_unadjusted, 3),
                   fdr_unadjusted = signif(fdr_unadjusted, 3),
                   HR_adjusted = round(HR_adjusted, 3), lo_adjusted = round(lo_adjusted, 3),
                   hi_adjusted = round(hi_adjusted, 3), p_adjusted = signif(p_adjusted, 3),
                   fdr_adjusted = signif(fdr_adjusted, 3), status)]
  save_tsv(out4, "23_technical_adjustment_with_rin.tsv")
  print(out4[, .(biotype, module, covariate_set, n, events, HR_unadjusted, HR_adjusted, status)])
  for (cs in unique(out4$covariate_set))
    msg("  [", cs, "] significant without adjustment: ", sum(out4[covariate_set == cs, fdr_unadjusted < FDR_ALPHA]),
        "; with: ", sum(out4[covariate_set == cs, fdr_adjusted < FDR_ALPHA]))
} else {
  msg("With/without contrast with RIN and necrosis skipped (RIN unavailable)")
}

# ---- 5. lncRNA network on the observed matrix ----
banner("5 | lncRNA network on the un-residualised matrix (production parameters)")
E_obs <- obs_expr(nets$lnc)                         # samples x genes
E_res <- nets$lnc$expr
stopifnot(identical(rownames(E_obs), rownames(E_res)),
          identical(colnames(E_obs), colnames(E_res)))
lnc_merge <- if (!is.null(nets$lnc$merge_cut)) nets$lnc$merge_cut else LNC_MERGE
lnc_split <- if (!is.null(nets$lnc$deep_split)) nets$lnc$deep_split else LNC_DEEPSPLIT
net_u_rds <- file.path(CACHE_DIR, "network_lncRNA_unadjusted.rds")
fp_u <- list(dim = dim(E_obs), seed = SEED, power = LNC_POWER, deep_split = lnc_split,
             merge_cut = lnc_merge, min_kme = LNC_MINKME, net = NETWORK_TYPE,
             tom = TOM_TYPE, minmod = MIN_MODULE_SIZE, block = MAX_BLOCK_SIZE)
net_u <- NULL
if (file.exists(net_u_rds)) {
  cached <- readRDS(net_u_rds)
  if (identical(cached$fingerprint, fp_u)) { net_u <- cached$net; msg("Unadjusted network from cache") }
}
if (is.null(net_u)) {
  t0 <- Sys.time()
  set.seed(SEED)
  bw <- blockwiseModules(
    E_obs, power = LNC_POWER, networkType = NETWORK_TYPE, TOMType = TOM_TYPE,
    minModuleSize = MIN_MODULE_SIZE, mergeCutHeight = lnc_merge, deepSplit = lnc_split,
    minKMEtoStay = LNC_MINKME, numericLabels = FALSE, pamRespectsDendro = FALSE,
    maxBlockSize = MAX_BLOCK_SIZE, saveTOMs = FALSE, verbose = 0)
  colors_u <- setNames(bw$colors, colnames(E_obs))
  MEs_u <- orderMEs(moduleEigengenes(E_obs, colors = colors_u)$eigengenes)
  rownames(MEs_u) <- rownames(E_obs)
  net_u <- list(tag = "lncRNA_unadjusted", power = LNC_POWER, deep_split = lnc_split,
                merge_cut = lnc_merge, colors = colors_u, MEs = MEs_u,
                samples = rownames(E_obs), n_modules = length(setdiff(unique(colors_u), "grey")))
  saveRDS(list(fingerprint = fp_u, net = net_u), net_u_rds)
  msg("Unadjusted network built in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}
colors_u <- net_u$colors
colors_a <- nets$lnc$colors[colnames(E_obs)]
msg("Unadjusted: ", net_u$n_modules, " modules, grey fraction ", round(mean(colors_u == "grey"), 3),
    "; production: ", nets$lnc$n_modules, " modules, grey fraction ", round(mean(colors_a == "grey"), 3))

# Colours do not correspond between networks. The Jaccard table maps them.
pc_obs <- pc1_of(E_obs)
pc_res_share <- pc1_of(E_res)$var_share
nf_net <- cohort$pct_noFeature[match(rownames(E_obs), cohort$sample_barcode)]
kwithin <- function(E, colors) {
  k <- intramodularConnectivity.fromExpr(E, colors, power = LNC_POWER,
                                         networkType = NETWORK_TYPE)
  setNames(k$kWithin, colnames(E))
}
msg("Computing intramodular connectivity (both networks) ...")
kW_u <- kwithin(E_obs, colors_u)
kW_a <- kwithin(E_res, colors_a)

network_rows <- function(label, colors, MEs, kW, pc_share) {
  sizes <- table(colors); ng <- length(colors)
  big <- names(sizes)[sizes == max(sizes[names(sizes) != "grey"])][1]
  rbindlist(lapply(names(sizes), function(m) {
    me <- if (paste0("ME", m) %in% colnames(MEs)) as.numeric(MEs[[paste0("ME", m)]]) else NA_real_
    r1 <- if (all(is.finite(me))) cor.test(me, pc_obs$pc) else NULL
    r2 <- if (all(is.finite(me))) cor.test(me, nf_net) else NULL
    data.table(network = label, module = m, n_genes = as.integer(sizes[[m]]),
               frac_genes = round(sizes[[m]] / ng, 4),
               r_PC1_observed = if (is.null(r1)) NA_real_ else round(unname(r1$estimate), 3),
               p_PC1_observed = if (is.null(r1)) NA_real_ else signif(r1$p.value, 3),
               r_noFeature = if (is.null(r2)) NA_real_ else round(unname(r2$estimate), 3),
               p_noFeature = if (is.null(r2)) NA_real_ else signif(r2$p.value, 3),
               mean_kWithin = round(mean(kW[colors == m]), 3),
               n_modules = length(setdiff(names(sizes), "grey")),
               largest_module = big, largest_module_frac = round(max(sizes[names(sizes) != "grey"]) / ng, 4),
               grey_frac = round(mean(colors == "grey"), 4),
               mean_kWithin_nongrey = round(mean(kW[colors != "grey"]), 3),
               pc1_var_share_matrix = round(pc_share, 4))
  }))
}
net_summary <- rbind(network_rows("unadjusted", colors_u, net_u$MEs, kW_u, pc_obs$var_share),
                     network_rows("adjusted",   colors_a, nets$lnc$MEs, kW_a, pc_res_share))
setorder(net_summary, network, -n_genes)
save_tsv(net_summary, "23_unadjusted_lncRNA_network_summary.tsv")
print(net_summary[, .(network, module, n_genes, frac_genes, r_PC1_observed, r_noFeature, mean_kWithin)])

ov <- CJ(module_unadjusted = sort(unique(colors_u)), module_adjusted = sort(unique(colors_a)))
ov[, c("n_unadjusted", "n_adjusted", "n_overlap") := {
  gu <- names(colors_u)[colors_u == module_unadjusted]
  ga <- names(colors_a)[colors_a == module_adjusted]
  list(length(gu), length(ga), length(intersect(gu, ga)))
}, by = .(module_unadjusted, module_adjusted)]
ov[, jaccard := round(n_overlap / (n_unadjusted + n_adjusted - n_overlap), 4)]
ov[, frac_of_unadjusted := round(n_overlap / n_unadjusted, 4)]
ov[, frac_of_adjusted   := round(n_overlap / n_adjusted, 4)]
ov[, best_match_for_unadjusted := jaccard == max(jaccard), by = module_unadjusted]
ov[, best_match_for_adjusted   := jaccard == max(jaccard), by = module_adjusted]
setorder(ov, module_unadjusted, -jaccard)
save_tsv(ov, "23_unadjusted_vs_adjusted_module_overlap.tsv")
print(ov[best_match_for_unadjusted == TRUE & module_unadjusted != "grey",
         .(module_unadjusted, n_unadjusted, module_adjusted, n_adjusted, n_overlap, jaccard)])

# ---- 6. observed PC1 against per-sample mean expression ----
banner("6 | PC1 versus per-sample mean expression and lncRNA FPKM share")
msg("cor(local observed PC1, cached lnc_axis) = ",
    round(stats::cor(pc_obs$pc, axis$lnc_axis[match(rownames(E_obs), axis$sample_barcode)],
                     use = "pairwise.complete.obs"), 4))
mean_expr <- rowMeans(E_obs)
lnc_share <- rep(NA_real_, nrow(E_obs))
raw_rds <- file.path(CACHE_DIR, "expr_raw.rds")
if (file.exists(raw_rds)) {
  msg("Reading expr_raw.rds for the lncRNA share of total FPKM ...")
  raw <- readRDS(raw_rds)
  # raw$fpkm has one column per file, so select the aliquot stage 01 kept by file_id.
  j <- match(cohort$file_id[match(rownames(E_obs), cohort$sample_barcode)],
             raw$sample_info$file_id)
  is_lnc <- !is.na(raw$gene_ann$gene_type) & raw$gene_ann$gene_type == "lncRNA"
  tot <- colSums(raw$fpkm[, j[!is.na(j)], drop = FALSE], na.rm = TRUE)
  lnc <- colSums(raw$fpkm[is_lnc, j[!is.na(j)], drop = FALSE], na.rm = TRUE)
  lnc_share[!is.na(j)] <- lnc / tot
  rm(raw); invisible(gc(verbose = FALSE))
} else msg("expr_raw.rds not found: lncRNA FPKM share unavailable")

E_mc  <- E_obs - apply(E_obs, 1, median)            # per-sample median centring
pc_mc <- pc1_of(E_mc)
pc_row <- function(label, pc, var_share) {
  s_mean <- spearman(pc, mean_expr); s_share <- spearman(pc, lnc_share)
  s_nf <- spearman(pc, nf_net); s_pc <- spearman(pc, pc_obs$pc)
  data.table(matrix = label, n = length(pc), pc1_var_share = round(var_share, 4),
             pearson_r_mean_expr = round(stats::cor(pc, mean_expr), 3),
             rho_mean_expr = round(s_mean[["rho"]], 3), p_mean_expr = signif(s_mean[["p"]], 3),
             rho_lnc_fpkm_share = round(s_share[["rho"]], 3), p_lnc_fpkm_share = signif(s_share[["p"]], 3),
             n_lnc_fpkm_share = s_share[["n"]],
             rho_noFeature = round(s_nf[["rho"]], 3), p_noFeature = signif(s_nf[["p"]], 3),
             rho_PC1_observed = round(s_pc[["rho"]], 3))
}
pc_check <- rbind(pc_row("observed", pc_obs$pc, pc_obs$var_share),
                  pc_row("sample_median_centred", pc_mc$pc, pc_mc$var_share))
save_tsv(pc_check, "23_pc1_normalisation_check.tsv"); print(pc_check, row.names = FALSE)

msg("Elapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
write_session_info("23_technical_axis_extended")
banner("23 | done")
