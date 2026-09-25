# 28_mutations_and_lncRNA_classes.R: driver mutations versus the non-feature fraction, and positional classes of network lncRNAs
# Part A: driver status (VHL, PBRM1, BAP1, SETD2, KDM5C, MTOR, TP53) and non-synonymous
#   count from GDC masked MAFs, tested against the metric and lncRNA axis and in Cox
#   models. MAF coverage is partial and tracks the metric, so 28_mutation_coverage.tsv
#   compares covered and uncovered patients.
# Part B: GENCODE v36 positional class of each network lncRNA (exonic sense, antisense,
#   intronic, sense overlapping, divergent, intergenic), compared on metric rho and PC1 loading.
# Inputs: caches and stage 07 results. MAF and GTF downloads are cached in data/derived/cache.
#   If GDC or GENCODE is unreachable, that part is skipped and 28_NOTE_*.txt is written.
# Outputs: 28_*.tsv tables.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(httr); library(jsonlite); library(survival)
  library(GenomicRanges); library(IRanges); library(S4Vectors)
})
banner("28 | Driver mutations vs the non-feature fraction; lncRNA positional classes")
set.seed(SEED)
t_start <- Sys.time()

ds     <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets   <- readRDS(file.path(CACHE_DIR, "networks.rds"))
axis   <- as.data.table(readRDS(file.path(CACHE_DIR, "lnc_global_axis.rds")))
cohort <- as.data.table(ds$cohort_full)
if (!"assigned_reads" %in% names(cohort)) cohort[, assigned_reads := libsize]
msg("Discovery cohort: ", nrow(cohort), " patients, ", sum(cohort$os_event), " deaths; ",
    "lncRNA axis for ", nrow(axis), " samples; ", nrow(nets$lnc$gene_tbl), " network lncRNAs")

DRIVERS <- c("VHL", "PBRM1", "BAP1", "SETD2", "KDM5C", "MTOR", "TP53")
NONSYN  <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
             "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation",
             "Translation_Start_Site")
GDC_API <- "https://api.gdc.cancer.gov"
CLIN_RHS <- "age + male + T_stage + N_pos + M1 + grade_num"

num1 <- function(x) if (is.null(x)) NA_real_ else suppressWarnings(as.numeric(x[[1]]))
chr1 <- function(x) if (is.null(x)) NA_character_ else as.character(x[[1]])
with_retry <- function(f, tries = 3, wait = 5) {
  out <- NULL
  for (k in seq_len(tries)) {
    out <- tryCatch(f(), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    msg("    attempt ", k, " failed: ", conditionMessage(out))
    if (k < tries) Sys.sleep(wait * k)
  }
  stop(conditionMessage(out))
}
gdc_post_json <- function(endpoint, body, timeout = 300) {
  r <- httr::POST(paste0(GDC_API, "/", endpoint),
                  body = jsonlite::toJSON(body, auto_unbox = TRUE),
                  httr::content_type_json(), httr::timeout(timeout))
  httr::stop_for_status(r)
  j <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)
  if (!is.null(j$data$pagination$total) && j$data$pagination$total > length(j$data$hits))
    warning("GDC returned ", length(j$data$hits), " of ", j$data$pagination$total,
            " records from /", endpoint)
  j$data$hits
}
ph_p <- function(fit, term) tryCatch(cox.zph(fit)$table[term, "p"], error = function(e) NA_real_)
cox_row <- function(rhs, term, dat, label, extra = NULL) {
  fit <- coxph(as.formula(paste("Surv(os_time, os_event) ~", rhs)), data = dat)
  s <- summary(fit)
  # n_terms counts every estimated coefficient, the exposure included.
  n_terms <- length(coef(fit))
  out <- data.table(model = label, term = term, n = s$n, events = s$nevent,
                    n_terms = n_terms, epv = round(s$nevent / n_terms, 2),
                    HR = round(s$conf.int[term, 1], 3), lo = round(s$conf.int[term, 3], 3),
                    hi = round(s$conf.int[term, 4], 3), p = signif(s$coefficients[term, 5], 3),
                    C = round(unname(s$concordance[1]), 3), ph_p = signif(ph_p(fit, term), 3))
  if (!is.null(extra) && extra %in% rownames(s$conf.int))
    out[, `:=`(HR_noFeature = round(s$conf.int[extra, 1], 3),
               lo_noFeature = round(s$conf.int[extra, 3], 3),
               hi_noFeature = round(s$conf.int[extra, 4], 3),
               p_noFeature  = signif(s$coefficients[extra, 5], 3))]
  out
}

# ---- Part A. Somatic mutation status from the GDC masked MAFs ----
part_a <- function() {
banner("A | Somatic mutation status (GDC open-access masked MAFs)")
maf_dir  <- file.path(CACHE_DIR, "28_maf")
list_rds <- file.path(CACHE_DIR, "28_maf_file_list.rds")
mut_rds  <- file.path(CACHE_DIR, "28_mutation_status.rds")
note_gdc <- file.path(RESULTS_DIR, "28_NOTE_gdc_unreachable.txt")
if (!dir.exists(maf_dir)) dir.create(maf_dir, recursive = TRUE)

# ---- A1. File list (one masked MAF per tumour aliquot) ----------------------
fetch_maf_list <- function() {
  filters <- list(op = "and", content = list(
    list(op = "in", content = list(field = "cases.project.project_id", value = list(PROJECT_ID))),
    list(op = "in", content = list(field = "data_category", value = list("Simple Nucleotide Variation"))),
    list(op = "in", content = list(field = "data_type", value = list("Masked Somatic Mutation"))),
    list(op = "in", content = list(field = "access", value = list("open")))))
  body <- list(filters = filters, format = "JSON", size = "5000",
               fields = paste("file_id", "file_name", "file_size", "md5sum",
                              "experimental_strategy", "analysis.workflow_type",
                              "cases.submitter_id", "cases.samples.submitter_id",
                              "cases.samples.sample_type", sep = ","))
  hits <- gdc_post_json("files", body)
  if (!length(hits)) stop("GDC /files returned no masked MAF records")
  rbindlist(lapply(hits, function(h) {
    cs  <- h$cases
    smp <- if (length(cs)) cs[[1]]$samples else list()
    st  <- vapply(smp, function(s) chr1(s$sample_type), character(1))
    sb  <- vapply(smp, function(s) chr1(s$submitter_id), character(1))
    tum <- !is.na(st) & !grepl("Normal", st)
    data.table(file_id = chr1(h$file_id), file_name = chr1(h$file_name),
               file_size = num1(h$file_size), md5sum = chr1(h$md5sum),
               strategy = chr1(h$experimental_strategy),
               workflow = chr1(h$analysis$workflow_type),
               patient = if (length(cs)) chr1(cs[[1]]$submitter_id) else NA_character_,
               n_cases = length(cs),
               tumour_sample = if (any(tum)) sb[tum][1] else NA_character_,
               tumour_sample_type = if (any(tum)) st[tum][1] else NA_character_)
  }))
}

if (file.exists(mut_rds)) {
  mut_cache <- readRDS(mut_rds)
  msg("Mutation status from cache (parsed ", format(mut_cache$parsed_at), "; ",
      nrow(mut_cache$status), " patients)")
} else {
  if (file.exists(list_rds)) {
    fl <- readRDS(list_rds); msg("MAF file list from cache: ", nrow(fl), " files")
  } else {
    msg("Querying GDC /files for TCGA-KIRC masked somatic MAFs ...")
    got <- tryCatch(with_retry(fetch_maf_list), error = function(e) e)
    if (inherits(got, "error")) {
      writeLines(c("28_mutations_and_lncRNA_classes.R: the GDC files endpoint could not be reached,",
                   "so somatic mutation calls were not obtained and Part A was skipped.",
                   paste0("Error: ", conditionMessage(got)), paste0("Time: ", format(Sys.time()))),
                 note_gdc)
      msg("GDC unreachable (", conditionMessage(got), "); Part A skipped")
      return(invisible(FALSE))
    }
    fl <- got; saveRDS(fl, list_rds)
  }
  msg(nrow(fl), " MAF files for ", uniqueN(fl$patient), " patients (",
      round(sum(fl$file_size) / 1e6, 1), " MB); workflow: ",
      paste(names(table(fl$workflow)), collapse = " | "))
  print(table(fl$tumour_sample_type, useNA = "ifany"))

  # ---- A2. Download in batches through POST /data (tar.gz), GET fallback ----
  fl[, path := file.path(maf_dir, file_id, file_name)]
  todo <- fl[!file.exists(path)]
  if (nrow(todo)) {
    msg("Downloading ", nrow(todo), " MAF files in batches of 100 ...")
    starts <- seq(1, nrow(todo), by = 100)
    for (s in starts) {
      ids <- todo$file_id[s:min(s + 99, nrow(todo))]
      ok <- tryCatch({
        tf <- tempfile(fileext = ".tar.gz")
        with_retry(function() {
          r <- httr::POST(paste0(GDC_API, "/data"), body = jsonlite::toJSON(list(ids = ids)),
                          httr::content_type_json(), httr::timeout(900),
                          httr::write_disk(tf, overwrite = TRUE))
          httr::stop_for_status(r); r
        })
        untar(tf, exdir = maf_dir); unlink(tf); TRUE
      }, error = function(e) { msg("  batch failed: ", conditionMessage(e)); FALSE })
      msg("  batch ", which(starts == s), " / ", length(starts), ": ",
          if (ok) "extracted" else "FAILED (per-file fallback follows)")
    }
    unlink(file.path(maf_dir, "MANIFEST.txt"))
    miss <- fl[!file.exists(path)]
    if (nrow(miss)) {
      msg("Per-file GET for ", nrow(miss), " files still missing ...")
      for (i in seq_len(nrow(miss))) {
        dir.create(dirname(miss$path[i]), showWarnings = FALSE, recursive = TRUE)
        tryCatch(with_retry(function() {
          r <- httr::GET(paste0(GDC_API, "/data/", miss$file_id[i]), httr::timeout(300),
                         httr::write_disk(miss$path[i], overwrite = TRUE))
          httr::stop_for_status(r); r
        }), error = function(e) msg("  ", miss$file_id[i], " failed: ", conditionMessage(e)))
      }
    }
  } else msg("All ", nrow(fl), " MAF files present in cache/28_maf")
  fl[, on_disk := file.exists(path)]
  fl[, md5_ok := NA]
  fl[on_disk == TRUE, md5_ok := unname(tools::md5sum(path)) == md5sum]
  if (any(fl$on_disk & !(fl$md5_ok %in% TRUE))) {
    bad <- fl[on_disk == TRUE & !md5_ok %in% TRUE]
    msg(nrow(bad), " files fail the md5 check; re-downloading once ...")
    for (i in seq_len(nrow(bad))) tryCatch({
      r <- httr::GET(paste0(GDC_API, "/data/", bad$file_id[i]), httr::timeout(300),
                     httr::write_disk(bad$path[i], overwrite = TRUE)); httr::stop_for_status(r)
    }, error = function(e) NULL)
    fl[, on_disk := file.exists(path)]
    fl[on_disk == TRUE, md5_ok := unname(tools::md5sum(path)) == md5sum]
  }
  msg("On disk: ", sum(fl$on_disk), " / ", nrow(fl), "; md5 verified: ", sum(fl$md5_ok %in% TRUE))
  if (sum(fl$on_disk) == 0) {
    writeLines(c("28_mutations_and_lncRNA_classes.R: the GDC data endpoint could not be reached,",
                 "so no MAF was downloaded and Part A was skipped.", paste0("Time: ", format(Sys.time()))),
               note_gdc)
    return(invisible(FALSE))
  }

  # ---- A3. Parse ---------------------------------------------------------------
  # fread() on text read through a gz connection is much faster than on the .gz file.
  MAF_SEL <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
               "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
               "Variant_Type", "Tumor_Sample_Barcode", "HGVSp_Short")
  parse_maf <- function(p) tryCatch({
    con <- gzfile(p, "rt"); on.exit(close(con), add = TRUE)
    ln  <- readLines(con, warn = FALSE)
    h   <- which(startsWith(ln, "Hugo_Symbol"))[1]
    if (is.na(h)) return(NULL)
    if (h == length(ln)) return(data.table())          # header only, no variants
    fread(text = ln[h:length(ln)], sep = "\t", quote = "", showProgress = FALSE,
          select = MAF_SEL)
  }, error = function(e) NULL)
  msg("Parsing ", sum(fl$on_disk), " MAF files ...")
  # This list must not be called `parsed`: `fl` gains a column of that name and
  # data.table would resolve the RHS to the column.
  maf_tabs <- lapply(seq_len(nrow(fl)), function(i) {
    if (!fl$on_disk[i]) return(NULL)
    d <- parse_maf(fl$path[i])
    if (is.null(d)) return(NULL)
    if (nrow(d)) d[, `:=`(file_id = fl$file_id[i], patient = fl$patient[i])]
    d
  })
  n_rows_v <- vapply(maf_tabs, function(d) if (is.null(d)) NA_integer_ else nrow(d), integer(1))
  ok_v     <- !vapply(maf_tabs, is.null, logical(1))
  fl[, parsed := ok_v]
  fl[, n_rows := n_rows_v]
  var_all <- rbindlist(maf_tabs[which(ok_v & n_rows_v > 0)])
  msg("  parsed ", sum(fl$parsed), " files (", sum(!fl$parsed & fl$on_disk), " unreadable); ",
      nrow(var_all), " variant rows; ", sum(fl$n_rows == 0, na.rm = TRUE), " files with no variants")
  # tumour aliquot from the MAF itself, to cross-check the file metadata
  maf_sample <- var_all[, .(maf_tumour_sample = substr(Tumor_Sample_Barcode[1], 1, 16)), by = file_id]
  fl <- merge(fl, maf_sample, by = "file_id", all.x = TRUE)
  fl[is.na(maf_tumour_sample), maf_tumour_sample := substr(tumour_sample, 1, 16)]

  # Primary-tumour aliquots only. A patient with two primary aliquots
  # contributes the union of their calls, and the same variant seen in both
  # aliquots is counted once (key = chromosome, start, end, ref, alt).
  fl_use <- fl[parsed == TRUE & tumour_sample_type == "Primary Tumor" &
               substr(maf_tumour_sample, 14, 15) == "01"]
  msg("  primary-tumour MAFs used: ", nrow(fl_use), " (", uniqueN(fl_use$patient), " patients); ",
      "excluded as non-primary: ", sum(fl$parsed) - nrow(fl_use))
  var_all <- var_all[file_id %in% fl_use$file_id]
  var_all[, vkey := paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2)]
  var_u <- unique(var_all, by = c("patient", "vkey"))
  var_u[, nonsyn := Variant_Classification %in% NONSYN]
  print(table(var_u$Variant_Classification))

  status <- fl_use[, .(n_maf_files = .N,
                       maf_tumour_samples = paste(sort(unique(maf_tumour_sample)), collapse = ";")),
                   by = patient]
  counts <- var_u[, .(n_variants_all = .N, n_nonsyn = sum(nonsyn)), by = patient]
  status <- merge(status, counts, by = "patient", all.x = TRUE)
  status[is.na(n_variants_all), `:=`(n_variants_all = 0L, n_nonsyn = 0L)]
  for (g in DRIVERS)
    set(status, j = g, value = as.integer(status$patient %in% var_u[nonsyn == TRUE & Hugo_Symbol == g, patient]))
  drv <- var_u[nonsyn == TRUE & Hugo_Symbol %in% DRIVERS,
               .(patient, file_id, Hugo_Symbol, Chromosome, Start_Position, End_Position,
                 Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type, HGVSp_Short)]
  mut_cache <- list(status = status, driver_variants = drv,
                    file_list = fl[, .(file_id, file_name, file_size, md5sum, md5_ok, patient,
                                       tumour_sample, tumour_sample_type, maf_tumour_sample,
                                       on_disk, parsed, n_rows,
                                       used = file_id %in% fl_use$file_id)],
                    parsed_at = Sys.time())
  saveRDS(mut_cache, mut_rds)
}
unlink(note_gdc)
status <- copy(mut_cache$status)

# ---- A4. Merge to the cohort, coverage and frequencies ---------------------
mut <- merge(cohort[, .(patient, sample_barcode)], status, by = "patient", all = TRUE)
mut[, in_cohort := !is.na(sample_barcode)]
mut[, has_maf := !is.na(n_maf_files)]
same_smp <- rep(NA, nrow(mut))
ii <- which(mut$has_maf & mut$in_cohort)
if (length(ii))
  same_smp[ii] <- unname(mapply(function(s, m) isTRUE(s %in% strsplit(m, ";", fixed = TRUE)[[1]]),
                                mut$sample_barcode[ii], mut$maf_tumour_samples[ii]))
mut[, same_sample_as_expression := same_smp]
setcolorder(mut, c("patient", "in_cohort", "sample_barcode", "has_maf", "n_maf_files",
                   "maf_tumour_samples", "same_sample_as_expression", "n_variants_all", "n_nonsyn", DRIVERS))
setorder(mut, -in_cohort, patient)
save_tsv(mut, "28_mutation_status.tsv")

cov_tbl <- data.table(quantity = c(
  "cohort patients", "cohort patients with a masked MAF", "pct of cohort covered",
  "cohort deaths among covered patients", "MAF patients in total", "MAF patients not in the cohort",
  "MAF files listed", "MAF files parsed", "MAF files used (primary tumour)",
  "covered patients whose MAF aliquot is the sequenced expression sample",
  "median non-synonymous variants per covered patient", "IQR lo", "IQR hi"),
  value = c(nrow(cohort), sum(mut$in_cohort & mut$has_maf),
            round(100 * mean(mut[in_cohort == TRUE, has_maf]), 1),
            sum(cohort$os_event[cohort$patient %in% mut[in_cohort == TRUE & has_maf == TRUE, patient]]),
            nrow(status), sum(!mut$in_cohort & mut$has_maf),
            nrow(mut_cache$file_list), sum(mut_cache$file_list$parsed), sum(mut_cache$file_list$used),
            sum(mut[in_cohort == TRUE & has_maf == TRUE, same_sample_as_expression], na.rm = TRUE),
            median(mut[in_cohort == TRUE & has_maf == TRUE, n_nonsyn]),
            quantile(mut[in_cohort == TRUE & has_maf == TRUE, n_nonsyn], 0.25),
            quantile(mut[in_cohort == TRUE & has_maf == TRUE, n_nonsyn], 0.75)),
  covered = NA_real_, uncovered = NA_real_, p = NA_real_, test = NA_character_)

# ---- A4b. Covered versus uncovered patients ------------------------------------
# Deaths, stage and non-feature fraction with and without a MAF. Part A results are
# conditional on the covered subset.
cv <- merge(cohort[, .(patient, os_time, os_event, T_stage, M1, pct_noFeature)],
            mut[in_cohort == TRUE, .(patient, has_maf)], by = "patient")
stopifnot(nrow(cv) == nrow(cohort), !anyNA(cv$has_maf))
p_chi <- function(a, b) suppressWarnings(chisq.test(table(a, b))$p.value)
sd_lr <- survdiff(Surv(os_time, os_event) ~ has_maf, data = cv)
p_lr  <- pchisq(sd_lr$chisq, length(sd_lr$n) - 1L, lower.tail = FALSE)
p_nf  <- wilcox.test(pct_noFeature ~ has_maf, data = cv, exact = FALSE)$p.value
dec   <- quantile(cv$pct_noFeature, c(0.1, 0.9))
sel <- function(q, cov_v, unc_v, pv = NA_real_, tst = NA_character_)
  data.table(quantity = q, value = NA_real_, covered = cov_v, uncovered = unc_v,
             p = pv, test = tst)
cov_cmp <- rbind(
  sel("MAF coverage | patients (n)", sum(cv$has_maf), sum(!cv$has_maf)),
  sel("MAF coverage | deaths (n)", sum(cv[has_maf == TRUE, os_event]),
      sum(cv[has_maf == FALSE, os_event])),
  sel("MAF coverage | death rate (%)", round(100 * mean(cv[has_maf == TRUE, os_event]), 1),
      round(100 * mean(cv[has_maf == FALSE, os_event]), 1),
      signif(p_chi(cv$has_maf, cv$os_event), 3), "chi-square, death by MAF availability"),
  sel("MAF coverage | overall survival", NA_real_, NA_real_, signif(p_lr, 3),
      "log-rank, overall survival by MAF availability"),
  sel("MAF coverage | M1 (%)", round(100 * mean(cv[has_maf == TRUE, M1]), 1),
      round(100 * mean(cv[has_maf == FALSE, M1]), 1),
      signif(p_chi(cv$has_maf, cv$M1), 3), "chi-square"),
  sel("MAF coverage | T3 or T4 (%)", round(100 * mean(cv[has_maf == TRUE, T_stage >= 3]), 1),
      round(100 * mean(cv[has_maf == FALSE, T_stage >= 3]), 1),
      signif(p_chi(cv$has_maf, cv$T_stage >= 3), 3), "chi-square"),
  sel("MAF coverage | median pct_noFeature (%)",
      round(median(cv[has_maf == TRUE, pct_noFeature]), 3),
      round(median(cv[has_maf == FALSE, pct_noFeature]), 3),
      signif(p_nf, 3), "Wilcoxon rank-sum: THE EXPOSURE ITSELF differs by MAF availability"),
  sel("MAF coverage | SD of pct_noFeature (percentage points)",
      round(sd(cv[has_maf == TRUE, pct_noFeature]), 3),
      round(sd(cv[has_maf == FALSE, pct_noFeature]), 3), NA_real_,
      "per-SD hazard ratios stay comparable because the SDs are close"),
  data.table(quantity = "MAF coverage | pct of the top decile of pct_noFeature with a MAF",
             value = round(100 * mean(cv[pct_noFeature >= dec[2], has_maf]), 1),
             covered = NA_real_, uncovered = NA_real_, p = NA_real_,
             test = paste0("top decile = pct_noFeature >= ", round(dec[2], 2), "%")),
  data.table(quantity = "MAF coverage | pct of the bottom decile of pct_noFeature with a MAF",
             value = round(100 * mean(cv[pct_noFeature <= dec[1], has_maf]), 1),
             covered = NA_real_, uncovered = NA_real_, p = NA_real_,
             test = paste0("bottom decile = pct_noFeature <= ", round(dec[1], 2), "%")))
cov_tbl <- rbind(cov_tbl, cov_cmp)
save_tsv(cov_tbl, "28_mutation_coverage.tsv"); print(cov_tbl)
msg("Coverage is NOT at random: death rate ", cov_cmp[quantity %like% "death rate", covered], "% covered vs ",
    cov_cmp[quantity %like% "death rate", uncovered], "% uncovered (p = ",
    cov_cmp[quantity %like% "death rate", p], "); median pct_noFeature ",
    cov_cmp[quantity %like% "median pct_noFeature", covered], "% vs ",
    cov_cmp[quantity %like% "median pct_noFeature", uncovered], "% (Wilcoxon p = ",
    cov_cmp[quantity %like% "median pct_noFeature", p],
    ") -- every Part A null is conditional on this subset")

mc <- mut[in_cohort == TRUE & has_maf == TRUE]
freq <- rbindlist(lapply(DRIVERS, function(g) data.table(
  gene = g, n_mutated = sum(mc[[g]]), n_tested = nrow(mc), pct = round(100 * mean(mc[[g]]), 1))))
save_tsv(freq, "28_mutation_frequency.tsv"); print(freq)

# ---- A5. Non-feature fraction and lncRNA axis by mutation status -----------
d <- merge(cohort[, .(patient, sample_barcode, os_time, os_event, age, sex, T_stage, N_pos, M1,
                      grade_num, pct_noFeature, pct_multimapping, assigned_reads)],
           mc[, c("patient", "n_nonsyn", DRIVERS), with = FALSE], by = "patient")
ax <- copy(axis)[, axis_z := as.numeric(scale(lnc_axis))]      # scaled over the network samples
d <- merge(d, ax[, .(sample_barcode, axis_z)], by = "sample_barcode", all.x = TRUE)
d[, `:=`(male = as.numeric(sex == "male"),
         nf = as.numeric(scale(pct_noFeature)),               # per SD over the covered set
         tmb_z = as.numeric(scale(log1p(n_nonsyn))))]
msg("Covered cohort for the tests: ", nrow(d), " patients (", sum(!is.na(d$axis_z)), " with the axis)")

wilcox_gene <- function(y, ylab, g) {
  ok <- is.finite(d[[y]]); x <- d[[y]][ok]; m <- d[[g]][ok]
  p <- if (sum(m == 1) >= 2 && sum(m == 0) >= 2)
         wilcox.test(x ~ factor(m, levels = c(0, 1)), exact = FALSE)$p.value else NA_real_
  data.table(outcome = ylab, variable = g, test = "Wilcoxon rank-sum", n = sum(ok),
             n_mutated = sum(m == 1), n_wildtype = sum(m == 0),
             median_mutated = round(median(x[m == 1]), 3), median_wildtype = round(median(x[m == 0]), 3),
             difference_of_medians = round(median(x[m == 1]) - median(x[m == 0]), 3),
             rho = NA_real_, p = signif(p, 3))
}
spearman_tmb <- function(y, ylab) {
  ok <- is.finite(d[[y]])
  ct <- cor.test(d[[y]][ok], d$n_nonsyn[ok], method = "spearman", exact = FALSE)
  data.table(outcome = ylab, variable = "n_nonsyn (total non-synonymous count)", test = "Spearman",
             n = sum(ok), n_mutated = NA_integer_, n_wildtype = NA_integer_,
             median_mutated = NA_real_, median_wildtype = NA_real_, difference_of_medians = NA_real_,
             rho = round(unname(ct$estimate), 3), p = signif(ct$p.value, 3))
}
outcomes <- c(pct_noFeature = "pct_noFeature (%)", axis_z = "lncRNA axis (observed PC1, z)")
nf_mut <- rbindlist(lapply(names(outcomes), function(y) {
  r <- rbind(rbindlist(lapply(DRIVERS, function(g) wilcox_gene(y, outcomes[[y]], g))),
             spearman_tmb(y, outcomes[[y]]))
  r[, fdr := signif(p.adjust(p, "BH"), 3)]                  # family = the tests on one outcome
  r
}))
save_tsv(nf_mut, "28_noFeature_vs_mutations.tsv")
print(nf_mut[, .(outcome, variable, n, n_mutated, median_mutated, median_wildtype, rho, p, fdr)])

# ---- A6. Cox models ----------------------------------------------------------
# One complete-case set for all models, so adding a term never changes the patients.
# Mutations are coded 0/1, the metric and mutation count are per SD.
need <- c("os_time", "os_event", "age", "male", "T_stage", "N_pos", "M1", "grade_num", "nf", DRIVERS)
dc <- d[complete.cases(d[, ..need])]
dc[, `:=`(nf = as.numeric(scale(pct_noFeature)), tmb_z = as.numeric(scale(log1p(n_nonsyn))))]
msg("Cox complete cases: ", nrow(dc), " patients, ", sum(dc$os_event), " events")

mut_models <- rbindlist(lapply(c(DRIVERS, "tmb_z"), function(g) {
  lab <- if (g == "tmb_z") "log1p non-synonymous count, per SD" else paste0(g, " non-synonymous mutation")
  r <- rbind(
    cox_row(g, g, dc, "mutation alone"),
    cox_row(paste(g, "+", CLIN_RHS), g, dc, "+ age, sex, T, N, M1, ordinal grade"),
    cox_row(paste(g, "+", CLIN_RHS, "+ nf"), g, dc, "+ clinical + non-feature fraction", extra = "nf"),
    fill = TRUE)
  # n_mutated_in_model counts the Cox complete cases, not the covered cohort of
  # 28_mutation_frequency.tsv.
  cbind(block = "mutation as exposure", exposure = lab,
        n_mutated_in_model = if (g == "tmb_z") NA_integer_ else sum(dc[[g]]), r)
}), fill = TRUE)
mut_models[block == "mutation as exposure" & exposure != "log1p non-synonymous count, per SD",
           fdr := signif(p.adjust(p, "BH"), 3), by = model]

nf_models <- list(
  `nf + clinical (covered subset)`          = paste("nf +", CLIN_RHS),
  `+ BAP1`                                  = paste("nf +", CLIN_RHS, "+ BAP1"),
  `+ PBRM1`                                 = paste("nf +", CLIN_RHS, "+ PBRM1"),
  `+ SETD2`                                 = paste("nf +", CLIN_RHS, "+ SETD2"),
  `+ BAP1 + PBRM1 + SETD2`                  = paste("nf +", CLIN_RHS, "+ BAP1 + PBRM1 + SETD2"),
  `+ all seven driver genes`                = paste("nf +", CLIN_RHS, "+", paste(DRIVERS, collapse = " + ")),
  `+ mutation count (log1p, per SD)`        = paste("nf +", CLIN_RHS, "+ tmb_z"),
  `+ all seven driver genes + mutation count` =
    paste("nf +", CLIN_RHS, "+", paste(DRIVERS, collapse = " + "), "+ tmb_z"))
nf_rows <- rbindlist(lapply(names(nf_models), function(nm) cox_row(nf_models[[nm]], "nf", dc, nm)))
nf_rows <- cbind(block = "non-feature fraction as exposure", exposure = "pct_noFeature, per SD",
                 n_mutated_in_model = NA_integer_, nf_rows)
surv_tbl <- rbind(mut_models, nf_rows, fill = TRUE)
setcolorder(surv_tbl, c("block", "exposure", "model", "term", "n", "events", "n_mutated_in_model",
                        "n_terms", "epv", "HR", "lo", "hi", "p", "fdr", "C", "ph_p",
                        "HR_noFeature", "lo_noFeature", "hi_noFeature", "p_noFeature"))
save_tsv(surv_tbl, "28_mutation_survival_models.tsv")
print(surv_tbl[, .(block, exposure, model, n, events, n_mutated_in_model, n_terms, epv,
                   HR, lo, hi, p, HR_noFeature)])
big <- surv_tbl[which.max(n_terms)]
msg("Largest model: ", big$model, " -- ", big$n_terms, " estimated terms, ", big$events,
    " events, EPV ", big$epv, " (MIN_EPV = ", MIN_EPV, ")")
if (big$epv < MIN_EPV) msg("  WARNING: the largest model is below MIN_EPV")
msg("Proportional hazards: VHL ph_p over its three models = ",
    paste(surv_tbl[exposure %like% "^VHL", ph_p], collapse = ", "))
invisible(TRUE)
}

# ---- Part B. Positional classes of the network lncRNAs ----
part_b <- function() {
banner("B | Positional classes of the 3,442 network lncRNAs (GENCODE v36)")
gtf_gz     <- file.path(CACHE_DIR, "28_gencode_v36.gtf.gz")
parsed_rds <- file.path(CACHE_DIR, "28_gencode_v36_parsed.rds")
note_gc    <- file.path(RESULTS_DIR, "28_NOTE_gencode_unreachable.txt")
GENCODE_URL <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz"
GENCODE_MIN_BYTES <- 40e6

# ---- B1. Annotation: download (cached) and parse (cached) --------------------
if (!file.exists(parsed_rds)) {
  if (!file.exists(gtf_gz) || file.size(gtf_gz) < GENCODE_MIN_BYTES) {
    msg("Downloading GENCODE v36 comprehensive annotation (about 45 MB) ...")
    tmp <- paste0(gtf_gz, ".part")
    got <- tryCatch({
      r <- httr::GET(GENCODE_URL, httr::write_disk(tmp, overwrite = TRUE), httr::timeout(1800))
      httr::stop_for_status(r)
      if (file.size(tmp) < GENCODE_MIN_BYTES) stop("truncated download (", file.size(tmp), " bytes)")
      file.rename(tmp, gtf_gz); TRUE
    }, error = function(e) e)
    if (inherits(got, "error")) {
      msg("  httr download failed (", conditionMessage(got), "); trying download.file")
      got <- tryCatch({
        old <- getOption("timeout"); options(timeout = 1800); on.exit(options(timeout = old), add = TRUE)
        download.file(GENCODE_URL, tmp, mode = "wb", quiet = TRUE)
        if (file.size(tmp) < GENCODE_MIN_BYTES) stop("truncated download")
        file.rename(tmp, gtf_gz); TRUE
      }, error = function(e) e)
    }
    if (inherits(got, "error")) {
      writeLines(c("28_mutations_and_lncRNA_classes.R: the GENCODE v36 GTF could not be downloaded,",
                   "so the positional classification (Part B) was skipped.",
                   paste0("URL: ", GENCODE_URL), paste0("Error: ", conditionMessage(got)),
                   paste0("Time: ", format(Sys.time()))), note_gc)
      msg("GENCODE unreachable; Part B skipped")
      return(invisible(FALSE))
    }
  }
  unlink(note_gc)
  msg("Parsing ", basename(gtf_gz), " (gene, transcript and exon records) ...")
  con <- gzfile(gtf_gz, "rt"); hdr <- readLines(con, 50); close(con)
  n_skip <- sum(startsWith(hdr, "#"))
  g <- fread(gtf_gz, sep = "\t", header = FALSE, skip = n_skip, quote = "", showProgress = FALSE,
             col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attr"))
  g <- g[feature %in% c("gene", "transcript", "exon")]
  get_attr <- function(a, key) {
    has <- grepl(paste0(key, ' "'), a, fixed = TRUE)
    out <- rep(NA_character_, length(a))
    out[has] <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", a[has], perl = TRUE)
    out
  }
  g[, gene_id := get_attr(attr, "gene_id")]
  g[, gene_type := get_attr(attr, "gene_type")]
  g[, gene_name := get_attr(attr, "gene_name")]
  g[, transcript_id := get_attr(attr, "transcript_id")]
  g[, transcript_type := get_attr(attr, "transcript_type")]
  g[, attr := NULL]
  g <- g[!grepl("_PAR_Y", gene_id, fixed = TRUE)]     # pseudo-autosomal duplicates on chrY
  ann <- list(
    genes = g[feature == "gene", .(chr, start, end, strand, gene_id, gene_type, gene_name)],
    tx    = g[feature == "transcript", .(chr, start, end, strand, gene_id, transcript_id, transcript_type)],
    exons = g[feature == "exon", .(chr, start, end, strand, gene_id, transcript_id)],
    source = basename(gtf_gz), parsed_at = Sys.time())
  rm(g); invisible(gc(verbose = FALSE))
  saveRDS(ann, parsed_rds)
} else {
  ann <- readRDS(parsed_rds)
  msg("GENCODE v36 tables from cache (", nrow(ann$genes), " genes, ", nrow(ann$tx),
      " transcripts, ", nrow(ann$exons), " exons)")
}
genes <- ann$genes; tx <- ann$tx; exons <- ann$exons
msg("GENCODE v36: ", nrow(genes), " genes (", sum(genes$gene_type == "protein_coding"),
    " protein-coding, ", sum(genes$gene_type == "lncRNA"), " lncRNA)")

# ---- B2. Gene sets --------------------------------------------------------------
gt <- as.data.table(nets$lnc$gene_tbl)
lnc_ids <- gt$gene_id
hit_v <- lnc_ids %in% genes$gene_id
if (!all(hit_v)) {
  # fall back to the unversioned identifier for any lncRNA whose version differs
  unv <- sub("\\..*$", "", genes$gene_id)
  alt <- genes$gene_id[match(sub("\\..*$", "", lnc_ids[!hit_v]), unv)]
  msg("  ", sum(!hit_v), " lncRNA identifiers not matched with version; ", sum(!is.na(alt)),
      " matched without version")
  lnc_ids[!hit_v] <- alt
}
id_map <- data.table(gene_id = gt$gene_id, gencode_id = lnc_ids)
genes_lnc <- genes[match(lnc_ids, gene_id)]
stopifnot(!anyNA(genes_lnc$gene_id))
genes_pc <- genes[gene_type == "protein_coding"]
tx_lnc <- tx[gene_id %in% genes_lnc$gene_id]
tx_pc  <- tx[gene_id %in% genes_pc$gene_id]
ex_lnc <- exons[gene_id %in% genes_lnc$gene_id]
ex_pc  <- exons[gene_id %in% genes_pc$gene_id]          # every transcript of a protein-coding gene
msg("lncRNA set: ", nrow(genes_lnc), " genes, ", nrow(tx_lnc), " transcripts, ", nrow(ex_lnc),
    " exons; protein-coding: ", nrow(genes_pc), " genes, ", nrow(tx_pc), " transcripts, ",
    nrow(ex_pc), " exons")

mkgr <- function(d) GRanges(d$chr, IRanges(d$start, d$end), strand = d$strand)
lnc_gr <- mkgr(genes_lnc); lnc_gr$gene_id <- genes_lnc$gene_id
pc_gr  <- mkgr(genes_pc);  pc_gr$gene_id <- genes_pc$gene_id; pc_gr$gene_name <- genes_pc$gene_name
lnc_ex <- mkgr(ex_lnc);    lnc_ex$gene_id <- ex_lnc$gene_id
pc_ex  <- mkgr(ex_pc);     pc_ex$gene_id <- ex_pc$gene_id; pc_ex$transcript_id <- ex_pc$transcript_id
pc_name <- setNames(genes_pc$gene_name, genes_pc$gene_id)

# Introns of every protein-coding transcript: the transcript span minus its exons.
msg("Building protein-coding transcript introns ...")
exl   <- split(pc_ex, pc_ex$transcript_id)
rng   <- range(exl)                       # one range per transcript (same chr/strand)
keep_tx <- lengths(rng) == 1L
if (!all(keep_tx)) {
  msg("  dropping ", sum(!keep_tx), " transcripts whose exons are not on one chromosome/strand")
  exl <- exl[keep_tx]; rng <- rng[keep_tx]
}
txr    <- unlist(rng, use.names = TRUE)
intr_l <- psetdiff(txr, exl)
names(intr_l) <- names(exl)
intr  <- unlist(intr_l, use.names = FALSE)
intr$transcript_id <- rep(names(intr_l), lengths(intr_l))
tx2g  <- setNames(tx_pc$gene_id, tx_pc$transcript_id)
intr$gene_id <- unname(tx2g[intr$transcript_id])
msg("  ", length(intr), " introns from ", length(exl), " transcripts")

# ---- B3. Overlap flags per lncRNA gene --------------------------------------------
strand_chr <- function(x) as.character(strand(x))
# exon-exon overlap with protein-coding exons, by strand relation
ov <- findOverlaps(lnc_ex, pc_ex, ignore.strand = TRUE)
ex_hits <- data.table(gene_id = lnc_ex$gene_id[queryHits(ov)],
                      pc_gene = pc_ex$gene_id[subjectHits(ov)],
                      same = strand_chr(lnc_ex)[queryHits(ov)] == strand_chr(pc_ex)[subjectHits(ov)])
ex_flag <- ex_hits[, .(exon_overlap_sense = any(same), exon_overlap_antisense = any(!same),
                       exon_partner_sense = paste(unique(pc_name[pc_gene[same]]), collapse = ";"),
                       exon_partner_antisense = paste(unique(pc_name[pc_gene[!same]]), collapse = ";")),
                   by = gene_id]
# gene span entirely within one protein-coding transcript intron
ov <- findOverlaps(lnc_gr, intr, type = "within", ignore.strand = TRUE)
in_hits <- data.table(gene_id = lnc_gr$gene_id[queryHits(ov)],
                      host = intr$gene_id[subjectHits(ov)],
                      same = strand_chr(lnc_gr)[queryHits(ov)] == strand_chr(intr)[subjectHits(ov)])
in_flag <- in_hits[, .(within_intron = TRUE, intron_host_same_strand = any(same),
                       intron_host = paste(unique(pc_name[host]), collapse = ";")), by = gene_id]
# gene-span overlap with any protein-coding gene span, by strand relation
ov <- findOverlaps(lnc_gr, pc_gr, ignore.strand = TRUE)
sp_hits <- data.table(gene_id = lnc_gr$gene_id[queryHits(ov)],
                      pc_gene = pc_gr$gene_id[subjectHits(ov)],
                      same = strand_chr(lnc_gr)[queryHits(ov)] == strand_chr(pc_gr)[subjectHits(ov)])
sp_flag <- sp_hits[, .(span_overlap_sense = any(same), span_overlap_antisense = any(!same),
                       span_partners = paste(unique(pc_name[pc_gene]), collapse = ";")), by = gene_id]
# divergent: a transcript TSS within 1 kb of a protein-coding transcript TSS on the
# opposite strand, in head-to-head geometry (the coding gene lies upstream)
tss_of <- function(d) ifelse(d$strand == "+", d$start, d$end)
lt <- tss_of(tx_lnc); pt <- tss_of(tx_pc)
lnc_tss <- GRanges(tx_lnc$chr, IRanges(lt, lt), strand = tx_lnc$strand); lnc_tss$gene_id <- tx_lnc$gene_id
pc_tss  <- GRanges(tx_pc$chr,  IRanges(pt, pt), strand = tx_pc$strand);  pc_tss$gene_id  <- tx_pc$gene_id
win <- GRanges(tx_lnc$chr, IRanges(pmax(1L, lt - 1000L), lt + 1000L))
ov <- findOverlaps(win, pc_tss, ignore.strand = TRUE)
q <- queryHits(ov); s <- subjectHits(ov)
sl <- strand_chr(lnc_tss)[q]; sp <- strand_chr(pc_tss)[s]
tl <- start(lnc_tss)[q]; tp <- start(pc_tss)[s]
h2h <- (sl == "+" & sp == "-" & tp <= tl) | (sl == "-" & sp == "+" & tp >= tl)
div_hits <- data.table(gene_id = lnc_tss$gene_id[q], partner = unname(pc_name[pc_tss$gene_id[s]]),
                       dist = abs(tp - tl), h2h = h2h)
div_hits <- div_hits[h2h == TRUE][, h2h := NULL]
setorder(div_hits, gene_id, dist)
div_flag <- div_hits[!duplicated(gene_id), .(gene_id, divergent_tss_partner = partner,
                                              divergent_tss_distance = dist)]
# distance to the nearest protein-coding gene (span, either strand)
dn <- distanceToNearest(lnc_gr, pc_gr, ignore.strand = TRUE)
near <- data.table(gene_id = lnc_gr$gene_id[queryHits(dn)],
                   nearest_pc_gene = pc_gr$gene_name[subjectHits(dn)],
                   nearest_pc_distance = mcols(dn)$distance,
                   nearest_pc_same_strand = strand_chr(lnc_gr)[queryHits(dn)] == strand_chr(pc_gr)[subjectHits(dn)])

# ---- B4. Gene features and the class call ----------------------------------------
n_tx <- tx_lnc[, .(n_transcripts = .N), by = gene_id]
ex_union <- lengths(reduce(split(lnc_ex, lnc_ex$gene_id)))
n_ex_tx <- ex_lnc[, .(n = .N), by = .(gene_id, transcript_id)][, .(max_exons_per_transcript = max(n)), by = gene_id]
cls <- data.table(gene_id = id_map$gene_id, gencode_id = id_map$gencode_id)
cls <- merge(cls, genes_lnc[, .(gencode_id = gene_id, gene_name, gencode_gene_type = gene_type,
                                chr, start, end, strand, gene_length = end - start + 1L)],
             by = "gencode_id", all.x = TRUE)
cls <- merge(cls, n_tx[, .(gencode_id = gene_id, n_transcripts)], by = "gencode_id", all.x = TRUE)
cls[, n_exons_union := unname(ex_union[gencode_id])]
cls <- merge(cls, n_ex_tx[, .(gencode_id = gene_id, max_exons_per_transcript)], by = "gencode_id", all.x = TRUE)
for (f in list(ex_flag, in_flag, sp_flag, div_flag, near))
  cls <- merge(cls, setnames(copy(f), "gene_id", "gencode_id"), by = "gencode_id", all.x = TRUE)
for (v in c("exon_overlap_sense", "exon_overlap_antisense", "within_intron", "intron_host_same_strand",
            "span_overlap_sense", "span_overlap_antisense"))
  set(cls, which(is.na(cls[[v]])), v, FALSE)
cls[, divergent := !is.na(divergent_tss_distance) & !span_overlap_sense & !span_overlap_antisense]
cls[, class_detail := fifelse(exon_overlap_sense, "exonic_sense",
                      fifelse(exon_overlap_antisense, "exonic_antisense",
                      fifelse(within_intron & intron_host_same_strand, "intronic_sense",
                      fifelse(within_intron, "intronic_antisense",
                      fifelse(span_overlap_antisense, "antisense_overlapping",
                      fifelse(span_overlap_sense, "sense_overlapping",
                      fifelse(divergent, "divergent", "intergenic")))))))]
CLASSES <- c("exonic_sense", "antisense", "intronic", "sense_overlapping", "divergent", "intergenic")
cls[, class := fifelse(class_detail == "exonic_sense", "exonic_sense",
              fifelse(class_detail %in% c("exonic_antisense", "antisense_overlapping"), "antisense",
              fifelse(class_detail %in% c("intronic_sense", "intronic_antisense"), "intronic",
                      class_detail)))]
cls[, class := factor(class, levels = CLASSES)]
# the name-based class used in stage 10, for comparison
name_class <- function(g) fifelse(grepl("-AS[0-9]*$", g), "antisense",
                          fifelse(grepl("-DT$", g), "divergent",
                          fifelse(grepl("-IT[0-9]*$", g), "intronic",
                          fifelse(grepl("^LINC", g), "lincRNA (named)",
                          fifelse(grepl("^(AC|AL|AP|Z|BX|CT|FP)[0-9]", g),
                                  "clone-based (unnamed)", "other")))))
cls[, name_class := name_class(gene_name)]

# ---- B5. Quality correlation, PC1 loading, expression, module -----------------
qc07 <- fread(file.path(RESULTS_DIR, "07_per_gene_quality_correlation.tsv"))[biotype == "lncRNA"]
X <- obs_expr(nets$lnc)
Xc <- scale(X, center = TRUE, scale = FALSE)
sv <- svd(Xc, nu = 1, nv = 1)
sc <- sv$u[, 1]; ld <- sv$v[, 1]
if (cor(sc, rowMeans(X)) < 0) { sc <- -sc; ld <- -ld }
r_axis <- cor(sc, axis$lnc_axis[match(rownames(X), axis$sample_barcode)])
# Per-gene quantities are measured on the lncRNA network samples, a subset of
# the cohort. N_LNC_SAMPLES records that denominator in the output tables.
N_LNC_SAMPLES <- nrow(X)
nf_sample <- cohort$pct_noFeature[match(rownames(X), cohort$sample_barcode)]
stopifnot(!anyNA(nf_sample))
ct_axis_nf <- cor.test(sc, nf_sample, method = "spearman", exact = FALSE)
msg("Observed lncRNA PC1: variance share ", round(sv$d[1]^2 / sum(Xc^2), 4),
    "; correlation of the recomputed score with cache/lnc_global_axis.rds = ", round(r_axis, 5),
    "; n = ", N_LNC_SAMPLES, " libraries")
msg("Sample-level Spearman(PC1 score, pct_noFeature) = ", round(unname(ct_axis_nf$estimate), 3),
    " -- the per-gene PC1 loading is therefore close to an algebraic restatement of the ",
    "per-gene correlation with the non-feature fraction, not independent corroboration of it")
stopifnot(r_axis > 0.999)
feat <- data.table(gene_id = colnames(X), pc1_loading = signif(ld, 5),
                   mean_log2_expr = round(colMeans(X), 4), sd_log2_expr = round(apply(X, 2, sd), 4),
                   frac_expressed = round(colMeans(X > log2(1 + LNC_MIN_FPKM)), 4))
cls <- merge(cls, gt[, .(gene_id, module, kME = round(kME, 4))], by = "gene_id", all.x = TRUE)
cls <- merge(cls, qc07[, .(gene_id, rho_noFeature)], by = "gene_id", all.x = TRUE)
cls <- merge(cls, feat, by = "gene_id", all.x = TRUE)
cls[, abs_rho := abs(rho_noFeature)]
stopifnot(nrow(cls) == nrow(gt), !anyNA(cls$rho_noFeature), !anyNA(cls$pc1_loading))
setcolorder(cls, c("gene_id", "gene_name", "gencode_id", "gencode_gene_type", "chr", "start", "end",
                   "strand", "gene_length", "n_transcripts", "n_exons_union", "max_exons_per_transcript",
                   "class", "class_detail", "name_class",
                   "exon_overlap_sense", "exon_overlap_antisense", "within_intron",
                   "intron_host_same_strand", "span_overlap_sense", "span_overlap_antisense",
                   "divergent", "exon_partner_sense", "exon_partner_antisense", "intron_host",
                   "span_partners", "divergent_tss_partner", "divergent_tss_distance",
                   "nearest_pc_gene", "nearest_pc_distance", "nearest_pc_same_strand",
                   "module", "kME", "rho_noFeature", "abs_rho", "pc1_loading",
                   "mean_log2_expr", "sd_log2_expr", "frac_expressed"))
setorder(cls, class, gene_name)
save_tsv(cls, "28_lncRNA_positional_classes.tsv")
msg("Class counts:"); print(table(cls$class)); print(table(cls$class_detail))

xt <- as.data.table(table(positional_class = cls$class, name_class = cls$name_class))
xt <- dcast(xt, positional_class ~ name_class, value.var = "N")
save_tsv(xt, "28_positional_vs_name_class.tsv"); print(xt)

# ---- B6. Per-class summary of the quality signal, length and expression ------
class_summary <- function(dt, grouping) {
  dt[, .(grouping = grouping, n_genes = .N, pct_of_set = round(100 * .N / nrow(cls), 1),
         n_samples = N_LNC_SAMPLES,
         median_abs_rho = round(median(abs_rho), 3),
         q25_abs_rho = round(quantile(abs_rho, 0.25), 3), q75_abs_rho = round(quantile(abs_rho, 0.75), 3),
         mean_rho = round(mean(rho_noFeature), 3),
         frac_abs_rho_gt_0.3 = round(mean(abs_rho > 0.3), 3),
         frac_abs_rho_gt_0.5 = round(mean(abs_rho > 0.5), 3),
         frac_rho_gt_0.3 = round(mean(rho_noFeature > 0.3), 3),
         frac_rho_lt_minus0.3 = round(mean(rho_noFeature < -0.3), 3),
         mean_PC1_loading = signif(mean(pc1_loading), 4),
         median_PC1_loading = signif(median(pc1_loading), 4),
         mean_abs_PC1_loading = signif(mean(abs(pc1_loading)), 4),
         frac_PC1_loading_positive = round(mean(pc1_loading > 0), 3),
         median_gene_length_kb = round(median(gene_length) / 1000, 2),
         q25_gene_length_kb = round(quantile(gene_length, 0.25) / 1000, 2),
         q75_gene_length_kb = round(quantile(gene_length, 0.75) / 1000, 2),
         # as.numeric(): median() of integers returns integer or double by group
         # size, which data.table rejects as inconsistent across by-groups.
         median_n_exons = as.numeric(median(n_exons_union)),
         median_n_transcripts = as.numeric(median(n_transcripts)),
         median_mean_log2_expr = round(median(mean_log2_expr), 3),
         q25_mean_log2_expr = round(quantile(mean_log2_expr, 0.25), 3),
         q75_mean_log2_expr = round(quantile(mean_log2_expr, 0.75), 3),
         median_frac_expressed = round(median(frac_expressed), 3),
         frac_grey = round(mean(module == "grey"), 3),
         median_nearest_pc_distance_kb =
           round(as.numeric(median(nearest_pc_distance, na.rm = TRUE)) / 1000, 1)),
     by = grp]
}
by_class <- rbind(
  class_summary(copy(cls)[, grp := "all"], "all"),
  class_summary(copy(cls)[, grp := as.character(class)], "class"),
  class_summary(copy(cls)[, grp := class_detail], "class_detail"))
setnames(by_class, "grp", "class")
setcolorder(by_class, c("grouping", "class"))
save_tsv(by_class, "28_quality_correlation_by_class.tsv")
print(by_class[, .(grouping, class, n_genes, n_samples, median_abs_rho, frac_abs_rho_gt_0.3,
                   mean_PC1_loading, median_gene_length_kb, median_n_exons, median_mean_log2_expr)])

# Tests across the six classes, and the class effect adjusted for expression
# level and gene length. n_samples is NA for annotation-only quantities.
kw <- function(v, label, ns = N_LNC_SAMPLES) {
  k <- kruskal.test(cls[[v]] ~ cls$class)
  data.table(quantity = label, test = "Kruskal-Wallis across the six positional classes",
             n_genes = nrow(cls), n_samples = ns, df = unname(k$parameter),
             statistic = round(unname(k$statistic), 2), p = signif(k$p.value, 3))
}
adj <- function(v, label) {
  f <- lm(cls[[v]] ~ mean_log2_expr + log10(gene_length) + class, data = cls)
  a <- anova(f)
  data.table(quantity = label,
             test = "class F-test adjusted for mean log2 expression and log10 gene length (sequential ANOVA, class last)",
             n_genes = nrow(cls), n_samples = N_LNC_SAMPLES, df = a["class", "Df"],
             statistic = round(a["class", "F value"], 2), p = signif(a["class", "Pr(>F)"], 3))
}
sp <- function(x, y, label, ns = N_LNC_SAMPLES) {
  ct <- cor.test(cls[[x]], cls[[y]], method = "spearman", exact = FALSE)
  data.table(quantity = label, test = "Spearman across all lncRNAs", n_genes = nrow(cls),
             n_samples = ns, df = NA_integer_,
             statistic = round(unname(ct$estimate), 3), p = signif(ct$p.value, 3))
}
tests <- rbind(
  kw("abs_rho", "|rho| with the non-feature fraction"),
  kw("rho_noFeature", "rho with the non-feature fraction (signed)"),
  kw("pc1_loading", "observed PC1 loading"),
  kw("gene_length", "gene length", NA_integer_), kw("n_exons_union", "exon count", NA_integer_),
  kw("n_transcripts", "transcript count", NA_integer_), kw("mean_log2_expr", "mean log2(FPKM+1)"),
  kw("frac_expressed", "fraction of samples expressed"),
  adj("abs_rho", "|rho| with the non-feature fraction"),
  adj("pc1_loading", "observed PC1 loading"),
  sp("abs_rho", "mean_log2_expr", "|rho| vs mean log2 expression (rho)"),
  sp("abs_rho", "gene_length", "|rho| vs gene length (rho)", NA_integer_),
  sp("abs_rho", "n_exons_union", "|rho| vs exon count (rho)", NA_integer_),
  sp("pc1_loading", "mean_log2_expr", "PC1 loading vs mean log2 expression (rho)"),
  sp("pc1_loading", "rho_noFeature", "PC1 loading vs rho with the non-feature fraction (rho)"),
  # The PC1 score is near-collinear with the metric at sample level, so the
  # loading is close to a rescaling of the per-gene rho.
  data.table(quantity = "observed PC1 sample score vs pct_noFeature (rho)",
             test = "Spearman across the discovery libraries in the lncRNA network (sample level, not gene level)",
             n_genes = NA_integer_, n_samples = N_LNC_SAMPLES, df = NA_integer_,
             statistic = round(unname(ct_axis_nf$estimate), 3),
             p = signif(ct_axis_nf$p.value, 3)))
save_tsv(tests, "28_class_comparison_tests.tsv"); print(tests)

# ---- B7. Class composition of every lncRNA module -----------------------------
tab <- table(cls$module, cls$class)
overall <- colSums(tab) / sum(tab)
mods <- c(setdiff(names(sort(table(cls$module), decreasing = TRUE)), "grey"), "grey")

# A simulated p cannot fall below 1/(B+1), so p values at that floor are flagged
# as upper bounds.
SIM_B <- 1e6L
at_floor <- function(p, B) isTRUE(p <= 1 / (B + 1) + 1e-12)
# Vectorised multinomial version of chisq.test()'s simulated goodness-of-fit
# test: same null and Pearson statistic, much faster at large B.
gof_sim <- function(x, prob, B) {
  x <- as.integer(x); prob <- as.numeric(prob) / sum(prob)
  n <- sum(x); E <- n * prob
  stat <- sum((x - E)^2 / E)
  null <- colSums((stats::rmultinom(B, n, prob) - E)^2 / E)
  list(statistic = stat, p.value = (1 + sum(null >= stat)) / (B + 1))
}
TEST_NOTE <- paste0("chisq_*: module vs all other lncRNAs, 2 x 6 chi-square; ",
                    "gof_*: module counts vs the class proportions of all ", nrow(cls),
                    " lncRNAs, goodness of fit; both with simulated p (B = ", format(SIM_B, scientific = FALSE),
                    "); *_p_at_floor TRUE means the p and its BH value are upper bounds of 1/(B+1), not estimates")
comp <- rbindlist(lapply(mods, function(m) {
  x <- tab[m, ]; rest <- colSums(tab) - x
  ct_tab <- rbind(x, rest); ct_tab <- ct_tab[, colSums(ct_tab) > 0, drop = FALSE]
  set.seed(SEED)
  ct <- suppressWarnings(chisq.test(ct_tab, simulate.p.value = TRUE, B = SIM_B))
  r <- data.table(module = m, n_genes = sum(x))
  for (k in CLASSES) set(r, j = paste0("n_", k), value = as.integer(x[[k]]))
  for (k in CLASSES) set(r, j = paste0("frac_", k), value = round(x[[k]] / sum(x), 3))
  for (k in CLASSES) set(r, j = paste0("enrichment_", k), value = round((x[[k]] / sum(x)) / overall[[k]], 2))
  # Two framings: module versus the other lncRNAs (2 x 6 chi-square), and module
  # counts versus whole-set class proportions (goodness of fit, as the
  # enrichment ratios describe).
  set.seed(SEED)
  gof <- gof_sim(x[CLASSES], overall[CLASSES], SIM_B)
  r[, `:=`(chisq_statistic = round(unname(ct$statistic), 2), chisq_p = signif(ct$p.value, 4),
           chisq_p_at_floor = at_floor(ct$p.value, SIM_B),
           gof_statistic = round(unname(gof$statistic), 2), gof_p = signif(gof$p.value, 4),
           gof_p_at_floor = at_floor(gof$p.value, SIM_B),
           sim_B = SIM_B, test = TEST_NOTE)]
  r
}))
comp[, chisq_fdr := signif(p.adjust(chisq_p, "BH"), 4)]
comp[, gof_fdr := signif(p.adjust(gof_p, "BH"), 4)]
all_row <- data.table(module = "all lncRNAs", n_genes = sum(tab))
for (k in CLASSES) set(all_row, j = paste0("n_", k), value = as.integer(colSums(tab)[[k]]))
for (k in CLASSES) set(all_row, j = paste0("frac_", k), value = round(overall[[k]], 3))
set.seed(SEED)
gtab <- tab[rownames(tab) != "grey", ]; gtab <- gtab[, colSums(gtab) > 0, drop = FALSE]
glob <- suppressWarnings(chisq.test(gtab, simulate.p.value = TRUE, B = SIM_B))
glob_row <- data.table(module = paste0("global (", nrow(gtab), " modules x ", ncol(gtab),
                                       " classes, grey excluded)"),
                       n_genes = sum(tab[rownames(tab) != "grey", ]),
                       chisq_statistic = round(unname(glob$statistic), 2), chisq_p = signif(glob$p.value, 4),
                       chisq_p_at_floor = at_floor(glob$p.value, SIM_B), sim_B = SIM_B,
                       test = paste0("module x class chi-square, simulated p (B = ",
                                     format(SIM_B, scientific = FALSE), ")"))
comp <- rbind(comp, all_row, glob_row, fill = TRUE)
save_tsv(comp, "28_module_composition_by_class.tsv")
print(comp[, c("module", "n_genes", paste0("frac_", CLASSES),
               "chisq_p", "chisq_p_at_floor", "chisq_fdr", "gof_p", "gof_p_at_floor", "gof_fdr"),
           with = FALSE])
msg(sum(comp$chisq_p_at_floor, na.rm = TRUE), " of the chi-square and ",
    sum(comp$gof_p_at_floor, na.rm = TRUE), " of the goodness-of-fit p-values are at the ",
    "simulation floor of ", signif(1 / (SIM_B + 1), 3), " and are upper bounds")

# Most enriched and depleted module per class. Grey (unassigned genes) is excluded.
nm <- comp[module %in% setdiff(mods, "grey")]
for (k in CLASSES) {
  hi <- nm[which.max(get(paste0("enrichment_", k)))]
  lo <- nm[which.min(get(paste0("enrichment_", k)))]
  msg("Class ", k, ": most enriched module = ", hi$module, " (",
      round(100 * hi[[paste0("frac_", k)]], 1), "%, x", hi[[paste0("enrichment_", k)]],
      "); most depleted = ", lo$module, " (", round(100 * lo[[paste0("frac_", k)]], 1),
      "%, x", lo[[paste0("enrichment_", k)]], ")")
}
for (m in intersect(c("blue", "greenyellow", "turquoise", "black", "lightcyan"), mods))
  msg("Module ", m, ": ", paste(sprintf("%s %.1f%% (x%.2f)", CLASSES,
      100 * unlist(comp[module == m, paste0("frac_", CLASSES), with = FALSE]),
      unlist(comp[module == m, paste0("enrichment_", CLASSES), with = FALSE])), collapse = ", "),
      "; chi-square p = ", comp[module == m, chisq_p],
      if (comp[module == m, chisq_p_at_floor]) " (at the simulation floor)" else "")
invisible(TRUE)
}

# ---- run both parts independently ----
ok_a <- tryCatch(part_a(), error = function(e) {
  msg("PART A FAILED: ", conditionMessage(e))
  writeLines(c("28_mutations_and_lncRNA_classes.R: Part A (mutations) failed.",
               paste0("Error: ", conditionMessage(e)), paste0("Time: ", format(Sys.time()))),
             file.path(RESULTS_DIR, "28_NOTE_partA_failed.txt"))
  FALSE
})
if (isTRUE(ok_a)) unlink(file.path(RESULTS_DIR, "28_NOTE_partA_failed.txt"))
ok_b <- tryCatch(part_b(), error = function(e) {
  msg("PART B FAILED: ", conditionMessage(e))
  writeLines(c("28_mutations_and_lncRNA_classes.R: Part B (positional classes) failed.",
               paste0("Error: ", conditionMessage(e)), paste0("Time: ", format(Sys.time()))),
             file.path(RESULTS_DIR, "28_NOTE_partB_failed.txt"))
  FALSE
})
if (isTRUE(ok_b)) unlink(file.path(RESULTS_DIR, "28_NOTE_partB_failed.txt"))

msg("Part A ", if (isTRUE(ok_a)) "completed" else "not completed", "; Part B ",
    if (isTRUE(ok_b)) "completed" else "not completed",
    "; elapsed ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
write_session_info("28_mutations_and_lncRNA_classes")
banner("28 | done")
