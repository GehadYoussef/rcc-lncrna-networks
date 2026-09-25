# 11_subtype_specificity.R: prognostic modules in papillary and chromophobe RCC.
#
# Downloads TCGA-KIRP and TCGA-KICH once (cached as subtype_<tag>.rds). Every
# cohort, TCGA-KIRC included, is scored the same way: log2(FPKM + 1) on the
# network transcripts, residualised on its own STAR metrics, and projected on
# the locked discovery loadings with cohort standardisation.
# Outputs: module effects under one covariate set (age, sex, T, M1, non-feature
# fraction), Cochran Q across cohorts that pass the events-per-variable rule,
# the leading axes in KIRP and KICH, lncRNA coupling to the catabolic module
# by cohort, and per-cohort STAR quality summaries. Run after 03, 08 and 09.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(httr)
  library(survival); library(ggplot2)
})
banner("11 | Subtype specificity: KIRP and KICH")
set.seed(SEED)

SUBTYPES <- c("TCGA-KIRP", "TCGA-KICH")
ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
surv <- readRDS(file.path(CACHE_DIR, "survival.rds"))
L    <- readRDS(file.path(CACHE_DIR, "locked_model.rds"))
stopifnot(identical(L$version, "v9"))
LOAD <- discovery_loadings(nets)

gdc <- function(endpoint, body) {
  r <- httr::POST(paste0("https://api.gdc.cancer.gov/", endpoint),
                  body = jsonlite::toJSON(body, auto_unbox = TRUE),
                  httr::content_type_json(), httr::timeout(300))
  httr::stop_for_status(r)
  jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"),
                     simplifyDataFrame = TRUE)$data$hits
}

# First non-missing string / largest numeric value across the nested records
# (diagnoses, follow_ups) the cases endpoint returns per patient.
dx1 <- function(l, c_) {
  if (is.null(l)) return(character(0))
  vapply(l, function(d) {
    if (is.null(d) || !is.data.frame(d) || !c_ %in% names(d) || nrow(d) == 0)
      return(NA_character_)
    v <- as.character(d[[c_]])
    v <- v[!is.na(v) & nzchar(v) & !tolower(v) %in% c("not reported", "unknown")]
    if (!length(v)) NA_character_ else v[1] }, character(1))
}
dxm <- function(l, c_) {
  if (is.null(l)) return(numeric(0))
  vapply(l, function(d) {
    if (is.null(d) || !is.data.frame(d) || !c_ %in% names(d) || nrow(d) == 0)
      return(NA_real_)
    v <- suppressWarnings(as.numeric(d[[c_]]))
    if (!any(is.finite(v))) NA_real_ else max(v, na.rm = TRUE) }, numeric(1))
}

parse_T <- function(x) {
  x <- toupper(trimws(as.character(x)))
  o <- rep(NA_integer_, length(x))
  o[grepl("^T1", x)] <- 1L; o[grepl("^T2", x)] <- 2L
  o[grepl("^T3", x)] <- 3L; o[grepl("^T4", x)] <- 4L
  o
}
parse_G <- function(x) {
  x <- toupper(trimws(as.character(x)))
  o <- rep(NA_integer_, length(x))
  o[x == "G1"] <- 1L; o[x == "G2"] <- 2L
  o[x == "G3"] <- 3L; o[x == "G4"] <- 4L
  o
}
# As in 01: sub-stages are collapsed to the Roman numeral.
parse_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "NOT REPORTED", "[NOT AVAILABLE]", "[UNKNOWN]")] <- NA
  rn <- sub("^STAGE\\s*", "", x)
  rn <- sub("[A-C]$", "", rn)
  out <- rep(NA_integer_, length(x))
  out[rn == "I"] <- 1L; out[rn == "II"] <- 2L
  out[rn == "III"] <- 3L; out[rn == "IV"] <- 4L
  out
}

# ---- build one subtype cohort (download, read, harmonise as in 01) ----
build_cohort <- function(proj) {
  tag <- sub("TCGA-", "", proj)
  dir <- file.path(DOWNLOAD_DIR, "validation_data", proj)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  cache_f <- file.path(CACHE_DIR, paste0("subtype_", tag, ".rds"))
  if (file.exists(cache_f)) { msg(proj, ": from cache"); return(readRDS(cache_f)) }

  filt <- list(op = "and", content = list(
    list(op = "in", content = list(field = "cases.project.project_id",
                                   value = list(proj))),
    list(op = "in", content = list(field = "data_category",
                                   value = list("Transcriptome Profiling"))),
    list(op = "in", content = list(field = "data_type",
                                   value = list("Gene Expression Quantification"))),
    list(op = "in", content = list(field = "analysis.workflow_type",
                                   value = list("STAR - Counts"))),
    list(op = "in", content = list(field = "access", value = list("open")))))
  hits <- gdc("files", list(filters = filt, format = "JSON", size = "5000",
                            fields = paste("file_id", "file_name",
                                           "cases.submitter_id",
                                           "cases.samples.submitter_id",
                                           "cases.samples.sample_type", sep = ",")))
  fm <- rbindlist(lapply(seq_len(nrow(hits)), function(i) {
    cs <- hits$cases[[i]]; sm <- cs$samples[[1]]
    data.table(file_id = hits$file_id[i], file_name = hits$file_name[i],
               patient = cs$submitter_id[1], sample_barcode = sm$submitter_id[1],
               sample_type = sm$sample_type[1])
  }))
  fm <- fm[sample_type == "Primary Tumor"]
  fm[, path := winlong(file.path(dir, file_id, file_name))]
  msg(proj, ": ", nrow(fm), " primary-tumour files")

  todo <- fm[!file.exists(path)]
  if (nrow(todo)) {
    ch <- split(todo$file_id, ceiling(seq_len(nrow(todo)) / 40))
    if (length(ch) > 1 && length(ch[[length(ch)]]) == 1) {
      ch[[length(ch) - 1]] <- c(ch[[length(ch) - 1]], ch[[length(ch)]])
      ch[[length(ch)]] <- NULL
    }
    for (i in seq_along(ch)) {
      msg("  downloading chunk ", i, "/", length(ch), " (", length(ch[[i]]), ")")
      tf <- tempfile(fileext = ".tar.gz"); ok <- FALSE
      for (a in 1:3) {
        r <- try(httr::POST("https://api.gdc.cancer.gov/data",
                            body = jsonlite::toJSON(list(ids = ch[[i]]),
                                                    auto_unbox = TRUE),
                            httr::content_type_json(),
                            httr::write_disk(tf, overwrite = TRUE),
                            httr::timeout(1800)), silent = TRUE)
        if (!inherits(r, "try-error") && httr::status_code(r) == 200) { ok <- TRUE; break }
      }
      if (!ok) stop("download failed for ", proj)
      if (length(ch[[i]]) == 1) {
        row <- fm[file_id == ch[[i]]]
        dir.create(file.path(dir, row$file_id), showWarnings = FALSE)
        file.copy(tf, winlong(file.path(dir, row$file_id, row$file_name)),
                  overwrite = TRUE)
      } else untar(tf, exdir = dir)
      unlink(tf)
    }
  }
  fm <- fm[file.exists(path)]

  msg(proj, ": reading ", nrow(fm), " count files ...")
  first <- fread(fm$path[1], skip = 1, showProgress = FALSE)
  kr  <- grepl("^ENSG", first$gene_id)
  ann <- data.table(gene_id = first$gene_id[kr], gene_name = first$gene_name[kr],
                    gene_type = first$gene_type[kr])
  fp  <- matrix(NA_real_, nrow(ann), nrow(fm),
                dimnames = list(ann$gene_id, fm$sample_barcode))
  qn  <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
  qc  <- matrix(NA_real_, nrow(fm), 4, dimnames = list(fm$sample_barcode, qn))
  lib <- numeric(nrow(fm))
  for (i in seq_len(nrow(fm))) {
    d <- fread(fm$path[i], skip = 1, showProgress = FALSE,
               select = c("gene_id", "unstranded", "fpkm_unstranded"))
    v <- setNames(as.numeric(d$unstranded), d$gene_id)
    qc[i, ] <- v[qn]
    d <- d[grepl("^ENSG", gene_id)]
    if (!identical(d$gene_id, ann$gene_id)) d <- d[match(ann$gene_id, d$gene_id)]
    fp[, i] <- d$fpkm_unstranded; lib[i] <- sum(d$unstranded, na.rm = TRUE)
    if (i %% 50 == 0) msg("    ", i, "/", nrow(fm))
  }
  tot <- lib + rowSums(qc, na.rm = TRUE)
  qcd <- data.table(sample_barcode = fm$sample_barcode,
                    pct_multimapping = 100 * qc[, "N_multimapping"] / tot,
                    pct_noFeature = 100 * qc[, "N_noFeature"] / tot,
                    assigned_reads = lib)

  h <- gdc("cases", list(
    filters = list(op = "in", content = list(field = "project.project_id",
                                             value = list(proj))),
    fields = "submitter_id", expand = "demographic,diagnoses,follow_ups",
    format = "JSON", size = "3000"))
  gt <- function(df, c_) if (!is.null(df) && c_ %in% names(df)) df[[c_]]
                         else rep(NA, nrow(h))
  g1 <- as.character(gt(h$demographic, "gender"))
  g2 <- as.character(gt(h$demographic, "sex_at_birth"))
  cl <- data.table(
    patient = h$submitter_id,
    sex_chr = ifelse(!is.na(g1) & nzchar(g1), g1, g2),
    vital   = as.character(gt(h$demographic, "vital_status")),
    dtd     = suppressWarnings(as.numeric(gt(h$demographic, "days_to_death"))),
    dtlf1   = dxm(h$diagnoses, "days_to_last_follow_up"),
    dtlf2   = dxm(h$follow_ups, "days_to_follow_up"),
    age_days= dxm(h$diagnoses, "age_at_diagnosis"),
    pT      = dx1(h$diagnoses, "ajcc_pathologic_t"),
    pN      = dx1(h$diagnoses, "ajcc_pathologic_n"),
    pM      = dx1(h$diagnoses, "ajcc_pathologic_m"),
    grade_raw = dx1(h$diagnoses, "tumor_grade"))
  cl[, `:=`(age = age_days / 365.25,
            sex = factor(tolower(sex_chr), levels = c("female", "male")),
            os_event = as.integer(vital == "Dead"),
            T_stage = parse_T(pT), grade_num = parse_G(grade_raw))]
  cl[, N_pos := ifelse(grepl("^N[1-9]", toupper(trimws(pN))), 1L, 0L)]
  cl[, M1    := ifelse(grepl("^M1", toupper(trimws(pM))), 1L, 0L)]
  cl[, os_time := ifelse(!is.na(dtd), dtd, pmax(dtlf1, dtlf2, na.rm = TRUE))]
  cl[!is.na(os_time) & os_time > OS_CENSOR_DAYS, os_event := 0L]
  cl[!is.na(os_time), os_time := pmin(os_time, OS_CENSOR_DAYS)]

  s <- data.table(patient = fm$patient, sample_barcode = fm$sample_barcode,
                  libsize = lib)
  setorder(s, patient, -libsize); s <- s[!duplicated(patient)]
  co <- merge(s, cl, by = "patient")
  co <- co[!is.na(os_time) & os_time > 0 & !is.na(os_event) & !is.na(age) &
           !is.na(sex)]
  out <- list(project = proj, fpkm = fp[, co$sample_barcode, drop = FALSE],
              gene_ann = ann, cohort = co,
              qc = qcd[match(co$sample_barcode, sample_barcode)])
  saveRDS(out, cache_f)
  out
}

# ---- AJCC stage and M reconciliation (as in 01) ----
# Stage is fetched separately (cached) so that MX at stage IV without T4 can be
# recoded M1. The reconciled flag is kept for the sensitivity fit.
fetch_stage <- function(obj) {
  tag <- sub("TCGA-", "", obj$project)
  f <- file.path(CACHE_DIR, paste0("subtype_stage_", tag, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  h <- tryCatch(gdc("cases", list(
         filters = list(op = "in", content = list(field = "project.project_id",
                                                  value = list(obj$project))),
         fields = "submitter_id,diagnoses.ajcc_pathologic_stage",
         format = "JSON", size = "3000")), error = function(e) NULL)
  if (is.null(h) || !length(h$submitter_id)) {
    msg(obj$project, ": GDC cases endpoint unreachable -- stage left missing, ",
        "no M reconciliation possible (result not cached; rerun to retry)")
    return(data.table(patient = obj$cohort$patient, stage_raw = NA_character_))
  }
  stage <- dx1(h$diagnoses, "ajcc_pathologic_stage")
  if (!length(stage)) stage <- rep(NA_character_, length(h$submitter_id))
  st <- data.table(patient = h$submitter_id, stage_raw = stage)
  saveRDS(st, f); st
}

add_stage <- function(obj) {
  st <- fetch_stage(obj)
  co <- copy(obj$cohort)
  co[, stage_raw := st$stage_raw[match(patient, st$patient)]]
  co[, stage_num := parse_stage(stage_raw)]
  co[, M1_as_coded := as.integer(M1)]
  co[, M1_imputed  := !grepl("^M[01]", toupper(trimws(pM)))]
  co[, M1_stage_reconciled := as.integer(M1_imputed & stage_num %in% 4L &
                                         T_stage %in% 1:3 & M1_as_coded == 0L)]
  co[M1_stage_reconciled == 1L, M1 := 1L]
  msg(obj$project, ": stage available for ", sum(!is.na(co$stage_num)), "/",
      nrow(co), ", stage IV n = ", sum(co$stage_num %in% 4L),
      ", M1 as coded n = ", sum(co$M1_as_coded == 1L),
      ", M1 reconciled from stage IV n = ", sum(co$M1_stage_reconciled == 1L),
      ", T missing n = ", sum(is.na(co$T_stage)))
  obj$cohort <- co
  obj
}

# ---- score a subtype cohort ----
# The observed matrix is kept for the axis analysis. The residualised matrix
# is projected on the discovery loadings, as CPTAC-3 is in 09.
prepare_subtype <- function(obj) {
  g_net <- unique(c(colnames(nets$mrna$expr), colnames(nets$lnc$expr)))
  g     <- intersect(rownames(obj$fpkm), g_net)
  E_obs <- t(log2(obj$fpkm[g, , drop = FALSE] + 1))
  q     <- obj$qc[match(rownames(E_obs), sample_barcode)]
  E     <- remove_technical(E_obs, tech_covariates(q))
  M     <- score_modules(E, LOAD)
  msg(obj$project, ": ", length(g), "/", length(g_net), " network transcripts present; ",
      ncol(M), "/", length(LOAD), " modules scored on ", nrow(M), " samples")
  list(E_obs = E_obs, q = q, M = M)
}

# ---- test the modules significant in ccRCC ----
target <- rbind(
  as.data.table(surv$mrna)[fdr_full < FDR_ALPHA,
    .(feature = paste0("mRNA_ME", module), biotype = "mRNA", module,
      HR_kirc = HR_full)],
  as.data.table(surv$lnc)[fdr_full < FDR_ALPHA,
    .(feature = paste0("lnc_ME", module), biotype = "lncRNA", module,
      HR_kirc = HR_full)])
msg("Testing ", nrow(target), " ccRCC-significant modules in other subtypes")

# Few events cause separation, so every cohort uses one small covariate set and
# cohorts below MIN_EPV events per variable are descriptive only.

fit_me <- function(y, ME, X) {
  d <- data.frame(ME = ME, X)
  fit <- tryCatch(coxph(y ~ ., data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  list(HR = s$conf.int["ME", "exp(coef)"],
       lo = s$conf.int["ME", "lower .95"], hi = s$conf.int["ME", "upper .95"],
       p  = s$coefficients["ME", "Pr(>|z|)"],
       logHR = s$coefficients["ME", "coef"], se = s$coefficients["ME", "se(coef)"],
       n = s$n, events = s$nevent)
}

test_cohort <- function(obj, M) {
  co <- obj$cohort[match(rownames(M), sample_barcode)]
  q  <- obj$qc[match(rownames(M), sample_barcode)]
  # Fixed across cohorts. Grade is omitted because it is not standard for
  # chromophobe tumours. Missing T drops the patient.
  X0 <- data.frame(age = as.numeric(co$age), male = as.numeric(co$sex == "male"),
                   T_stage = as.numeric(co$T_stage), M1 = as.numeric(co$M1),
                   noFeat = as.numeric(q$pct_noFeature))
  cc <- complete.cases(X0) & is.finite(co$os_time) & !is.na(co$os_event) &
        !is.na(co$M1_as_coded)
  n_T_missing <- sum(is.na(co$T_stage))
  X <- X0[cc, , drop = FALSE]
  X$age <- as.numeric(scale(X$age)); X$noFeat <- as.numeric(scale(X$noFeat))
  keep <- vapply(X, function(z) length(unique(z)) > 1, logical(1))
  X <- X[, keep, drop = FALSE]
  # Sensitivity: M1 as coded (MX -> 0), ignoring the stage-IV reconciliation.
  Xc <- X
  if ("M1" %in% names(Xc)) {
    Xc$M1 <- as.numeric(co$M1_as_coded[cc])
    if (length(unique(Xc$M1)) < 2) Xc$M1 <- NULL
  }
  y <- Surv(co$os_time[cc], co$os_event[cc])
  n_ev <- sum(co$os_event[cc])
  epv  <- n_ev / (ncol(X) + 1)
  adequate <- epv >= MIN_EPV
  msg(obj$project, ": ", nrow(M), " scored, ", sum(cc), " complete cases (T missing ",
      n_T_missing, "), ", n_ev, " events, ", ncol(X) + 1, " parameters, EPV = ",
      round(epv, 1), if (adequate) "  [adequate]" else
      "  [UNDERPOWERED -- descriptive only, excluded from inference]")
  summ <- data.table(cohort = obj$project, n_scored = nrow(M), n_complete = sum(cc),
                     events = n_ev, n_parameters = ncol(X) + 1L,
                     epv = round(epv, 1), adequate = adequate,
                     n_T_missing = n_T_missing,
                     n_stage_available = sum(!is.na(co$stage_num)),
                     n_stage_iv = sum(co$stage_num %in% 4L),
                     n_M1_as_coded = sum(co$M1_as_coded == 1L, na.rm = TRUE),
                     n_M1_reconciled = sum(co$M1_stage_reconciled == 1L, na.rm = TRUE))
  eff <- rbindlist(lapply(target$feature, function(f) {
    if (!f %in% colnames(M)) return(NULL)
    # HR per SD of the score within the analysed (complete-case) patients.
    me <- as.numeric(scale(M[cc, f]))
    a <- fit_me(y, me, X); if (is.null(a)) return(NULL)
    b <- fit_me(y, me, Xc)
    data.table(cohort = obj$project, feature = f, adequate = adequate,
               epv = round(epv, 1), n = a$n, events = a$events,
               n_scored = nrow(M), n_T_missing = n_T_missing,
               HR = a$HR, lo = a$lo, hi = a$hi, p = a$p, logHR = a$logHR, se = a$se,
               HR_M1coded = if (is.null(b)) NA_real_ else b$HR,
               lo_M1coded = if (is.null(b)) NA_real_ else b$lo,
               hi_M1coded = if (is.null(b)) NA_real_ else b$hi,
               p_M1coded  = if (is.null(b)) NA_real_ else b$p)
  }))
  list(effects = eff, summary = summ)
}

# TCGA-KIRC through the same path: residualised matrices, samples in both
# networks, same loadings and covariate set.
kirc_cl <- as.data.table(ds$cohort_full)
for (v in c("M1_as_coded", "M1_stage_reconciled", "stage_num")) {
  if (!v %in% names(kirc_cl)) {
    msg("dataset.rds lacks ", v, " (written by an older 01_build_data.R); treating as ",
        if (v == "M1_as_coded") "M1" else "missing")
    if (v == "M1_as_coded") kirc_cl[, M1_as_coded := as.integer(M1)]
    else if (v == "M1_stage_reconciled") kirc_cl[, M1_stage_reconciled := NA_integer_]
    else kirc_cl[, stage_num := NA_integer_]
  }
}
kirc_obj <- list(project = "TCGA-KIRC", cohort = kirc_cl,
                 qc = kirc_cl[, .(sample_barcode, pct_noFeature, pct_multimapping,
                                  assigned_reads = libsize)])
.ks <- Reduce(intersect, list(rownames(nets$mrna$expr), rownames(nets$lnc$expr),
                              kirc_cl$sample_barcode))
kirc_M <- cbind(score_modules(nets$mrna$expr[.ks, , drop = FALSE],
                              LOAD[grep("^mRNA_", names(LOAD))]),
                score_modules(nets$lnc$expr[.ks, , drop = FALSE],
                              LOAD[grep("^lnc_",  names(LOAD))]))
msg("TCGA-KIRC: ", nrow(kirc_M), " samples in both networks, ", ncol(kirc_M),
    " module scores")
kirc_res <- test_cohort(kirc_obj, kirc_M)
all_eff  <- list("TCGA-KIRC" = kirc_res$effects)
all_summ <- list("TCGA-KIRC" = kirc_res$summary)
scores   <- list("TCGA-KIRC" = kirc_M)
prepared <- list()

for (proj in SUBTYPES) {
  obj <- add_stage(build_cohort(proj))
  msg(proj, ": n = ", nrow(obj$cohort), ", events = ", sum(obj$cohort$os_event))
  pr  <- prepare_subtype(obj)
  prepared[[proj]] <- pr
  scores[[proj]]   <- pr$M
  r <- test_cohort(obj, pr$M)
  all_eff[[proj]]  <- r$effects
  all_summ[[proj]] <- r$summary
}
save_tsv(rbindlist(all_summ), "11_subtype_cohort_summary.tsv")

res <- rbindlist(all_eff)
res <- merge(res, target[, .(feature, biotype, module)], by = "feature")

# Flag separation or non-convergence.
res[, implausible := !is.finite(HR) | HR > 20 | HR < 0.05 | se > 3]
if (any(res$implausible))
  msg("Flagged ", sum(res$implausible),
      " implausible fits (separation / non-convergence); excluded from inference")
setorder(res, biotype, module, cohort)
save_tsv(res[, .(biotype, module, cohort, n, events, epv, adequate,
                 implausible, HR = round(HR, 3), lo = round(lo, 3),
                 hi = round(hi, 3), p = signif(p, 3),
                 n_scored, n_T_missing,
                 HR_M1coded = round(HR_M1coded, 3), lo_M1coded = round(lo_M1coded, 3),
                 hi_M1coded = round(hi_M1coded, 3), p_M1coded = signif(p_M1coded, 3))],
         "11_subtype_module_effects.tsv")
print(res[, .(module, cohort, n, events, epv, adequate, HR = round(HR, 3),
              p = signif(p, 3))])

# ---- formal heterogeneity across subtypes (Cochran Q) -----------------------
usable <- res[adequate == TRUE & implausible == FALSE]
msg("Cohorts contributing to heterogeneity tests: ",
    paste(sort(unique(usable$cohort)), collapse = ", "))
het <- usable[is.finite(se) & se > 0,
  {
    w <- 1 / se^2
    mu <- sum(w * logHR) / sum(w)
    Q  <- sum(w * (logHR - mu)^2)
    df <- .N - 1
    .(k = .N, cohorts = paste(sort(cohort), collapse = ";"),
      pooled_HR = round(exp(mu), 3), Q = round(Q, 2), df = df,
      p_heterogeneity = signif(pchisq(Q, df, lower.tail = FALSE), 3),
      I2 = round(max(0, (Q - df) / Q) * 100))
  }, by = .(biotype, module)]
# Add the per-cohort estimates that entered Q.
est_of <- function(coh, suffix) {
  usable[cohort == coh, setNames(.(biotype, module, round(HR, 3), round(lo, 3),
                                   round(hi, 3), signif(p, 3)),
                                 c("biotype", "module", paste0(c("HR", "lo", "hi", "p"),
                                                              "_", suffix)))]
}
het <- merge(het, est_of("TCGA-KIRC", "kirc"), by = c("biotype", "module"), all.x = TRUE)
het <- merge(het, est_of("TCGA-KIRP", "kirp"), by = c("biotype", "module"), all.x = TRUE)
save_tsv(het, "11_subtype_heterogeneity.tsv")
print(het)
msg("Modules with significant between-subtype heterogeneity (p < 0.05): ",
    paste(het[p_heterogeneity < 0.05, module], collapse = ", "))

g <- ggplot(res[implausible == FALSE], aes(HR, module, colour = cohort)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y",
                position = position_dodge(0.6)) +
  geom_point(size = 2, position = position_dodge(0.6)) +
  scale_x_log10() + facet_grid(biotype ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("TCGA-KIRC" = "#B2182B",
                                 "TCGA-KIRP" = "#2166AC",
                                 "TCGA-KICH" = "#1B7837")) +
  labs(title = "Is the prognostic axis clear-cell specific?",
       subtitle = "Module membership and loadings locked in TCGA-KIRC; age, sex, T, M1, non-feature fraction",
       x = "Hazard ratio per 1 SD (95% CI)", y = NULL, colour = NULL) +
  theme_bw() + theme(legend.position = "bottom")
save_fig(g, "11_subtype_specificity", 8, 5.5)

# ---- leading axis of the observed matrices in KIRP and KICH ----
# PC1 as in 07, with its variance share and Spearman correlations with the STAR
# metrics and mean expression. The protein-coding matrix is the control.
banner("Axis replication in KIRP and KICH")
axis_stats <- function(E, q, cohort, mat) {
  X  <- scale(E, center = TRUE, scale = FALSE)
  sv <- svd(X, nu = 0, nv = 1)
  vs <- sv$d^2 / sum(sv$d^2)
  pc1 <- as.numeric(X %*% sv$v[, 1]); me <- rowMeans(E)
  if (cor(pc1, me) < 0) pc1 <- -pc1
  sp <- function(z) {
    ct <- suppressWarnings(cor.test(pc1, z, method = "spearman"))
    c(unname(ct$estimate), ct$p.value)
  }
  a <- sp(q$pct_noFeature); b <- sp(q$pct_multimapping)
  d <- sp(log10(q$assigned_reads)); e <- sp(me)
  data.table(cohort = cohort, matrix = mat, n = nrow(E), n_genes = ncol(E),
             pc1_var_share = round(vs[1], 4), pc2_var_share = round(vs[2], 4),
             pc3_var_share = round(vs[3], 4), pc4_var_share = round(vs[4], 4),
             pc5_var_share = round(vs[5], 4),
             rho_noFeature = round(a[1], 3), p_noFeature = signif(a[2], 3),
             rho_multimap  = round(b[1], 3), p_multimap  = signif(b[2], 3),
             rho_log_depth = round(d[1], 3), p_log_depth = signif(d[2], 3),
             rho_mean_expr = round(e[1], 3), p_mean_expr = signif(e[2], 3))
}
axis_tbl <- rbindlist(lapply(SUBTYPES, function(proj) {
  pr <- prepared[[proj]]
  gl <- intersect(colnames(pr$E_obs), colnames(nets$lnc$expr))
  gm <- intersect(colnames(pr$E_obs), colnames(nets$mrna$expr))
  rbind(axis_stats(pr$E_obs[, gl, drop = FALSE], pr$q, proj, "lncRNA"),
        axis_stats(pr$E_obs[, gm, drop = FALSE], pr$q, proj, "protein_coding"))
}))
save_tsv(axis_tbl, "11_axis_replication_subtypes.tsv")
print(axis_tbl[, .(cohort, matrix, n, n_genes, pc1_var_share, rho_noFeature,
                   rho_multimap, rho_log_depth, rho_mean_expr)])

# ---- coupling of prognostic lncRNA scores with the catabolic module ----
# The coupling that stage 19 measures in KIRC, repeated in every cohort.
# CPTAC-3 scores are the locked model's Ev.
banner("lncRNA coupling to the catabolic protein-coding module, by cohort")
scores[["CPTAC-3"]] <- L$Ev
sig_lnc <- as.data.table(surv$lnc)[fdr_full < FDR_ALPHA, module]
# The catabolic module is found by its enrichment, since colours can change
# when the network is rebuilt.
ref_mod <- catabolic_reference_module(); ref_col <- paste0("mRNA_ME", ref_mod)
coupling <- rbindlist(lapply(names(scores), function(coh) {
  M <- scores[[coh]]
  rbindlist(lapply(sig_lnc, function(m) {
    a <- paste0("lnc_ME", m)
    if (!all(c(a, ref_col) %in% colnames(M))) return(NULL)
    ct <- cor.test(M[, a], M[, ref_col], method = "pearson")
    cs <- suppressWarnings(cor(M[, a], M[, ref_col], method = "spearman"))
    data.table(cohort = coh, module = m, reference_module = ref_mod, n = nrow(M),
               pearson_r = round(unname(ct$estimate), 3),
               CI_lo = round(ct$conf.int[1], 3), CI_hi = round(ct$conf.int[2], 3),
               p = signif(ct$p.value, 3), spearman = round(cs, 3))
  }))
}))
coupling[, cohort := factor(cohort, levels = c("TCGA-KIRC", "CPTAC-3", SUBTYPES))]
setorder(coupling, module, cohort)
coupling[, cohort := as.character(cohort)]
save_tsv(coupling, "11_lnc_black_coupling_by_cohort.tsv")
print(coupling)

# ---- per-cohort STAR quality summaries ----
# Library protocols differ (polyA in TCGA, ribo-depleted total RNA in CPTAC-3),
# so metrics are summarised per cohort.
banner("Quality metrics by cohort")
vd <- readRDS(file.path(CACHE_DIR, "validation_dataset.rds"))
vq <- vd$qc
if (is.null(vq)) {
  msg("validation_dataset.rds carries no qc table (written by an older 08_validation_data.R); ",
      "reading valid_star_qc.rds")
  vq <- as.data.table(readRDS(file.path(CACHE_DIR, "valid_star_qc.rds")))
  vq <- vq[sample_barcode %in% vd$cohort$sample_barcode]
}
qsum <- function(q, cohort, protocol) {
  qi <- function(z, nm) {
    v <- quantile(as.numeric(z), c(0.5, 0.25, 0.75), na.rm = TRUE)
    setNames(as.list(round(v, 3)), paste0(nm, c("_median", "_q25", "_q75")))
  }
  as.data.table(c(list(cohort = cohort, protocol = protocol, n = nrow(q)),
                  qi(q$pct_noFeature, "noFeature"), qi(q$pct_multimapping, "multimap"),
                  qi(q$assigned_reads, "assigned_reads")))
}
qual <- rbind(
  qsum(kirc_obj$qc, "TCGA-KIRC", "polyA"),
  qsum(as.data.table(vq), "CPTAC-3", "ribo-depleted total RNA"),
  rbindlist(lapply(SUBTYPES, function(proj) qsum(prepared[[proj]]$q, proj, "polyA"))))
save_tsv(qual, "11_cohort_quality_metrics.tsv")
print(qual)

write_session_info("11_subtype_specificity")
banner("11 | done")
