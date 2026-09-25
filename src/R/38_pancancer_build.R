# 38_pancancer_build.R: build all 33 TCGA project cohorts through one pipeline.
# Per project: fetch the open STAR-Counts manifest, download primary tumours
# (peripheral blood for TCGA-LAML) and solid-tissue normals with md5 checks,
# read protein-coding and lncRNA FPKM and the STAR summary rows, apply the
# discovery library-failure rule, keep the deepest library per patient and
# group as in 01, and take plate and tissue source site from the aliquot barcode.
# Builds and describes cohorts only. Requires results/38_decision_rules_lock.json.
# Outputs: cache/pan_<project>.rds, results 38_pancancer_file_manifest.tsv
#   (one row per file, retained or reason excluded) and 38_pancancer_cohorts.tsv.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(httr) })
banner("38 | Pan-cancer build: 33 TCGA projects")
set.seed(SEED)
stopifnot(file.exists(file.path(RESULTS_DIR, "38_decision_rules_lock.json")))
dir.create(PAN_DATA_DIR, showWarnings = FALSE, recursive = TRUE)

with_retry <- function(f, tries = 4, wait = 10) {
  for (i in seq_len(tries)) {
    out <- tryCatch(f(), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    msg("  retry ", i, "/", tries, ": ", conditionMessage(out)); Sys.sleep(wait * i)
  }
  stop(out)
}
gdc_hits <- function(endpoint, body) with_retry(function() {
  body$format <- "JSON"      # without it (and accept_json) the API may answer in XML
  r <- httr::POST(paste0("https://api.gdc.cancer.gov/", endpoint),
                  body = jsonlite::toJSON(body, auto_unbox = TRUE),
                  httr::content_type_json(), httr::accept_json(), httr::timeout(600))
  httr::stop_for_status(r)
  jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)$data$hits
})

# ---- 1. projects and manifests ---------------------------------------------
proj_hits <- gdc_hits("projects", list(
  filters = list(op = "=", content = list(field = "program.name", value = "TCGA")),
  fields = "project_id", size = "100"))
PROJECTS <- sort(vapply(proj_hits, function(h) h$project_id, ""))
msg(length(PROJECTS), " TCGA projects")
stopifnot(length(PROJECTS) == 33)
# PAN_ONLY=TCGA-X,TCGA-Y restricts the run, and tables then get a _partial suffix.
ONLY <- Sys.getenv("PAN_ONLY")
if (nzchar(ONLY)) PROJECTS <- intersect(PROJECTS, strsplit(ONLY, ",")[[1]])
SUFFIX <- if (nzchar(ONLY)) "_partial" else ""

manifest_of <- function(proj) {
  f <- file.path(CACHE_DIR, paste0("pan_manifest_", proj, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  hits <- gdc_hits("files", list(
    filters = list(op = "and", content = list(
      list(op = "=", content = list(field = "cases.project.project_id", value = proj)),
      list(op = "=", content = list(field = "data_type", value = "Gene Expression Quantification")),
      list(op = "=", content = list(field = "analysis.workflow_type", value = "STAR - Counts")),
      list(op = "=", content = list(field = "access", value = "open")))),
    fields = paste("file_id,file_name,md5sum,file_size,cases.submitter_id,",
                   "cases.samples.submitter_id,cases.samples.sample_type,",
                   "associated_entities.entity_submitter_id", sep = ""),
    size = "5000"))
  m <- rbindlist(lapply(hits, function(h) {
    cs <- h$cases[[1]]; sm <- cs$samples[[1]]
    data.table(project = proj, file_id = h$file_id, file_name = h$file_name,
               md5 = h$md5sum, file_size = h$file_size, patient = cs$submitter_id,
               sample_barcode = sm$submitter_id, sample_type = sm$sample_type,
               aliquot = h$associated_entities[[1]]$entity_submitter_id)
  }))
  saveRDS(m, f); m
}

# Files already downloaded by earlier stages are used in place.
local_path <- function(proj, file_id, file_name) {
  cand <- c(file.path(PROJECT_ROOT, "GDCdata", proj, "Transcriptome_Profiling",
                      "Gene_Expression_Quantification", file_id, file_name),
            file.path(ANALYSIS_DIR, "validation_data", proj, file_id, file_name),
            file.path(PAN_DATA_DIR, proj, file_id, file_name))
  hit <- cand[file.exists(winlong(cand))]
  if (length(hit)) hit[1] else cand[3]
}

download_missing <- function(m) {
  # Remove truncated files from an interrupted run (md5 mismatch).
  mine <- m[startsWith(path, PAN_DATA_DIR) & file.exists(winlong(path))]
  if (nrow(mine)) {
    bad <- mine[unname(tools::md5sum(winlong(path))) != md5]
    if (nrow(bad)) { msg("  ", nrow(bad), " partial file(s) from an earlier run removed"); unlink(winlong(bad$path)) }
  }
  todo <- m[!file.exists(winlong(path))]
  if (!nrow(todo)) return(invisible(0L))
  dir <- file.path(PAN_DATA_DIR, todo$project[1]); dir.create(dir, showWarnings = FALSE)
  ch <- split(todo$file_id, ceiling(seq_len(nrow(todo)) / 40))
  # A transfer can return HTTP 200 with a truncated archive, so each chunk is
  # md5-checked and failed files are refetched, up to four attempts.
  for (i in seq_along(ch)) {
    ids <- ch[[i]]
    for (attempt in 1:4) {
      msg("  ", todo$project[1], ": downloading chunk ", i, "/", length(ch), " (", length(ids),
          " files)", if (attempt > 1) paste0(", attempt ", attempt) else "")
      tf <- tempfile(fileext = ".tar.gz")
      with_retry(function() {
        r <- httr::POST("https://api.gdc.cancer.gov/data",
                        body = jsonlite::toJSON(list(ids = ids), auto_unbox = FALSE),
                        httr::content_type_json(), httr::write_disk(tf, overwrite = TRUE),
                        httr::timeout(3600))
        if (httr::status_code(r) != 200) stop("HTTP ", httr::status_code(r))
      })
      if (length(ids) == 1) {
        row <- todo[file_id == ids]
        dir.create(file.path(dir, row$file_id), showWarnings = FALSE)
        file.copy(tf, winlong(file.path(dir, row$file_id, row$file_name)), overwrite = TRUE)
      } else suppressWarnings(tryCatch(untar(tf, exdir = dir), error = function(e) NULL))
      unlink(tf)
      rows <- todo[file_id %in% ids]
      pth <- winlong(file.path(dir, rows$file_id, rows$file_name))
      good <- file.exists(pth)
      good[good] <- unname(tools::md5sum(pth[good])) == rows$md5[good]
      if (any(!good)) unlink(pth[!good])
      ids <- rows$file_id[!good]
      if (!length(ids)) break
      msg("    ", length(ids), " file(s) missing or failing md5; fetching again")
    }
    if (length(ids)) stop(todo$project[1], ": ", length(ids), " files could not be downloaded intact")
  }
  invisible(nrow(todo))
}

QN <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
read_project <- function(m) {
  first <- fread(winlong(m$path[1]), skip = 1, showProgress = FALSE)
  keep  <- first$gene_type %in% c("protein_coding", "lncRNA")
  ann   <- data.table(gene_id = first$gene_id[keep], gene_name = first$gene_name[keep],
                      gene_type = first$gene_type[keep])
  fp <- matrix(NA_real_, nrow(ann), nrow(m), dimnames = list(ann$gene_id, m$file_id))
  qc <- matrix(NA_real_, nrow(m), 4, dimnames = list(m$file_id, QN))
  assigned <- numeric(nrow(m))
  for (i in seq_len(nrow(m))) {
    d <- fread(winlong(m$path[i]), skip = 1, showProgress = FALSE,
               select = c("gene_id", "gene_type", "unstranded", "fpkm_unstranded"))
    v <- setNames(as.numeric(d$unstranded), d$gene_id)
    qc[i, ] <- v[QN]
    g <- d[grepl("^ENSG", gene_id)]
    assigned[i] <- sum(as.numeric(g$unstranded), na.rm = TRUE)
    g <- g[gene_type %in% c("protein_coding", "lncRNA")]
    if (!identical(g$gene_id, ann$gene_id)) g <- g[match(ann$gene_id, g$gene_id)]
    fp[, i] <- g$fpkm_unstranded
    if (i %% 200 == 0) msg("    read ", i, "/", nrow(m))
  }
  # Percent of all reads: four STAR summary rows plus reads assigned to ENSG genes.
  tot <- assigned + rowSums(qc)
  m[, `:=`(assigned_reads = assigned,
           pct_noFeature = 100 * qc[, "N_noFeature"] / tot,
           pct_multimapping = 100 * qc[, "N_multimapping"] / tot)]
  list(fpkm = fp, ann = ann, samples = m)
}

# ---- 2. per project ----------------------------------------------------------
fm_all <- list(); coh <- list()
for (proj in PROJECTS) {
  cache_f <- file.path(CACHE_DIR, paste0("pan_", proj, ".rds"))
  m <- manifest_of(proj)
  m[, group := fifelse(sample_type %in% PAN_TUMOUR_TYPES, "tumour",
                fifelse(sample_type == PAN_NORMAL_TYPE, "normal", NA_character_))]
  mu <- m[!is.na(group)]
  if (file.exists(cache_f)) {
    obj <- readRDS(cache_f); msg(proj, ": from cache")
  } else {
    mu[, path := mapply(local_path, project, file_id, file_name)]
    n_dl <- download_missing(mu)
    mu[, path := mapply(local_path, project, file_id, file_name)]
    stopifnot(all(file.exists(winlong(mu$path))))
    # md5 of every file used, downloaded or pre-existing
    mu[, md5_local := unname(tools::md5sum(winlong(path)))]
    bad <- mu[md5_local != md5]
    if (nrow(bad)) stop(proj, ": md5 mismatch for ", nrow(bad), " files, e.g. ", bad$file_id[1])
    msg(proj, ": ", nrow(mu), " files (", n_dl, " downloaded), md5 verified; reading ...")
    obj <- read_project(mu)
    s <- obj$samples
    s[, fail := assigned_reads < 10e6 | pct_noFeature > 30]
    # Aliquot barcode: tissue source site is field 2, plate is field 6.
    s[, `:=`(plate = substr(aliquot, 22, 25), tss = substr(aliquot, 6, 7))]
    s[, keep := FALSE]
    ok <- s[fail == FALSE][order(patient, group, -assigned_reads)]
    ok <- ok[!duplicated(ok[, .(patient, group)])]
    s[file_id %in% ok$file_id, keep := TRUE]
    s[, reason := fifelse(fail, "failed library",
                  fifelse(!keep, "duplicate library of a retained patient", "retained"))]
    obj$samples <- s
    obj$fpkm <- obj$fpkm[, s[keep == TRUE, file_id], drop = FALSE]
    obj$project <- proj
    saveRDS(obj, cache_f)
  }
  s <- obj$samples
  fm_all[[proj]] <- merge(m[, .(project, file_id, file_name, md5, sample_type, group)],
                          s[, .(file_id, reason)], by = "file_id", all.x = TRUE)
  coh[[proj]] <- s[, .(n_files = .N, n_failed = sum(fail), n_retained = sum(keep),
                       median_noFeature = median(pct_noFeature[keep]),
                       q25_noFeature = quantile(pct_noFeature[keep], 0.25),
                       q75_noFeature = quantile(pct_noFeature[keep], 0.75),
                       n_plates = uniqueN(plate[keep]), n_tss = uniqueN(tss[keep])),
                   by = group][, project := proj]
  rm(obj); gc(verbose = FALSE)
}

fm <- rbindlist(fm_all)
fm[is.na(reason), reason := "not a primary tumour or solid-tissue normal"]
save_tsv(fm[order(project, sample_type, file_id)], paste0("38_pancancer_file_manifest", SUFFIX, ".tsv"))
co <- rbindlist(coh)
co[, eligible := fifelse(group == "tumour",
                         n_retained >= PAN_MIN_TUMOURS & !project %in% PAN_EXCLUDE_INFERENCE,
                         n_retained >= PAN_MIN_NORMALS)]
setcolorder(co, c("project", "group"))
save_tsv(co[order(project, -rank(group))], paste0("38_pancancer_cohorts", SUFFIX, ".tsv"))
msg("eligible tumour projects: ", co[group == "tumour" & eligible == TRUE, .N],
    "; projects with eligible normals: ", co[group == "normal" & eligible == TRUE, .N])
write_session_info("38_pancancer_build")
msg("38 done")
