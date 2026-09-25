# 01_validate_manifest.R: manifest validation, input acquisition and checksums.
# Validates the dataset manifest and sample sheet in data/manifests, lists every
# input file (GEO matrices, author annotations, TISCH2 cell tables, extra
# cohorts), downloads missing files into data/raw/singlecell from public URLs
# (GEO files via the NCBI HTTPS tree) and records SHA-256 for each local input.
# Outputs: results/singlecell/01_acquisition_log.tsv,
#          01_manifest_resolved.tsv, 01_manifest_validation.tsv
# Flag:    --no-download   audit local files only

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
banner("01  manifest validation and acquisition")
start_log("01_validate_manifest")

NO_DOWNLOAD <- "--no-download" %in% commandArgs(trailingOnly = TRUE)
m <- read_manifest(); man <- m$manifest; smp <- m$samples

problems <- validate_manifest(man, smp)
save_tsv(data.table(check = c("manifest_schema", "n_datasets", "n_samples", "n_samples_included"),
                    result = c(if (length(problems)) paste(problems, collapse = " | ") else "ok",
                               nrow(man), nrow(smp), sum(smp$include))),
         "01_manifest_validation.tsv")
if (length(problems)) stop("manifest validation failed:\n  ", paste(problems, collapse = "\n  "))
msg("manifest valid: ", nrow(man), " datasets, ", nrow(smp), " sample rows")

# ---- build the file list -----------------------------------------------------
files <- list()
geo_rows <- smp[format %in% c("10x_h5", "10x_mtx_v2", "10x_mtx_v3", "dense_txt")]
for (i in seq_len(nrow(geo_rows))) {
  r <- geo_rows[i]
  for (fn in sample_files(r$file_prefix, r$format))
    files[[length(files) + 1]] <- data.table(
      dataset_name = r$dataset_name, sample_id = r$sample_id, role = "expression_raw_counts",
      filename = fn, url = geo_sample_url(r$gsm, fn),
      dest = file.path(RAW_DIR, r$dataset_name, "geo", fn), fetch = TRUE)
}
tisch <- man[tisch_id != ""]
for (i in seq_len(nrow(tisch))) {
  ds <- tisch$dataset_name[i]
  base <- sprintf("https://tisch.compbio.cn/static/data/%s/%s", ds, ds)
  files[[length(files) + 1]] <- data.table(
    dataset_name = ds, sample_id = "", role = "tisch2_cell_metainfo",
    filename = paste0(ds, "_CellMetainfo_table.tsv"), url = paste0(base, "_CellMetainfo_table.tsv"),
    dest = file.path(RAW_DIR, ds, "tisch2", paste0(ds, "_CellMetainfo_table.tsv")), fetch = TRUE)
  files[[length(files) + 1]] <- data.table(
    dataset_name = ds, sample_id = "", role = "tisch2_expression_lognorm",
    filename = paste0(ds, "_Expression.zip"), url = paste0(base, "_Expression.zip"),
    dest = file.path(RAW_DIR, ds, "tisch2", paste0(ds, "_Expression.zip")),
    fetch = isTRUE(CFG$download$fetch_tisch2_expression))
}
extra_f <- EXTRA_FILES
if (file.exists(extra_f)) {
  ex <- fread(extra_f, colClasses = "character")
  files[[length(files) + 1]] <- ex[, .(dataset_name, sample_id = "", role, filename, url,
                                       dest = file.path(RAW_DIR, dataset_name, "author", filename),
                                       fetch = toupper(fetch_required) == "TRUE")]
}
files <- rbindlist(files)
files[, dest := normalizePath(dest, winslash = "/", mustWork = FALSE)]
msg(nrow(files), " input files registered (", sum(files$fetch), " required)")

# ---- download ------------------------------------------------------------------
files[, status := NA_character_]
files[, retrieved_utc := NA_character_]
for (i in seq_len(nrow(files))) {
  f <- files[i]
  if (!f$fetch) { files[i, status := "deferred"]; next }
  if (NO_DOWNLOAD) {
    files[i, status := if (file.exists(f$dest)) "present" else "missing"]
    next
  }
  msg("[", i, "/", nrow(files), "] ", f$dataset_name, " ", f$filename)
  res <- download_resumable(f$url, f$dest, ua = CFG$download$user_agent,
                            retries = CFG$download$retries,
                            wait = CFG$download$retry_wait_seconds,
                            timeout = CFG$download$timeout_seconds)
  files[i, status := res$status]
  if (res$status == "downloaded") files[i, retrieved_utc := format(Sys.time(), tz = "UTC", usetz = TRUE)]
  if (res$status == "failed") msg("WARNING: download failed for ", f$url)
}

# ---- archive extraction ----------------------------------------------------------
# Bundled cohorts are unpacked into <dataset>/extracted/ unless a completion
# marker exists. Checksums refer to the untouched archives.
extract_bundle <- function(f) {
  out <- file.path(dirname(dirname(f)), "extracted")
  marker <- file.path(out, paste0(".extracted_", basename(f)))
  if (file.exists(marker)) return(invisible("present"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  msg("extracting ", basename(f))
  if (grepl("\\.tar(\\.gz)?$", f)) utils::untar(f, exdir = out)
  for (z in list.files(out, pattern = "[.]zip$", full.names = TRUE)) { utils::unzip(z, exdir = out); unlink(z) }
  if (grepl("[.]rds[.]gz$", f)) {
    # GSE222703 deposits a gzip of an already-gzipped .rds
    con_in <- gzfile(f, "rb"); con_out <- file(file.path(out, sub("[.]gz$", "", basename(f))), "wb")
    while (length(b <- readBin(con_in, "raw", 1e7))) writeBin(b, con_out)
    close(con_in); close(con_out)
  }
  writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), marker)
  invisible("extracted")
}
for (i in which(files$status %in% c("present", "downloaded") & grepl("[.](tar|tar[.]gz|rds[.]gz)$", files$filename)))
  extract_bundle(files$dest[i])

# ---- checksums -------------------------------------------------------------------
# Hashes are cached against size and modification time so large inputs are not
# re-hashed on every run.
hash_cache_f <- file.path(RAW_DIR, ".sha256_cache.tsv")
hc <- if (file.exists(hash_cache_f)) fread(hash_cache_f, colClasses = "character") else
  data.table(dest = character(), bytes = character(), mtime = character(), sha256 = character())
files[, `:=`(bytes = NA_real_, sha256 = NA_character_)]
for (i in which(file.exists(files$dest))) {
  d <- files$dest[i]; b <- as.character(file.size(d)); mt <- as.character(as.numeric(file.mtime(d)))
  hit <- hc[dest == d & bytes == b & mtime == mt]
  sh <- if (nrow(hit)) hit$sha256[1] else {
    msg("hashing ", basename(d)); s <- sha256_file(d)
    hc <- rbind(hc[dest != d], data.table(dest = d, bytes = b, mtime = mt, sha256 = s)); s }
  files[i, `:=`(bytes = as.numeric(b), sha256 = sh)]
  if (is.na(files$retrieved_utc[i]))
    files[i, retrieved_utc := format(file.mtime(d), tz = "UTC", usetz = TRUE)]
}
fwrite(hc, hash_cache_f, sep = "\t")

log_out <- files[, .(dataset_name, sample_id, role, filename, source_url = url,
                     local_path = sub(paste0("^", SUB_ROOT, "/"), "", dest),
                     status, retrieved_utc, bytes, sha256)]
save_tsv(log_out, "01_acquisition_log.tsv")

# ---- resolved manifest ------------------------------------------------------------
res_man <- copy(man)
for (i in seq_len(nrow(res_man))) {
  ds <- res_man$dataset_name[i]
  fx <- files[dataset_name == ds & fetch == TRUE]
  expr <- fx[role %in% c("expression_raw_counts", "expression_and_annotation", "expression")]
  meta <- fx[role %in% c("tisch2_cell_metainfo", "author_annotation", "expression_and_annotation")]
  rel <- function(x) paste(unique(dirname(sub(paste0("^", SUB_ROOT, "/"), "", x))), collapse = ";")
  if (nrow(expr)) res_man$local_expression_path[i] <- rel(expr$dest)
  if (nrow(meta)) res_man$local_metadata_path[i] <- rel(meta$dest)
  # dataset-level digest: SHA-256 over the sorted per-file digests
  if (nrow(fx) && all(!is.na(fx$sha256)))
    res_man$sha256[i] <- digest::digest(paste(sort(fx$sha256), collapse = ""), algo = "sha256", serialize = FALSE)
  # TISCH2 expression summaries and documentation are audit-only and not required
  req <- fx[!role %in% c("tisch2_expression_lognorm", "documentation")]
  res_man[i, acquisition_status := if (!nrow(req)) "deferred" else
    if (all(req$status %in% c("present", "downloaded"))) "complete" else "incomplete"]
}
save_tsv(res_man, "01_manifest_resolved.tsv")

summ <- files[, .(n_files = .N, n_ok = sum(status %in% c("present", "downloaded")),
                  n_missing = sum(status %in% c("missing", "failed")),
                  n_deferred = sum(status == "deferred"), GB = round(sum(bytes, na.rm = TRUE) / 1e9, 3)),
              by = dataset_name]
print(summ)
msg("acquisition log written; ", sum(files$status %in% c("missing", "failed")), " required files missing")
write_session_info("01_validate_manifest")
