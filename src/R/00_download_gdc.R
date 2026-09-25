# 00_download_gdc.R: fetch the TCGA-KIRC files that 01 reads
# Downloads the discovery cohort from the Genomic Data Commons into the local
# GDCdata/ layout that 01_build_data.R reads:
#   GDCdata/TCGA-KIRC/Transcriptome_Profiling/Gene_Expression_Quantification/
#       <file_id>/<file_name>        every open-access STAR-Counts file
#   GDCdata/TCGA-KIRC/Clinical/Clinical_Supplement/
#       <file_id>/<file_name>        every BCR XML clinical supplement
# Files already present are checked against their GDC md5. Only missing or
# corrupt files are fetched, in chunks of 40. Other cohorts are fetched by the
# stages that use them (08 CPTAC-3, 11 TCGA-KIRP and TCGA-KICH, 38 all TCGA).
# The analysis used GDC Data Release 46.0. The GDC serves only its current
# release, which is recorded in results/00_gdc_download_manifest.tsv.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(httr) })
banner("00 | Download the TCGA-KIRC files from the GDC")

with_retry <- function(f, tries = 4, wait = 10) {
  for (i in seq_len(tries)) {
    out <- tryCatch(f(), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    msg("  retry ", i, "/", tries, ": ", conditionMessage(out)); Sys.sleep(wait * i)
  }
  stop(out)
}
gdc_post <- function(endpoint, body) with_retry(function() {
  body$format <- "JSON"
  r <- httr::POST(paste0("https://api.gdc.cancer.gov/", endpoint),
                  body = jsonlite::toJSON(body, auto_unbox = TRUE),
                  httr::content_type_json(), httr::accept_json(), httr::timeout(600))
  httr::stop_for_status(r)
  jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), simplifyVector = TRUE)
})

release <- with_retry(function() jsonlite::fromJSON(httr::content(
  httr::GET("https://api.gdc.cancer.gov/status", httr::accept_json()), "text", encoding = "UTF-8")))$data_release
msg("GDC ", release)

manifest <- function(filters) {
  h <- gdc_post("files", list(filters = list(op = "and", content = filters),
                              fields = "file_id,file_name,md5sum,file_size", size = "5000"))$data$hits
  as.data.table(h)[, .(file_id, file_name, md5 = md5sum, file_size)]
}
f_proj <- list(op = "=", content = list(field = "cases.project.project_id", value = PROJECT_ID))
f_open <- list(op = "=", content = list(field = "access", value = "open"))
sets <- list(
  expression = list(dir = EXPR_DIR, m = manifest(list(f_proj, f_open,
    list(op = "=", content = list(field = "data_type", value = "Gene Expression Quantification")),
    list(op = "=", content = list(field = "analysis.workflow_type", value = "STAR - Counts"))))),
  clinical = list(dir = CLIN_XML_DIR, m = manifest(list(f_proj, f_open,
    list(op = "=", content = list(field = "data_type", value = "Clinical Supplement")),
    list(op = "=", content = list(field = "data_format", value = "BCR XML"))))))

fetch <- function(ids, dir) {
  tf <- tempfile(fileext = ".tar.gz")
  with_retry(function() {
    r <- httr::POST("https://api.gdc.cancer.gov/data",
                    body = jsonlite::toJSON(list(ids = ids), auto_unbox = FALSE),
                    httr::content_type_json(), httr::write_disk(tf, overwrite = TRUE), httr::timeout(3600))
    if (httr::status_code(r) != 200) stop("HTTP ", httr::status_code(r))
  })
  tf
}
out <- list()
for (nm in names(sets)) {
  S <- sets[[nm]]; m <- copy(S$m); dir.create(S$dir, recursive = TRUE, showWarnings = FALSE)
  m[, path := winlong(file.path(S$dir, file_id, file_name))]
  ok <- file.exists(m$path)
  ok[ok] <- unname(tools::md5sum(m$path[ok])) == m$md5[ok]
  todo <- m[!ok]
  msg(nm, ": ", nrow(m), " files at the GDC, ", sum(ok), " already present and verified, ", nrow(todo), " to fetch")
  if (nrow(todo)) for (ch in split(todo$file_id, ceiling(seq_len(nrow(todo)) / 40))) {
    ids <- ch
    for (attempt in 1:4) {
      tf <- fetch(ids, S$dir)
      if (length(ids) == 1) {
        row <- m[file_id == ids]; dir.create(file.path(S$dir, row$file_id), showWarnings = FALSE)
        file.copy(tf, winlong(file.path(S$dir, row$file_id, row$file_name)), overwrite = TRUE)
      } else suppressWarnings(tryCatch(untar(tf, exdir = S$dir), error = function(e) NULL))
      unlink(tf)
      r <- m[file_id %in% ids]
      good <- file.exists(r$path); good[good] <- unname(tools::md5sum(r$path[good])) == r$md5[good]
      if (any(!good)) unlink(r$path[!good])
      ids <- r$file_id[!good]
      if (!length(ids)) break
      msg("  ", length(ids), " file(s) missing or failing md5; fetching again")
    }
    if (length(ids)) stop(nm, ": ", length(ids), " files could not be downloaded intact")
  }
  out[[nm]] <- m[, .(set = nm, file_id, file_name, md5, file_size)]
}
man <- rbindlist(out)[, gdc_release := release]
save_tsv(man, "00_gdc_download_manifest.tsv")
msg("done: ", nrow(man), " files verified against the GDC md5")
write_session_info("00_download_gdc")
