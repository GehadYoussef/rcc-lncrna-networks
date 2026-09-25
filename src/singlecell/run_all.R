# run_all.R: runs the single-cell pipeline stages in order.
#
#   Rscript src/singlecell/run_all.R                  full run (downloads if needed)
#   Rscript src/singlecell/run_all.R --no-download    use local inputs only
#   Rscript src/singlecell/run_all.R --new-lock       replace an existing candidate or bulk lock
#   Rscript src/singlecell/run_all.R --force          run 06 even if the data-quality gate fails
#   Rscript src/singlecell/run_all.R --report-only    rebuild the report (12) only
#   Rscript src/singlecell/run_all.R --phase6-only    pseudobulk and within-dataset DE (06-07)
#   Rscript src/singlecell/run_all.R --phase7-only    meta-analysis (08)
#   Rscript src/singlecell/run_all.R --phase8-only    candidate lock (09)
#   Rscript src/singlecell/run_all.R --phase9-only    bulk handoff and prognostic model (10-11)
#
# Any failing stage stops the run. The report (12) runs before 06 because 06
# checks its data-quality gate.

args <- commandArgs(trailingOnly = TRUE)
here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
rscript <- file.path(R.home("bin"), "Rscript")

stages <- c("01_validate_manifest.R", "02_import_dataset.R", "03_harmonize_metadata.R",
            "04_gene_annotation.R", "05_cell_qc_and_annotations.R", "12_render_report.R",
            "06_build_pseudobulk.R", "07_within_dataset_de.R", "08_meta_analysis.R", "09_candidate_lock.R",
            "10_bulk_handoff.R", "11_bulk_prognostic.R")
if ("--report-only" %in% args) stages <- "12_render_report.R"
if ("--phase6-only" %in% args) stages <- c("06_build_pseudobulk.R", "07_within_dataset_de.R")
if ("--phase7-only" %in% args) stages <- "08_meta_analysis.R"
if ("--phase8-only" %in% args) stages <- "09_candidate_lock.R"
if ("--phase9-only" %in% args) stages <- c("10_bulk_handoff.R", "11_bulk_prognostic.R")

pass <- intersect(args, c("--no-download", "--new-lock", "--force"))
t0 <- Sys.time()
for (s in stages) {
  message(format(Sys.time(), "[%H:%M:%S] "), "running ", s)
  st <- system2(rscript, c(shQuote(file.path(here, s)), pass))
  if (st != 0) stop("stage ", s, " failed with status ", st)
}
message("single-cell stages complete in ", format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
