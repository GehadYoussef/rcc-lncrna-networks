# run_all.R: run the bulk pipeline in order
#   Rscript src/R/run_all.R
# Stages cache their output in data/derived/cache and skip work already done.
# Delete an .rds to force its stage to re-run:
#   expr_raw.rds           re-read the 614 STAR count files
#   gdc_file_map.rds        re-query the GDC file/barcode mapping
#   gdc_clinical.rds        re-query GDC harmonised clinical data
#   grade_xml.rds           re-parse the local BCR clinical XML files
#   dataset.rds             rebuild the analysis cohort and expression matrices
#   networks.rds            rebuild both WGCNA networks (slow)
#   module_loadings.rds     refit the discovery module rotations (03)
#   survival.rds            re-run the Cox models
#   estimate_scores.rds     re-run ESTIMATE on the discovery cohort (07)
#   star_qc.rds             rebuild the discovery STAR metric table (07)
#   lnc_global_axis.rds     recompute the leading lncRNA axis (07)
#   validation_dataset.rds  rebuild the CPTAC-3 cohort (08)
#   valid_star_qc.rds       re-read the CPTAC-3 STAR summaries (08)
#   locked_model.rds        re-lock the validation models (09)
#   subtype_*.rds           rebuild the KIRP/KICH cohorts (11)
#   biospecimen_kirc.rds    re-query GDC biospecimen data (23)
#   network_lncRNA_unadjusted.rds   rebuild the un-residualised lncRNA network (23)
#   module_preservation_*.rds       re-run modulePreservation (24)
# The network stages (02, 12 section D, 23, 24) and 38 are the slow ones.

R_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1], winslash = "/")) else getwd()
})

# 02b_lncrna_tuning.R is not run here. It prints candidate lncRNA network
# parameters, which are set by hand in 00_config.R (LNC_POWER, LNC_DEEPSPLIT,
# LNC_MERGE) so the values used are visible in the config.
# The distance-correlation stages (04, 04b) are optional. RUN_DCOR <- TRUE runs
# them after 03. Tables and figures run near the end because they read results
# from every earlier stage.
RUN_DCOR <- FALSE

# The single-cell pipeline (src/singlecell) runs between stages 35 and 36. Its
# stages 10-11 read the bulk caches from 01-08, and bulk stages 36, 42 and 44
# read its results. With RUN_SINGLECELL <- FALSE it is skipped and those three
# stages stop for lack of input.
RUN_SINGLECELL <- TRUE

steps <- c(# 00 skips files already present and verified.
           "00_download_gdc.R",
           "01_build_data.R", "02_wgcna.R", "07_purity_qc.R", "03_survival.R",
           if (RUN_DCOR) c("04_dcor.R", "04b_dcor_diagnostics.R"),
           "05_ml.R", "06_enrichment.R",
           "08_validation_data.R", "09_validate.R", "10_paper_analyses.R",
           "11_subtype_specificity.R", "12_model_diagnostics.R",
           "15_sensitivity_tissue_source_site.R", "18_signsplit_enrichment.R",
           "19_eigengene_correlation.R", "20_turquoise_units_and_locked_model.R",
           "22_locked_model_export.R", "23_technical_axis_extended.R",
           "24_module_preservation.R", "25_matched_normal_control.R",
           "26_published_signatures.R", "27_metric_biology.R",
           "28_mutations_and_lncRNA_classes.R", "29_normalisation_check.R",
           "30_endpoint_sensitivity.R", "33_sex_stratified.R",
           "34_proximal_tubule_reading.R", "35_axis_projection_subtypes.R",
           "SINGLECELL",
           # 36 reads the single-cell results, 37 the unadjusted-network cache written by 23.
           "36_singlecell_module_localisation.R", "37_network_rewiring.R",
           # Pan-cancer stages. 38 downloads every TCGA project (tens of GB) and
           # requires results/38_decision_rules_lock.json. 39-41 check the lock.
           # 42 reads the single-cell results.
           "38_pancancer_build.R", "39_pancancer_axis.R", "40_pancancer_network.R",
           "41_pancancer_metric_survival.R", "42_singlecell_truth_correction.R",
           "43_pancancer_dose_response.R",
           # 44 compares the residualisation with RUV-III (ruv package from CRAN).
           "44_ruv3_prps_comparison.R",
           "14_tables.R", "13_figures.R", "32_supplementary_figures.R",
           # 31 indexes every results file, so it runs last.
           "31_supplementary_data_and_index.R")

# Looks for the single-cell driver in src/singlecell for either layout.
run_singlecell <- function() {
  if (!RUN_SINGLECELL) { cat("RUN_SINGLECELL is FALSE: single-cell pipeline skipped\n"); return(invisible()) }
  cand <- c(file.path(dirname(R_DIR), "singlecell"),
            file.path(PROJECT_ROOT, "submission", "src", "singlecell"))
  sc <- cand[file.exists(file.path(cand, "run_all.R"))][1]
  if (is.na(sc)) stop("single-cell pipeline not found in ", paste(cand, collapse = " or "))
  sc_results <- file.path(dirname(dirname(sc)), "results", "singlecell")
  if (file.exists(file.path(sc_results, "bulk_model_lock.json"))) {
    cat("single-cell outputs present in ", sc_results, ": not re-run\n", sep = ""); return(invisible())
  }
  # The candidate lock is immutable: if it exists, only the bulk handoff stages
  # (10-11) run. Otherwise the whole single-cell pipeline runs.
  args <- if (file.exists(file.path(sc_results, "candidate_lock.json"))) "--phase9-only" else character(0)
  st <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(sc, "run_all.R")), args))
  if (st != 0) stop("single-cell pipeline failed with status ", st)
}

t_start <- Sys.time()
for (s in steps) {
  cat("\n\n########## ", s, " ##########\n\n", sep = "")
  t0 <- Sys.time()
  if (identical(s, "SINGLECELL")) run_singlecell() else source(file.path(R_DIR, s), echo = FALSE)
  cat("\n---- ", s, " finished in ",
      round(difftime(Sys.time(), t0, units = "mins"), 1), " min ----\n", sep = "")
}
cat("\nPipeline complete in ",
    round(difftime(Sys.time(), t_start, units = "mins"), 1), " min\n", sep = "")
