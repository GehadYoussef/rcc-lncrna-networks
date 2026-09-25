# 12_render_report.R: data-audit report and data-quality gate (run after 01-05).
# Evaluates the feasibility thresholds in singlecell.yml and writes them to
# results/singlecell/12_go_no_go_criteria.tsv, which 06 checks before running.
# The report (reports/singlecell_data_audit.md) covers the dataset manifest and
# provenance, the schema and QC audit, lncRNA coverage per dataset and the full
# configuration used.

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
banner("12  data audit report")
start_log("12_render_report")

FE <- CFG$feasibility; PB <- CFG$pseudobulk
rd <- function(f) if (file.exists(file.path(RESULTS_DIR, f))) read_tsv(f) else data.table()
md_table <- function(dt, digits = 3) {
  if (!nrow(dt)) return("_(no rows)_\n")
  dt <- copy(as.data.table(dt))
  for (c in names(dt)) if (is.numeric(dt[[c]])) dt[[c]] <- format(round(dt[[c]], digits), big.mark = ",", trim = TRUE)
  dt[is.na(dt)] <- ""
  esc <- function(x) gsub("\\|", "\\\\|", as.character(x))
  paste0("| ", paste(names(dt), collapse = " | "), " |\n",
         "|", paste(rep("---", ncol(dt)), collapse = "|"), "|\n",
         paste0(apply(dt, 1, function(r) paste0("| ", paste(esc(r), collapse = " | "), " |")), collapse = "\n"), "\n")
}

man  <- rd("01_manifest_resolved.tsv")
acq  <- rd("01_acquisition_log.tsv")
sch  <- rd("02_schema_audit.tsv")
join <- rd("02_label_join_audit.tsv")
conc <- rd("03_label_concordance.tsv")
flg  <- rd("03_sample_flags.tsv")
dup  <- rd("03_duplicate_library_check.tsv")
lma  <- rd("04_lncRNA_matrix_audit.tsv")
ret  <- rd("04_tisch2_gene_retention.tsv")
ref  <- rd("04_reference_provenance.tsv")
exc  <- rd("05_cell_exclusions.tsv")
mk   <- rd("05_marker_check.tsv")
dcc  <- rd("05_donor_compartment_cells.tsv")
cov  <- rd("05_lncRNA_donor_coverage.tsv")
det  <- rd("05_lncRNA_detection_by_compartment.tsv")
mdis <- rd("05_marker_discordance.tsv")
lnc_pass <- function(ds) { f <- file.path(DERIVED_DIR, paste0(ds, ".lnc_malignant_passing.rds"))
  if (file.exists(f)) readRDS(f)[passing == TRUE, unique(gene_key)] else character() }

# ---- data-quality gate ----
ev <- function(ds, comp) if (nrow(dcc)) dcc[dataset == ds & compartment == comp & evaluable == TRUE, uniqueN(donor_id)] else 0L
disc <- man[analysis_role == "discovery", dataset_name]
disc <- if (length(disc)) disc[1] else NA_character_
repl_ds <- man[analysis_role %in% c("replication", "replication_candidate"), dataset_name]
repl_mal <- if (nrow(dcc)) dcc[dataset %in% repl_ds & compartment == "malignant" & evaluable == TRUE,
                               .(dataset, donor_id)] else data.table(dataset = character(), donor_id = character())
lnc_disc <- if (nrow(cov) && !is.na(disc)) cov[dataset == disc, lnc_passing_donor_filter] else 0L
if (!length(lnc_disc)) lnc_disc <- 0L
disc_keys <- if (!is.na(disc)) lnc_pass(disc) else character()
overlap <- rbindlist(lapply(repl_ds, function(ds) data.table(dataset = ds, replication_passing = length(lnc_pass(ds)),
                                                             overlap_with_discovery = length(intersect(disc_keys, lnc_pass(ds))))))
overlap_any <- length(intersect(disc_keys, unique(unlist(lapply(repl_ds, lnc_pass)))))
save_tsv(overlap, "12_lncRNA_discovery_replication_overlap.tsv")

crit <- data.table(
  criterion = c(
    sprintf("Discovery dataset (%s) imported with integer counts", disc),
    sprintf("Discovery: evaluable malignant donors >= %d", FE$min_donors_malignant_discovery),
    sprintf("Discovery: evaluable normal-epithelial donors >= %d", FE$min_donors_normal_epithelial_discovery),
    sprintf("Discovery: lncRNAs passing donor-level filter >= %d", FE$min_lnc_passing_discovery),
    sprintf("Replication: independent evaluable malignant donors >= %d", FE$min_replication_malignant_donors),
    "Replication: >= 2 datasets with evaluable malignant donors",
    sprintf("Replication: discovery-passing lncRNAs also passing in >= 1 replication dataset >= %d", FE$min_lnc_overlap_replication)),
  observed = c(
    if (nrow(sch)) paste(sch[dataset == disc, paste(status, integer_counts)], collapse = "") else "not run",
    ev(disc, "malignant"), ev(disc, "normal_epithelial"), lnc_disc,
    nrow(repl_mal), uniqueN(repl_mal$dataset), overlap_any))
crit[, pass := c(
  nrow(sch) > 0 && isTRUE(sch[dataset == disc, status == "imported" & integer_counts == TRUE]),
  ev(disc, "malignant") >= FE$min_donors_malignant_discovery,
  ev(disc, "normal_epithelial") >= FE$min_donors_normal_epithelial_discovery,
  lnc_disc >= FE$min_lnc_passing_discovery,
  nrow(repl_mal) >= FE$min_replication_malignant_donors,
  uniqueN(repl_mal$dataset) >= 2,
  overlap_any >= FE$min_lnc_overlap_replication)]
discovery_go <- all(crit$pass[1:4])
replication_go <- all(crit$pass[5:7])
tisch2_only_mal <- if (nrow(dcc)) dcc[dataset %in% man[tisch_id != "", dataset_name] & dataset %in% repl_ds &
                                      compartment == "malignant" & evaluable == TRUE, uniqueN(paste(dataset, donor_id))] else 0L
decision <- if (discovery_go && replication_go) "GO" else if (discovery_go) "CONDITIONAL GO (discovery only)" else "NO-GO"
save_tsv(crit, "12_go_no_go_criteria.tsv")

# ---- report -----------------------------------------------------------------------------
L <- character()
add <- function(...) L <<- c(L, paste0(...))
add("# Single-cell lncRNA discovery: data audit\n")
add(sprintf("Generated %s from `src/singlecell` (config sha256 `%s`).\n",
            format(Sys.time(), "%Y-%m-%d %H:%M %Z"), substr(CFG_HASH, 1, 16)))
add("No differential expression has been run. No survival, stage, grade or response field was read by any step.\n")

add("## Recommendation: **", decision, "**\n")
add(md_table(crit))
add("\n")
add(sprintf("- Tumour specificity (malignant vs normal renal epithelium) is evaluable only where normal-epithelial donors exist: %s.",
            if (nrow(dcc)) paste(dcc[compartment == "normal_epithelial" & evaluable == TRUE, .(n = uniqueN(donor_id)), by = dataset][
              , paste0(dataset, " (", n, ")")], collapse = ", ") else "none"))
add(sprintf("- Of the TISCH2 datasets other than the discovery set, %d donor(s) have an evaluable malignant compartment.", tisch2_only_mal))
add(if (!discovery_go) "- Discovery criteria are not met: candidate discovery must not proceed on these files. Options are raw-read re-quantification with a lncRNA-aware reference, or a full-length/total-RNA single-cell dataset." else
      "- Discovery criteria are met in the discovery dataset.")
add(if (!replication_go) "- Replication criteria are not met: any candidate list would rest on a single cohort. Add the deferred cohorts (manifest rows marked `replication_candidate`) before phase 6, or restrict claims to single-cohort discovery with leave-one-donor-out stability." else
      "- Replication criteria are met.")
add("\n")

add("## 1. Dataset manifest and provenance\n")
if (nrow(man)) add(md_table(man[, .(dataset_name, geo_accession, analysis_role, treatment, primary_metastatic,
                                     expression_type, integer_counts_available, normal_tissue_available,
                                     acquisition_status, sha256 = substr(sha256, 1, 12))]))
if (nrow(acq)) {
  add("\nAcquisition summary (full per-file log with URLs, retrieval time and SHA-256 in `results/singlecell/01_acquisition_log.tsv`):\n")
  add(md_table(acq[, .(files = .N, retrieved = sum(status %in% c("present", "downloaded")),
                       deferred = sum(status == "deferred"), missing = sum(status %in% c("missing", "failed")),
                       GB = sum(bytes, na.rm = TRUE) / 1e9), by = dataset_name]))
}
if (length(CFG$exclude_donors)) {
  add("\nExcluded donors (configured in singlecell.yml):\n")
  for (ds in names(CFG$exclude_donors)) for (dn in names(CFG$exclude_donors[[ds]]))
    add("- ", ds, " ", dn, ": ", CFG$exclude_donors[[ds]][[dn]])
}
add("\nDataset notes:\n")
if (nrow(man)) for (i in seq_len(nrow(man))) add("- **", man$dataset_name[i], "**: ", man$notes[i])
add("\n")

add("## 2. Schema and QC audit\n")
add(md_table(sch))
if (nrow(join)) { add("\nCell-label joins (annotated cells found in the raw matrices):\n"); add(md_table(join)) }
if (nrow(conc)) {
  add("\nAuthor vs TISCH2 compartment agreement (cells labelled by both):\n")
  add(md_table(conc[, .(cells = sum(N), agree = sum(N[agree == TRUE]), pct_agree = 100 * sum(N[agree == TRUE]) / sum(N)), by = dataset]))
  add("\nLargest disagreements:\n")
  add(md_table(head(conc[as.character(agree) %in% c("FALSE", "false")][order(-N)], 12)))
}
rann <- rd("03_reference_annotation_summary.tsv")
if (nrow(rann)) {
  add("\nReference-based labels for datasets deposited without annotation (SingleR, protein-coding genes only):\n")
  add(md_table(rann))
}
if (nrow(flg)) { add("\nSample flags:\n"); add(md_table(flg[flags != ""])) }
add(sprintf("\nIdentical libraries across datasets/samples (sorted barcode + UMI signature): %d.\n", nrow(dup)))
if (nrow(exc)) {
  add("\nCell exclusions by reason:\n")
  add(md_table(dcast(exc, dataset ~ reason, value.var = "N", fun.aggregate = sum)))
}
if (nrow(dcc)) {
  add(sprintf("\nEvaluable donors per compartment (>= %d cells after QC):\n", PB$min_cells_per_donor_compartment))
  add(md_table(dcast(dcc[evaluable == TRUE], dataset + analysis_role ~ compartment, value.var = "donor_id",
                     fun.aggregate = function(x) length(unique(x)))))
}
if (nrow(mdis)) {
  add("\nMarker discordance (mean detection of the ccRCC panel by assigned compartment). Flagged labels probably hide malignant cells; they are reported, not relabelled:\n")
  add(md_table(mdis[order(dataset, -ccrcc_panel_mean_detection)], digits = 2))
}
if (nrow(mk)) {
  add("\nMarker check (fraction of cells detecting each marker, by assigned compartment; ccRCC markers shown):\n")
  add(md_table(mk[panel %in% c("malignant_ccrcc", "renal_epithelial")], digits = 2))
}
add("\n")

add("## 3. lncRNA coverage\n")
if (nrow(ref)) add("Reference: ", paste(ref[field %in% c("genome_build", "gencode_release", "n_lncRNA_genes", "sha256"),
                                             paste0(field, " = ", value)], collapse = "; "), "\n")
add("\nMatrix-level audit:\n"); add(md_table(lma))
add("\nDonor-level coverage in malignant cells:\n"); add(md_table(cov))
if (nrow(overlap)) { add("\nOverlap of discovery-passing lncRNAs with each replication dataset:\n"); add(md_table(overlap)) }
if (nrow(det)) { add("\nDetection prevalence by compartment:\n"); add(md_table(det)) }
if (nrow(ret)) { add("\nGenes retained in TISCH2 re-processed matrices:\n"); add(md_table(ret)) }
add("\nFigures: `figures/singlecell/05_composition.svg`, `figures/singlecell/05_lncRNA_coverage.svg`.\n")

add("## 4. Configuration used\n")
add("```yaml\n", paste(readLines(CFG_FILE), collapse = "\n"), "\n```\n")

out <- file.path(REPORT_DIR, "singlecell_data_audit.md")
writeLines(L, out)
msg("report written: ", out, " | decision: ", decision)
write_session_info("12_render_report")
