# 09_candidate_lock.R: immutable candidate lock (run after 08, before 10-11).
# Freezes the outcome-blind selection from 08. Every input table is first
# scanned for outcome-like columns and the lock is refused if any is found.
# lock_id is a SHA-256 over candidates, evidence, rules and input digests,
# without timestamps, so identical inputs reproduce it. A rerun with the same
# id only records a verification. A changed selection stops unless --new-lock
# is given, which archives the old lock under results/singlecell/lock_archive/.
# Outputs: results/singlecell/locked_lncRNA_candidates.tsv,
#          candidate_evidence_long.tsv, candidate_exclusions.tsv (first failed
#          rule per lncRNA), candidate_lock.json (provenance and lock id)

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
banner("09  candidate lock")
start_log("09_candidate_lock")
NEW_LOCK <- "--new-lock" %in% commandArgs(TRUE)
# --refresh-provenance: with an unchanged lock_id, re-stamps git and code provenance
REFRESH <- "--refresh-provenance" %in% commandArgs(TRUE)

man <- read_manifest()$manifest
disc_ds <- man[analysis_role == "discovery", dataset_name][1]
E  <- read_tsv("08_candidate_evaluation.tsv")
M  <- read_tsv("08_meta_lncRNA.tsv.gz")
DE <- read_tsv("07_de_lncRNA.tsv.gz")
acq <- read_tsv("01_acquisition_log.tsv")
ref <- read_tsv("04_reference_provenance.tsv")
pbs <- read_tsv("06_pseudobulk_samples.tsv")
resc <- if (file.exists(file.path(RESULTS_DIR, "05_malignant_rescue.tsv"))) read_tsv("05_malignant_rescue.tsv") else data.table()

# ---- outcome blindness check ---------------------------------------------------------
fed <- list(`08_candidate_evaluation` = E, `08_meta_lncRNA` = M, `07_de_lncRNA` = DE,
            `06_pseudobulk_samples` = pbs, `01_acquisition_log` = acq)
outcome_hits <- rbindlist(lapply(names(fed), function(n) {
  h <- grep(OUTCOME_PATTERN, names(fed[[n]]), ignore.case = TRUE, value = TRUE)
  if (length(h)) data.table(table = n, column = h) else NULL
}))
if (nrow(outcome_hits)) stop("outcome-like columns found in selection inputs: ",
                             paste(outcome_hits[, paste(table, column, sep = "$")], collapse = ", "))
cells_md <- list.files(DERIVED_DIR, pattern = "[.]coldata_qc[.]rds$", full.names = TRUE)
for (f in cells_md) assert_no_outcome_fields(readRDS(f), basename(f))
removed <- if (file.exists(file.path(RESULTS_DIR, "02_removed_outcome_fields.tsv"))) read_tsv("02_removed_outcome_fields.tsv") else data.table()

# ---- candidates ------------------------------------------------------------------------
sel <- E[tier %in% c("tier1_tumour_specific", "tier2_malignant_compartment_enriched")]
if (!nrow(sel)) {
  stop("no candidates selected in 08; nothing to lock")
}
sel[, tier_order := match(tier, c("tier1_tumour_specific", "tier2_malignant_compartment_enriched"))]
setorder(sel, tier_order, provisional_rank)
sel[, lock_rank := seq_len(.N)]
bulk <- read_tsv("04_bulk_lncRNA_features.tsv")
ids <- bulk[, .(tcga_gene_id = gene_id[cohort == "TCGA-KIRC"][1], cptac_gene_id = gene_id[cohort == "CPTAC-3"][1]), by = key]
cand <- merge(sel, ids, by.x = "gene_key", by.y = "key", all.x = TRUE, sort = FALSE)
setorder(cand, lock_rank)
cand_out <- cand[, .(lock_rank, tier, gene_key, symbol, tcga_gene_id, cptac_gene_id, provisional_score,
                     tme_pooled_log2FC = B_pooled_log2FC, tme_pooled_FDR = B_pooled_FDR, tme_I2 = B_I2, tme_k = B_k,
                     tme_repl_same_direction = B_repl_same_direction, tme_repl_tested = B_repl_tested,
                     normal_pooled_log2FC = A_pooled_log2FC, normal_pooled_FDR = A_pooled_FDR, normal_I2 = A_I2, normal_k = A_k,
                     normal_repl_same_direction = A_repl_same_direction, normal_repl_tested = A_repl_tested,
                     heterogeneity_warning = (B_heterogeneity_flag %in% TRUE) | (A_heterogeneity_flag %in% TRUE),
                     discovery_malignant_donor_detection = B_discovery_detection_malignant)]

# evidence: every per-dataset estimate (inferential and descriptive) plus pooled rows
ev_ds <- DE[gene_key %in% cand$gene_key & contrast %in% unique(M$contrast),
            .(gene_key, symbol, contrast, source = dataset, analysis_role, inference, paired, log2FC, SE, CI_low, CI_high,
              p = PValue, FDR = FDR_lncRNA, n_donors_ref, n_donors_test, donor_detection_ref, donor_detection_test,
              lodo_direction_flip)]
ev_pool <- M[gene_key %in% cand$gene_key,
             .(gene_key, symbol, contrast, source = "pooled_random_effects", analysis_role = "meta", inference = meta_method,
               paired = NA, log2FC = pooled_log2FC, SE = pooled_SE, CI_low = pooled_CI_low, CI_high = pooled_CI_high,
               p = pooled_p, FDR = pooled_FDR, n_donors_ref = NA_integer_, n_donors_test = NA_integer_,
               donor_detection_ref = NA_real_, donor_detection_test = NA_real_, lodo_direction_flip = !(loo_direction_stable %in% c(TRUE, NA)))]
evidence <- rbind(ev_ds, ev_pool)
setorder(evidence, gene_key, contrast, source)

exclusions <- E[tier == "not_selected", .(gene_key, symbol, first_failed_rule, tme_pooled_log2FC = B_pooled_log2FC,
                                          tme_pooled_FDR = B_pooled_FDR, normal_pooled_log2FC = A_pooled_log2FC,
                                          normal_pooled_FDR = A_pooled_FDR, discovery_donor_filter, bulk_expressed_both)]

# ---- rules and input digests ------------------------------------------------------------
rules <- list(pseudobulk = CFG$pseudobulk, pseudobulk_build = CFG$pseudobulk_build, cell_qc = CFG$cell_qc,
              exclude_donors = CFG$exclude_donors, rescue_malignant = CFG$rescue_malignant,
              reference_annotation = CFG$reference_annotation, feasibility = CFG$feasibility, de = CFG$de,
              replication = CFG$replication, meta_analysis = CFG$meta_analysis, candidate_rules = CFG$candidate_rules,
              candidate_ranking_weights = CFG$candidate_ranking_weights,
              lncrna_biotypes = CFG$annotation$lncrna_biotypes)
input_digests <- c(
  setNames(as.list(acq[!is.na(sha256), sha256]), acq[!is.na(sha256), paste0("raw:", dataset_name, "/", filename)]),
  list(`reference:gencode_v36_gtf` = ref[field == "sha256", value],
       `reference:gencode_legacy_gtf` = ref[field == "legacy_sha256", value][1] %||% NA_character_,
       `table:06_pseudobulk_samples` = table_digest(pbs),
       `table:05_malignant_rescue` = table_digest(resc),
       `table:07_de_lncRNA` = table_digest(DE),
       `table:08_meta_lncRNA` = table_digest(M),
       `table:08_candidate_evaluation` = table_digest(E)))
lock_id <- lock_digest(cand_out, evidence, rules, input_digests)
msg("lock_id ", lock_id)

# ---- code provenance ----------------------------------------------------------------------
git <- function(...) tryCatch(system2("git", c("-C", shQuote(SUB_ROOT), ...), stdout = TRUE, stderr = FALSE),
                              error = function(e) character(), warning = function(w) character())
head_commit <- git("rev-parse", "HEAD")
dirty <- git("status", "--porcelain", "--", "src/singlecell", "src/config")
code_files <- c(list.files(SC_DIR, pattern = "[.]R$", full.names = TRUE, recursive = TRUE),
                list.files(file.path(SUB_ROOT, "src", "config"), full.names = TRUE))
code_sha <- setNames(vapply(code_files, sha256_file, character(1)), sub(paste0("^", SUB_ROOT, "/"), "", normalizePath(code_files, winslash = "/")))

lock <- list(
  lock_id = lock_id,
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stage = "single-cell outcome-blind candidate selection (phase 8 lock)",
  git = list(head_commit = if (length(head_commit)) head_commit else NA,
             single_cell_code_uncommitted = length(dirty) > 0,
             uncommitted_paths = dirty,
             note = if (length(dirty)) "code/config not committed at lock time; per-file SHA-256 below identify the exact code" else ""),
  code_sha256 = as.list(code_sha),
  config_sha256 = CFG_HASH,
  annotation = list(genome_build = CFG$annotation$genome_build, gencode_release = CFG$annotation$gencode_release,
                    gencode_sha256 = ref[field == "sha256", value],
                    legacy_gencode_release = CFG$annotation$legacy_gencode_release),
  discovery_dataset = disc_ds,
  datasets = man[, .(dataset_name, analysis_role, geo_accession, sha256)],
  input_digests = input_digests,
  rules = rules,
  n_candidates = nrow(cand_out),
  n_tier1_tumour_specific = cand_out[tier == "tier1_tumour_specific", .N],
  n_tier2_malignant_compartment_enriched = cand_out[tier == "tier2_malignant_compartment_enriched", .N],
  candidates_ordered = cand_out[, .(lock_rank, tier, gene_key, symbol)],
  outcome_blindness = list(
    confirmed = TRUE,
    statement = paste("No survival, stage, grade, treatment-response or other outcome field was read by any step",
                      "that produced this selection. Outcome-like metadata columns were removed at import (listed in",
                      "removed_at_import) and every selection input table and per-cell metadata file was scanned",
                      "for outcome-like columns before locking; none were found. Bulk cohorts contributed only",
                      "expression-based feature presence (04_bulk_lncRNA_features)."),
    scanned_tables = names(fed), scanned_cell_metadata = basename(cells_md),
    removed_at_import = if (nrow(removed)) removed else list()),
  output_files = list(candidates = "results/singlecell/locked_lncRNA_candidates.tsv",
                      evidence = "results/singlecell/candidate_evidence_long.tsv",
                      exclusions = "results/singlecell/candidate_exclusions.tsv"),
  session = list(R = R.version.string, platform = R.version$platform,
                 packages = lapply(c(edgeR = "edgeR", metafor = "metafor", SingleR = "SingleR", data.table = "data.table",
                                     Matrix = "Matrix", SingleCellExperiment = "SingleCellExperiment"),
                                   function(p) as.character(utils::packageVersion(p))))
)

# ---- immutability ---------------------------------------------------------------------------
lock_f <- file.path(RESULTS_DIR, "candidate_lock.json")
outs <- c(cand = "locked_lncRNA_candidates.tsv", ev = "candidate_evidence_long.tsv", ex = "candidate_exclusions.tsv")
write_lock <- function(previous = NULL) {
  if (!is.null(previous)) lock$previous_lock_id <<- previous
  save_tsv(cand_out, outs[["cand"]]); save_tsv(evidence, outs[["ev"]]); save_tsv(exclusions, outs[["ex"]])
  lock$output_sha256 <<- as.list(setNames(vapply(file.path(RESULTS_DIR, outs), sha256_file, character(1)), outs))
  jsonlite::write_json(lock, lock_f, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
  for (f in c(lock_f, file.path(RESULTS_DIR, outs))) Sys.chmod(f, mode = "0444")
}
if (file.exists(lock_f)) {
  old <- jsonlite::read_json(lock_f)
  if (identical(old$lock_id, lock_id) && REFRESH) {
    arch <- file.path(RESULTS_DIR, "lock_archive", lock_id, paste0("provenance_", format(Sys.time(), "%Y%m%dT%H%M%S")))
    dir.create(arch, recursive = TRUE, showWarnings = FALSE)
    for (f in c(lock_f, file.path(RESULTS_DIR, outs))) if (file.exists(f)) { Sys.chmod(f, "0644"); file.copy(f, arch) }
    lock$provenance_refreshed_from <- old$created_utc
    write_lock()
    msg("lock ", substr(lock_id, 1, 16), " reproduced; provenance re-stamped (previous JSON archived in ", arch, ")")
  } else if (identical(old$lock_id, lock_id)) {
    msg("existing lock reproduced exactly (", substr(lock_id, 1, 16), "); files left untouched")
    ver <- data.table(verified_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), lock_id = lock_id, result = "reproduced",
                      git_head = if (length(head_commit)) head_commit else NA_character_,
                      single_cell_code_uncommitted = length(dirty) > 0)
    vf <- file.path(RESULTS_DIR, "candidate_lock_verifications.tsv")
    fwrite(ver, vf, sep = "\t", append = file.exists(vf))
  } else if (!NEW_LOCK) {
    stop("a different candidate lock already exists (", substr(old$lock_id, 1, 16), " vs ", substr(lock_id, 1, 16),
         "). The selection or its inputs changed. Rerun with --new-lock to archive the old lock and write a new one.")
  } else {
    arch <- file.path(RESULTS_DIR, "lock_archive", old$lock_id)
    dir.create(arch, recursive = TRUE, showWarnings = FALSE)
    for (f in c(lock_f, file.path(RESULTS_DIR, outs))) if (file.exists(f)) { Sys.chmod(f, "0644"); file.copy(f, arch, overwrite = FALSE); unlink(f) }
    msg("archived previous lock ", old$lock_id, " -> ", arch)
    write_lock(previous = old$lock_id)
  }
} else {
  write_lock()
}

msg("locked ", nrow(cand_out), " candidates (tier1 ", lock$n_tier1_tumour_specific, ", tier2 ",
    lock$n_tier2_malignant_compartment_enriched, "); code uncommitted: ", lock$git$single_cell_code_uncommitted)
write_session_info("09_candidate_lock")
