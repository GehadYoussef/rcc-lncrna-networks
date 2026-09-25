# 10_bulk_handoff.R: carries locked candidates to TCGA-KIRC and CPTAC-3.
# Run after 09 and after the bulk pipeline (src/R) has built its caches in
# data/derived/cache. Verifies the single-cell lock digests, then, without reading any
# survival field, matches candidates to both bulk matrices and measures detection
# (FPKM >= detection_fpkm) and association (R2, Spearman rho) with the three STAR
# library-quality metrics, the technical axis that confounds lncRNA prognostic
# associations in bulk. Candidates failing the configured rules are removed and a
# second lock, chained to the single-cell lock, is written.
# Outputs: results/singlecell/10_*.tsv, bulk_retained_candidates.tsv and
#          bulk_candidate_lock.json (read-only), data/derived/singlecell/
#          bulk_candidate_expression.rds

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
banner("10  bulk handoff (no survival data)")
start_log("10_bulk_handoff")

BH <- CFG$bulk_handoff
pth <- function(x) normalizePath(file.path(SUB_ROOT, x), winslash = "/", mustWork = TRUE)

# ---- verify the single-cell lock ------------------------------------------------------
lock <- jsonlite::read_json(file.path(RESULTS_DIR, "candidate_lock.json"))
for (nm in names(lock$output_sha256)) {
  if (!identical(sha256_file(file.path(RESULTS_DIR, nm)), lock$output_sha256[[nm]]))
    stop("locked file ", nm, " does not match its digest in candidate_lock.json")
}
cand <- read_tsv("locked_lncRNA_candidates.tsv")
msg("single-cell lock ", substr(lock$lock_id, 1, 16), " verified: ", nrow(cand), " candidates")

# ---- bulk cohorts, outcome fields removed on load ------------------------------------------
strip <- function(dt) { dt <- as.data.table(dt); dt[, grep("^os_|vital|dtd|dtlf|death", names(dt), value = TRUE) := NULL]; dt }
d <- readRDS(bulk_path(BH$tcga_dataset))
tcl <- strip(d$cohort_adj)[, .(patient, sample_barcode, file_id)]
rm(d)
xr <- readRDS(bulk_path(BH$tcga_expr))
tq <- as.data.table(readRDS(bulk_path(BH$tcga_star_qc)))
v <- readRDS(bulk_path(BH$cptac_dataset))
vcl_all <- strip(v$cohort)
vcl <- vcl_all[complete.cases(vcl_all[, .(age, sex, T_stage, N_pos, M1, grade_num)]), .(patient, sample_barcode, file_id)]
vq <- as.data.table(v$qc)
for (t in list(tcl, tq, vcl, vq)) assert_no_outcome_fields(t, "bulk handoff inputs")

tcga <- list(name = "TCGA-KIRC", fpkm = xr$fpkm[, match(tcl$file_id, xr$sample_info$file_id)], ann = as.data.table(xr$gene_ann),
             samples = tcl$sample_barcode, qc = tq[match(tcl$sample_barcode, sample_barcode)])
colnames(tcga$fpkm) <- tcl$sample_barcode
cptac <- list(name = "CPTAC-3", fpkm = v$fpkm[, match(vcl$sample_barcode, colnames(v$fpkm))], ann = as.data.table(v$gene_ann),
              samples = vcl$sample_barcode, qc = vq[match(vcl$sample_barcode, sample_barcode)])
rm(xr, v); invisible(gc(verbose = FALSE))
for (co in list(tcga, cptac)) {
  if (anyNA(co$fpkm)) stop(co$name, ": expression columns missing for some patients")
  if (anyNA(co$qc$assigned_reads)) stop(co$name, ": STAR metrics missing for some patients")
}
msg("TCGA-KIRC n=", length(tcga$samples), "; CPTAC-3 n=", length(cptac$samples))

# ---- 1-3: availability, reliability, technical association ---------------------------------
H <- new.env(); H$R_DIR <- dirname(bulk_path(BH$helpers))
sys.source(bulk_path(BH$helpers), envir = H)
helpers_sha <- sha256_file(bulk_path(BH$helpers))

avail <- list(); tech <- list(); expr_out <- list()
for (co in list(tcga, cptac)) {
  key <- strip_ensembl_version(rownames(co$fpkm))
  idx <- lapply(cand$gene_key, function(k) which(key == k))
  n_match <- lengths(idx)
  a <- data.table(cohort = co$name, gene_key = cand$gene_key, symbol = cand$symbol, n_rows_matched = n_match,
                  status = fifelse(n_match == 0, "missing", fifelse(n_match > 1, "ambiguous", "present")))
  ok <- which(n_match == 1)
  E <- log2(t(co$fpkm[unlist(idx[ok]), , drop = FALSE]) + 1)             # samples x genes
  colnames(E) <- cand$gene_key[ok]
  F <- t(co$fpkm[unlist(idx[ok]), , drop = FALSE])
  a[ok, `:=`(detection_fraction = colMeans(F >= BH$detection_fpkm), median_log2fpkm = apply(E, 2, median),
             sd_log2fpkm = apply(E, 2, sd))]
  cov <- H$tech_covariates(co$qc)
  r2 <- technical_r2(E, cov)
  rho <- sapply(colnames(cov), function(m) apply(E, 2, function(g) suppressWarnings(cor(g, cov[, m], method = "spearman"))))
  tech[[co$name]] <- data.table(cohort = co$name, gene_key = colnames(E), technical_r2 = r2,
                                rho_noFeature = rho[, "pct_noFeature"], rho_multimapping = rho[, "pct_multimapping"],
                                rho_log_depth = rho[, "log_depth"])
  avail[[co$name]] <- a
  expr_out[[co$name]] <- list(log2fpkm = E, star = cov, samples = co$samples)
}
A <- rbindlist(avail); TA <- rbindlist(tech)
save_tsv(A, "10_bulk_feature_availability.tsv")
save_tsv(merge(TA, unique(A[, .(gene_key, symbol)]), by = "gene_key")[order(cohort, -technical_r2)], "10_technical_association.tsv")

# ---- 4: pre-specified removal and second lock --------------------------------------------------
w <- dcast(merge(A, TA, by = c("cohort", "gene_key"), all.x = TRUE),
           gene_key ~ cohort, value.var = c("status", "detection_fraction", "technical_r2"))
setnames(w, gsub("-", "_", names(w)))
dec <- merge(cand[, .(lock_rank, tier, gene_key, symbol)], w, by = "gene_key")[order(lock_rank)]
dec[, decision_reason := fcase(
  status_TCGA_KIRC != "present" | status_CPTAC_3 != "present", "not_uniquely_present_in_both_cohorts",
  detection_fraction_TCGA_KIRC < BH$min_detection_fraction | detection_fraction_CPTAC_3 < BH$min_detection_fraction,
    "unreliable_quantification",
  technical_r2_TCGA_KIRC >= BH$technical_r2_max | technical_r2_CPTAC_3 >= BH$technical_r2_max, "technically_sensitive",
  default = "retained")]
dec[, retained := decision_reason == "retained"]
save_tsv(dec, "10_bulk_candidate_decisions.tsv")
kept <- dec[retained == TRUE]
msg("retained ", nrow(kept), " of ", nrow(dec), " (tier1 ", kept[tier == "tier1_tumour_specific", .N], ", tier2 ",
    kept[tier == "tier2_malignant_compartment_enriched", .N], ")")
print(dec[, .N, by = .(tier, decision_reason)])

saveRDS(lapply(expr_out, function(e) { e$log2fpkm <- e$log2fpkm[, intersect(kept$gene_key, colnames(e$log2fpkm)), drop = FALSE]; e }),
        file.path(DERIVED_DIR, "bulk_candidate_expression.rds"))

rules <- CFG$bulk_handoff
bulk_lock_id <- digest::digest(paste(lock$lock_id, table_digest(kept[, .(lock_rank, tier, gene_key)]),
                                     digest::digest(jsonlite::toJSON(rules, auto_unbox = TRUE, digits = NA), algo = "sha256", serialize = FALSE),
                                     table_digest(TA), table_digest(A), sep = "|"), algo = "sha256", serialize = FALSE)
bl <- list(bulk_lock_id = bulk_lock_id, parent_single_cell_lock_id = lock$lock_id,
           created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
           rules = rules, helpers_sha256 = helpers_sha,
           cohorts = list(TCGA_KIRC = list(n = length(tcga$samples), patient_set = "analysis/cache/dataset.rds $cohort_adj"),
                          CPTAC_3 = list(n = length(cptac$samples), patient_set = "validation_dataset $cohort, complete CLIN_TERMS")),
           n_locked = nrow(dec), n_retained = nrow(kept),
           removed = dec[retained == FALSE, .(lock_rank, gene_key, symbol, decision_reason)],
           retained = kept[, .(lock_rank, tier, gene_key, symbol)],
           outcome_blindness = "No survival or outcome field was read: survival columns were dropped on load and all tables were checked by assert_no_outcome_fields().")
bf <- file.path(RESULTS_DIR, "bulk_candidate_lock.json"); rf <- file.path(RESULTS_DIR, "bulk_retained_candidates.tsv")
if (file.exists(bf)) {
  old <- jsonlite::read_json(bf)
  if (identical(old$bulk_lock_id, bulk_lock_id)) {
    msg("bulk candidate lock reproduced exactly (", substr(bulk_lock_id, 1, 16), ")")
  } else if (!"--new-lock" %in% commandArgs(TRUE)) {
    stop("a different bulk candidate lock exists; rerun with --new-lock to archive it and replace")
  } else {
    arch <- file.path(RESULTS_DIR, "lock_archive", paste0("bulk_", old$bulk_lock_id)); dir.create(arch, recursive = TRUE, showWarnings = FALSE)
    for (f in c(bf, rf)) { Sys.chmod(f, "0644"); file.copy(f, arch); unlink(f) }
    bl$previous_bulk_lock_id <- old$bulk_lock_id
  }
}
if (!file.exists(bf)) {
  save_tsv(kept[, .(lock_rank, tier, gene_key, symbol)], basename(rf))
  bl$retained_sha256 <- sha256_file(rf)
  jsonlite::write_json(bl, bf, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
  Sys.chmod(c(bf, rf), "0444")
  msg("bulk candidate lock written: ", substr(bulk_lock_id, 1, 16))
}
write_session_info("10_bulk_handoff")
