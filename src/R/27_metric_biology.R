# 27_metric_biology.R: biological correlates of the non-feature read fraction and what accounts for its hazard
# Tests whether the non-feature fraction tracks a tissue-intrinsic RNA state
# (hypoxia, proliferation, composition) or pre-analytical variation.
# The observed protein-coding matrix carries the metric itself, so scores
# computed on it (observed eigengenes, principal components, random gene sets)
# are upper bounds, not biological shares. Negative controls and residualised
# counterparts are reported beside them.
# Sections: axis checks and scores, correlations, Cox attenuation with
# collinearity diagnostics, the converse, reverse check, variance partition.
# Inputs: caches only (no network access) plus results of stages 06, 07 and 23.
# Outputs: 27_*.tsv tables.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(survival) })
banner("27 | Biology of the non-feature fraction and attenuation of its hazard")
set.seed(SEED)
t_start <- Sys.time()

RANDOM_SET_DRAWS <- 20L   # random gene sets for the variance-partition control
SIM_DRAWS        <- 300L  # synthetic collinear proxies for the attenuation null

ds     <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets   <- readRDS(file.path(CACHE_DIR, "networks.rds"))
est    <- as.data.table(readRDS(file.path(CACHE_DIR, "estimate_scores.rds")))
axis   <- as.data.table(readRDS(file.path(CACHE_DIR, "lnc_global_axis.rds")))
qc     <- as.data.table(readRDS(file.path(CACHE_DIR, "star_qc.rds")))
cohort <- as.data.table(ds$cohort_full)
if (!"tss" %in% names(cohort)) cohort[, tss := tstrsplit(sample_barcode, "-", keep = 2)[[1]]]
msg("Discovery cohort: ", nrow(cohort), " patients, ", sum(cohort$os_event), " deaths")

# Biospecimen measurements from stage 23. Plate and batch are read as character
# so that leading zeros ("0864") survive.
bio_f <- file.path(RESULTS_DIR, "23_biospecimen_kirc.tsv")
stopifnot(file.exists(bio_f))
bio <- fread(bio_f, colClasses = list(character = c("plate", "batch", "aliquot_id")))
stopifnot(all(c("rin", "plate", "batch", "pct_tumor_nuclei_mean", "pct_stromal_mean",
                "pct_necrosis_mean") %in% names(bio)))

# ---- helpers ----------------------------------------------------------------
# Levels with fewer than min_n samples are pooled into "other".
pool_levels <- function(x, min_n = 10) {
  x <- as.character(x); tb <- table(x)
  x[!is.na(x) & x %in% names(tb)[tb < min_n]] <- "other"
  factor(x)
}
spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(c(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}
tertile <- function(x) {
  q <- quantile(x, c(1 / 3, 2 / 3), na.rm = TRUE)
  cut(x, c(-Inf, q, Inf), labels = c("T1 (low)", "T2", "T3 (high)"))
}
zc <- function(x) as.numeric(scale(x))
# Leading PC as stage 07 defines it: column-centred, unscaled, sign-aligned to
# per-sample mean expression.
pc1_07 <- function(E) {
  s  <- svd(scale(E, center = TRUE, scale = FALSE), nu = 1, nv = 1)
  pc <- s$u[, 1]
  if (stats::cor(pc, rowMeans(E)) < 0) pc <- -pc
  list(pc = setNames(pc, rownames(E)), var_share = s$d[1]^2 / sum(s$d^2))
}
# Leading k PCs from the sample Gram matrix (n << p, so cheaper than a full SVD
# and equivalent up to sign).
leading_pcs <- function(E, k) {
  Ec <- scale(E, center = TRUE, scale = FALSE)
  ev <- eigen(tcrossprod(Ec), symmetric = TRUE)
  U  <- ev$vectors[, seq_len(k), drop = FALSE]
  dimnames(U) <- list(rownames(E), paste0("PC", seq_len(k)))
  U
}
# Unit-variance PC1 score of one gene set, with the fit_module_loadings()
# convention, computed from the Gram matrix so that many random sets are affordable.
set_eigengene <- function(E, genes) {
  X   <- E[, genes, drop = FALSE]
  mu  <- colMeans(X); sdv <- apply(X, 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
  Z   <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
  ev  <- eigen(tcrossprod(Z), symmetric = TRUE)
  s   <- ev$vectors[, 1]
  if (stats::cor(s, rowMeans(Z)) < 0) s <- -s
  setNames(zc(s), rownames(E))
}
# Disjoint random gene sets of the given sizes, mirroring the disjoint modules.
random_sets <- function(sizes, universe) {
  perm <- sample(universe); idx <- c(0L, cumsum(sizes))
  lapply(seq_along(sizes), function(i) perm[(idx[i] + 1L):idx[i + 1L]])
}
# Variance inflation factor of `target` given `others` (car::vif does not apply
# to a stratified coxph).
vif_of <- function(target, others, dat) {
  others <- others[nzchar(others)]
  if (!length(others)) return(1)
  r2 <- summary(lm(as.formula(paste(target, "~", paste(others, collapse = " + "))),
                   data = dat))$r.squared
  1 / (1 - r2)
}
nf_of <- function(ids) qc$pct_noFeature[match(ids, qc$sample_barcode)]

# ---- 1. Axes, biological scores and the per-sample table ----
banner("1 | Axes, biological scores and the per-sample table")
E_obs <- obs_expr(nets$mrna); E_res <- nets$mrna$expr
stopifnot(identical(rownames(E_obs), rownames(E_res)),
          identical(colnames(E_obs), colnames(E_res)))
msg("Protein-coding matrices: ", nrow(E_obs), " samples x ", ncol(E_obs), " genes")

pc_mrna <- pc1_07(E_obs)
pc_lnc  <- pc1_07(obs_expr(nets$lnc))
s_pc    <- spearman(pc_mrna$pc, nf_of(names(pc_mrna$pc)))

# The recomputed lncRNA PC1 is compared with the cached axis value by value,
# because a correlation is blind to scale and location.
lnc_cached <- axis$lnc_axis[match(names(pc_lnc$pc), axis$sample_barcode)]
stopifnot(!anyNA(lnc_cached))
lnc_maxdiff <- max(abs(pc_lnc$pc - lnc_cached))
lnc_alleq   <- all.equal(as.numeric(pc_lnc$pc), as.numeric(lnc_cached))
lnc_alleq_s <- if (isTRUE(lnc_alleq)) {
  "TRUE (all.equal on the values reports no difference)"
} else {
  paste(lnc_alleq, collapse = "; ")
}

# The protein-coding PC1 is checked against the values reported by stage 07.
r07a <- fread(file.path(RESULTS_DIR, "07_pc_variance_explained.tsv"))[
  matrix == "protein_coding_observed" & pc == 1]
r07b <- fread(file.path(RESULTS_DIR, "07_axis_vs_purity_qc_correlations.tsv"))[
  matrix == "protein_coding" & variable == "pct_noFeature"]
axis_check <- data.table(
  quantity = c("protein_coding_PC1_var_explained", "protein_coding_PC1_rho_noFeature",
               "protein_coding_PC1_n", "lncRNA_PC1_vs_cached_lnc_global_axis"),
  check = c(rep("absolute difference from the value reported in 07", 3),
            sprintf("element-wise identity of the %d axis values (all.equal + maximum absolute difference)",
                    length(lnc_cached))),
  recomputed_27 = c(round(pc_mrna$var_share, 4), round(s_pc[["rho"]], 3), s_pc[["n"]], NA_real_),
  reported_07   = c(r07a$var_explained, r07b$spearman_rho, r07b$n, NA_real_),
  n = c(s_pc[["n"]], s_pc[["n"]], s_pc[["n"]], length(lnc_cached)),
  max_abs_diff = c(abs(round(pc_mrna$var_share, 4) - r07a$var_explained),
                   abs(round(s_pc[["rho"]], 3) - r07b$spearman_rho),
                   abs(s_pc[["n"]] - r07b$n), lnc_maxdiff),
  tolerance = c(5e-4, 5e-3, 0, 1e-8))
axis_check[, agree := max_abs_diff <= tolerance]
axis_check[, note := c(rep("07_pc_variance_explained.tsv / 07_axis_vs_purity_qc_correlations.tsv", 3),
                       paste("cache/lnc_global_axis.rds; all.equal result:", lnc_alleq_s))]
save_tsv(axis_check, "27_axis_check.tsv"); print(axis_check)
if (!all(axis_check$agree)) warning("27: recomputed axis differs from 07 (see 27_axis_check.tsv)")

# ---- hypoxia and proliferation scores, lncRNA FPKM share (all genes) --------
# Scores use the full FPKM matrix, so no signature gene is lost to the network
# MAD filter. Score = mean per-gene z-score of log2(FPKM + 1) across the cohort.
# Hypoxia: Buffa et al., Br J Cancer 2010 (PMID 20087356), via MSigDB M34030,
# which maps 50 of 51 source probes. ALDOA is added from the publication.
# Probe 212639_x_at is annotated TUBA1B by Affymetrix but TUBA1A by MSigDB.
# TUBA1B is primary and TUBA1A a sensitivity variant.
HYPOXIA_MSIGDB <- c(
  "ACOT7", "ADM", "AK4", "ANKRD37", "ANLN", "BNIP3", "CA9", "CDKN3", "CHCHD2",
  "CORO1C", "CTSV", "DDIT4", "ENO1", "ESRP1", "GAPDH", "GPI", "HILPDA", "HK2",
  "KIF20A", "KIF4A", "LDHA", "LRRC42", "MAD2L2", "MAP7D1", "MCTS1", "MIF",
  "MRGBP", "MRPL13", "MRPL15", "MRPS17", "NDRG1", "P4HA1", "PFKP", "PGAM1",
  "PGK1", "PNP", "PSMA7", "PSRC1", "SEC61G", "SHCBP1", "SLC16A1", "SLC25A32",
  "SLC2A1", "TPI1", "TUBA1B", "TUBA1C", "TUBB6", "UTP11", "VEGFA", "YKT6")
HYPOXIA_ALT_TUBULIN <- "TUBA1A"          # MSigDB's mapping of probe 212639_x_at
HYPOXIA_PROV <- list(
  source_database    = "MSigDB (Human Molecular Signatures Database), collection C2:CGP",
  source_accession   = "M34030 / BUFFA_HYPOXIA_METAGENE",
  source_exact       = "Supplementary Table S5, common hypoxia signature ranked by common connectivity score",
  source_publication = "Buffa FM, Harris AL, West CM, Miller CJ. Br J Cancer 2010;102(2):428-435",
  source_pmid        = "20087356",
  source_platform    = "AFFY_HG_U133 (51 source probe sets, 50 mapped to symbols)",
  retrieval          = "gene-symbol export from download_geneset.jsp, retrieved 2026-09-10; hard-coded here from that export and from the publication")
hyp_symbols <- c(HYPOXIA_MSIGDB, "ALDOA")
hyp_mapnote <- c(
  ifelse(HYPOXIA_MSIGDB == "TUBA1B",
         "source probe 212639_x_at; MSigDB maps it to TUBA1A (NCBI 7846), the Affymetrix HG-U133 annotation and published re-implementations give TUBA1B (NCBI 10376); TUBA1B used, TUBA1A reported as a sensitivity variant",
         "MSigDB symbol export, one-to-one probe-to-symbol mapping"),
  "not in the MSigDB symbol export (one source probe, 238996_x_at, is unmapped); added from the publication")
stopifnot(length(hyp_symbols) == 51L, !anyDuplicated(hyp_symbols),
          length(hyp_mapnote) == 51L)
PROLIF_SYMBOLS <- c("MKI67", "TOP2A", "PCNA", "MCM2", "BUB1", "CCNB1", "CDK1",
                    "TYMS", "RRM2", "AURKA")
PROLIF_PROV <- list(
  source_database = "none (specified for this analysis)", source_accession = "",
  source_exact = "ten canonical proliferation markers", source_publication = "",
  source_pmid = "", source_platform = "",
  retrieval = "hard-coded for this analysis; score = mean z of log2(FPKM+1)")

msg("Reading expr_raw.rds for the hypoxia, proliferation and lncRNA-share scores ...")
raw <- readRDS(file.path(CACHE_DIR, "expr_raw.rds"))
# One column per file, so select the aliquot kept in stage 01 by file_id.
j <- match(cohort$file_id, raw$sample_info$file_id)
stopifnot(!anyNA(j), !anyDuplicated(j))
ann <- as.data.table(raw$gene_ann)
is_lnc <- !is.na(ann$gene_type) & ann$gene_type == "lncRNA"
tot_fpkm  <- colSums(raw$fpkm[, j, drop = FALSE], na.rm = TRUE)
lnc_share <- setNames(colSums(raw$fpkm[is_lnc, j, drop = FALSE], na.rm = TRUE) / tot_fpkm,
                      cohort$sample_barcode)

# Symbol -> one Ensembl gene: protein-coding entries preferred, then the
# highest median FPKM across the cohort (a few symbols carry a PAR_Y copy).
pick_genes <- function(symbols) rbindlist(lapply(symbols, function(s) {
  i <- which(!is.na(ann$gene_name) & ann$gene_name == s)
  if (!length(i)) return(data.table(symbol = s, gene_id = NA_character_, gene_type = NA_character_,
                                    n_ids_for_symbol = 0L, median_fpkm = NA_real_))
  med <- vapply(i, function(k) median(raw$fpkm[k, j]), numeric(1))
  o <- order(ann$gene_type[i] != "protein_coding", -med)[1]
  data.table(symbol = s, gene_id = ann$gene_id[i[o]], gene_type = ann$gene_type[i[o]],
             n_ids_for_symbol = length(i), median_fpkm = round(med[o], 3))
}))
score_set <- function(gt, label) {
  g <- gt[!is.na(gene_id)]
  M <- log2(t(raw$fpkm[match(g$gene_id, ann$gene_id), j, drop = FALSE]) + 1)
  rownames(M) <- cohort$sample_barcode; colnames(M) <- g$symbol
  sdv <- apply(M, 2, sd); keep <- is.finite(sdv) & sdv > 0
  msg(label, ": ", nrow(gt), " symbols, ", nrow(g), " found in the GDC annotation, ",
      sum(keep), " with non-zero variance; score = mean z of log2(FPKM+1)")
  list(score = setNames(rowMeans(scale(M[, keep, drop = FALSE])), rownames(M)),
       used = g$symbol[keep], M = M)
}
hyp_genes <- pick_genes(hyp_symbols); hyp <- score_set(hyp_genes, "Hypoxia (Buffa 2010)")
pro_genes <- pick_genes(PROLIF_SYMBOLS); pro <- score_set(pro_genes, "Proliferation")

# Hypoxia sensitivity variants: (a) genes at median FPKM >= 1 only, since a
# near-zero gene can proxy the metric. (b) TUBA1A in place of TUBA1B, and both.
sub_score <- function(sc, syms) {
  syms <- intersect(syms, colnames(sc$M))
  setNames(rowMeans(scale(sc$M[, syms, drop = FALSE])), rownames(sc$M))
}
hyp_hi_sym  <- hyp_genes[symbol %in% hyp$used & !is.na(median_fpkm) & median_fpkm >= 1, symbol]
hyp_lo_sym  <- setdiff(hyp$used, hyp_hi_sym)
hyp_hi      <- sub_score(hyp, hyp_hi_sym)
alt_genes   <- pick_genes(HYPOXIA_ALT_TUBULIN)
alt_M       <- log2(t(raw$fpkm[match(alt_genes$gene_id, ann$gene_id), j, drop = FALSE]) + 1)
rownames(alt_M) <- cohort$sample_barcode; colnames(alt_M) <- alt_genes$symbol
hyp_altM    <- cbind(hyp$M[, setdiff(hyp$used, "TUBA1B"), drop = FALSE], alt_M)
hyp_alt     <- setNames(rowMeans(scale(hyp_altM)), rownames(hyp_altM))
hyp_bothM   <- cbind(hyp$M[, hyp$used, drop = FALSE], alt_M)
hyp_both    <- setNames(rowMeans(scale(hyp_bothM)), rownames(hyp_bothM))
rm(raw); invisible(gc(verbose = FALSE))
msg("Hypoxia sensitivity: ", length(hyp_hi_sym), " of ", length(hyp$used),
    " genes at median FPKM >= 1 (excluded: ", paste(sort(hyp_lo_sym), collapse = ", "), ")")
msg("Pearson r of the full and median-FPKM >= 1 hypoxia scores: ",
    round(stats::cor(hyp$score, hyp_hi[names(hyp$score)]), 3),
    "; Spearman rho: ", round(spearman(hyp$score, hyp_hi[names(hyp$score)])[["rho"]], 3))
msg("Pearson r of the TUBA1B (primary) and TUBA1A (alternative mapping) hypoxia scores: ",
    round(stats::cor(hyp$score, hyp_alt[names(hyp$score)]), 4))

gene_table <- function(gt, sc, set_name, mapnote, prov, extra_flags = NULL) {
  nf <- nf_of(cohort$sample_barcode)
  ax <- axis$lnc_axis[match(cohort$sample_barcode, axis$sample_barcode)]
  gtab <- as.data.table(nets$mrna$gene_tbl)
  out <- rbindlist(lapply(seq_len(nrow(gt)), function(k) {
    gid  <- gt$gene_id[k]; sym <- gt$symbol[k]
    used <- !is.na(gid) && sym %in% sc$used
    v    <- if (used) sc$M[, sym] else rep(NA_real_, nrow(cohort))
    s1 <- spearman(v, nf); s2 <- spearman(v, ax)
    in_net <- !is.na(gid) && gid %in% gtab$gene_id
    data.table(gene_set = set_name, symbol = sym, gene_id = gid, gene_type = gt$gene_type[k],
               n_ids_for_symbol = gt$n_ids_for_symbol[k], median_fpkm = gt$median_fpkm[k],
               used_in_score = used,
               in_median_fpkm_ge1_subset = used && !is.na(gt$median_fpkm[k]) && gt$median_fpkm[k] >= 1,
               in_mrna_network = in_net,
               module = if (in_net) gtab$module[match(gid, gtab$gene_id)] else NA_character_,
               n = s1[["n"]], rho_noFeature = round(s1[["rho"]], 3), p_noFeature = signif(s1[["p"]], 3),
               rho_lnc_axis = round(s2[["rho"]], 3), p_lnc_axis = signif(s2[["p"]], 3),
               source_database = prov$source_database, source_accession = prov$source_accession,
               source_exact_source = prov$source_exact, source_publication = prov$source_publication,
               source_pmid = prov$source_pmid, source_platform = prov$source_platform,
               retrieval = prov$retrieval, mapping_note = mapnote[k])
  }))
  if (!is.null(extra_flags)) out <- cbind(out, extra_flags)
  out
}
hyp_tbl <- gene_table(hyp_genes, hyp, "Buffa_2010_hypoxia_metagene", hyp_mapnote, HYPOXIA_PROV)
# ALDOA comes from the publication, not the MSigDB export.
hyp_tbl[symbol == "ALDOA",
        `:=`(source_database = "publication (not in the MSigDB symbol export)",
             source_accession = "", source_platform = "")]
pro_tbl <- gene_table(pro_genes, pro, "proliferation_10_genes",
                      rep("hard-coded for this analysis", nrow(pro_genes)), PROLIF_PROV)
save_tsv(hyp_tbl, "27_hypoxia_gene_set.tsv")
save_tsv(pro_tbl, "27_proliferation_gene_set.tsv")
msg("Hypoxia genes used: ", sum(hyp_tbl$used_in_score), " of ", nrow(hyp_tbl),
    " (missing: ", paste(hyp_tbl[used_in_score == FALSE, symbol], collapse = ", "), ")")
msg("Genes in the hypoxia set correlating with the metric above |rho| 0.5: ",
    paste(hyp_tbl[abs(rho_noFeature) > 0.5,
                  sprintf("%s (rho %.2f, median FPKM %.3f)", symbol, rho_noFeature, median_fpkm)],
          collapse = "; "))
msg("Spearman rho hypoxia vs proliferation score: ",
    round(spearman(hyp$score, pro$score[names(hyp$score)])[["rho"]], 3))

# ---- protein-coding module eigengenes: residualised and observed -----------
# Residualised = the discovery scores. Observed = the same modules with the
# rotation refitted on the observed matrix, as in stages 12 and 23.
LOAD  <- discovery_loadings(nets)
S_res <- discovery_scores(nets, LOAD)$mrna
L_obs <- fit_module_loadings(E_obs, nets$mrna$gene_tbl, "mRNA_ME")
S_obs <- score_modules(E_obs, L_obs)
mods  <- intersect(colnames(S_res), colnames(S_obs))
stopifnot(length(mods) == nets$mrna$n_modules, identical(rownames(S_res), rownames(S_obs)))
cat_mod <- catabolic_reference_module(); cat_col <- paste0("mRNA_ME", cat_mod)
stopifnot(cat_col %in% mods)
# The part of an observed eigengene orthogonal to its residualised counterpart
# is, for a metric-tracking module, close to a copy of the metric.
nf_506 <- nf_of(rownames(S_obs))
eig_cmp <- rbindlist(lapply(mods, function(m) {
  orth <- resid(lm(S_obs[, m] ~ S_res[, m]))
  data.table(
    module = sub("^mRNA_ME", "", m), n_genes = length(LOAD[[m]]$genes), n = nrow(S_obs),
    var_explained_residualised = round(LOAD[[m]]$var_explained, 4),
    var_explained_observed     = round(L_obs[[m]]$var_explained, 4),
    pearson_r_observed_vs_residualised = round(stats::cor(S_obs[, m], S_res[, m]), 3),
    pearson_r_observed_vs_metric       = round(stats::cor(S_obs[, m], nf_506), 3),
    pearson_r_residualised_vs_metric   = round(stats::cor(S_res[, m], nf_506), 3),
    r_orthogonal_part_vs_metric        = round(stats::cor(orth, nf_506), 3),
    is_catabolic_reference = m == cat_col)
}))
save_tsv(eig_cmp, "27_eigengene_observed_vs_residualised.tsv"); print(eig_cmp)

# ---- per-sample table -------------------------------------------------------
ps <- cohort[, .(sample_barcode, patient, tss, plate, os_time, os_event)]
if (all(is.na(ps$plate))) ps[, plate := bio$plate[match(sample_barcode, bio$sample_barcode)]]
ps[, lnc_axis := axis$lnc_axis[match(sample_barcode, axis$sample_barcode)]]
ps[, pc_axis  := unname(pc_mrna$pc[sample_barcode])]
ps[, `:=`(lnc_axis_z = zc(lnc_axis), pc_axis_z = zc(pc_axis))]
ps <- merge(ps, qc[, .(sample_barcode, pct_noFeature, pct_multimapping,
                       log10_assigned_reads = log10(assigned_reads))],
            by = "sample_barcode", all.x = TRUE)
ps <- merge(ps, bio[, .(sample_barcode, rin, batch, pct_tumor_nuclei_mean, pct_stromal_mean,
                        pct_necrosis_mean)], by = "sample_barcode", all.x = TRUE)
ps <- merge(ps, est[, .(sample_barcode, StromalScore, ImmuneScore, ESTIMATEScore)],
            by = "sample_barcode", all.x = TRUE)
ps[, hypoxia_score           := unname(hyp$score[sample_barcode])]
ps[, hypoxia_score_fpkm_ge1  := unname(hyp_hi[sample_barcode])]
ps[, hypoxia_score_TUBA1A    := unname(hyp_alt[sample_barcode])]
ps[, hypoxia_score_both_TUBA := unname(hyp_both[sample_barcode])]
ps[, prolif_score            := unname(pro$score[sample_barcode])]
ps[, lnc_fpkm_share          := unname(lnc_share[sample_barcode])]
for (m in mods) {
  set(ps, j = paste0("obs_", m), value = S_obs[match(ps$sample_barcode, rownames(S_obs)), m])
  set(ps, j = paste0("res_", m), value = S_res[match(ps$sample_barcode, rownames(S_res)), m])
}
setcolorder(ps, c("sample_barcode", "lnc_axis", "pc_axis", "pct_noFeature", "pct_multimapping",
                  "log10_assigned_reads", "rin", "plate", "tss", "lnc_axis_z", "pc_axis_z",
                  "batch", "patient", "os_time", "os_event"))
save_tsv(ps, "27_axis_scores_per_sample.tsv")
msg("Per-sample table: ", nrow(ps), " rows, ", ncol(ps), " columns; lncRNA axis for ",
    sum(!is.na(ps$lnc_axis)), ", protein-coding PC1 for ", sum(!is.na(ps$pc_axis)),
    ", RIN for ", sum(!is.na(ps$rin)), ", plate for ", sum(!is.na(ps$plate) & nzchar(ps$plate)))

# ---- 2. Correlations of the metric and the axis with the biological candidates ----
banner("2 | Spearman correlations with biological candidates")
cands <- rbind(
  data.table(candidate = mods, column = paste0("res_", mods),
             candidate_class = "protein_coding_module_eigengene", matrix = "residualised",
             sensitivity_variant = FALSE),
  data.table(candidate = mods, column = paste0("obs_", mods),
             candidate_class = "protein_coding_module_eigengene", matrix = "observed",
             sensitivity_variant = FALSE),
  data.table(candidate = c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
             column = c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
             candidate_class = "ESTIMATE", matrix = NA_character_, sensitivity_variant = FALSE),
  data.table(candidate = c("slide_tumour_nuclei_pct", "slide_stromal_pct"),
             column = c("pct_tumor_nuclei_mean", "pct_stromal_mean"),
             candidate_class = "slide_histology", matrix = NA_character_, sensitivity_variant = FALSE),
  data.table(candidate = "hypoxia_score_Buffa2010", column = "hypoxia_score",
             candidate_class = "hypoxia_metagene", matrix = "observed", sensitivity_variant = FALSE),
  data.table(candidate = "proliferation_score_10genes", column = "prolif_score",
             candidate_class = "proliferation", matrix = "observed", sensitivity_variant = FALSE),
  data.table(candidate = "lnc_fpkm_share", column = "lnc_fpkm_share",
             candidate_class = "intronic_signal_proxy", matrix = "observed", sensitivity_variant = FALSE),
  # Sensitivity variants: reported but excluded from the FDR families.
  data.table(candidate = "hypoxia_score_Buffa2010_median_FPKM_ge1",
             column = "hypoxia_score_fpkm_ge1", candidate_class = "hypoxia_metagene",
             matrix = "observed", sensitivity_variant = TRUE),
  data.table(candidate = "hypoxia_score_Buffa2010_TUBA1A_mapping",
             column = "hypoxia_score_TUBA1A", candidate_class = "hypoxia_metagene",
             matrix = "observed", sensitivity_variant = TRUE),
  data.table(candidate = "hypoxia_score_Buffa2010_both_tubulins",
             column = "hypoxia_score_both_TUBA", candidate_class = "hypoxia_metagene",
             matrix = "observed", sensitivity_variant = TRUE))
cands[, is_catabolic_reference := candidate == cat_col]
# BH families are split by matrix: each module appears once per matrix, and n
# differs between families. The pooled correction is also reported.
cands[, fdr_family := fifelse(!is.na(matrix) & matrix == "observed" &
                                candidate_class == "protein_coding_module_eigengene",
                              "observed protein-coding module eigengenes (n = 506)",
                     fifelse(!is.na(matrix) & matrix == "residualised",
                              "residualised protein-coding module eigengenes (n = 506)",
                              "non-eigengene candidates (ESTIMATE, slide histology, hypoxia, proliferation, lncRNA FPKM share; n = 524-528)"))]
FAMILY_NOTE <- paste(
  "BH-FDR family = one exposure x one matrix family (observed eigengenes, residualised eigengenes,",
  "non-eigengene candidates); families are split by matrix because each module appears twice, once per",
  "matrix, and because n differs between families (506 for eigengenes, 524-528 otherwise).",
  "fdr_pooled_all_candidates is the single-family correction over all primary candidates of one exposure.",
  "Sensitivity variants are reported but excluded from every family.")
EXPOSURES <- c(pct_noFeature = "pct_noFeature", lnc_axis = "lnc_axis")
cor_tbl <- rbindlist(lapply(names(EXPOSURES), function(ex) {
  out <- rbindlist(lapply(seq_len(nrow(cands)), function(k) {
    s <- spearman(ps[[EXPOSURES[[ex]]]], ps[[cands$column[k]]])
    data.table(exposure = ex, candidate = cands$candidate[k],
               candidate_class = cands$candidate_class[k], matrix = cands$matrix[k],
               sensitivity_variant = cands$sensitivity_variant[k],
               is_catabolic_reference = cands$is_catabolic_reference[k],
               fdr_family = cands$fdr_family[k],
               n = s[["n"]], spearman_rho = round(s[["rho"]], 3), p = signif(s[["p"]], 3))
  }))
  out[sensitivity_variant == FALSE, fdr := signif(p.adjust(p, "BH"), 3), by = fdr_family]
  out[sensitivity_variant == FALSE, n_tests_in_family := .N, by = fdr_family]
  out[sensitivity_variant == FALSE,
      fdr_pooled_all_candidates := signif(p.adjust(p, "BH"), 3)]
  out[sensitivity_variant == FALSE, n_tests_pooled := .N]
  out[sensitivity_variant == TRUE, fdr_family := "sensitivity variant (excluded from every family)"]
  out[order(sensitivity_variant, -abs(spearman_rho))]
}))

# Gene-level profile of each protein-coding module (stages 06 and 07): whether
# the whole module moves with the metric, and its top GO term.
pg <- fread(file.path(RESULTS_DIR, "07_per_gene_quality_correlation.tsv"))[biotype == "protein_coding"]
go_top10 <- fread(file.path(RESULTS_DIR, "06_GO_enrichment_top10_per_module.tsv"))
go_all   <- fread(file.path(RESULTS_DIR, "06_GO_enrichment_all_modules.tsv"))
# Modules missing from the top-ten file fall back to the all-modules file.
go_best <- rbind(
  go_top10[order(p.adjust)][!duplicated(module),
    .(module, top_GO_term = Description, top_GO_padj = signif(p.adjust, 3),
      top_GO_count = Count, top_GO_source = "06_GO_enrichment_top10_per_module.tsv")],
  go_all[order(p.adjust)][!duplicated(module),
    .(module, top_GO_term = Description, top_GO_padj = signif(p.adjust, 3),
      top_GO_count = Count, top_GO_source = "06_GO_enrichment_all_modules.tsv")])[
  !duplicated(module)]
gtab <- as.data.table(nets$mrna$gene_tbl)
hubs <- gtab[module != "grey"][order(module, -kME)][
  , .(top10_kME_genes = paste(head(gene_name, 10), collapse = ";")), by = module]
rho_e <- function(ex, mat) cor_tbl[exposure == ex &
                                     candidate_class == "protein_coding_module_eigengene" &
                                     matrix == mat]
mod_tbl <- rbindlist(lapply(mods, function(m) {
  col <- sub("^mRNA_ME", "", m); g <- pg[module == col]; a <- abs(g$rho_noFeature)
  data.table(module = col, n_genes = nrow(g),
             gene_level_median_rho_noFeature = round(median(g$rho_noFeature), 3),
             gene_level_median_abs_rho = round(median(a), 3),
             gene_level_q25_abs_rho = round(unname(quantile(a, 0.25)), 3),
             gene_level_q75_abs_rho = round(unname(quantile(a, 0.75)), 3),
             gene_level_frac_abs_rho_gt_0.3 = round(mean(a > 0.3), 3),
             gene_level_frac_abs_rho_gt_0.5 = round(mean(a > 0.5), 3),
             gene_level_frac_rho_gt_0.3  = round(mean(g$rho_noFeature > 0.3), 3),
             gene_level_frac_rho_lt_m0.3 = round(mean(g$rho_noFeature < -0.3), 3),
             eigengene_rho_noFeature_observed     = rho_e("pct_noFeature", "observed")[candidate == m, spearman_rho],
             eigengene_rho_noFeature_residualised = rho_e("pct_noFeature", "residualised")[candidate == m, spearman_rho],
             eigengene_rho_lnc_axis_observed      = rho_e("lnc_axis", "observed")[candidate == m, spearman_rho],
             pearson_r_observed_vs_residualised   = eig_cmp[module == col, pearson_r_observed_vs_residualised],
             r_orthogonal_part_vs_metric          = eig_cmp[module == col, r_orthogonal_part_vs_metric],
             is_catabolic_reference = m == cat_col,
             top_GO_term  = go_best[module == col, top_GO_term][1],
             top_GO_padj  = go_best[module == col, top_GO_padj][1],
             top_GO_count = go_best[module == col, top_GO_count][1],
             top_GO_source = go_best[module == col, top_GO_source][1],
             top10_kME_genes = hubs[module == col, top10_kME_genes][1],
             note = paste0("per-gene Spearman rho with pct_noFeature on observed log2(FPKM+1) ",
                           "(07_per_gene_quality_correlation.tsv, n = ", g$n[1], " samples)"))
}))
mod_tbl <- mod_tbl[order(-abs(eigengene_rho_noFeature_observed))]
# Modules without an enriched GO term are labelled explicitly.
mod_tbl[is.na(top_GO_term),
        `:=`(top_GO_source = "none: no enriched term for this module in 06",
             note = paste(note, "No Gene Ontology term is reported for this module in",
                          "06_GO_enrichment_top10_per_module.tsv or 06_GO_enrichment_all_modules.tsv."))]
save_tsv(mod_tbl, "27_module_gene_level_metric_correlation.tsv")
print(mod_tbl[, .(module, n_genes, gene_level_median_rho_noFeature, gene_level_median_abs_rho,
                  gene_level_frac_abs_rho_gt_0.3, eigengene_rho_noFeature_observed,
                  eigengene_rho_noFeature_residualised, top_GO_term)])

# Carry the gene-level profile onto the eigengene correlation rows.
cor_tbl <- merge(cor_tbl,
  mod_tbl[, .(candidate = paste0("mRNA_ME", module),
              module_n_genes = n_genes,
              module_gene_level_median_rho_noFeature = gene_level_median_rho_noFeature,
              module_gene_level_median_abs_rho = gene_level_median_abs_rho,
              module_gene_level_frac_abs_rho_gt_0.3 = gene_level_frac_abs_rho_gt_0.3,
              module_top_GO_term = top_GO_term)],
  by = "candidate", all.x = TRUE, sort = FALSE)
cor_tbl[, interpretable_as_biology := fifelse(
  candidate_class == "protein_coding_module_eigengene" & matrix == "observed", FALSE,
  fifelse(candidate == "lnc_fpkm_share", FALSE, TRUE))]
cor_tbl[, note := FAMILY_NOTE]
cor_tbl[candidate_class == "protein_coding_module_eigengene" & matrix == "observed",
        note := paste("Eigengene computed on the OBSERVED matrix, which carries the metric, so this",
                      "correlation is not identifiable as biology; the module-level gene columns are",
                      "the interpretable part.", FAMILY_NOTE)]
cor_tbl[candidate_class == "protein_coding_module_eigengene" & matrix == "residualised",
        note := paste("Eigengene computed on the technically residualised matrix, from which the",
                      "metric has been removed by construction; a near-zero correlation here is",
                      "expected and is not evidence that the module carries no biology.", FAMILY_NOTE)]
cor_tbl[candidate == "lnc_fpkm_share",
        note := paste("The lncRNA share of total FPKM is a compositional restatement of unassigned",
                      "read content, not an independent biological candidate.", FAMILY_NOTE)]
cor_tbl[sensitivity_variant == TRUE,
        note := paste("Sensitivity variant of the hypoxia metagene; excluded from every FDR family.",
                      FAMILY_NOTE)]
setcolorder(cor_tbl, c("exposure", "candidate", "candidate_class", "matrix", "sensitivity_variant",
                       "is_catabolic_reference", "fdr_family", "n", "spearman_rho", "p", "fdr",
                       "n_tests_in_family", "fdr_pooled_all_candidates", "n_tests_pooled",
                       "interpretable_as_biology"))
cor_tbl <- cor_tbl[order(exposure, sensitivity_variant, -abs(spearman_rho))]
save_tsv(cor_tbl, "27_metric_biology_correlations.tsv")
print(cor_tbl[, .(exposure, candidate, matrix, n, spearman_rho, p, fdr, fdr_pooled_all_candidates)],
      nrows = 200)
msg("Hypoxia vs the metric: full score rho ",
    cor_tbl[exposure == "pct_noFeature" & candidate == "hypoxia_score_Buffa2010", spearman_rho],
    "; median-FPKM >= 1 subset rho ",
    cor_tbl[exposure == "pct_noFeature" & candidate == "hypoxia_score_Buffa2010_median_FPKM_ge1", spearman_rho],
    "; TUBA1A mapping rho ",
    cor_tbl[exposure == "pct_noFeature" & candidate == "hypoxia_score_Buffa2010_TUBA1A_mapping", spearman_rho])

# The five observed eigengenes most correlated with the metric. They are
# selected on the exposure, so their attenuation is an upper bound.
top5 <- cor_tbl[exposure == "pct_noFeature" & candidate_class == "protein_coding_module_eigengene" &
                matrix == "observed"][order(-abs(spearman_rho))][1:5]
top5_mods <- top5$candidate
msg("Five observed protein-coding eigengenes most correlated with the metric (selected ON the exposure): ",
    paste(sprintf("%s (rho %.2f)", sub("^mRNA_ME", "", top5_mods), top5$spearman_rho), collapse = ", "))

# ---- 3. Attenuation of the non-feature hazard by each candidate, and the converse ----
banner("3 | Attenuation of the non-feature hazard by each candidate")
d <- merge(ps, cohort[, .(sample_barcode, age, sex, T_stage, N_pos, M1, grade_num)],
           by = "sample_barcode")
d[, male := as.numeric(sex == "male")]
obs_cols <- paste0("obs_", mods); res_cols <- paste0("res_", mods)
# One complete-case set for every model, so that adding a candidate changes the
# model and not the patients. The lncRNA-axis rows use their own subset.
NEED <- c("os_time", "os_event", "age", "male", "T_stage", "N_pos", "M1", "grade_num",
          "pct_noFeature", "hypoxia_score", "prolif_score", obs_cols, res_cols, "StromalScore",
          "ImmuneScore", "pct_tumor_nuclei_mean", "rin", "plate", "lnc_fpkm_share")
dc <- d[complete.cases(d[, ..NEED]) & !is.na(plate) & nzchar(plate)]
add_z <- function(dat) {
  dat <- copy(dat)
  dat[, `:=`(nf_z = zc(pct_noFeature), nf_log_z = zc(log10(pct_noFeature)),
             axis_z = zc(lnc_axis), hypoxia_z = zc(hypoxia_score), prolif_z = zc(prolif_score),
             hypoxia_alt_z = zc(hypoxia_score_TUBA1A), hypoxia_hi_z = zc(hypoxia_score_fpkm_ge1),
             stromal_z = zc(StromalScore), immune_z = zc(ImmuneScore),
             nuclei_z = zc(pct_tumor_nuclei_mean), rin_z = zc(rin), lncshare_z = zc(lnc_fpkm_share),
             plate_f = pool_levels(plate), tss_f = pool_levels(tss))]
  for (m in mods) {
    set(dat, j = paste0("z_", m),    value = zc(dat[[paste0("obs_", m)]]))
    set(dat, j = paste0("zres_", m), value = zc(dat[[paste0("res_", m)]]))
  }
  dat
}
dc    <- add_z(dc)
dc_ax <- add_z(dc[!is.na(lnc_axis)])
msg("Common complete-case set: ", nrow(dc), " patients, ", sum(dc$os_event), " events (",
    nlevels(dc$plate_f), " plate strata); with the lncRNA axis: ", nrow(dc_ax), " patients, ",
    sum(dc_ax$os_event), " events")

CLIN      <- "age + male + T_stage + N_pos + M1 + grade_num"
CLIN_V    <- c("age", "male", "T_stage", "N_pos", "M1", "grade_num")
top5_z    <- paste0("z_", top5_mods)
top5_zres <- paste0("zres_", top5_mods)
all_z     <- paste0("z_", mods)
all_zres  <- paste0("zres_", mods)
top5_lab  <- sub("^mRNA_ME", "", top5_mods)

# Each attenuation model carries its own interpretation flags.
AM <- function(label, terms, matrix, interpretable, selected_on_exposure, note)
  list(label = label, terms = terms, matrix = matrix, interpretable = interpretable,
       selected_on_exposure = selected_on_exposure, note = note)
N_OBS <- paste("Candidate computed on the OBSERVED protein-coding matrix, which carries the metric.",
               "Any apparent attenuation here is collinearity, not mediation: see",
               "27_collinearity_diagnostics.tsv for the correlation with the exposure, the variance",
               "inflation, the standard-error inflation and the synthetic-proxy null.")
N_RES <- paste("Candidate computed on the technically residualised matrix, from which the metric has",
               "been removed by construction; this row is the interpretable counterpart of the",
               "observed-eigengene row for the same module.")
att_models <- c(
  list(AM("base: exposure + age, sex, T, N, M1, ordinal grade", character(0), NA_character_,
          TRUE, FALSE, "Reference model for every row of this exposure and scale."),
       AM("+ hypoxia metagene (Buffa 2010)", "hypoxia_z", "observed", TRUE, FALSE,
          paste("A-priori annotated gene set, not selected on the exposure. Computed on observed",
                "expression, so see the median-FPKM >= 1 sensitivity score in",
                "27_metric_biology_correlations.tsv.")),
       AM("+ proliferation score (10 genes)", "prolif_z", "observed", TRUE, FALSE,
          "A-priori annotated gene set, not selected on the exposure.")),
  setNames(lapply(seq_along(top5_z), function(i)
    AM(sprintf("+ %s observed eigengene", top5_lab[i]), top5_z[i], "observed", FALSE, TRUE,
       paste(N_OBS, "Module selected on its correlation with the exposure."))),
    sprintf("+ %s observed eigengene", top5_lab)),
  setNames(lapply(seq_along(top5_zres), function(i)
    AM(sprintf("+ %s residualised eigengene", top5_lab[i]), top5_zres[i], "residualised", TRUE, TRUE,
       paste(N_RES, "Module selected on its correlation with the exposure on the observed matrix."))),
    sprintf("+ %s residualised eigengene", top5_lab)),
  list(AM("+ five observed eigengenes most correlated with the metric, jointly", top5_z,
          "observed", FALSE, TRUE,
          paste(N_OBS, "The block is SELECTED ON THE EXPOSURE and is therefore guaranteed to be the",
                "block that attenuates most; it is an upper bound, not an unbiased estimate.")),
       AM("+ the same five modules' residualised eigengenes, jointly", top5_zres,
          "residualised", TRUE, TRUE, N_RES),
       AM("+ all 13 observed protein-coding eigengenes, jointly", all_z, "observed", FALSE, FALSE,
          paste(N_OBS, "Not selected on the exposure, but still computed on the observed matrix.")),
       AM("+ all 13 residualised protein-coding eigengenes, jointly", all_zres,
          "residualised", TRUE, FALSE, N_RES),
       AM("+ ESTIMATE stromal and immune scores", c("stromal_z", "immune_z"), "observed", TRUE, FALSE,
          "Composition scores from expression, not selected on the exposure."),
       AM("+ slide tumour nuclei (%)", "nuclei_z", NA_character_, TRUE, FALSE,
          "Independent measurement (pathology slide review, 23)."),
       AM("+ RIN", "rin_z", NA_character_, TRUE, FALSE,
          "Independent measurement (RNA integrity number of the sequenced analyte, 23)."),
       AM("+ plate stratum", "strata(plate_f)", NA_character_, TRUE, FALSE,
          "Independent measurement (sequencing plate of the RNA aliquot, 23)."),
       AM("+ lncRNA share of total FPKM", "lncshare_z", "observed", FALSE, FALSE,
          paste("The lncRNA share of total FPKM is a compositional restatement of unassigned read",
                "content (Spearman rho about 0.73 with the metric), not an independent candidate.")),
       AM(paste("+ all candidates jointly, observed eigengenes (hypoxia, proliferation, five",
                "observed eigengenes, ESTIMATE, tumour nuclei, RIN, plate stratum)"),
          c("hypoxia_z", "prolif_z", top5_z, "stromal_z", "immune_z", "nuclei_z", "rin_z",
            "strata(plate_f)"), "observed", FALSE, TRUE, N_OBS),
       AM(paste("+ all candidates jointly, residualised eigengenes (hypoxia, proliferation, five",
                "residualised eigengenes, ESTIMATE, tumour nuclei, RIN, plate stratum)"),
          c("hypoxia_z", "prolif_z", top5_zres, "stromal_z", "immune_z", "nuclei_z", "rin_z",
            "strata(plate_f)"), "residualised", TRUE, TRUE, N_RES)))
names(att_models) <- vapply(att_models, function(a) a$label, character(1))
cox_fit <- function(rhs, dat) coxph(as.formula(paste("Surv(os_time, os_event) ~", rhs)), data = dat)

fit_attenuation <- function(term, exposure_label, scale_label, dat, sd_raw) {
  base_rhs <- paste(term, "+", CLIN)
  f0 <- cox_fit(base_rhs, dat); b0 <- coef(f0)[[term]]
  se0 <- sqrt(diag(vcov(f0)))[[term]]
  rbindlist(lapply(att_models, function(a) {
    extra <- a$terms
    rhs   <- if (length(extra)) paste(base_rhs, "+", paste(extra, collapse = " + ")) else base_rhs
    fit   <- cox_fit(rhs, dat); s <- summary(fit); b <- coef(fit)[[term]]
    se    <- sqrt(diag(vcov(fit)))[[term]]
    cov_terms  <- extra[!grepl("^strata", extra)]
    has_strata <- length(cov_terms) < length(extra)
    single     <- length(cov_terms) == 1 && !has_strata
    # LRT of the candidate block against the base model (undefined if strata change).
    lrt <- if (length(cov_terms) && !has_strata) 2 * (fit$loglik[2] - f0$loglik[2]) else NA_real_
    data.table(exposure = exposure_label, exposure_scale = scale_label,
               exposure_sd_raw = round(sd_raw, 4), model = a$label,
               candidate_matrix = a$matrix, interpretable = a$interpretable,
               selected_on_exposure = a$selected_on_exposure,
               candidate_terms = paste(extra, collapse = " + "), n = s$n, events = s$nevent,
               HR = round(exp(b), 3), lo = round(s$conf.int[term, 3], 3),
               hi = round(s$conf.int[term, 4], 3), p = signif(s$coefficients[term, 5], 3),
               se_logHR = round(se, 4), se_inflation_vs_base = round(se / se0, 3),
               HR_base = round(exp(b0), 3), delta_HR = round(exp(b) - exp(b0), 3),
               delta_logHR = round(b - b0, 4),
               HR_ratio_vs_base = round(exp(b - b0), 3),
               # Percentage change is unstable near zero, so it is suppressed below |b0| = 0.05.
               pct_change_logHR = if (abs(b0) >= 0.05) round(100 * (b - b0) / b0, 1) else NA_real_,
               C = round(unname(s$concordance[1]), 3),
               ph_p = signif(tryCatch(cox.zph(fit)$table[term, "p"], error = function(e) NA_real_), 3),
               candidate_HR = if (single) round(s$conf.int[cov_terms, 1], 3) else NA_real_,
               candidate_lo = if (single) round(s$conf.int[cov_terms, 3], 3) else NA_real_,
               candidate_hi = if (single) round(s$conf.int[cov_terms, 4], 3) else NA_real_,
               candidate_p  = if (single) signif(s$coefficients[cov_terms, 5], 3) else NA_real_,
               candidate_block_df = length(cov_terms),
               candidate_block_lrt_p = if (is.finite(lrt))
                 signif(pchisq(lrt, df = length(cov_terms), lower.tail = FALSE), 3) else NA_real_,
               note = a$note)
  }))
}
att <- rbind(
  fit_attenuation("nf_z",     "pct_noFeature", "linear, per SD", dc,    sd(dc$pct_noFeature)),
  fit_attenuation("nf_log_z", "pct_noFeature", "log10, per SD",  dc,    sd(log10(dc$pct_noFeature))),
  fit_attenuation("axis_z",   "lnc_axis",      "per SD",         dc_ax, sd(dc_ax$lnc_axis)))
save_tsv(att, "27_metric_attenuation.tsv")
print(att[, .(exposure, exposure_scale, model, candidate_matrix, interpretable, n, events,
              HR, lo, hi, p, HR_base, se_inflation_vs_base)], nrows = 100)

# ---- collinearity diagnostics ----------------------------------------------
# For each single candidate: correlation with the exposure, VIF, SE inflation,
# and the exposure HR when the candidate is replaced by a synthetic covariate
# with the same correlation and no survival information (collinearity null).
banner("3b | Collinearity diagnostics for the attenuation table")
sim_base <- as.data.frame(dc[, c("os_time", "os_event", CLIN_V, "nf_z"), with = FALSE])
simulate_proxy <- function(target_r, B = SIM_DRAWS, seed = SEED) {
  set.seed(seed)
  x  <- (sim_base$nf_z - mean(sim_base$nf_z)) / sd(sim_base$nf_z)
  n  <- length(x); r <- min(abs(target_r), 0.999) * sign(target_r)
  fm <- as.formula(paste("Surv(os_time, os_event) ~ nf_z +", CLIN, "+ proxy_z"))
  hr <- vapply(seq_len(B), function(b) {
    e <- rnorm(n); e <- resid(lm(e ~ x)); e <- e / sd(e)
    w <- r * x + sqrt(1 - r^2) * e
    dd <- sim_base; dd$proxy_z <- (w - mean(w)) / sd(w)
    exp(unname(coef(coxph(fm, data = dd))[["nf_z"]]))
  }, numeric(1))
  c(median = median(hr), lo = unname(quantile(hr, 0.025)), hi = unname(quantile(hr, 0.975)))
}
coll_cands <- rbind(
  data.table(candidate = "hypoxia_score_Buffa2010", term = "hypoxia_z", candidate_matrix = "observed"),
  data.table(candidate = "proliferation_score_10genes", term = "prolif_z", candidate_matrix = "observed"),
  data.table(candidate = paste0(top5_lab, " observed eigengene"), term = top5_z,
             candidate_matrix = "observed"),
  data.table(candidate = paste0(top5_lab, " residualised eigengene"), term = top5_zres,
             candidate_matrix = "residualised"),
  data.table(candidate = "StromalScore", term = "stromal_z", candidate_matrix = NA_character_),
  data.table(candidate = "ImmuneScore", term = "immune_z", candidate_matrix = NA_character_),
  data.table(candidate = "slide_tumour_nuclei_pct", term = "nuclei_z", candidate_matrix = NA_character_),
  data.table(candidate = "RIN", term = "rin_z", candidate_matrix = NA_character_),
  data.table(candidate = "lnc_fpkm_share", term = "lncshare_z", candidate_matrix = "observed"))
f_base   <- cox_fit(paste("nf_z +", CLIN), dc)
b_base   <- coef(f_base)[["nf_z"]]; se_base <- sqrt(diag(vcov(f_base)))[["nf_z"]]
SIM_MIN_R <- 0.3   # simulate the null only where collinearity is material
coll <- rbindlist(lapply(seq_len(nrow(coll_cands)), function(k) {
  tm <- coll_cands$term[k]
  f1 <- cox_fit(paste("nf_z +", CLIN, "+", tm), dc)
  b1 <- coef(f1)[["nf_z"]]; se1 <- sqrt(diag(vcov(f1)))[["nf_z"]]
  r  <- stats::cor(dc$nf_z, dc[[tm]]); rs <- spearman(dc$nf_z, dc[[tm]])[["rho"]]
  # Correlation of the part orthogonal to the residualised eigengene.
  orth_r <- NA_real_
  if (grepl("^z_mRNA_ME", tm)) {
    rt <- sub("^z_", "zres_", tm)
    orth_r <- stats::cor(resid(lm(dc[[tm]] ~ dc[[rt]])), dc$nf_z)
  }
  sim <- if (abs(r) >= SIM_MIN_R) simulate_proxy(r) else c(median = NA, lo = NA, hi = NA)
  data.table(exposure = "pct_noFeature", exposure_scale = "linear, per SD",
             candidate = coll_cands$candidate[k], candidate_matrix = coll_cands$candidate_matrix[k],
             term = tm, n = nrow(dc), events = sum(dc$os_event),
             pearson_r_with_exposure = round(r, 3), spearman_rho_with_exposure = round(rs, 3),
             vif_exposure_full_design = round(vif_of("nf_z", c(CLIN_V, tm), dc), 3),
             vif_pairwise = round(1 / (1 - r^2), 3),
             r_orthogonal_part_with_exposure = round(orth_r, 3),
             HR_exposure_base = round(exp(b_base), 3), HR_exposure_with_candidate = round(exp(b1), 3),
             se_logHR_base = round(se_base, 4), se_logHR_with_candidate = round(se1, 4),
             se_inflation_ratio = round(se1 / se_base, 3),
             sim_n_draws = if (is.na(sim[["median"]])) NA_integer_ else SIM_DRAWS,
             sim_target_r = if (is.na(sim[["median"]])) NA_real_ else round(r, 3),
             sim_median_HR = round(unname(sim[["median"]]), 3),
             sim_lo_2.5 = round(unname(sim[["lo"]]), 3),
             sim_hi_97.5 = round(unname(sim[["hi"]]), 3),
             sim_covers_observed_HR = if (is.na(sim[["median"]])) NA
                                      else exp(b1) >= sim[["lo"]] & exp(b1) <= sim[["hi"]],
             interpretable = !(coll_cands$candidate_matrix[k] %in% "observed" &
                                 grepl("eigengene|lnc_fpkm_share", coll_cands$candidate[k])),
             note = if (abs(r) >= SIM_MIN_R)
               paste0("Synthetic null: ", SIM_DRAWS, " draws of a covariate with Pearson r ",
                      round(r, 3), " with the exposure and no survival information, added to the ",
                      "base model; the interval is the 2.5 to 97.5 per cent range of the exposure's ",
                      "hazard ratio under that null. If the observed attenuated hazard ratio lies ",
                      "inside it, the attenuation is explained by collinearity alone.")
             else paste0("Correlation with the exposure below ", SIM_MIN_R,
                         "; the synthetic-proxy null was not simulated."))
}))
save_tsv(coll, "27_collinearity_diagnostics.tsv")
print(coll[, .(candidate, candidate_matrix, pearson_r_with_exposure, vif_exposure_full_design,
                HR_exposure_with_candidate, se_inflation_ratio, sim_median_HR, sim_lo_2.5,
                sim_hi_97.5, sim_covers_observed_HR)], nrows = 50)

# ---- the converse: each candidate with and without the metric ---------------
banner("4 | Each candidate's hazard with and without the non-feature fraction")
conv <- rbind(
  data.table(candidate = "hypoxia_score_Buffa2010", term = "hypoxia_z",
             candidate_class = "hypoxia_metagene", candidate_matrix = "observed", interpretable = TRUE),
  data.table(candidate = "proliferation_score_10genes", term = "prolif_z",
             candidate_class = "proliferation", candidate_matrix = "observed", interpretable = TRUE),
  data.table(candidate = mods, term = all_z,
             candidate_class = "protein_coding_module_eigengene",
             candidate_matrix = "observed", interpretable = FALSE),
  data.table(candidate = mods, term = all_zres,
             candidate_class = "protein_coding_module_eigengene",
             candidate_matrix = "residualised", interpretable = TRUE),
  data.table(candidate = c("StromalScore", "ImmuneScore"), term = c("stromal_z", "immune_z"),
             candidate_class = "ESTIMATE", candidate_matrix = NA_character_, interpretable = TRUE),
  data.table(candidate = "slide_tumour_nuclei_pct", term = "nuclei_z",
             candidate_class = "slide_histology", candidate_matrix = NA_character_, interpretable = TRUE),
  data.table(candidate = "RIN", term = "rin_z", candidate_class = "RNA_integrity",
             candidate_matrix = NA_character_, interpretable = TRUE),
  data.table(candidate = "lnc_fpkm_share", term = "lncshare_z",
             candidate_class = "intronic_signal_proxy", candidate_matrix = "observed",
             interpretable = FALSE))
conv[, top5_by_metric_correlation := term %in% c(top5_z, top5_zres)]
conv[, is_catabolic_reference := candidate == cat_col]
conv_tbl <- rbindlist(lapply(seq_len(nrow(conv)), function(k) {
  tm <- conv$term[k]
  f1 <- cox_fit(paste(tm, "+", CLIN), dc);           s1 <- summary(f1); b1 <- coef(f1)[[tm]]
  f2 <- cox_fit(paste(tm, "+", CLIN, "+ nf_z"), dc); s2 <- summary(f2); b2 <- coef(f2)[[tm]]
  data.table(candidate = conv$candidate[k], candidate_class = conv$candidate_class[k],
             candidate_matrix = conv$candidate_matrix[k], interpretable = conv$interpretable[k],
             top5_by_metric_correlation = conv$top5_by_metric_correlation[k],
             is_catabolic_reference = conv$is_catabolic_reference[k],
             n = s1$n, events = s1$nevent,
             HR_without_metric = round(exp(b1), 3), lo_without = round(s1$conf.int[tm, 3], 3),
             hi_without = round(s1$conf.int[tm, 4], 3), p_without = signif(s1$coefficients[tm, 5], 3),
             HR_with_metric = round(exp(b2), 3), lo_with = round(s2$conf.int[tm, 3], 3),
             hi_with = round(s2$conf.int[tm, 4], 3), p_with = signif(s2$coefficients[tm, 5], 3),
             # Percentage change only where |base log HR| >= 0.05.
             delta_logHR = round(b2 - b1, 3),
             HR_ratio_with_over_without = round(exp(b2 - b1), 3),
             pct_change_logHR = if (abs(b1) >= 0.05) round(100 * (b2 - b1) / b1, 1) else NA_real_,
             base_logHR = round(b1, 4),
             metric_HR_in_joint_model = round(s2$conf.int["nf_z", 1], 3),
             metric_p_in_joint_model  = signif(s2$coefficients["nf_z", 5], 3),
             note = if (identical(conv$candidate_matrix[k], "observed") &&
                        conv$candidate_class[k] == "protein_coding_module_eigengene")
                      N_OBS
                    else if (identical(conv$candidate_matrix[k], "residualised")) N_RES
                    else if (conv$candidate[k] == "lnc_fpkm_share")
                      "Compositional restatement of unassigned read content, not an independent candidate."
                    else "Candidate measured independently of, or a priori with respect to, the exposure.")
}))
conv_tbl[, `:=`(fdr_without = signif(p.adjust(p_without, "BH"), 3),
                fdr_with    = signif(p.adjust(p_with, "BH"), 3))]
save_tsv(conv_tbl, "27_candidate_hr_with_without_metric.tsv")
print(conv_tbl[, .(candidate, candidate_matrix, n, events, HR_without_metric, p_without,
                   HR_with_metric, p_with, delta_logHR, HR_ratio_with_over_without,
                   metric_HR_in_joint_model)], nrows = 60)

# ---- 5. Reverse check: the non-feature fraction by clinical and technical groups ----
banner("5 | Reverse check: the non-feature fraction across groups")
# Computed on all discovery patients, not the complete-case set.
REV_NOTE <- paste0("Computed on all ", nrow(d), " discovery patients with the metric; the attenuation ",
                   "and variance-partition tables use the ", nrow(dc), "-patient complete-case set.")
nf07 <- fread(file.path(RESULTS_DIR, "07_noFeature_vs_clinical.tsv"))
reuse <- nf07[, .(variable, test, n,
                  n_groups  = fifelse(test == "Spearman", NA_integer_, 2L),
                  statistic = fifelse(test == "Spearman", as.numeric(rho), NA_real_),
                  p,
                  group_medians = fifelse(test == "Spearman", NA_character_,
                    sprintf("%s=%.2f (n=%d); %s=%.2f (n=%d)", group1, median_group1, n_group1,
                            group2, median_group2, n_group2)),
                  excluded_levels = NA_character_,
                  source = "07_noFeature_vs_clinical.tsv (re-used, not recomputed)",
                  note = REV_NOTE)]
# Levels smaller than drop_n are dropped: a Kruskal-Wallis group of one carries
# no information.
kw_row <- function(variable, grp, dat, drop_n = 3L) {
  ok <- !is.na(grp) & is.finite(dat$pct_noFeature)
  g  <- droplevels(factor(grp[ok])); y <- dat$pct_noFeature[ok]
  tb <- table(g); small <- names(tb)[tb < drop_n]
  if (length(small)) {
    keep <- !(as.character(g) %in% small)
    g <- droplevels(g[keep]); y <- y[keep]
  }
  k  <- kruskal.test(y ~ g); med <- tapply(y, g, median); ns <- table(g)
  data.table(variable = variable, test = "Kruskal-Wallis", n = length(y), n_groups = nlevels(g),
             statistic = round(unname(k$statistic), 2), p = signif(k$p.value, 3),
             group_medians = paste(sprintf("%s=%.2f (n=%d)", names(med), med,
                                           as.integer(ns[names(med)])), collapse = "; "),
             excluded_levels = if (length(small))
               paste0(paste(sprintf("%s (n=%d)", small, as.integer(table(droplevels(factor(grp[ok])))[small])),
                            collapse = "; "), " dropped: fewer than ", drop_n, " samples")
               else NA_character_,
             source = "27_metric_biology.R (computed)", note = REV_NOTE)
}
rev_tbl <- rbind(reuse, rbindlist(list(
  kw_row("sex", d$sex, d),
  kw_row("age_tertile", tertile(d$age), d),
  kw_row("RIN_tertile", tertile(d$rin), d),
  kw_row("plate (levels with >= 10 samples, others pooled)", pool_levels(d$plate), d),
  kw_row("batch (levels with >= 10 samples, others pooled)", pool_levels(d$batch), d))))
save_tsv(rev_tbl, "27_noFeature_reverse_check.tsv")
print(rev_tbl[, .(variable, test, n, n_groups, statistic, p, excluded_levels, source)])

# ---- 6. Variance partition with negative controls ----
banner("6 | Variance of the non-feature fraction: what each block can support")
# Any score computed on the observed matrix partly restates the exposure.
# Negative controls (random gene sets of the same sizes, leading PCs) show the
# R2 such blocks reach without biology. Only the annotated-biology block is a
# biological share, and only the technical block a technical share.
annot_terms <- c("hypoxia_z", "prolif_z", "stromal_z", "immune_z", "nuclei_z")
tech_terms  <- c("rin_z", "plate_f", "tss_f")
top1_mod    <- top5_mods[1]; top1_lab <- sub("^mRNA_ME", "", top1_mod)

# Negative controls: the first draw gives the table row, all RANDOM_SET_DRAWS
# draws give the reported range.
mod_sizes <- vapply(mods, function(m) length(L_obs[[m]]$genes), integer(1))
universe  <- colnames(E_obs)
stopifnot(sum(mod_sizes) <= length(universe))
msg("Negative controls: ", length(mod_sizes), " random gene sets of sizes ",
    paste(range(mod_sizes), collapse = "-"), " drawn from ", length(universe),
    " protein-coding genes, ", RANDOM_SET_DRAWS, " draws per matrix")
set.seed(SEED)
rnd_draws_obs <- lapply(seq_len(RANDOM_SET_DRAWS), function(b) {
  sets <- random_sets(mod_sizes, universe)
  M <- vapply(sets, function(g) set_eigengene(E_obs, g), numeric(nrow(E_obs)))
  dimnames(M) <- list(rownames(E_obs), paste0("rndobs", seq_along(sets))); M
})
rnd_draws_res <- lapply(seq_len(RANDOM_SET_DRAWS), function(b) {
  sets <- random_sets(mod_sizes, universe)
  M <- vapply(sets, function(g) set_eigengene(E_res, g), numeric(nrow(E_res)))
  dimnames(M) <- list(rownames(E_res), paste0("rndres", seq_along(sets))); M
})
PC_obs <- leading_pcs(E_obs, length(mods))
colnames(PC_obs) <- paste0("pcobs", seq_len(ncol(PC_obs)))
msg("Negative controls built in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1),
    " min of elapsed time so far")

# Attach the controls to the complete-case set (defined for every network
# sample, so the set is unchanged).
attach_cols <- function(dat, M) {
  i <- match(dat$sample_barcode, rownames(M))
  for (cn in colnames(M)) set(dat, j = cn, value = zc(M[i, cn]))
  dat
}
dc <- attach_cols(dc, rnd_draws_obs[[1]])
dc <- attach_cols(dc, rnd_draws_res[[1]])
dc <- attach_cols(dc, PC_obs)
stopifnot(nrow(dc) == sum(complete.cases(dc[, c(colnames(rnd_draws_obs[[1]]),
                                                colnames(PC_obs)), with = FALSE])))
rnd_obs_terms <- colnames(rnd_draws_obs[[1]]); rnd_res_terms <- colnames(rnd_draws_res[[1]])
pc_terms      <- colnames(PC_obs)

# Blocks. `interpretable_as_biology` is TRUE only where no predictor derives
# from the observed expression matrix.
BLK <- function(label, terms, kind, interpretable, note)
  list(label = label, terms = terms, kind = kind, interpretable = interpretable, note = note)
N_CIRC <- paste("Contains scores computed on the OBSERVED protein-coding matrix, which carries the",
                "metric; regressing the metric on quantities that contain it explains most of it by",
                "construction. Upper bound on the observed expression space, NOT a biological share.")
blocks <- list(
  BLK("annotated biology only (hypoxia, proliferation, ESTIMATE stromal and immune, slide tumour nuclei)",
      annot_terms, "biological", TRUE,
      paste("The interpretable biological estimate: a-priori annotated quantities, none of them",
            "selected on the exposure. This is the number to quote as the biological share.")),
  BLK("sensitivity: annotated biology, hypoxia under the alternative TUBA1A mapping of probe 212639_x_at",
      c("hypoxia_alt_z", annot_terms[-1]), "biological_sensitivity", TRUE,
      paste("As the annotated-biology block, but with the MSigDB mapping of source probe 212639_x_at",
            "(TUBA1A) in place of the Affymetrix mapping (TUBA1B) used in the primary score. One gene",
            "of 51; the two hypoxia scores correlate at r 0.998.")),
  BLK("sensitivity: annotated biology, hypoxia restricted to genes at median FPKM >= 1",
      c("hypoxia_hi_z", annot_terms[-1]), "biological_sensitivity", TRUE,
      paste("As the annotated-biology block, but with the hypoxia score recomputed on the 45 genes at",
            "median FPKM >= 1, so that a near-zero gene cannot act as a proxy for the metric inside a",
            "biological candidate. The two hypoxia scores correlate at r 0.99.")),
  BLK("observed expression space, all (annotated biology + 13 observed protein-coding eigengenes) [UPPER BOUND, NOT IDENTIFIABLE AS BIOLOGY]",
      c(annot_terms, all_z), "observed_expression_space", FALSE,
      paste(N_CIRC, "It is not a biological share.")),
  BLK("13 observed protein-coding eigengenes alone", all_z, "observed_expression_space", FALSE, N_CIRC),
  BLK(sprintf("observed %s eigengene alone (the module most correlated with the metric)", top1_lab),
      paste0("z_", top1_mod), "observed_expression_space", FALSE,
      paste(N_CIRC, "One observed eigengene reproduces most of the 18-degree-of-freedom block.")),
  BLK("negative control: 13 eigengenes of RANDOM gene sets of the same sizes, observed matrix",
      rnd_obs_terms, "negative_control", FALSE,
      paste0(N_CIRC, " Negative control: the sets carry no biological grouping at all. Row = the ",
             "single draw under set.seed(SEED); RANDOM_SET_DRAWS = ", RANDOM_SET_DRAWS,
             " draws in total, mean and range given in random_control_summary.")),
  BLK("negative control: 13 leading principal components of the observed matrix",
      pc_terms, "negative_control", FALSE,
      paste(N_CIRC, "Negative control: 13 unsupervised components of the same matrix, no biology,",
            "and a HIGHER R2 than the 'biological' block.")),
  BLK("13 residualised protein-coding eigengenes", all_zres, "residualised_expression_space", FALSE,
      paste("The same 13 modules on the technically residualised matrix, from which the metric has",
            "been removed by construction; the near-zero R2 is guaranteed by that construction and",
            "is not evidence that the modules carry no biology.")),
  BLK("negative control: 13 random gene-set eigengenes, residualised matrix",
      rnd_res_terms, "negative_control", FALSE,
      paste0("Companion to the observed random-set control on the residualised matrix. Row = the ",
             "single draw under set.seed(SEED); ", RANDOM_SET_DRAWS, " draws in total.")),
  BLK(sprintf("annotated biology + observed eigengenes excluding the %s eigengene", top1_lab),
      c(annot_terms, paste0("z_", setdiff(mods, top1_mod))), "observed_expression_space", FALSE,
      paste(N_CIRC, "Retained to show that dropping the single most metric-correlated module does",
            "not rescue the block.")),
  BLK("technical (RIN, plate, site)", tech_terms, "technical", FALSE,
      paste("Measured technical variables, none derived from expression: the interpretable TECHNICAL",
            "share. Not a biological block, hence interpretable_as_biology = FALSE.")),
  BLK("annotated biology + technical (the two interpretable blocks)",
      c(annot_terms, tech_terms), "combined", FALSE,
      "Union of the two blocks that can be interpreted; not itself a biological share."),
  BLK("observed expression space + technical (all terms)", c(annot_terms, all_z, tech_terms),
      "combined", FALSE, N_CIRC))

partition <- function(yexpr, label) {
  f    <- function(terms) lm(as.formula(paste(yexpr, "~", paste(terms, collapse = " + "))), data = dc)
  r2   <- function(m) summary(m)$r.squared
  ar2  <- function(m) summary(m)$adj.r.squared
  pR2  <- function(red, full) (deviance(red) - deviance(full)) / deviance(red)
  fp   <- function(m) { fs <- summary(m)$fstatistic; unname(pf(fs[1], fs[2], fs[3], lower.tail = FALSE)) }
  m_t  <- f(tech_terms); m_annot <- f(annot_terms)
  # Random-set control R2 across all draws.
  rnd_summary <- function(draws) vapply(seq_along(draws), function(b) {
    dd <- copy(dc); M <- draws[[b]]; i <- match(dd$sample_barcode, rownames(M))
    for (cn in colnames(M)) set(dd, j = cn, value = zc(M[i, cn]))
    summary(lm(as.formula(paste(yexpr, "~", paste(colnames(M), collapse = " + "))), data = dd))$r.squared
  }, numeric(1))
  rnd_v <- list(obs = rnd_summary(rnd_draws_obs), res = rnd_summary(rnd_draws_res))
  rnd_txt <- function(v) sprintf(" Mean R2 over %d draws %.4f, range %.4f to %.4f.",
                                 length(v), mean(v), min(v), max(v))
  rbindlist(lapply(blocks, function(bk) {
    m_b <- f(bk$terms)
    # Partial R2: non-technical blocks given the technical block, and the
    # technical block given annotated biology.
    if (bk$kind == "technical") {
      cond_lab <- "annotated biology only"
      m_red <- m_annot; m_full <- f(c(annot_terms, tech_terms))
    } else if (bk$kind == "combined") {
      cond_lab <- NA_character_; m_red <- NULL; m_full <- NULL
    } else {
      cond_lab <- "technical (RIN, plate, site)"
      m_red <- m_t; m_full <- f(c(bk$terms, tech_terms))
    }
    note <- bk$note
    v <- if (identical(bk$terms, rnd_obs_terms)) rnd_v$obs
         else if (identical(bk$terms, rnd_res_terms)) rnd_v$res else NULL
    if (!is.null(v)) note <- paste0(note, rnd_txt(v))
    data.table(outcome = yexpr, outcome_scale = label, n = nrow(dc),
               block = bk$label, block_kind = bk$kind,
               interpretable_as_biology = bk$interpretable,
               df = m_b$rank - 1L, R2 = round(r2(m_b), 4), adj_R2 = round(ar2(m_b), 4),
               n_control_draws  = if (is.null(v)) NA_integer_ else length(v),
               R2_mean_over_draws = if (is.null(v)) NA_real_ else round(mean(v), 4),
               R2_min_over_draws  = if (is.null(v)) NA_real_ else round(min(v), 4),
               R2_max_over_draws  = if (is.null(v)) NA_real_ else round(max(v), 4),
               partial_R2_conditioned_on = cond_lab,
               partial_R2_given_other_block = if (is.null(m_red)) NA_real_ else round(pR2(m_red, m_full), 4),
               p = signif(fp(m_b), 3),
               p_partial = if (is.null(m_red)) NA_real_ else signif(anova(m_red, m_full)[["Pr(>F)"]][2], 3),
               note = note)
  }))
}
vp <- rbind(partition("pct_noFeature", "linear"), partition("log10(pct_noFeature)", "log10"))
save_tsv(vp, "27_noFeature_variance_partition.tsv")
print(vp[, .(outcome_scale, block, df, R2, adj_R2, interpretable_as_biology,
             partial_R2_given_other_block)], nrows = 40)
msg("Interpretable biological share (annotated biology, linear scale): R2 ",
    vp[outcome_scale == "linear" & block_kind == "biological", R2],
    "; interpretable technical share: R2 ",
    vp[outcome_scale == "linear" & block_kind == "technical", R2],
    "; observed-expression-space upper bound: R2 ",
    vp[outcome_scale == "linear" & grepl("^observed expression space, all", block), R2])

msg("Elapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
write_session_info("27")
write_session_info("27_metric_biology")
banner("27 | done")
