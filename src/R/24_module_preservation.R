# 24_module_preservation.R: preservation of the TCGA-KIRC modules in external cohorts.
#
# WGCNA modulePreservation() (Langfelder et al. 2011) with the residualised
# TCGA-KIRC network as reference. Zsummary < 2 means no preservation, 2 to 10
# weak to moderate, > 10 strong. medianRank does not favour large modules.
# Test cohorts: CPTAC-3, TCGA-KIRP and TCGA-KICH, as log2(FPKM + 1) residualised
# within cohort on their own STAR metrics, as in 09 and 11. Run after 08 and 11.
# Inputs: caches networks.rds, validation_dataset.rds, subtype_KIRP.rds, subtype_KICH.rds.
# Outputs: results 24_module_preservation_{lncRNA,mRNA,summary}.tsv, figure
#   24_module_preservation_Zsummary, cached modulePreservation() objects.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(WGCNA); library(ggplot2)
})
banner("24 | Module preservation in CPTAC-3, TCGA-KIRP and TCGA-KICH")
set.seed(SEED)
enableWGCNAThreads(N_THREADS)

PRESERVATION_N_PERM    <- 100    # permutations, lncRNA network (primary)
MRNA_N_PERM            <- 30     # fewer: its largest modules exceed 2,000 genes
PRESERVATION_MAX_MODULE <- 1000   # larger modules are sub-sampled by WGCNA
RUN_MRNA_PRESERVATION  <- TRUE

nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
val  <- readRDS(file.path(CACHE_DIR, "validation_dataset.rds"))
vco  <- as.data.table(val$cohort)
# STAR metrics from validation_dataset.rds, else valid_star_qc.rds (same table).
vqc <- if (!is.null(val$qc)) as.data.table(val$qc) else {
  msg("validation_dataset.rds carries no qc table; reading valid_star_qc.rds")
  as.data.table(readRDS(file.path(CACHE_DIR, "valid_star_qc.rds")))
}
kirp <- readRDS(file.path(CACHE_DIR, "subtype_KIRP.rds"))
kich <- readRDS(file.path(CACHE_DIR, "subtype_KICH.rds"))
msg("Test cohorts: CPTAC-3 n = ", nrow(vco), ", KIRP n = ", nrow(kirp$cohort),
    ", KICH n = ", nrow(kich$cohort))

# -----------------------------------------------------------------------------
# 1. Test matrices: same genes, same transform, same within-cohort residualisation
# -----------------------------------------------------------------------------
# Genes missing or constant in a test cohort are dropped from that cohort only.
build_test_matrix <- function(fpkm, samples, gene_ids, qc, label) {
  samples <- intersect(samples, colnames(fpkm))
  gi <- intersect(gene_ids, rownames(fpkm))
  E  <- t(log2(fpkm[gi, samples, drop = FALSE] + 1))
  q  <- qc[match(rownames(E), sample_barcode)]
  cov <- tech_covariates(q)
  ok  <- stats::complete.cases(cov)
  if (!all(ok)) msg(label, ": ", sum(!ok), " samples without STAR metrics dropped")
  E <- E[ok, , drop = FALSE]; cov <- cov[ok, , drop = FALSE]
  E <- E[, colSums(!is.finite(E)) == 0, drop = FALSE]
  E <- remove_technical(E, cov)
  sdv <- apply(E, 2, sd)
  E <- E[, is.finite(sdv) & sdv > 0, drop = FALSE]
  msg(label, ": ", nrow(E), " samples x ", ncol(E), "/", length(gene_ids),
      " network transcripts after residualisation")
  E
}

test_sets_for <- function(net) {
  g <- colnames(net$expr)
  list(CPTAC3 = build_test_matrix(val$fpkm,  vco$sample_barcode,          g, vqc,
                                  paste0(net$tag, " CPTAC-3")),
       KIRP   = build_test_matrix(kirp$fpkm, kirp$cohort$sample_barcode, g,
                                  as.data.table(kirp$qc), paste0(net$tag, " KIRP")),
       KICH   = build_test_matrix(kich$fpkm, kich$cohort$sample_barcode, g,
                                  as.data.table(kich$qc), paste0(net$tag, " KICH")))
}
TEST_LABELS <- c(CPTAC3 = "CPTAC-3", KIRP = "TCGA-KIRP", KICH = "TCGA-KICH")

# -----------------------------------------------------------------------------
# 2. modulePreservation() with the discovery colours as reference
# -----------------------------------------------------------------------------
reference_colours <- function(net) {
  cols <- net$colors
  if (is.null(names(cols))) names(cols) <- net$gene_tbl$gene_id
  cols <- as.character(cols[colnames(net$expr)]); names(cols) <- colnames(net$expr)
  stopifnot(!anyNA(cols))
  cols
}

run_preservation <- function(net, tests, n_perm, tag) {
  cache_f <- file.path(CACHE_DIR, paste0("module_preservation_", tag, ".rds"))
  n_test  <- vapply(tests, nrow, integer(1))
  if (file.exists(cache_f)) {
    mp <- readRDS(cache_f)
    # Valid only for the same permutation count and test sample sizes.
    if (identical(mp$n_perm, n_perm) &&
        identical(unname(mp$test_samples), unname(n_test))) {
      msg(tag, ": modulePreservation() result from cache (", basename(cache_f), ")")
      return(mp)
    }
    msg(tag, ": cached result is stale (", mp$n_perm, " permutations, test n = ",
        paste(mp$test_samples, collapse = "/"), "; requested ", n_perm,
        " permutations, test n = ", paste(n_test, collapse = "/"), "); recomputing")
  }
  cols <- reference_colours(net)
  multiData  <- c(list(KIRC = list(data = net$expr)),
                  lapply(tests, function(E) list(data = E)))
  # Test sets get the reference colours restricted to their own genes.
  multiColor <- c(list(KIRC = cols),
                  lapply(tests, function(E) cols[colnames(E)]))
  msg(tag, ": reference ", nrow(net$expr), " x ", ncol(net$expr), ", ",
      sum(cols != "grey"), " assigned genes in ",
      length(setdiff(unique(cols), "grey")), " modules; ", n_perm, " permutations")
  t0 <- proc.time()[["elapsed"]]
  mp <- modulePreservation(multiData, multiColor,
                           dataIsExpr        = TRUE,
                           networkType       = "signed",
                           referenceNetworks = 1,
                           nPermutations     = n_perm,
                           randomSeed        = SEED,
                           quickCor          = 0,
                           maxModuleSize     = PRESERVATION_MAX_MODULE,
                           maxGoldModuleSize = PRESERVATION_MAX_MODULE,
                           savePermutedStatistics = FALSE,
                           verbose           = 1)
  el <- proc.time()[["elapsed"]] - t0
  msg(sprintf("%s: modulePreservation() finished in %.1f min", tag, el / 60))
  mp$elapsed_sec  <- el
  mp$n_perm       <- n_perm
  mp$test_samples <- n_test
  mp$test_genes   <- vapply(tests, ncol, integer(1))
  saveRDS(mp, cache_f)
  mp
}

# -----------------------------------------------------------------------------
# 3. Tidy the nested result into one row per module per test cohort
# -----------------------------------------------------------------------------
# Nested lists are named "ref.<set>" and "inColumnsAlsoPresentIn.<set>" and
# are matched by set name, not position.
pick_set <- function(lst, set) {
  i <- grep(paste0("\\.", set, "$"), names(lst))
  if (length(i) != 1) return(NULL)
  lst[[i]]
}
# WGCNA returns the Z tables as data frames but the observed tables as plain
# matrices, so columns are looked up by colnames() and read with [, nm].
col_or_na <- function(df, nm) {
  if (is.null(df) || !nm %in% colnames(df)) return(NA_real_)
  as.numeric(df[, nm])
}

tidy_preservation <- function(mp, network, ref = "KIRC") {
  rbindlist(lapply(names(TEST_LABELS), function(set) {
    Zp <- pick_set(pick_set(mp$preservation$Z,        ref), set)
    Op <- pick_set(pick_set(mp$preservation$observed, ref), set)
    Zq <- pick_set(pick_set(mp$quality$Z,             ref), set)
    Oq <- pick_set(pick_set(mp$quality$observed,      ref), set)
    if (is.null(Zp)) { msg("  no preservation table for ", set); return(NULL) }
    mods <- rownames(Zp)
    al <- function(df, nm) { v <- col_or_na(df, nm)
                             if (length(v) == 1 && is.na(v)) rep(NA_real_, length(mods))
                             else v[match(mods, rownames(df))] }
    data.table(
      network        = network,
      test_cohort    = TEST_LABELS[[set]],
      module         = mods,
      moduleSize     = al(Zp, "moduleSize"),
      Zsummary       = al(Zp, "Zsummary.pres"),
      Z.density      = al(Zp, "Z.density.pres"),
      Z.connectivity = al(Zp, "Z.connectivity.pres"),
      medianRank     = al(Op, "medianRank.pres"),
      medianRank.density      = al(Op, "medianRankDensity.pres"),
      medianRank.connectivity = al(Op, "medianRankConnectivity.pres"),
      Zsummary.qual  = al(Zq, "Zsummary.qual"),
      medianRank.qual= al(Oq, "medianRank.qual"),
      n_samples_test = unname(mp$test_samples[[set]]),
      n_genes_test   = unname(mp$test_genes[[set]]),
      n_permutations = mp$n_perm)
  }))
}
classify_Z <- function(z) fifelse(is.na(z), NA_character_,
                          fifelse(z < 2, "none", fifelse(z <= 10, "weak_moderate", "strong")))

# -----------------------------------------------------------------------------
# 4. lncRNA network (primary)
# -----------------------------------------------------------------------------
banner("24 | lncRNA network")
t_lnc <- proc.time()[["elapsed"]]
lnc_tests <- test_sets_for(nets$lnc)
mp_lnc  <- run_preservation(nets$lnc, lnc_tests, PRESERVATION_N_PERM, "lncRNA")
res_lnc <- tidy_preservation(mp_lnc, "lncRNA")
res_lnc[, preservation := classify_Z(Zsummary)]
setorder(res_lnc, test_cohort, -Zsummary, na.last = TRUE)
save_tsv(res_lnc[, .(network, test_cohort, module, moduleSize,
                     Zsummary = round(Zsummary, 2), Z.density = round(Z.density, 2),
                     Z.connectivity = round(Z.connectivity, 2), medianRank,
                     medianRank.density, medianRank.connectivity,
                     Zsummary.qual = round(Zsummary.qual, 2), medianRank.qual,
                     preservation, n_samples_test, n_genes_test, n_permutations)],
         "24_module_preservation_lncRNA.tsv")
print(res_lnc[!module %in% c("gold", "grey"),
              .(test_cohort, module, moduleSize, Zsummary = round(Zsummary, 1),
                medianRank, preservation)])
msg(sprintf("lncRNA network: %.1f min including test-matrix construction",
            (proc.time()[["elapsed"]] - t_lnc) / 60))

# -----------------------------------------------------------------------------
# 5. Protein-coding network (switchable)
# -----------------------------------------------------------------------------
res_mrna <- NULL
if (RUN_MRNA_PRESERVATION) {
  banner("24 | protein-coding network")
  t_mrna <- proc.time()[["elapsed"]]
  mrna_tests <- test_sets_for(nets$mrna)
  mp_mrna  <- run_preservation(nets$mrna, mrna_tests, MRNA_N_PERM, "mRNA")
  res_mrna <- tidy_preservation(mp_mrna, "mRNA")
  res_mrna[, preservation := classify_Z(Zsummary)]
  setorder(res_mrna, test_cohort, -Zsummary, na.last = TRUE)
  save_tsv(res_mrna[, .(network, test_cohort, module, moduleSize,
                        Zsummary = round(Zsummary, 2), Z.density = round(Z.density, 2),
                        Z.connectivity = round(Z.connectivity, 2), medianRank,
                        medianRank.density, medianRank.connectivity,
                        Zsummary.qual = round(Zsummary.qual, 2), medianRank.qual,
                        preservation, n_samples_test, n_genes_test, n_permutations)],
           "24_module_preservation_mRNA.tsv")
  print(res_mrna[!module %in% c("gold", "grey"),
                 .(test_cohort, module, moduleSize, Zsummary = round(Zsummary, 1),
                   medianRank, preservation)])
  msg(sprintf("protein-coding network: %.1f min including test-matrix construction",
              (proc.time()[["elapsed"]] - t_mrna) / 60))
} else {
  msg("RUN_MRNA_PRESERVATION = FALSE: protein-coding network skipped")
}

# -----------------------------------------------------------------------------
# 6. Summary across networks and cohorts, and figure
# -----------------------------------------------------------------------------
res_all <- rbind(res_lnc, res_mrna)
summ <- res_all[!module %in% c("gold", "grey"),
  .(n_modules = .N,
    n_strong = sum(preservation == "strong", na.rm = TRUE),
    n_weak_moderate = sum(preservation == "weak_moderate", na.rm = TRUE),
    n_none = sum(preservation == "none", na.rm = TRUE),
    median_Zsummary = round(median(Zsummary, na.rm = TRUE), 2),
    min_Zsummary = round(min(Zsummary, na.rm = TRUE), 2),
    max_Zsummary = round(max(Zsummary, na.rm = TRUE), 2),
    n_samples_test = n_samples_test[1], n_genes_test = n_genes_test[1],
    n_permutations = n_permutations[1]),
  by = .(network, test_cohort)]
save_tsv(summ, "24_module_preservation_summary.tsv")
print(summ)

pd <- res_all[!module %in% c("gold", "grey") & is.finite(Zsummary)]
if (nrow(pd)) {
  pd[, cohort := factor(test_cohort, levels = unname(TEST_LABELS))]
  g <- ggplot(pd, aes(moduleSize, Zsummary, colour = module)) +
    geom_hline(yintercept = c(2, 10), linetype = 2, colour = "grey50") +
    geom_point(size = 2.2) +
    scale_colour_identity() + scale_x_log10() +
    facet_grid(network ~ cohort, scales = "free_y") +
    labs(title = "Preservation of TCGA-KIRC modules in external cohorts",
         subtitle = "WGCNA modulePreservation(); dashed lines at Zsummary = 2 and 10",
         x = "Module size (genes)", y = "Zsummary") +
    theme_bw()
  save_fig(g, "24_module_preservation_Zsummary", 9, 5.5)
}

write_session_info("24_module_preservation")
banner("24 | done")
