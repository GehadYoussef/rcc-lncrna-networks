# 02_wgcna.R: weighted co-expression networks for mRNA and lncRNA (Langfelder and Horvath 2008).
#
# Signed networks in a single block, seeded for reproducible module colours.
# Soft power is the lowest reaching the scale-free fit target unless fixed in
# 00_config.R. With ADJUST_TECHNICAL, library-quality covariates are regressed
# out before module detection.
# Module-trait correlations use clinicopathological traits only, since
# correlating with os_time ignores censoring. Survival is modelled in 03.
# Inputs: cache dataset.rds. Outputs: caches networks.rds and network_<tag>.rds,
#   results 02_<tag>_{soft_threshold,module_genes,module_sizes,module_trait}.tsv
#   and 02_network_summary.tsv, plus soft-threshold, module-trait and dendrogram figures.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(WGCNA); library(ggplot2)
})
banner("02 | WGCNA")

enableWGCNAThreads(N_THREADS)
ds     <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
cohort <- as.data.table(ds$cohort_full)

# -----------------------------------------------------------------------------
# Network construction for one biotype
# -----------------------------------------------------------------------------
build_network <- function(obj, tag, power_fixed = NA,
                          deep_split = DEEP_SPLIT,
                          merge_cut = MERGE_CUT_HEIGHT,
                          min_kme = MIN_KME_TO_STAY,
                          variant = "raw") {
  banner(paste("Network:", tag))

  # Per-biotype cache keyed on every parameter that can change the modules.
  fp <- list(tag = tag, dim = dim(obj$expr), seed = SEED,
             power_fixed = power_fixed, deep_split = deep_split,
             merge_cut = merge_cut, min_kme = min_kme,
             variant = variant, tech_adj = ADJUST_TECHNICAL,
             net = NETWORK_TYPE,
             tom = TOM_TYPE, minmod = MIN_MODULE_SIZE,
             block = MAX_BLOCK_SIZE, rsq = RSQ_CUT)
  cache_f <- file.path(CACHE_DIR, paste0("network_", tag, ".rds"))
  if (file.exists(cache_f)) {
    cached <- readRDS(cache_f)
    if (identical(cached$fingerprint, fp)) {
      msg(tag, ": unchanged parameters -- reusing cached network")
      return(cached$net)
    }
    msg(tag, ": parameters changed -- rebuilding network")
  }

  datExpr <- t(obj$expr)                       # samples x genes
  msg(tag, ": ", nrow(datExpr), " samples x ", ncol(datExpr), " genes")

  # --- gene/sample QC ---
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    msg("Dropping ", sum(!gsg$goodGenes), " genes and ",
        sum(!gsg$goodSamples), " samples failing goodSamplesGenes()")
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  }

  # --- sample outliers by standardised connectivity (Z.k < -2.5) -------------
  A   <- adjacency(t(datExpr), type = "distance")
  k   <- colSums(A) - 1
  Z.k <- scale(k)[, 1]
  outl <- which(Z.k < -2.5)
  msg(tag, ": ", length(outl), " sample outlier(s) removed (Z.k < -2.5)")
  if (length(outl)) datExpr <- datExpr[-outl, , drop = FALSE]
  rm(A); gc(verbose = FALSE)

  # Observed expression, kept before any adjustment. Modules are defined on the
  # adjusted matrix, while dCor (04) and the predictive model (05) use observed levels.
  datExpr_obs <- datExpr

  # --- technical (library-quality) adjustment before module detection -------
  # Stromal and immune composition is left in as tumour biology.
  if (ADJUST_TECHNICAL) {
    cq <- cohort[match(rownames(datExpr), sample_barcode)]
    cov <- cbind(pct_noFeature    = cq$pct_noFeature,
                 pct_multimapping = cq$pct_multimapping,
                 log_depth        = log10(cq$libsize))
    ok <- stats::complete.cases(cov)
    if (all(ok)) {
      before <- abs(cor(svd(scale(datExpr, TRUE, FALSE), nu = 1)$u[, 1],
                        cq$pct_noFeature))
      datExpr <- remove_technical(datExpr, cov)
      after <- abs(cor(svd(scale(datExpr, TRUE, FALSE), nu = 1)$u[, 1],
                       cq$pct_noFeature))
      msg(tag, ": technical adjustment -- |cor(PC1, pct_noFeature)| ",
          round(before, 3), " -> ", round(after, 3))
    } else {
      msg(tag, ": technical covariates incomplete for ", sum(!ok),
          " samples; adjustment SKIPPED")
    }
  }

  # Optional variant: remove leading PCs before network construction.
  if (!identical(variant, "raw")) {
    datExpr <- remove_leading_pcs(datExpr, variant)
    msg(tag, ": network built on residuals after removing ",
        if (variant == "resid_PC1") "PC1" else "PC1-PC2")
  }

  # --- soft-thresholding power ---------------------------------------------
  powers <- c(1:10, seq(12, 30, 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0,
                           networkType = NETWORK_TYPE, blockSize = 5000)
  fit <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  ok  <- which(fit >= RSQ_CUT)
  auto_power <- if (length(ok)) sft$fitIndices$Power[min(ok)] else DEFAULT_POWER
  power <- if (!is.na(power_fixed)) power_fixed else auto_power
  r2_at <- fit[match(power, sft$fitIndices$Power)]
  msg(tag, ": soft power = ", power, " (scale-free R2 = ", round(r2_at, 3),
      if (!is.na(power_fixed) && power_fixed != auto_power)
        paste0(" -- set by tuning; automatic rule would give ", auto_power)
      else if (!length(ok)) " -- fit target not reached, fallback used" else "",
      ")  deepSplit = ", deep_split, ", mergeCutHeight = ", merge_cut)

  sft_dt <- as.data.table(sft$fitIndices)[, .(Power, SFT_R2 = round(fit, 3),
                                              mean_k = round(mean.k., 1))]
  save_tsv(sft_dt, paste0("02_", tag, "_soft_threshold.tsv"))

  p <- ggplot(sft_dt, aes(Power, SFT_R2)) +
    geom_hline(yintercept = RSQ_CUT, linetype = 2, colour = "red") +
    geom_vline(xintercept = power, linetype = 3) +
    geom_line() + geom_point() +
    labs(title = paste0(tag, ": scale-free topology fit (", NETWORK_TYPE, " network)"),
         subtitle = paste0("chosen power = ", power),
         x = "Soft-threshold power", y = expression(Signed~R^2)) +
    theme_bw()
  save_fig(p, paste0("02_", tag, "_soft_threshold"), 6, 4)

  # --- module detection ------------------------------------------------------
  set.seed(SEED)                                   # stable module colours
  net <- blockwiseModules(
    datExpr,
    power            = power,
    networkType      = NETWORK_TYPE,
    TOMType          = TOM_TYPE,
    minModuleSize    = MIN_MODULE_SIZE,
    mergeCutHeight   = merge_cut,
    deepSplit        = deep_split,
    minKMEtoStay     = min_kme,
    numericLabels    = FALSE,
    pamRespectsDendro= FALSE,
    maxBlockSize     = MAX_BLOCK_SIZE,
    saveTOMs         = FALSE,
    verbose          = 0
  )
  colors <- net$colors
  MEs <- orderMEs(moduleEigengenes(datExpr, colors = colors)$eigengenes)
  rownames(MEs) <- rownames(datExpr)

  msg(tag, ": ", length(unique(colors)) - ("grey" %in% colors),
      " modules (excluding grey), ", sum(colors == "grey"), " unassigned genes")

  # --- module membership (kME) ----------------------------------------------
  kME <- as.data.frame(cor(datExpr, MEs, use = "p"))
  ann <- as.data.table(obj$ann)[match(colnames(datExpr), gene_id)]
  gene_tbl <- data.table(
    gene_id   = colnames(datExpr),
    gene_name = ann$gene_name,
    module    = colors,
    kME       = vapply(seq_along(colors), function(i)
                  kME[i, paste0("ME", colors[i])], numeric(1))
  )
  setorder(gene_tbl, module, -kME)
  save_tsv(gene_tbl, paste0("02_", tag, "_module_genes.tsv"))

  mod_sizes <- gene_tbl[, .(n_genes = .N), by = module][order(-n_genes)]
  save_tsv(mod_sizes, paste0("02_", tag, "_module_sizes.tsv"))
  print(mod_sizes)

  # --- module-trait correlations (clinicopathological traits only) ----------
  cl <- cohort[match(rownames(MEs), sample_barcode)]
  traits <- data.frame(
    age        = cl$age,
    male       = as.integer(cl$sex == "male"),
    stage      = cl$stage_num,
    grade      = cl$grade_num
  )
  rownames(traits) <- rownames(MEs)
  mt_cor <- cor(MEs, traits, use = "pairwise.complete.obs")
  mt_p   <- corPvalueStudent(mt_cor, nrow(MEs))
  mt <- data.table(module = rownames(mt_cor))
  for (tr in colnames(mt_cor)) {
    mt[[paste0("r_", tr)]] <- round(mt_cor[, tr], 3)
    mt[[paste0("p_", tr)]] <- signif(mt_p[, tr], 3)
  }
  save_tsv(mt, paste0("02_", tag, "_module_trait.tsv"))

  save_fig_base(paste0("02_", tag, "_module_trait_heatmap"), 7.4, 9.5, {
  par(mar = c(6, 9, 3, 2))
  labeledHeatmap(Matrix = mt_cor, xLabels = colnames(mt_cor),
                 yLabels = rownames(mt_cor), ySymbols = rownames(mt_cor),
                 colorLabels = FALSE, colors = blueWhiteRed(50),
                 textMatrix = paste0(signif(mt_cor, 2), "\n(",
                                     signif(mt_p, 1), ")"),
                 setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1, 1),
                 main = paste(tag, ": module-trait relationships"))
  })

  save_fig_base(paste0("02_", tag, "_dendrogram"), 10.6, 5.9, {
  plotDendroAndColors(net$dendrograms[[1]],
                      colors[net$blockGenes[[1]]],
                      "Module", dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste(tag, ": gene dendrogram and modules"))
  })

  out <- list(tag = tag, power = power, deep_split = deep_split,
              merge_cut = merge_cut, variant = variant,
              colors = colors, MEs = MEs,
              gene_tbl = gene_tbl, samples = rownames(datExpr),
              expr = datExpr, expr_obs = datExpr_obs,
              n_modules = length(setdiff(unique(colors), "grey")))
  saveRDS(list(fingerprint = fp, net = out), cache_f)
  out
}

net_mrna <- build_network(ds$mrna, "mRNA")
net_lnc  <- build_network(
  ds$lnc, "lncRNA",
  power_fixed = LNC_POWER,
  deep_split  = if (is.na(LNC_DEEPSPLIT)) DEEP_SPLIT else LNC_DEEPSPLIT,
  merge_cut   = if (is.na(LNC_MERGE)) MERGE_CUT_HEIGHT else LNC_MERGE,
  min_kme     = if (is.na(LNC_MINKME)) MIN_KME_TO_STAY else LNC_MINKME,
  variant     = LNC_VARIANT)

saveRDS(list(mrna = net_mrna, lnc = net_lnc),
        file.path(CACHE_DIR, "networks.rds"))

summ <- data.table(
  biotype   = c("mRNA", "lncRNA"),
  n_genes   = c(ncol(net_mrna$expr), ncol(net_lnc$expr)),
  n_samples = c(nrow(net_mrna$expr), nrow(net_lnc$expr)),
  soft_power= c(net_mrna$power, net_lnc$power),
  n_modules = c(net_mrna$n_modules, net_lnc$n_modules),
  n_grey    = c(sum(net_mrna$colors == "grey"), sum(net_lnc$colors == "grey"))
)
save_tsv(summ, "02_network_summary.tsv")
print(summ)

write_session_info("02_wgcna")
banner("02 | done")
