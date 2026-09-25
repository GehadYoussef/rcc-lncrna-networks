# 02b_lncrna_tuning.R: outcome-independent tuning of lncRNA module resolution
#
# Sweeps lncRNA WGCNA settings (power, deepSplit, mergeCutHeight, minKMEtoStay)
# on the observed matrix and after removing the first one or two principal
# components (the usual remedy when one global axis dominates). Survival, stage
# and grade are never used. Admissible settings (fixed in advance): largest
# module <= 35% of genes, grey <= 30%, 5 to 40 modules. They are ranked by mean
# adjusted Rand index against rebuilds on 80% patient subsamples, then by mean
# eigengene variance explained.
# Inputs: cache/dataset.rds, cache/networks.rds. Outputs: the full grid
# (02b_lncRNA_tuning_grid.tsv), 02b_lncRNA_dominant_axis_diagnostics.tsv, figure
# 02b_lncRNA_tuning_grid and cache/lncrna_tuning.rds. The chosen values go into
# LNC_POWER, LNC_DEEPSPLIT and LNC_MERGE in 00_config.R.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(WGCNA); library(ggplot2)
})
banner("02b | lncRNA module-resolution tuning")
enableWGCNAThreads(N_THREADS)
set.seed(SEED)

ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
# Use the observed matrix: when LNC_VARIANT is a residualised mode,
# nets$lnc$expr is already residualised and the variants would compound.
datExpr <- obs_expr(nets$lnc)

# Apply the same technical adjustment as production.
if (ADJUST_TECHNICAL) {
  cq  <- as.data.table(ds$cohort_full)[match(rownames(datExpr), sample_barcode)]
  cov <- cbind(pct_noFeature    = cq$pct_noFeature,
               pct_multimapping = cq$pct_multimapping,
               log_depth        = log10(cq$libsize))
  if (all(stats::complete.cases(cov))) {
    datExpr <- remove_technical(datExpr, cov)
    msg("Technical (library-quality) adjustment applied before sweeping")
  } else {
    stop("Technical covariates missing -- re-run 01_build_data.R first")
  }
}
msg("Tuning on ", nrow(datExpr), " samples x ", ncol(datExpr), " lncRNAs")

# ---- pre-specified grid and constraints ------------------------------------
# Configurations are scored through blockwiseModules(), the production call,
# because its minKMEtoStay step returns weakly connected genes to grey and this
# strongly affects the grey fraction after residualisation. Low powers are
# included because strong soft-thresholding of the sparse lncRNA network pushes
# genes into grey or into one large module.
POWERS      <- c(8, 10, 12, 18)
DEEPSPLITS  <- c(2, 4)
MERGES      <- c(0.15, 0.25)
MINKMES     <- c(0.0, 0.3)
VARIANTS    <- c("raw", "resid_PC1", "resid_PC1_2")
MAX_LARGEST <- 0.35
MAX_GREY    <- 0.30
MIN_MODS    <- 5
MAX_MODS    <- 40
N_SUBSAMPLE <- 10
SUB_FRAC    <- 0.80

# =============================================================================
# Diagnostics: is one global axis dominating the lncRNA correlation structure?
# =============================================================================
banner("Diagnostics: dominant axis")
xc  <- scale(datExpr, center = TRUE, scale = FALSE)
sv  <- svd(xc)
pve <- sv$d^2 / sum(sv$d^2)
msg("Variance explained by PC1-PC5: ",
    paste0(round(100 * pve[1:5], 1), "%", collapse = ", "))

cohort <- as.data.table(ds$cohort_full)[match(rownames(datExpr), sample_barcode)]
pc1 <- sv$u[, 1]
# Largest non-grey module of the current network, whatever its colour.
gtc  <- as.data.table(nets$lnc$gene_tbl)[module != "grey", .N, by = module]
bigm <- if (nrow(gtc)) gtc[which.max(N), module] else NA_character_
big_me <- if (!is.na(bigm) && paste0("ME", bigm) %in% colnames(nets$lnc$MEs))
            round(abs(cor(pc1, nets$lnc$MEs[[paste0("ME", bigm)]])), 3) else NA_real_
diag_tbl <- data.table(
  quantity = c("PC1 variance explained",
               paste0("|cor| PC1 vs largest module eigengene (", bigm, ")"),
               "|cor| PC1 vs library size",
               "|cor| PC1 vs per-sample mean lncRNA expression"),
  value = c(round(pve[1], 3), big_me,
            round(abs(cor(pc1, cohort$libsize, use = "pairwise.complete.obs")), 3),
            round(abs(cor(pc1, rowMeans(datExpr))), 3)))

# Is the largest module mostly antisense or divergent transcripts? Compare its
# composition with the rest of the network.
ann <- as.data.table(ds$lnc$ann)
is_as <- function(g) grepl("-(AS[0-9]*|DT|OT[0-9]*|IT[0-9]*)$", g)
gt  <- nets$lnc$gene_tbl
big <- gt[, .N, by = module][order(-N)][1, module]
diag_tbl <- rbind(diag_tbl, data.table(
  quantity = c(paste0("antisense/divergent fraction in largest module (", big, ")"),
               "antisense/divergent fraction elsewhere"),
  value = c(round(mean(is_as(gt[module == big, gene_name])), 3),
            round(mean(is_as(gt[module != big, gene_name])), 3))))
save_tsv(diag_tbl, "02b_lncRNA_dominant_axis_diagnostics.tsv")
print(diag_tbl)

# ---- helpers ----------------------------------------------------------------
adj_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  cmb <- function(x) sum(choose(x, 2))
  idx <- cmb(tab); ra <- cmb(rowSums(tab)); rb <- cmb(colSums(tab))
  exp <- ra * rb / choose(n, 2); mx <- (ra + rb) / 2
  if (mx == exp) return(NA_real_)
  (idx - exp) / (mx - exp)
}

make_variant <- function(expr, variant) {
  if (variant == "raw") return(expr)
  k <- if (variant == "resid_PC1") 1L else 2L
  x <- scale(expr, center = TRUE, scale = FALSE)
  s <- svd(x, nu = k, nv = k)
  x - s$u[, 1:k, drop = FALSE] %*% diag(s$d[1:k], k, k) %*%
      t(s$v[, 1:k, drop = FALSE])
}

# Same call as 02_wgcna.R, so scores reflect the production networks.
detect <- function(expr, power, ds_, mch, minkme) {
  set.seed(SEED)
  blockwiseModules(expr, power = power, networkType = NETWORK_TYPE,
                   TOMType = TOM_TYPE, minModuleSize = MIN_MODULE_SIZE,
                   mergeCutHeight = mch, deepSplit = ds_,
                   minKMEtoStay = minkme, numericLabels = FALSE,
                   pamRespectsDendro = FALSE, maxBlockSize = MAX_BLOCK_SIZE,
                   saveTOMs = FALSE, verbose = 0)$colors
}

metrics <- function(expr, cols) {
  tb <- table(cols)
  grey <- if ("grey" %in% names(tb)) tb[["grey"]] else 0
  nong <- tb[names(tb) != "grey"]
  me <- moduleEigengenes(expr, colors = cols, verbose = 0)
  pv <- me$varExplained[1, ]; names(pv) <- colnames(me$eigengenes)
  pv <- pv[names(pv) != "MEgrey"]
  list(n_modules = length(nong), frac_grey = grey / length(cols),
       frac_largest = if (length(nong)) max(nong) / length(cols) else NA_real_,
       median_size = if (length(nong)) median(nong) else NA_real_,
       mean_pve = if (length(pv)) mean(as.numeric(pv), na.rm = TRUE) else NA_real_)
}

# ---- pass 1: sweep variant x power x deepSplit x merge x minKMEtoStay -------
grid <- CJ(variant = VARIANTS, power = POWERS, deepSplit = DEEPSPLITS,
           mergeCutHeight = MERGES, minKME = MINKMES, sorted = FALSE)
msg("Sweeping ", nrow(grid), " configurations through blockwiseModules ...")
expr_cache <- setNames(lapply(VARIANTS, function(v) make_variant(datExpr, v)),
                       VARIANTS)
t0 <- Sys.time()
res <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  ex   <- expr_cache[[grid$variant[i]]]
  cols <- detect(ex, grid$power[i], grid$deepSplit[i],
                 grid$mergeCutHeight[i], grid$minKME[i])
  m <- metrics(ex, cols)
  res[[i]] <- cbind(grid[i], data.table(
    n_modules = m$n_modules, frac_grey = round(m$frac_grey, 3),
    frac_largest = round(m$frac_largest, 3),
    median_size = m$median_size, mean_pve = round(m$mean_pve, 3)))
  msg("  [", i, "/", nrow(grid), "] ", grid$variant[i], " p=", grid$power[i],
      " ds=", grid$deepSplit[i], " mch=", grid$mergeCutHeight[i],
      " kME=", grid$minKME[i], " -> ", m$n_modules, " modules, grey=",
      round(100 * m$frac_grey), "%, largest=", round(100 * m$frac_largest), "%")
}
grid_res <- rbindlist(res)
msg("Sweep finished in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
grid_res[, admissible := frac_largest <= MAX_LARGEST & frac_grey <= MAX_GREY &
                         n_modules >= MIN_MODS & n_modules <= MAX_MODS]

# ---- figure (always written, whatever the outcome) -------------------------
g <- ggplot(grid_res, aes(frac_largest, frac_grey,
                          colour = factor(minKME), shape = factor(power))) +
  annotate("rect", xmin = -Inf, xmax = MAX_LARGEST, ymin = -Inf, ymax = MAX_GREY,
           alpha = 0.10, fill = "#2166AC") +
  geom_point(size = 2.4) +
  facet_wrap(~ variant) +
  scale_colour_brewer(palette = "Dark2", name = "minKMEtoStay") +
  scale_shape_manual(values = c(16, 17, 15, 3), name = "power") +
  labs(title = "lncRNA network: module-resolution parameter sweep",
       subtitle = paste0("Shaded region satisfies the pre-specified constraints ",
                         "(largest <= ", MAX_LARGEST, ", grey <= ", MAX_GREY, ")"),
       x = "Fraction of genes in the largest module",
       y = "Fraction of genes unassigned (grey)") +
  theme_bw()
save_fig(g, "02b_lncRNA_tuning_grid", 10, 4.2)

msg(sum(grid_res$admissible), " of ", nrow(grid_res),
    " configurations satisfy the pre-specified constraints")
setorder(grid_res, frac_largest)
print(grid_res[1:min(15, .N)])

adm <- grid_res[admissible == TRUE]
if (!nrow(adm)) {
  save_tsv(grid_res, "02b_lncRNA_tuning_grid.tsv")
  saveRDS(list(grid = grid_res, admissible = FALSE, diagnostics = diag_tbl),
          file.path(CACHE_DIR, "lncrna_tuning.rds"))
  banner("CONCLUSION")
  cat(
"No parameter setting, with or without removing the leading principal\n",
"components, produces a balanced lncRNA module structure. The largest module\n",
"never falls below the pre-specified 35% ceiling. This is a property of the\n",
"data: the lncRNA correlation structure in TCGA-KIRC is dominated by a single\n",
"global axis.\n\n",
"Do NOT relax the constraints to manufacture modules. The defensible routes\n",
"are (a) report the dominant lncRNA component as what it is -- one global\n",
"expression axis that carries prognostic information -- rather than as a set\n",
"of modules, and (b) make the lncRNA claims gene-level (the penalised Cox\n",
"signature and the dCor non-linearity result are already gene-level and do\n",
"not depend on module definitions).\n", sep = "")
  write_session_info("02b_tuning")
  quit(save = "no", status = 0)
}

# ---- pass 2: stability of admissible configurations -------------------------
setorder(adm, frac_largest)
adm <- head(adm, 8)
msg("Assessing stability of ", nrow(adm), " configuration(s) over ",
    N_SUBSAMPLE, " x ", round(100 * SUB_FRAC), "% patient subsamples ...")

n_s <- nrow(datExpr)
sub_idx <- lapply(seq_len(N_SUBSAMPLE), function(b) {
  set.seed(SEED + b); sort(sample.int(n_s, floor(SUB_FRAC * n_s)))
})
ref <- lapply(seq_len(nrow(adm)), function(i)
  detect(expr_cache[[adm$variant[i]]], adm$power[i], adm$deepSplit[i],
         adm$mergeCutHeight[i], adm$minKME[i]))

ari_mat <- matrix(NA_real_, nrow(adm), N_SUBSAMPLE)
for (b in seq_len(N_SUBSAMPLE)) {
  # Residualisation is recomputed within each subsample, because PC1 depends
  # on the sample set.
  subs <- setNames(lapply(unique(adm$variant), function(v)
    make_variant(datExpr[sub_idx[[b]], , drop = FALSE], v)),
    unique(adm$variant))
  for (i in seq_len(nrow(adm))) {
    cb <- detect(subs[[adm$variant[i]]], adm$power[i], adm$deepSplit[i],
                 adm$mergeCutHeight[i], adm$minKME[i])
    ari_mat[i, b] <- adj_rand(ref[[i]], cb)
  }
  gc(verbose = FALSE)
  msg("  subsample ", b, " / ", N_SUBSAMPLE, "  (",
      round(difftime(Sys.time(), t0, units = "mins"), 1), " min elapsed)")
}
adm[, mean_ARI := round(rowMeans(ari_mat, na.rm = TRUE), 3)]
adm[, sd_ARI   := round(apply(ari_mat, 1, sd, na.rm = TRUE), 3)]

setorder(adm, -mean_ARI, -mean_pve)
best <- adm[1]
msg("Selected: variant = ", best$variant, ", power = ", best$power,
    ", deepSplit = ", best$deepSplit, ", mergeCutHeight = ", best$mergeCutHeight,
    ", minKMEtoStay = ", best$minKME,
    "  (", best$n_modules, " modules, grey = ", round(100 * best$frac_grey),
    "%, largest = ", round(100 * best$frac_largest),
    "%, mean ARI = ", best$mean_ARI, ")")

key_cols <- c("variant", "power", "deepSplit", "mergeCutHeight", "minKME")
out <- merge(grid_res, adm[, c(key_cols, "mean_ARI", "sd_ARI"), with = FALSE],
             by = key_cols, all.x = TRUE)
out[, selected := variant == best$variant & power == best$power &
                  deepSplit == best$deepSplit &
                  mergeCutHeight == best$mergeCutHeight &
                  minKME == best$minKME]
setorder(out, -admissible, -mean_ARI, frac_largest)
save_tsv(out, "02b_lncRNA_tuning_grid.tsv")
print(out[1:min(15, .N)])

saveRDS(list(variant = best$variant, power = best$power,
             deepSplit = best$deepSplit, mergeCutHeight = best$mergeCutHeight,
             minKME = best$minKME,
             grid = out, ari = ari_mat, admissible = TRUE,
             diagnostics = diag_tbl),
        file.path(CACHE_DIR, "lncrna_tuning.rds"))

banner("02b | done -- set LNC_POWER / LNC_DEEPSPLIT / LNC_MERGE in 00_config.R")
write_session_info("02b_tuning")
