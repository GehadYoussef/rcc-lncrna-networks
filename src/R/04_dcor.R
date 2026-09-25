# 04_dcor.R: distance correlation between mRNA and lncRNA hub genes.
#
# Hub genes of survival-associated modules on each side are paired. dCor is
# calibrated against a permutation null (BH FDR), the top pairs get an exact
# per-pair permutation test, and pairs are also ranked by dCor in excess of the
# stronger of |Pearson| and |Spearman|. Stops if either side has no significant module.
# Inputs: caches networks.rds, survival.rds.
# Outputs: cache dcor_matrix.rds, results 04_top_dcor_pairs.tsv,
#   04_all_significant_dcor_pairs.tsv, 04_top_nonlinear_pairs.tsv and
#   04_nonlinear_lncRNA_concentration.tsv, figures 04_dcor_null_vs_observed and
#   04_dcor_heatmap_top40.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(energy); library(ggplot2)
})
banner("04 | mRNA-lncRNA distance correlation")
set.seed(SEED)

nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
surv <- readRDS(file.path(CACHE_DIR, "survival.rds"))

# ---- choose modules ---------------------------------------------------------
pick_modules <- function(s, tag) {
  m <- s[fdr_adj < FDR_ALPHA, module]
  crit <- "FDR-significant after adjustment"
  if (!length(m)) { m <- s[fdr_uni < FDR_ALPHA, module]
                    crit <- "FDR-significant unadjusted (no module survived adjustment)" }
  if (!length(m)) stop("No significant ", tag, " modules -- dCor cannot proceed. ",
                       "This is a real negative result, not a file to be written empty.")
  msg(tag, ": ", length(m), " module(s) carried forward [", crit, "]: ",
      paste(m, collapse = ", "))
  m
}
mods_m <- pick_modules(surv$mrna, "mRNA")
mods_l <- pick_modules(surv$lnc,  "lncRNA")

# ---- hub genes --------------------------------------------------------------
hubs <- function(net, mods, cap) {
  g <- net$gene_tbl[module %in% mods]
  setorder(g, module, -kME)
  g <- g[, head(.SD, DCOR_TOP_HUBS), by = module]
  setorder(g, -kME)
  head(g, cap)
}
hub_m <- hubs(nets$mrna, mods_m, DCOR_MAX_PER_SIDE)
hub_l <- hubs(nets$lnc,  mods_l, DCOR_MAX_PER_SIDE)
msg("Hub genes: ", nrow(hub_m), " mRNA x ", nrow(hub_l), " lncRNA = ",
    format(nrow(hub_m) * nrow(hub_l), big.mark = ","), " pairs")

# The two networks may have dropped different sample outliers, so intersect.
common <- intersect(rownames(nets$mrna$expr), rownames(nets$lnc$expr))
msg("Samples common to both networks: ", length(common))

Xm <- obs_expr(nets$mrna)[common, hub_m$gene_id, drop = FALSE]
Xl <- obs_expr(nets$lnc)[common,  hub_l$gene_id, drop = FALSE]
stopifnot(!is.null(colnames(Xm)), !is.null(colnames(Xl)))

# ---- observed dCor matrix (cached) ------------------------------------------
dm_rds <- file.path(CACHE_DIR, "dcor_matrix.rds")
dm <- if (file.exists(dm_rds)) readRDS(dm_rds) else NULL
if (!is.null(dm) &&
    identical(rownames(dm), colnames(Xm)) &&
    identical(colnames(dm), colnames(Xl))) {
  msg("dCor matrix loaded from cache")
} else {
  msg("Computing dCor ...")
  t0 <- Sys.time()
  dm <- matrix(NA_real_, nrow = ncol(Xm), ncol = ncol(Xl),
               dimnames = list(colnames(Xm), colnames(Xl)))
  for (i in seq_len(ncol(Xm))) {
    xi <- Xm[, i]
    for (j in seq_len(ncol(Xl))) dm[i, j] <- energy::dcor(xi, Xl[, j])
    if (i %% 25 == 0) msg("  ", i, " / ", ncol(Xm), "  (",
                          round(difftime(Sys.time(), t0, units = "mins"), 1), " min)")
  }
  saveRDS(dm, dm_rds)
}
stopifnot(!is.null(dimnames(dm)), !anyNA(dm))
msg("dCor matrix: ", nrow(dm), " x ", ncol(dm),
    "  range ", round(min(dm), 3), " - ", round(max(dm), 3))

# ---- global permutation null ------------------------------------------------
# Screening null: random gene pairs with one member's samples permuted. The
# exact null depends on each pair's marginals, so top hits are re-tested below.
msg("Building null distribution (", DCOR_N_PERM, " draws) ...")
null_vals <- vapply(seq_len(DCOR_N_PERM), function(b) {
  i <- sample.int(ncol(Xm), 1); j <- sample.int(ncol(Xl), 1)
  energy::dcor(Xm[, i], Xl[sample.int(nrow(Xl)), j])
}, numeric(1))

obs <- as.vector(dm)
p_emp <- (1 + vapply(obs, function(o) sum(null_vals >= o), numeric(1))) /
         (1 + DCOR_N_PERM)
fdr   <- p.adjust(p_emp, "BH")

pairs <- data.table(
  mRNA_gene_id   = rep(rownames(dm), times = ncol(dm)),
  lncRNA_gene_id = rep(colnames(dm), each  = nrow(dm)),
  dCor = obs, p_screen = p_emp, fdr_screen = fdr
)
pairs[hub_m, on = .(mRNA_gene_id = gene_id),
      `:=`(mRNA_gene = i.gene_name, mRNA_module = i.module)]
pairs[hub_l, on = .(lncRNA_gene_id = gene_id),
      `:=`(lncRNA_gene = i.gene_name, lncRNA_module = i.module)]

# ---- dCor in excess of monotone correlation ---------------------------------
# Hub genes are linearly co-expressed by construction, so raw dCor mostly
# reflects correlation. The excess is taken over max(|Pearson|, |Spearman|)
# so that monotone dependence in skewed transcripts is not counted as non-linear.
pm <- cor(Xm, Xl, method = "pearson")
sm <- cor(Xm, Xl, method = "spearman")
pairs[, pearson     := as.vector(pm)]
pairs[, spearman    := as.vector(sm)]
pairs[, monotone_max := pmax(abs(pearson), abs(spearman))]
pairs[, dcor_excess := dCor - monotone_max]

# Flag antisense and divergent lncRNAs paired with their own locus gene,
# which correlate through shared promoters or read-through.
stem <- sub("-(AS[0-9]*|DT|OT[0-9]*|IT[0-9]*)$", "", pairs$lncRNA_gene)
pairs[, same_locus := !is.na(stem) & stem != lncRNA_gene & stem == mRNA_gene]

setorder(pairs, -dCor)

n_sig <- sum(pairs$fdr_screen < FDR_ALPHA)
msg("Pairs with screening FDR < ", FDR_ALPHA, ": ",
    format(n_sig, big.mark = ","), " / ", format(nrow(pairs), big.mark = ","),
    " (", round(100 * n_sig / nrow(pairs)), "%)")
if (n_sig / nrow(pairs) > 0.5)
  msg("NOTE: the screening null permutes away ALL structure, so almost every ",
      "co-expressed pair clears it. Treat the screening FDR as a floor, not ",
      "as evidence of a specific relationship; rank on dcor_excess instead.")
msg("Median |Pearson| among the top 200 by dCor: ",
    round(median(abs(head(pairs, 200)$pearson)), 3))
msg("Median |Spearman| among the top 200 by dCor: ",
    round(median(abs(head(pairs, 200)$spearman)), 3))
msg("Pairs where dCor exceeds every monotone measure by > 0.10: ",
    sum(pairs$dcor_excess > 0.10), " / ", nrow(pairs))
msg("Same-locus (antisense/divergent) pairs in the top 200 by dCor: ",
    sum(head(pairs, 200)$same_locus))

# ---- exact per-pair permutation test on the strongest pairs -----------------
top <- head(pairs, DCOR_CONFIRM_TOP)
msg("Exact permutation test (R = 999) on the top ", nrow(top), " pairs ...")
top[, p_exact := vapply(seq_len(.N), function(k)
  energy::dcor.test(Xm[, mRNA_gene_id[k]], Xl[, lncRNA_gene_id[k]],
                    R = 999)$p.value, numeric(1))]
top[, fdr_exact := p.adjust(p_exact, "BH")]

save_tsv(top[, .(mRNA_gene, lncRNA_gene, mRNA_module, lncRNA_module,
                 dCor = round(dCor, 4), pearson = round(pearson, 3),
                 dcor_excess = round(dcor_excess, 3), same_locus,
                 p_screen = signif(p_screen, 3),
                 fdr_screen = signif(fdr_screen, 3),
                 p_exact = signif(p_exact, 3), fdr_exact = signif(fdr_exact, 3),
                 mRNA_gene_id, lncRNA_gene_id)],
         "04_top_dcor_pairs.tsv")
save_tsv(pairs[fdr_screen < FDR_ALPHA,
               .(mRNA_gene, lncRNA_gene, mRNA_module, lncRNA_module,
                 dCor = round(dCor, 4), pearson = round(pearson, 3),
                 dcor_excess = round(dcor_excess, 3), same_locus,
                 fdr_screen = signif(fdr_screen, 3))],
         "04_all_significant_dcor_pairs.tsv")
cat("\n-- Top 20 by raw dCor (dominated by linear co-expression) --\n")
print(head(top[, .(mRNA_gene, lncRNA_gene, dCor = round(dCor, 3),
                   pearson = round(pearson, 3), same_locus,
                   fdr_exact = signif(fdr_exact, 3))], 20))

# ---- pairs ranked by non-linear excess --------------------------------------
nl <- pairs[fdr_screen < FDR_ALPHA & same_locus == FALSE]
setorder(nl, -dcor_excess)
nl_top <- head(nl, DCOR_CONFIRM_TOP)
msg("Exact permutation test on the top ", nrow(nl_top),
    " pairs ranked by non-linear excess ...")
nl_top[, p_exact := vapply(seq_len(.N), function(k)
  energy::dcor.test(Xm[, mRNA_gene_id[k]], Xl[, lncRNA_gene_id[k]],
                    R = 999)$p.value, numeric(1))]
nl_top[, fdr_exact := p.adjust(p_exact, "BH")]
save_tsv(nl_top[, .(mRNA_gene, lncRNA_gene, mRNA_module, lncRNA_module,
                    dCor = round(dCor, 4), pearson = round(pearson, 3),
                    spearman = round(spearman, 3),
                    monotone_max = round(monotone_max, 3),
                    dcor_excess = round(dcor_excess, 3),
                    p_exact = signif(p_exact, 3),
                    fdr_exact = signif(fdr_exact, 3),
                    mRNA_gene_id, lncRNA_gene_id)],
         "04_top_nonlinear_pairs.tsv")
cat("\n-- Top 20 by excess over the STRONGER of Pearson/Spearman --\n")
print(head(nl_top[, .(mRNA_gene, lncRNA_gene, dCor = round(dCor, 3),
                      pearson = round(pearson, 3),
                      spearman = round(spearman, 3),
                      excess = round(dcor_excess, 3),
                      fdr_exact = signif(fdr_exact, 3))], 20))

# Concentration: pairs per lncRNA, since one lncRNA against a co-expressed
# block counts as one relationship.
conc <- nl_top[, .N, by = lncRNA_gene][order(-N)]
msg("Distinct lncRNAs in the list: ", nrow(conc),
    "; largest share held by one lncRNA: ",
    round(100 * conc$N[1] / sum(conc$N)), "%")
save_tsv(conc, "04_nonlinear_lncRNA_concentration.tsv")

# ---- figures ----------------------------------------------------------------
hist_dt <- rbind(data.table(dcor = obs,       set = "observed"),
                 data.table(dcor = null_vals, set = "permuted null"))
g <- ggplot(hist_dt, aes(dcor, fill = set)) +
  geom_density(alpha = 0.55, colour = NA) +
  scale_fill_manual(values = c(observed = "#B2182B", `permuted null` = "grey60"),
                    name = NULL) +
  labs(title = "Distance correlation: observed mRNA-lncRNA pairs vs permutation null",
       x = "dCor", y = "Density") +
  theme_bw()
save_fig(g, "04_dcor_null_vs_observed", 7, 4)

sub_m <- head(rownames(dm)[order(rowMeans(dm), decreasing = TRUE)], 40)
sub_l <- head(colnames(dm)[order(colMeans(dm), decreasing = TRUE)], 40)
hm <- dm[sub_m, sub_l]
rownames(hm) <- hub_m$gene_name[match(sub_m, hub_m$gene_id)]
colnames(hm) <- hub_l$gene_name[match(sub_l, hub_l$gene_id)]
save_fig_base("04_dcor_heatmap_top40", 10, 9.4,
  pheatmap::pheatmap(hm, main = "dCor: top 40 mRNA vs top 40 lncRNA hub genes",
                     fontsize_row = 6, fontsize_col = 6, silent = FALSE))

write_session_info("04_dcor")
banner("04 | done")
