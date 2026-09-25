# 04b_dcor_diagnostics.R: artefact checks on the top non-linear (dCor) pairs
#
# Tests whether the non-linear ranking of stage 04 reflects real dependence or
# the behaviour of a few transcripts:
#   1. concentration of the list on single lncRNAs and on one mRNA module
#   2. marginal distributions, since zero inflation or bimodality can give high
#      dCor with near-zero Pearson r
#   3. dCor after trimming 2.5% tails and after dropping samples at the lncRNA
#      detection floor
#   4. co-expression among the partners of the leading lncRNA
# Inputs: cache/networks.rds, results/04_top_nonlinear_pairs.tsv.
# Outputs: 04b_nonlinear_lncRNA_distributions.tsv, 04b_nonlinear_robustness.tsv,
# figure 04b_nonlinear_scatter.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(energy); library(ggplot2)
})
banner("04b | dCor artefact diagnostics")
set.seed(SEED)

nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
nl   <- fread(file.path(RESULTS_DIR, "04_top_nonlinear_pairs.tsv"))
Xm   <- obs_expr(nets$mrna); Xl <- obs_expr(nets$lnc)
common <- intersect(rownames(Xm), rownames(Xl))
Xm <- Xm[common, ]; Xl <- Xl[common, ]

# ---- 1. Concentration of the non-linear list on single transcripts ----
banner("Concentration of the non-linear list")
lnc_tab <- nl[, .N, by = lncRNA_gene][order(-N)]
msg("Top 200 non-linear pairs involve ", nrow(lnc_tab), " distinct lncRNAs")
print(head(lnc_tab, 10))
msg("Share of the list held by the single commonest lncRNA: ",
    round(100 * lnc_tab$N[1] / sum(lnc_tab$N)), "%")

mod_tab <- nl[, .N, by = .(lncRNA_gene, mRNA_module)][order(-N)]
msg("\nWhich mRNA modules do the partners come from?")
print(head(mod_tab, 10))

# ---- 2. Marginal distributions of the implicated lncRNAs ----
banner("Marginal distributions")
focus <- head(lnc_tab$lncRNA_gene, 6)
ann_l <- as.data.table(nets$lnc$gene_tbl)
ids   <- ann_l[match(focus, gene_name), gene_id]

dist_tbl <- rbindlist(lapply(seq_along(focus), function(i) {
  v <- Xl[, ids[i]]
  data.table(gene = focus[i],
             n_pairs_in_top200 = lnc_tab$N[i],
             frac_at_floor = round(mean(v <= min(v) + 1e-8), 3),
             frac_below_1  = round(mean(v < 1), 3),
             skewness = round(mean((v - mean(v))^3) / sd(v)^3, 2),
             kurtosis = round(mean((v - mean(v))^4) / sd(v)^4, 2),
             p01 = round(quantile(v, 0.01), 2), median = round(median(v), 2),
             p99 = round(quantile(v, 0.99), 2), max = round(max(v), 2))
}))
print(dist_tbl)
save_tsv(dist_tbl, "04b_nonlinear_lncRNA_distributions.tsv")

# Cohort-wide reference over all lncRNAs.
allv <- apply(Xl, 2, function(v) c(mean(v <= min(v) + 1e-8),
                                   mean((v - mean(v))^4) / sd(v)^4))
msg("\nCohort-wide lncRNA reference:")
msg("  median fraction-at-floor = ", round(median(allv[1, ]), 3),
    "   (90th pct = ", round(quantile(allv[1, ], 0.9), 3), ")")
msg("  median kurtosis          = ", round(median(allv[2, ]), 2),
    "   (90th pct = ", round(quantile(allv[2, ], 0.9), 2), ")")

# ---- 3. Robustness of the top pairs: trimming and floor removal ----
banner("Robustness of the top pairs")
top <- head(nl, 25)
rob <- rbindlist(lapply(seq_len(nrow(top)), function(k) {
  x <- Xm[, top$mRNA_gene_id[k]]; y <- Xl[, top$lncRNA_gene_id[k]]
  # (a) trim the most extreme 2.5% of each variable
  keep_t <- x > quantile(x, 0.025) & x < quantile(x, 0.975) &
            y > quantile(y, 0.025) & y < quantile(y, 0.975)
  # (b) drop samples where the lncRNA is at its detection floor
  keep_f <- y > min(y) + 1e-8
  # (c) Spearman rank correlation, for comparison
  data.table(mRNA = top$mRNA_gene[k], lncRNA = top$lncRNA_gene[k],
             dCor_all = round(top$dCor[k], 3),
             dCor_trimmed = round(dcor(x[keep_t], y[keep_t]), 3),
             n_trimmed = sum(keep_t),
             dCor_no_floor = round(dcor(x[keep_f], y[keep_f]), 3),
             n_no_floor = sum(keep_f),
             spearman = round(cor(x, y, method = "spearman"), 3))
}))
rob[, retained_pct_trim := round(100 * dCor_trimmed / dCor_all)]
print(rob)
save_tsv(rob, "04b_nonlinear_robustness.tsv")
msg("\nMedian dCor retained after 2.5% trimming: ",
    round(median(rob$retained_pct_trim)), "%")

# ---- 4. Co-expression among the mRNA partners of the leading lncRNA ----
banner("Independence of the mRNA partners")
lead <- lnc_tab$lncRNA_gene[1]
pids <- unique(nl[lncRNA_gene == lead, mRNA_gene_id])
if (length(pids) > 2) {
  cm <- cor(Xm[, pids])
  msg(lead, " has ", length(pids), " partners in the top 200")
  msg("  median pairwise |r| among those partners = ",
      round(median(abs(cm[upper.tri(cm)])), 3))
  msg("  -> if this is high, they are one co-expressed block and the ",
      length(pids), " 'hits' are ONE relationship, not ", length(pids), ".")
  ev <- eigen(cor(Xm[, pids]), only.values = TRUE)$values
  msg("  first eigenvalue explains ", round(100 * ev[1] / sum(ev)),
      "% of the variance among partners")
}

# ---- 5. Scatter plots of the top pairs ----
plt <- rbindlist(lapply(seq_len(min(9, nrow(top))), function(k)
  data.table(pair = paste0(top$mRNA_gene[k], " ~ ", top$lncRNA_gene[k],
                           "\ndCor=", round(top$dCor[k], 2),
                           "  r=", round(top$pearson[k], 2)),
             x = Xm[, top$mRNA_gene_id[k]], y = Xl[, top$lncRNA_gene_id[k]])))
g <- ggplot(plt, aes(x, y)) +
  geom_point(alpha = 0.35, size = 0.9) +
  geom_smooth(method = "loess", se = FALSE, colour = "#B2182B", linewidth = 0.6) +
  facet_wrap(~ pair, scales = "free") +
  labs(title = "Top non-linear dCor pairs: the actual scatter",
       subtitle = "A real non-linear dependence should show visible structure, not a few outliers or a zero-inflated stripe",
       x = "mRNA log2(FPKM+1)", y = "lncRNA log2(FPKM+1)") +
  theme_bw(base_size = 8)
save_fig(g, "04b_nonlinear_scatter", 9, 8)
msg("Wrote figures/04b_nonlinear_scatter.png -- inspect before believing anything")

write_session_info("04b_dcor_diagnostics")
banner("04b | done")
