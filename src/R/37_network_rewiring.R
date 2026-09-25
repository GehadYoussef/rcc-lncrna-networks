# 37_network_rewiring.R: how the library-quality axis reorganises the lncRNA network
# Compares the unadjusted lncRNA network with the network built on expression
# residualised on the three STAR metrics (same 511 libraries and parameters):
#   1. adjusted Rand index between the partitions (all genes, and genes
#      assigned in both)
#   2. for each adjusted module, the share of its genes placed in the
#      unadjusted dominant module, another module or grey
#   3. GENCODE v36 positional class of genes in the unadjusted dominant module.
# Inputs: cache network_lncRNA_unadjusted.rds (23), network_lncRNA.rds (02),
#         results/28_lncRNA_positional_classes.tsv (28, optional).
# Outputs: 37_network_partition_agreement.tsv, 37_adjusted_module_fate.tsv,
#          37_dominant_module_positional_class.tsv

source(file.path(if (exists("R_DIR")) R_DIR else "analysis/R", "00_config.R"))
suppressPackageStartupMessages(library(data.table))
banner("37  network rewiring by the library-quality axis")

u <- readRDS(file.path(CACHE_DIR, "network_lncRNA_unadjusted.rds"))$net$colors
a <- readRDS(file.path(CACHE_DIR, "network_lncRNA.rds"))$net$colors
g <- intersect(names(u), names(a))
stopifnot(length(g) == 3442)
u <- u[g]; a <- a[g]

ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  s_ij <- sum(choose(tab, 2)); s_i <- sum(choose(rowSums(tab), 2)); s_j <- sum(choose(colSums(tab), 2))
  e <- s_i * s_j / choose(n, 2)
  (s_ij - e) / ((s_i + s_j) / 2 - e)
}
both <- u != "grey" & a != "grey"
agree <- data.table(
  comparison = c("all genes, grey as a class", "genes assigned in both networks"),
  n_genes = c(length(g), sum(both)),
  adjusted_rand_index = round(c(ari(u, a), ari(u[both], a[both])), 4),
  n_modules_unadjusted = length(setdiff(unique(u), "grey")),
  n_modules_adjusted = length(setdiff(unique(a), "grey")))
save_tsv(agree, "37_network_partition_agreement.tsv"); print(agree)

dom <- names(sort(table(u[u != "grey"]), decreasing = TRUE))[1]
fate <- data.table(gene = g, adjusted = a, unadjusted = u)[adjusted != "grey",
  .(n_genes = .N,
    pct_in_unadjusted_dominant = round(100 * mean(unadjusted == dom), 1),
    pct_in_other_unadjusted_module = round(100 * mean(unadjusted != dom & unadjusted != "grey"), 1),
    pct_unadjusted_grey = round(100 * mean(unadjusted == "grey"), 1)), by = adjusted]
fate[, unadjusted_dominant_module := dom]
setorder(fate, -pct_in_unadjusted_dominant)
save_tsv(fate, "37_adjusted_module_fate.tsv"); print(fate)

pc_file <- file.path(RESULTS_DIR, "28_lncRNA_positional_classes.tsv")
if (file.exists(pc_file)) {
  pc <- fread(pc_file)
  idcol <- intersect(c("gene_id", "gene"), names(pc))[1]
  clcol <- intersect(c("positional_class", "class"), names(pc))[1]
  pc <- pc[, .(gene = get(idcol), class = get(clcol))]
  d <- merge(data.table(gene = g, in_dominant = u == dom), pc, by = "gene")
  cls <- d[, .(n_network = .N, n_in_dominant = sum(in_dominant),
               pct_in_dominant = round(100 * mean(in_dominant), 1)), by = class]
  cls[, overall_pct_in_dominant := round(100 * mean(d$in_dominant), 1)]
  cls[, fisher_p := vapply(class, function(k) fisher.test(table(d$class == k, d$in_dominant))$p.value, 1)]
  cls[, fdr := p.adjust(fisher_p, "BH")]
  cls[, fdr_family := "across positional classes"]
  setorder(cls, -pct_in_dominant)
  save_tsv(cls, "37_dominant_module_positional_class.tsv"); print(cls)
} else msg("28 positional class table not found; section 3 skipped")

write_session_info("37_network_rewiring")
msg("37 done")
