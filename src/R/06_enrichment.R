# 06_enrichment.R: functional enrichment and hub genes of the prognostic modules
#
# Writes hub-gene tables for prognostic mRNA and lncRNA modules and runs GO
# biological process over-representation for each prognostic mRNA module.
# Inputs: cache networks.rds (02) and survival.rds (03).
# Outputs: 06_<biotype>_hub_genes.tsv, 06_GO_enrichment_*.tsv and one GO dot
# plot per module.
# The GO universe is the set of genes that entered the network. A whole-genome
# background would inflate significance for an expression-filtered gene set.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
  library(ggplot2)
})
banner("06 | Functional enrichment and hub genes")

nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
surv <- readRDS(file.path(CACHE_DIR, "survival.rds"))

strip_ver <- function(x) sub("\\..*$", "", x)

# ---- hub gene tables --------------------------------------------------------
hub_table <- function(net, s, tag, n_top = 25) {
  sig <- s[fdr_adj < FDR_ALPHA, module]
  if (!length(sig)) sig <- s[fdr_uni < FDR_ALPHA, module]
  if (!length(sig)) { msg("No prognostic ", tag, " modules; skipping"); return(NULL) }
  g <- net$gene_tbl[module %in% sig]
  setorder(g, module, -kME)
  top <- g[, head(.SD, n_top), by = module]
  top <- merge(top, s[, .(module, HR_adj, fdr_adj)], by = "module")
  setorder(top, fdr_adj, module, -kME)
  save_tsv(top[, .(module, HR_adj = round(HR_adj, 3),
                   fdr_adj = signif(fdr_adj, 3),
                   gene_name, gene_id, kME = round(kME, 3))],
           paste0("06_", tag, "_hub_genes.tsv"))
  msg(tag, ": wrote top ", n_top, " hub genes for ", length(sig), " module(s)")
  print(head(top[, .(module, gene_name, kME = round(kME, 3))], 25))
  top
}
hub_m <- hub_table(nets$mrna, surv$mrna, "mRNA")
hub_l <- hub_table(nets$lnc,  surv$lnc,  "lncRNA")

# ---- GO over-representation for prognostic mRNA modules ---------------------
sig_m <- surv$mrna[fdr_adj < FDR_ALPHA, module]
if (!length(sig_m)) sig_m <- surv$mrna[fdr_uni < FDR_ALPHA, module]

universe <- strip_ver(nets$mrna$gene_tbl$gene_id)
all_go <- list()

for (m in sig_m) {
  genes <- strip_ver(nets$mrna$gene_tbl[module == m, gene_id])
  msg("GO enrichment for mRNA module ", m, " (", length(genes), " genes) ...")
  eg <- tryCatch(
    enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
             keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
             pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE),
    error = function(e) { msg("  failed: ", conditionMessage(e)); NULL })
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) {
    msg("  no enriched GO BP terms at q < 0.10"); next
  }
  df <- as.data.table(as.data.frame(eg))
  df[, module := m]
  all_go[[m]] <- df
  msg("  ", nrow(df), " enriched terms; top: ", df$Description[1])

  p <- dotplot(eg, showCategory = 15) +
    ggtitle(paste0("mRNA module ", m, ": GO biological process")) +
    theme(axis.text.y = element_text(size = 7))
  save_fig(p, paste0("06_GO_mRNA_", m), 7.5, 6)
}

if (length(all_go)) {
  go_all <- rbindlist(all_go, fill = TRUE)
  setcolorder(go_all, c("module", "ID", "Description", "GeneRatio",
                        "BgRatio", "pvalue", "p.adjust", "qvalue", "Count"))
  save_tsv(go_all[, .(module, ID, Description, GeneRatio, BgRatio,
                      pvalue = signif(pvalue, 3), p.adjust = signif(p.adjust, 3),
                      Count, geneID)],
           "06_GO_enrichment_all_modules.tsv")
  save_tsv(go_all[, head(.SD, 10), by = module][
             , .(module, Description, GeneRatio, p.adjust = signif(p.adjust, 3), Count)],
           "06_GO_enrichment_top10_per_module.tsv")
  msg("GO results written for ", length(all_go), " module(s)")
} else {
  msg("No GO enrichment results to write")
}

write_session_info("06_enrichment")
banner("06 | done")
