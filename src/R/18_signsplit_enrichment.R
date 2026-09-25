# 18_signsplit_enrichment.R: sign-stratified guilt-by-association enrichment
# For each prognostic lncRNA module, protein-coding partners are split by the
# sign of their correlation with the module eigengene. The top 300 partners of
# each sign are tested for GO biological-process enrichment.
# Reads: cache networks.rds and survival.rds.
# Writes: results/gba_direction_summary.tsv and
#   results/10_lncRNA_guilt_by_association_GO_signed.tsv.
if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))

source(file.path(R_DIR, "ora_hypergeometric.R"))
CACHE <- CACHE_DIR; OUT <- RESULTS_DIR
obs_expr <- function(net) if (!is.null(net$expr_obs)) net$expr_obs else net$expr
nets <- readRDS(file.path(CACHE,"networks.rds")); surv <- readRDS(file.path(CACHE,"survival.rds"))
sig_l <- as.data.table(surv$lnc)[fdr_full < 0.05, module]
hr    <- as.data.table(surv$lnc)[fdr_full < 0.05, .(module, HR_full)]
Em <- obs_expr(nets$mrna); MEs_l <- nets$lnc$MEs
common <- intersect(rownames(Em), rownames(MEs_l)); universe <- sub("\\..*$","", colnames(Em))
gomap <- build_go_map(universe)

# Protein-coding module carrying the catabolic programme ("black"). It is
# located by enrichment because module colours change on rebuild.
gt_m <- as.data.table(nets$mrna$gene_tbl)
ref_mod <- catabolic_reference_module()
black <- gt_m[module == ref_mod, gene_id]
cat("reference catabolic protein-coding module:", ref_mod, "with", length(black), "genes\n\n")

all_res <- list(); summ <- list()
for (m in sig_l) {
  me <- MEs_l[[paste0("ME", m)]][match(common, rownames(MEs_l))]
  r  <- cor(Em[common, ], me, use="pairwise.complete.obs")[,1]
  for (sgn in c("positive","negative")) {
    rs <- if (sgn=="positive") r[r > 0] else r[r < 0]
    if (length(rs) < 50) next
    top <- names(sort(abs(rs), decreasing=TRUE))[1:min(300, length(rs))]
    res <- ora(sub("\\..*$","",top), universe, gomap)
    res <- res[p.adjust < 0.05]
    nblack <- sum(top %in% black)                      # overlap with the catabolic module
    summ[[paste(m,sgn)]] <- data.table(
      module=m, HR=round(hr[module==m, HR_full],3), sign=sgn,
      n_partners=length(rs), n_absr_gt_.5=sum(abs(rs)>0.5),
      mean_r_top=round(mean(r[top]),3),
      pct_of_top300_in_black=round(100*nblack/length(top),1),
      n_GO=nrow(res),
      top_term=if(nrow(res)) res$Description[1] else NA_character_,
      top_padj=if(nrow(res)) signif(res$p.adjust[1],3) else NA_real_)
    if (nrow(res)) { d <- res[1:min(10,.N)][, `:=`(module=m, sign=sgn)]; all_res[[paste(m,sgn)]] <- d }
  }
}
S <- rbindlist(summ, fill=TRUE); setorder(S, module, sign)
fwrite(S, file.path(OUT,"gba_direction_summary.tsv"), sep="\t")
cat("================ DIRECTION SUMMARY ================\n"); print(S)
G <- rbindlist(all_res, fill=TRUE)
fwrite(G[, .(module, sign, Description, GeneRatio, BgRatio, Count,
             p.adjust=signif(p.adjust,3))],
       file.path(OUT,"10_lncRNA_guilt_by_association_GO_signed.tsv"), sep="\t")
cat("\n================ TOP 5 TERMS PER MODULE x SIGN ================\n")
for (m in sig_l) for (s in c("positive","negative")) {
  d <- G[module==m & sign==s]
  if (!nrow(d)) next
  cat("\n--", m, "/", s, "--\n")
  print(d[1:min(5,.N), .(Description, GeneRatio, p.adjust=signif(p.adjust,3))])
}
cat("\nSIGNSPLIT_DONE\n")
