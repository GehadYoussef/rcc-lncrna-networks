# 19_eigengene_correlation.R: lncRNA eigengene coupling to the catabolic module
# Annotation-free test: Pearson and Spearman correlation of each prognostic
# lncRNA module eigengene with the eigengene of the protein-coding catabolic
# ("black") module.
# Reads: cache networks.rds and survival.rds.
# Writes: results/lnc_vs_black_eigengene_correlation.tsv.
if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))

suppressPackageStartupMessages(library(data.table))
CACHE <- CACHE_DIR; OUT <- RESULTS_DIR
nets <- readRDS(file.path(CACHE,"networks.rds")); surv <- readRDS(file.path(CACHE,"survival.rds"))
MEl <- nets$lnc$MEs; MEm <- nets$mrna$MEs
cm <- intersect(rownames(MEl), rownames(MEm))
cat("samples with both eigengenes:", length(cm), "\n\n")
sig_l <- as.data.table(surv$lnc)[fdr_full<0.05, .(module, HR=round(HR_full,3))]
ref_mod <- catabolic_reference_module(); blk <- MEm[cm, paste0("ME", ref_mod)]   # catabolic module, located by enrichment
cat("=== correlation of each prognostic lncRNA eigengene with the protein-coding", toupper(ref_mod), "(catabolic) eigengene ===
")
out <- rbindlist(lapply(seq_len(nrow(sig_l)), function(i){
  m <- sig_l$module[i]; v <- MEl[cm, paste0("ME", m)]
  ct <- cor.test(v, blk, method="pearson"); cs <- cor(v, blk, method="spearman")
  data.table(module=m, HR=sig_l$HR[i], reference_module=ref_mod, pearson_r=round(unname(ct$estimate),3),
             CI_lo=round(ct$conf.int[1],3), CI_hi=round(ct$conf.int[2],3),
             p=signif(ct$p.value,3), spearman=round(cs,3))}))
print(out)
fwrite(out, file.path(OUT, "lnc_vs_black_eigengene_correlation.tsv"), sep="\t")
cat("\nreference module HR (protein-coding", ref_mod, "):", as.data.table(surv$mrna)[module==ref_mod, round(HR_full,3)], "\n")
cat("\nannotation package versions used for GO:\n")
for (p in c("org.Hs.eg.db","GO.db","AnnotationDbi")) if (requireNamespace(p, quietly = TRUE)) cat(" ", p, as.character(packageVersion(p)), "\n")
