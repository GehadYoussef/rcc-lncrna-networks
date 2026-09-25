# 36_singlecell_module_localisation.R: which cells express each bulk co-expression module
# For each of the 29 discovery modules (13 protein-coding, 16 lncRNA), tests in
# independent ccRCC single-cell cohorts whether member genes shift towards malignant
# cells, using the donor-level pseudobulk contrasts of src/singlecell stage 07
# (malignant vs normal epithelium, malignant vs pooled tumour microenvironment),
# inferential rows only. Statistic: median log2 fold change of tested members against
# 2,000 random same-class gene sets matched on abundance deciles. z scores are combined
# by unweighted Stouffer, with BH within contrast across the 29 modules.
# The single-cell cohorts carry no survival data.
# Inputs: results/02_*_module_genes.tsv, results/singlecell/07_de_all_genes.tsv.gz
#   and 08_meta_lncRNA.tsv.gz. Run after the single-cell pipeline.
# Outputs: results/36_sc_module_*.tsv,
#   figures/SupplementaryFigureS10_sc_module_localisation.

source(file.path(if (exists("R_DIR")) R_DIR else "analysis/R", "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

banner("36  single-cell localisation of the 29 discovery modules")
set.seed(SEED)

# Single-cell results are in results/singlecell, or under the repository when
# run from a standalone analysis directory beside it.
SC_DIR      <- Filter(dir.exists, c(file.path(PROJECT_ROOT, "submission", "results", "singlecell"),
                                    file.path(PROJECT_ROOT, "results", "singlecell")))[1]
if (is.na(SC_DIR)) stop("single-cell results directory not found")
N_NULL      <- 2000
N_BINS      <- 10
CONTRASTS   <- c("malignant_vs_normal_epithelial", "malignant_vs_pooled_tme")
PROGNOSTIC  <- data.table(network = c("protein-coding", "protein-coding", "lncRNA", "lncRNA", "lncRNA"),
                          module  = c("green", "purple", "blue", "greenyellow", "turquoise"))

# ---- 1. module membership and gene universes ----
read_mod <- function(f, network) {
  x <- fread(file.path(RESULTS_DIR, f))
  x[, `:=`(gene_key = sub("\\..*$", "", gene_id), network = network)]
  x[, .(network, module, gene_key, gene_name)]
}
mods <- rbind(read_mod("02_mRNA_module_genes.tsv", "protein-coding"),
              read_mod("02_lncRNA_module_genes.tsv", "lncRNA"))
stopifnot(mods[network == "protein-coding", .N] == 12000,
          mods[network == "lncRNA", .N] == 3442,
          !anyDuplicated(mods$gene_key))
n_mod <- mods[module != "grey", uniqueN(paste(network, module))]
stopifnot(n_mod == 29)
msg("universe: 12,000 protein-coding + 3,442 lncRNA genes; ", n_mod, " modules")

# ---- 2. single-cell differential expression ----
de <- fread(file.path(SC_DIR, "07_de_all_genes.tsv.gz"))
de <- de[inference == "inferential" & contrast %in% CONTRASTS & is.finite(log2FC)]
de[, abund := log10((mean_cpm_ref + mean_cpm_test) / 2 + 1e-3)]
de <- merge(de, mods, by = "gene_key", suffixes = c("_sc", ""))
msg("inferential dataset x contrast pairs: ",
    paste(unique(de[, paste(dataset, contrast)]), collapse = "; "))

# ---- 3. per-dataset module statistic against an abundance-matched null ----
# Abundance matching is needed because single-cell detection depends steeply
# on abundance, and lncRNAs are low in abundance.
one_cell <- function(d) {
  d[, bin := cut(abund, breaks = unique(quantile(abund, seq(0, 1, length.out = N_BINS + 1))),
                 include.lowest = TRUE, labels = FALSE), by = network]
  rbindlist(lapply(split(d[module != "grey"], by = c("network", "module")), function(m) {
    net  <- m$network[1]
    pool <- d[network == net]
    need <- m[, .N, by = bin]
    by_bin <- split(pool$log2FC, pool$bin)
    null <- vapply(seq_len(N_NULL), function(k) {
      median(unlist(lapply(seq_len(nrow(need)), function(i)
        sample(by_bin[[as.character(need$bin[i])]], need$N[i]))))
    }, numeric(1))
    obs <- median(m$log2FC)
    z   <- (obs - mean(null)) / sd(null)
    p   <- min(1, 2 * min(mean(null >= obs), mean(null <= obs)) + 1 / N_NULL)
    data.table(network = net, module = m$module[1], n_members_tested = nrow(m),
               median_log2FC = obs, null_mean = mean(null), null_sd = sd(null),
               z = z, p_empirical = p,
               pct_sig_up = 100 * mean(m$FDR_all_genes < 0.05 & m$log2FC > 0, na.rm = TRUE),
               pct_sig_down = 100 * mean(m$FDR_all_genes < 0.05 & m$log2FC < 0, na.rm = TRUE))
  }))
}
per_ds <- de[, one_cell(copy(.SD)), by = .(contrast, dataset)]
per_ds[, fdr_family := "within dataset x contrast, across 29 modules"]
per_ds[, fdr := p.adjust(p_empirical, "BH"), by = .(contrast, dataset)]
setorder(per_ds, contrast, network, module, dataset)
save_tsv(per_ds, "36_sc_module_localisation_per_dataset.tsv")

# ---- 4. combine across datasets ----
comb <- per_ds[, .(n_datasets = .N,
                   datasets = paste(dataset, collapse = ";"),
                   median_members_tested = as.numeric(median(n_members_tested)),
                   mean_median_log2FC = mean(median_log2FC),
                   mean_shift_vs_null = mean(median_log2FC - null_mean),
                   stouffer_z = sum(z) / sqrt(.N),
                   n_datasets_positive = sum(z > 0),
                   n_datasets_negative = sum(z < 0)),
               by = .(contrast, network, module)]
comb[, p_stouffer := 2 * pnorm(-abs(stouffer_z))]
comb[, fdr := p.adjust(p_stouffer, "BH"), by = contrast]
comb[, fdr_family := "within contrast, across 29 modules"]
comb[, direction := fifelse(fdr < 0.05, fifelse(stouffer_z > 0, "malignant-enriched", "malignant-depleted"),
                            "not localised")]
comb <- merge(comb, PROGNOSTIC[, .(network, module, prognostic = TRUE)],
              by = c("network", "module"), all.x = TRUE)
comb[is.na(prognostic), prognostic := FALSE]
setorder(comb, contrast, -prognostic, network, -stouffer_z)
save_tsv(comb, "36_sc_module_localisation_combined.tsv")
print(comb[prognostic == TRUE, .(contrast, network, module, n_datasets, median_members_tested,
                                 mean_median_log2FC = round(mean_median_log2FC, 2),
                                 stouffer_z = round(stouffer_z, 2), fdr = signif(fdr, 2),
                                 n_datasets_positive, direction)])

# ---- 5. coverage ----
cov <- merge(mods[module != "grey", .(module_size = .N), by = .(network, module)],
             de[module != "grey", .(tested = uniqueN(gene_key)), by = .(contrast, network, module)],
             by = c("network", "module"), all.x = TRUE)
cov[is.na(tested), tested := 0L]
cov[, pct_tested_any_dataset := round(100 * tested / module_size, 1)]
save_tsv(cov, "36_sc_module_gene_coverage.tsv")

# ---- 6. descriptive pooled meta-analysis summary for the lncRNA modules ----
meta <- fread(file.path(SC_DIR, "08_meta_lncRNA.tsv.gz"))[contrast %in% CONTRASTS]
meta <- merge(meta[, .(contrast, gene_key, pooled_log2FC, pooled_FDR)],
              mods[network == "lncRNA"], by = "gene_key")
msum <- meta[module != "grey", .(n_tested = .N, median_pooled_log2FC = median(pooled_log2FC),
                                 n_sig_up = sum(pooled_FDR < 0.05 & pooled_log2FC > 0),
                                 n_sig_down = sum(pooled_FDR < 0.05 & pooled_log2FC < 0),
                                 sig_up_members = paste(gene_name[pooled_FDR < 0.05 & pooled_log2FC > 0],
                                                        collapse = ";")),
             by = .(contrast, module)]
setorder(msum, contrast, -median_pooled_log2FC)
save_tsv(msum, "36_sc_module_lncRNA_pooled_meta_summary.tsv")

# ---- 7. figure ----
fd <- copy(comb)
fd[, label := paste0(fifelse(network == "lncRNA", "lncRNA ", "PC "), module)]
fd[, contrast_lab := fifelse(contrast == CONTRASTS[1], "vs normal epithelium", "vs tumour microenvironment")]
ord <- fd[contrast == CONTRASTS[2]][order(mean_shift_vs_null), label]
fd[, label := factor(label, levels = ord)]
p <- ggplot(fd, aes(mean_shift_vs_null, label, colour = prognostic, shape = fdr < 0.05)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 2.2) +
  facet_wrap(~contrast_lab, nrow = 1) +
  scale_colour_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "grey55"),
                      labels = c(`TRUE` = "prognostic module", `FALSE` = "other module"), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR >= 0.05"), name = NULL) +
  labs(x = "Median log2 fold change above the abundance-matched null, mean over cohorts", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
save_fig(p, "SupplementaryFigureS10_sc_module_localisation", width = 7.2, height = 5.2)

write_session_info("36_singlecell_module_localisation")
msg("36 done")
