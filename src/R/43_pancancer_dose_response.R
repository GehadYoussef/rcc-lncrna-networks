# 43_pancancer_dose_response.R: does axis strength scale with the spread of library quality?
# Post hoc analysis, outside the decision rules in
# results/38_decision_rules_lock.json. Every output row is labelled as such.
# Four spread measures of the non-feature fraction per eligible tumour project (SD, IQR,
# coefficient of variation, MAD of its log), each correlated across projects (Spearman)
# with the lncRNA PC1 |rho|, the lncRNA minus protein-coding |rho| margin and
# plate-explained variance (39), and with the largest-module share and ARI (40).
# Run after 39 and 40.
# Outputs: 43_pancancer_dose_response_projects.tsv (one row per project),
#          43_pancancer_dose_response.tsv (one row per spread measure x outcome).

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages(library(data.table))
banner("43 | Pan-cancer dose-response (post hoc)")

ax <- fread(file.path(RESULTS_DIR, "39_pancancer_axis.tsv"))[group == "tumour" & eligible == TRUE]
nw <- fread(file.path(RESULTS_DIR, "40_pancancer_network.tsv"))
sp <- rbindlist(lapply(ax$project, function(p) {
  s <- readRDS(file.path(CACHE_DIR, paste0("pan_", p, ".rds")))$samples[keep == TRUE & group == "tumour"]
  x <- s$pct_noFeature
  data.table(project = p, n = length(x), sd_metric = sd(x), iqr_metric = IQR(x),
             cv_metric = sd(x) / mean(x), mad_log_metric = mad(log(x)))
}))
d <- merge(merge(ax[, .(project, abs_rho_lnc = abs(rho_lnc_pc1_noFeature),
                        margin = abs(rho_lnc_pc1_noFeature) - abs(rho_pc_pc1_noFeature),
                        r2_plate_metric)], sp, by = "project"),
           nw[, .(project, largest_frac_observed, ari_all_genes, rule_C_pass)], by = "project")
d[, analysis := "post hoc (not a locked rule)"]
save_tsv(d[order(-mad_log_metric)], "43_pancancer_dose_response_projects.tsv")

spread <- c("mad_log_metric", "cv_metric", "iqr_metric", "sd_metric")
outc <- c(abs_rho_lnc = "|rho| of the lncRNA PC1 with the metric",
          margin = "lncRNA minus protein-coding |rho|",
          largest_frac_observed = "largest observed module, share of genes",
          ari_all_genes = "adjusted Rand index, observed vs adjusted",
          r2_plate_metric = "metric variance explained by plate")
res <- rbindlist(lapply(spread, function(v) rbindlist(lapply(names(outc), function(o) {
  ok <- is.finite(d[[v]]) & is.finite(d[[o]])
  ct <- cor.test(d[[v]][ok], d[[o]][ok], method = "spearman", exact = FALSE)
  data.table(spread_measure = v, outcome = outc[[o]], n_projects = sum(ok),
             spearman_rho = unname(ct$estimate), p = ct$p.value)
}))))
res[, fdr := p.adjust(p, "BH")]
res[, fdr_family := "all spread measures x outcomes (20 tests)"]
res[, analysis := "post hoc (not a locked rule)"]
save_tsv(res, "43_pancancer_dose_response.tsv"); print(res[spread_measure == "mad_log_metric"])
write_session_info("43_pancancer_dose_response")
msg("43 done")
