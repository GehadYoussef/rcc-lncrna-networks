# 08_meta_analysis.R: cross-dataset replication and meta-analysis (run after 07).
# For each contrast and lncRNA: random-effects meta-analysis (metafor REML, DL fallback
# recorded) of inferential estimates from datasets in replication$pooled_roles, with BH
# FDR, tau^2, Q and I^2, leave-one-dataset-out range, replication direction counts
# outside the discovery dataset, and direction agreement in descriptive datasets (not
# pooled). No outcome data are read. The candidate rules are then evaluated and a
# provisional ranking written. The candidate lock is written by 09.
# Outputs: results/singlecell/08_meta_lncRNA.tsv.gz, 08_replication_summary.tsv,
#          08_candidate_evaluation.tsv, 08_candidate_rule_attrition.tsv,
#          figures/singlecell/08_forest_<contrast>

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
suppressPackageStartupMessages({ library(ggplot2) })
if (!requireNamespace("metafor", quietly = TRUE)) stop("metafor is required for 08")
banner("08  cross-dataset replication and meta-analysis")
start_log("08_meta_analysis")

MA <- CFG$meta_analysis; RP <- CFG$replication; CR <- CFG$candidate_rules; W <- CFG$candidate_ranking_weights
man <- read_manifest()$manifest
disc_ds <- man[analysis_role == "discovery", dataset_name][1]
de <- read_tsv("07_de_lncRNA.tsv.gz")
assert_no_outcome_fields(de, "meta-analysis input")
msg("discovery dataset: ", disc_ds, "; methods: ", MA$method)

meta_rows <- list(); rep_summ <- list()
for (cn in names(CFG$de$contrasts)) {
  inf <- de[contrast == cn & inference == "inferential" & analysis_role %in% unlist(RP$pooled_roles)]
  dsc <- de[contrast == cn & inference == "descriptive" & analysis_role %in% unlist(RP$descriptive_roles)]
  if (!nrow(inf)) { msg(cn, ": no inferential datasets"); next }
  k_by <- inf[, .N, by = gene_key]
  keys <- k_by[N >= MA$min_datasets, gene_key]
  datasets_used <- sort(unique(inf$dataset))
  msg(cn, ": ", length(datasets_used), " inferential datasets (", paste(datasets_used, collapse = ", "), "); ",
      length(keys), " lncRNAs with k >= ", MA$min_datasets)
  if (!length(keys)) next
  do_loo <- cn %in% unlist(MA$leave_one_out_contrasts)
  split_inf <- split(inf[gene_key %in% keys], by = "gene_key")
  split_dsc <- split(dsc, by = "gene_key")

  res <- rbindlist(lapply(names(split_inf), function(k) {
    x <- split_inf[[k]]
    m <- meta_one(x$log2FC, x$SE, MA$method)
    loo <- if (do_loo) meta_loo(x$log2FC, x$SE, MA$method) else list(min_est = NA_real_, max_est = NA_real_, max_p = NA_real_)
    dsc_x <- split_dsc[[k]]
    d <- x[dataset == disc_ds]
    ref_sign <- if (nrow(d)) sign(d$log2FC) else sign(m$estimate)
    r <- x[dataset != disc_ds]
    data.table(contrast = cn, gene_key = k, symbol = x$symbol[1], k = m$k, datasets = paste(x$dataset, collapse = ";"),
               pooled_log2FC = m$estimate, pooled_SE = m$se, pooled_CI_low = m$ci_low, pooled_CI_high = m$ci_high,
               pooled_p = m$pval, tau2 = m$tau2, Q = m$Q, Q_p = m$Q_pval, I2 = m$I2, meta_method = m$method,
               loo_min_log2FC = loo$min_est, loo_max_log2FC = loo$max_est, loo_max_p = loo$max_p,
               discovery_tested = nrow(d) > 0,
               discovery_log2FC = if (nrow(d)) d$log2FC else NA_real_,
               discovery_FDR = if (nrow(d)) d$FDR_lncRNA else NA_real_,
               discovery_detection_malignant = if (nrow(d)) d$donor_detection_test else NA_real_,
               direction_reference = if (nrow(d)) "discovery" else "pooled",
               repl_tested = nrow(r), repl_same_direction = sum(sign(r$log2FC) == ref_sign),
               repl_same_direction_fdr = sum(sign(r$log2FC) == ref_sign & r$FDR_lncRNA < CFG$de$fdr, na.rm = TRUE),
               descriptive_tested = if (is.null(dsc_x)) 0L else nrow(dsc_x),
               descriptive_same_direction = if (is.null(dsc_x)) 0L else sum(sign(dsc_x$log2FC) == ref_sign))
  }))
  res[, pooled_FDR := p.adjust(pooled_p, "BH")]
  res[, repl_fraction := fifelse(repl_tested > 0, repl_same_direction / repl_tested, NA_real_)]
  res[, heterogeneity_flag := !is.na(I2) & I2 > MA$i2_flag]
  res[, loo_direction_stable := fifelse(k >= 3, sign(loo_min_log2FC) == sign(loo_max_log2FC) &
                                          sign(loo_min_log2FC) == sign(pooled_log2FC), NA)]
  meta_rows[[cn]] <- res
  rep_summ[[cn]] <- res[, .(contrast = cn, datasets_pooled = paste(datasets_used, collapse = ";"), lnc_meta_tested = .N,
                            up_fdr_lfc = sum(pooled_FDR < MA$fdr & pooled_log2FC >= MA$min_abs_log2fc, na.rm = TRUE),
                            down_fdr_lfc = sum(pooled_FDR < MA$fdr & pooled_log2FC <= -MA$min_abs_log2fc, na.rm = TRUE),
                            up_replicated = sum(pooled_FDR < MA$fdr & pooled_log2FC >= MA$min_abs_log2fc &
                                                  repl_same_direction >= RP$min_datasets_concordant &
                                                  repl_fraction >= RP$min_fraction_concordant, na.rm = TRUE),
                            heterogeneity_flagged = sum(heterogeneity_flag),
                            reml_fallbacks = sum(meta_method == "DL_fallback", na.rm = TRUE))]
}
M <- rbindlist(meta_rows)
save_tsv(M, "08_meta_lncRNA.tsv.gz")
S <- rbindlist(rep_summ)
save_tsv(S, "08_replication_summary.tsv")
print(S)

# ---- candidate rule evaluation (no lock) ---------------------------------------------
passes <- function(x) {
  x[, .(gene_key,
        pass_effect = !is.na(pooled_FDR) & pooled_FDR < MA$fdr & pooled_log2FC >= MA$min_abs_log2fc,
        pass_replication = !is.na(repl_fraction) & repl_same_direction >= RP$min_datasets_concordant &
          repl_fraction >= RP$min_fraction_concordant & direction_reference == "discovery" & discovery_log2FC > 0,
        pass_loo = if (isTRUE(CR$require_loo_direction_stable)) is.na(loo_direction_stable) | loo_direction_stable else TRUE,
        pooled_log2FC, pooled_FDR, I2, heterogeneity_flag, k, repl_fraction, repl_same_direction, repl_tested,
        discovery_log2FC, discovery_FDR, discovery_detection_malignant)]
}
A <- passes(M[contrast == CR$tumour_specificity_contrast]); setnames(A, names(A)[-1], paste0("A_", names(A)[-1]))
B <- passes(M[contrast == CR$compartment_contrast]);       setnames(B, names(B)[-1], paste0("B_", names(B)[-1]))
E <- merge(B, A, by = "gene_key", all = TRUE)
sym <- unique(M[, .(gene_key, symbol)], by = "gene_key")
E <- merge(sym, E, by = "gene_key", all.y = TRUE)

disc_pass_f <- file.path(DERIVED_DIR, paste0(disc_ds, ".lnc_malignant_passing.rds"))
disc_pass <- if (file.exists(disc_pass_f)) readRDS(disc_pass_f)[passing == TRUE, unique(gene_key)] else character()
bulk <- if (file.exists(file.path(RESULTS_DIR, "04_bulk_lncRNA_features.tsv"))) read_tsv("04_bulk_lncRNA_features.tsv") else data.table()
bulk_both <- if (nrow(bulk)) bulk[bulk_expressed == TRUE, .N, by = key][N == 2, key] else character()

E[, discovery_donor_filter := gene_key %in% disc_pass]
E[, bulk_expressed_both := gene_key %in% bulk_both]
E[, A_testable := !is.na(A_k) & A_k >= MA$min_datasets]
# Where the normal-epithelium contrast was tested in fewer than min_datasets
# datasets, a single inferential estimate significantly lower in malignant
# cells excludes the gene from tier 2.
A_single <- de[contrast == CR$tumour_specificity_contrast & inference == "inferential" &
                 analysis_role %in% unlist(RP$pooled_roles)]
A_contrary <- A_single[FDR_lncRNA < CFG$de$fdr & log2FC <= -CFG$de$min_abs_log2fc, unique(gene_key)]
E[, A_single_dataset_contrary := !A_testable & gene_key %in% A_contrary]
E[, pass_B := B_pass_effect %in% TRUE & B_pass_replication %in% TRUE & B_pass_loo %in% TRUE]
E[, pass_A := A_pass_effect %in% TRUE & A_pass_replication %in% TRUE & A_pass_loo %in% TRUE]
E[, pass_common := (!isTRUE(CR$require_discovery_donor_filter) | discovery_donor_filter) &
                   (!isTRUE(CR$require_bulk_expressed_both) | bulk_expressed_both)]
E[, tier := fifelse(pass_B & pass_common & A_testable & pass_A, "tier1_tumour_specific",
            fifelse(pass_B & pass_common & !A_testable & !A_single_dataset_contrary,
                    "tier2_malignant_compartment_enriched", "not_selected"))]
E[, first_failed_rule := fifelse(tier != "not_selected", "",
    fifelse(!(B_pass_effect %in% TRUE), "B_pooled_effect_or_FDR",
    fifelse(!(B_pass_replication %in% TRUE), "B_replication_direction",
    fifelse(!(B_pass_loo %in% TRUE), "B_leave_one_dataset_out",
    fifelse(isTRUE(CR$require_discovery_donor_filter) & !discovery_donor_filter, "discovery_donor_expression",
    fifelse(isTRUE(CR$require_bulk_expressed_both) & !bulk_expressed_both, "bulk_expression_TCGA_CPTAC",
    fifelse(A_single_dataset_contrary, "A_single_dataset_lower_than_normal_epithelium",
    fifelse(A_testable & !(A_pass_effect %in% TRUE), "A_tested_not_enriched_vs_normal_epithelium",
    fifelse(A_testable & !(A_pass_replication %in% TRUE), "A_replication_direction",
    fifelse(A_testable & !(A_pass_loo %in% TRUE), "A_leave_one_dataset_out", "other"))))))))))]

# provisional ranking within selected lncRNAs (weights from config)
sel <- E$tier != "not_selected"
s01 <- function(x) { r <- range(x, na.rm = TRUE); if (!is.finite(diff(r)) || diff(r) == 0) rep(1, length(x)) else (x - r[1]) / diff(r) }
E[, provisional_score := NA_real_]
if (any(sel)) E[sel, provisional_score := W$meta_effect * s01(B_pooled_log2FC) + W$replication_fraction * B_repl_fraction +
                  W$malignant_detection * fifelse(is.na(B_discovery_detection_malignant), 0, B_discovery_detection_malignant)]
setorder(E, tier, -provisional_score, na.last = TRUE)
E[, provisional_rank := if (tier[1] != "not_selected") seq_len(.N) else NA_integer_, by = tier]
save_tsv(E, "08_candidate_evaluation.tsv")

att <- E[, .N, by = .(tier, first_failed_rule)][order(tier, -N)]
save_tsv(att, "08_candidate_rule_attrition.tsv")
print(att)
msg("tier1 (tumour-specific): ", E[tier == "tier1_tumour_specific", .N],
    "; tier2 (malignant-compartment enriched, no normal comparator): ", E[tier == "tier2_malignant_compartment_enriched", .N])

# ---- forest plots for top provisional candidates ------------------------------------
de_all <- read_tsv("07_de_lncRNA.tsv.gz")
for (cn in unlist(MA$leave_one_out_contrasts)) {
  top <- head(E[tier != "not_selected"][order(tier, provisional_rank), gene_key], 15)
  if (!length(top)) next
  per <- de_all[contrast == cn & gene_key %in% top & inference == "inferential",
                .(gene_key, symbol, dataset, log2FC, lo = CI_low, hi = CI_high, kind = "dataset")]
  pool <- M[contrast == cn & gene_key %in% top, .(gene_key, symbol, dataset = "Pooled (REML)", log2FC = pooled_log2FC,
                                                   lo = pooled_CI_low, hi = pooled_CI_high, kind = "pooled")]
  pd <- rbind(per, pool)
  if (!nrow(pd)) next
  pd[, symbol := factor(symbol, levels = rev(unique(E[gene_key %in% top][match(top, gene_key), symbol])))]
  p <- ggplot(pd, aes(x = log2FC, y = dataset, colour = kind)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y") + geom_point(aes(shape = kind), size = 1.6) +
    facet_wrap(~ symbol, ncol = 3) + scale_colour_manual(values = c(dataset = "grey30", pooled = "#b2182b")) +
    labs(x = "log2 fold change (malignant vs comparator), 95% CI", y = NULL, colour = NULL, shape = NULL,
         title = paste("Donor-level effects and random-effects pooling:", cn),
         subtitle = "Top provisional candidates (not locked); each dataset estimate uses independent donors") +
    theme_bw(base_size = 7) + theme(legend.position = "bottom")
  save_fig(p, paste0("08_forest_", cn), width = 9, height = 10)
}
write_session_info("08_meta_analysis")
