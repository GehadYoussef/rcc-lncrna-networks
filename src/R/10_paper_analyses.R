# 10_paper_analyses.R: validation summaries and lncRNA module characterisation.
#
#   A. Bootstrap CI on the CPTAC-3 delta C-index, tabulated from 09.
#   B. Time-dependent AUC in CPTAC-3, comparator alone and with module eigengenes.
#   D. Kaplan-Meier in CPTAC-3 with risk tertiles cut in TCGA.
#   E. lncRNA module members, transcript classes and guilt-by-association GO.
# Uses the locked model from 09 without refitting.
# Inputs: caches locked_model.rds, networks.rds, survival.rds, dataset.rds and
#   results 09_delta_cindex_validation.tsv.
# Outputs: results 10_*.tsv and figures 10_validation_KM_tertiles_<comparator>.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(survminer)
  library(timeROC); library(ggplot2)
  library(clusterProfiler); library(org.Hs.eg.db)
})
banner("10 | Paper-grade analyses")
set.seed(SEED)

L    <- readRDS(file.path(CACHE_DIR, "locked_model.rds"))
stopifnot(identical(L$version, "v9"))
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
surv <- readRDS(file.path(CACHE_DIR, "survival.rds"))
ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))

COMPARATORS <- c("clinical", "augmented")
lpv <- function(X, b) as.numeric(X[, names(b), drop = FALSE] %*% b)

# Locked linear predictors per comparator: comparator alone and comparator +
# module eigengenes, in discovery (t) and CPTAC-3 (v).
LP <- lapply(setNames(COMPARATORS, COMPARATORS), function(cmp) {
  X <- L$X[[cmp]]; b <- L$models[[cmp]]
  list(clin_v = lpv(X$Xv, b$b_clin),
       full_v = lpv(cbind(X$Xv, L$Ev), b$b_full),
       clin_t = lpv(X$Xt, b$b_clin),
       full_t = lpv(cbind(X$Xt, L$Et), b$b_full))
})
nv <- nrow(L$vcl); ev_v <- sum(L$vcl$os_event)
msg("CPTAC-3 evaluated set: n = ", nv, ", events = ", ev_v)

# -----------------------------------------------------------------------------
# A. Bootstrap CI on the validation increment (from 09)
# -----------------------------------------------------------------------------
banner("A | Bootstrap delta C in the validation cohort (from 09)")
d09 <- fread(file.path(RESULTS_DIR, "09_delta_cindex_validation.tsv"))
sel <- d09[(cohort %like% "validation" | cohort %like% "CPTAC") &
           standardisation == "cohort"]
stopifnot(nrow(sel) >= 1, all(COMPARATORS %in% sel$comparator))
boot_tbl <- sel[, .(comparator,
                    delta_C_validation = round(delta_C, 4),
                    lo = round(lo, 4), hi = round(hi, 4),
                    p_bootstrap = signif(p_boot, 3),
                    prop_resamples_positive = round(prop_positive, 3),
                    n_resamples = if ("n_boot" %in% names(d09)) n_boot else BOOT_B)]
boot_tbl <- boot_tbl[order(match(comparator, COMPARATORS))]
save_tsv(boot_tbl, "10_bootstrap_delta_cindex.tsv")
print(boot_tbl)

# -----------------------------------------------------------------------------
# B. Time-dependent AUC in validation
# -----------------------------------------------------------------------------
banner("B | Time-dependent AUC, validation cohort")
times <- EVAL_TIMES_YRS * 365.25
auc <- rbindlist(lapply(COMPARATORS, function(cmp) {
  rbindlist(lapply(c("comparator only", "+ module eigengenes"), function(nm) {
    lp <- if (nm == "comparator only") LP[[cmp]]$clin_v else LP[[cmp]]$full_v
    tr <- timeROC(T = L$vcl$os_time, delta = L$vcl$os_event, marker = lp,
                  cause = 1, times = times, iid = TRUE)
    data.table(comparator = cmp, model = nm, years = EVAL_TIMES_YRS,
               n = nv, events = ev_v,
               AUC = as.numeric(tr$AUC),
               se = as.numeric(tr$inference$vect_sd_1))
  }))
}))
auc[, `:=`(lo = AUC - 1.96 * se, hi = AUC + 1.96 * se)]
save_tsv(auc, "10_validation_time_auc.tsv")
print(auc)

# -----------------------------------------------------------------------------
# D. Kaplan-Meier in validation, cut-points fixed in TCGA
# -----------------------------------------------------------------------------
banner("D | Kaplan-Meier, validation, TCGA-defined tertiles")
km_tbl <- rbindlist(lapply(COMPARATORS, function(cmp) {
  cuts <- quantile(LP[[cmp]]$full_t, c(1/3, 2/3))
  grp  <- cut(LP[[cmp]]$full_v, breaks = c(-Inf, cuts, Inf),
              labels = c("Low", "Intermediate", "High"))
  kmd <- data.frame(os_time = L$vcl$os_time, os_event = L$vcl$os_event, Risk = grp)
  print(table(kmd$Risk))
  fit <- survfit(Surv(os_time, os_event) ~ Risk, data = kmd)
  lr  <- survdiff(Surv(os_time, os_event) ~ Risk, data = kmd)
  p_lr <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  msg(cmp, " comparator: validation log-rank across TCGA-defined risk tertiles: p = ",
      signif(p_lr, 3))
  save_fig_base(paste0("10_validation_KM_tertiles_", cmp), 6.3, 6.7,
    print(ggsurvplot(fit, data = kmd, pval = TRUE, risk.table = TRUE,
                     conf.int = TRUE,
                     palette = c("#2166AC", "#999999", "#B2182B"),
                     xlab = "Days from diagnosis", ylab = "Overall survival",
                     legend.title = "Risk group (TCGA cut-points)",
                     title = paste0("CPTAC-3 validation: locked risk score (",
                                    cmp, " comparator + modules)"))))
  data.table(comparator = cmp, group = names(table(kmd$Risk)),
             n = as.integer(table(kmd$Risk)),
             events = as.integer(tapply(kmd$os_event, kmd$Risk, sum)),
             logrank_p = signif(p_lr, 3))
}))
save_tsv(km_tbl, "10_validation_km_tertiles.tsv")

# -----------------------------------------------------------------------------
# E. Characterisation of the prognostic lncRNA modules
# -----------------------------------------------------------------------------
banner("E | lncRNA module characterisation")
sig_l <- as.data.table(surv$lnc)[fdr_full < FDR_ALPHA, module]
gt_l  <- as.data.table(nets$lnc$gene_tbl)
ann_m <- as.data.table(ds$mrna$ann)

# Transcript class from GENCODE gene-name conventions.
classify <- function(g) fifelse(grepl("-AS[0-9]*$", g), "antisense",
                        fifelse(grepl("-DT$", g), "divergent",
                        fifelse(grepl("-IT[0-9]*$", g), "intronic",
                        fifelse(grepl("^LINC", g), "lincRNA (named)",
                        fifelse(grepl("^(AC|AL|AP|Z|BX|CT|FP)[0-9]", g),
                                "clone-based (unnamed)", "other")))))
memb <- gt_l[module %in% sig_l]
memb[, class := classify(gene_name)]
save_tsv(memb[order(module, -kME),
              .(module, gene_name, gene_id, kME = round(kME, 3), class)],
         "10_lncRNA_module_members.tsv")
cls <- dcast(memb[, .N, by = .(module, class)], module ~ class,
             value.var = "N", fill = 0)
save_tsv(cls, "10_lncRNA_module_composition.tsv")
print(cls)

# Guilt-by-association, since lncRNAs have little GO annotation: GO BP
# enrichment of the protein-coding genes most correlated with each discovery
# module score. Partners are ranked on the observed and on the residualised
# coding matrix, at 100, 300 and 500 partners. The universe is the analysed
# protein-coding set.
S_l  <- discovery_scores(nets)$lnc
mats <- list(observed = obs_expr(nets$mrna), residualised = nets$mrna$expr)
N_PARTNERS <- c(100L, 300L, 500L)
gba <- list(); partners <- list(); gba_summary <- list()
for (mt in names(mats)) {
  Em <- mats[[mt]]
  common   <- intersect(rownames(Em), rownames(S_l))
  universe <- sub("\\..*$", "", colnames(Em))
  for (m in sig_l) {
    me  <- S_l[common, paste0("lnc_ME", m)]
    r   <- cor(Em[common, ], me, use = "pairwise.complete.obs")[, 1]
    ord <- names(sort(abs(r), decreasing = TRUE))
    kmax <- max(N_PARTNERS); top_all <- ord[seq_len(kmax)]
    partners[[paste(mt, m)]] <- data.table(
      module = m, matrix = mt, rank = seq_len(kmax), gene_id = top_all,
      gene_name = ann_m$gene_name[match(top_all, ann_m$gene_id)],
      r = round(r[top_all], 3))
    msg("lncRNA module ", m, " [", mt, "]: n = ", length(common), " samples; ",
        sum(abs(r) > 0.5), " coding genes with |r| > 0.5, ",
        sum(abs(r) > 0.3), " with |r| > 0.3")
    for (k in N_PARTNERS) {
      top <- ord[seq_len(k)]
      eg <- tryCatch(enrichGO(gene = sub("\\..*$", "", top), universe = universe,
                              OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
                              pAdjustMethod = "BH", pvalueCutoff = 0.05,
                              qvalueCutoff = 0.10, readable = TRUE),
                     error = function(e) NULL)
      egd <- if (is.null(eg)) data.frame() else as.data.frame(eg)
      gba_summary[[paste(mt, m, k)]] <- data.table(
        module = m, matrix = mt, n_partners = k, n_samples = length(common),
        min_abs_r = round(min(abs(r[top])), 3),
        n_go_terms = nrow(egd),
        top_term = if (nrow(egd)) egd$Description[1] else NA_character_,
        top_p_adjust = if (nrow(egd)) signif(egd$p.adjust[1], 3) else NA_real_)
      msg("  ", k, " partners (min |r| = ", round(min(abs(r[top])), 3), "): ",
          if (nrow(egd) == 0) "no enriched GO BP terms"
          else paste0(nrow(egd), " GO BP terms, top: ", egd$Description[1]))
      if (nrow(egd) > 0) {
        d <- as.data.table(egd)[1:min(10, .N)]
        d[, `:=`(module = m, matrix = mt, n_partners = k)]
        gba[[paste(mt, m, k)]] <- d
      }
    }
  }
}
save_tsv(rbindlist(partners), "10_lncRNA_module_partners.tsv")
save_tsv(rbindlist(gba_summary), "10_lncRNA_guilt_by_association_summary.tsv")
if (length(gba)) {
  g_all <- rbindlist(gba, fill = TRUE)
  save_tsv(g_all[, .(module, matrix, n_partners, Description, GeneRatio,
                     p.adjust = signif(p.adjust, 3), Count)],
           "10_lncRNA_guilt_by_association_GO.tsv")
  print(g_all[, .(module, matrix, n_partners, Description,
                  p.adjust = signif(p.adjust, 3))])
}

# Summary table: discovery effect and CPTAC-3 replication under the augmented
# and parsimonious models from 09.
rep_cols <- intersect(c("module", "HR_cptac", "lo", "hi", "p_cptac", "fdr_cptac",
                        "same_direction", "n", "events"),
                      names(L$replication))
hdr <- merge(as.data.table(surv$lnc)[fdr_full < FDR_ALPHA,
               .(module, n_genes, HR_tcga = round(HR_full, 3),
                 fdr_tcga = signif(fdr_full, 3))],
             if (!is.null(L$replication))
               as.data.table(L$replication)[biotype == "lncRNA", ..rep_cols]
             else data.table(module = character()),
             by = "module", all.x = TRUE)
if (!is.null(L$replication_parsimonious)) {
  rp <- as.data.table(L$replication_parsimonious)[biotype == "lncRNA"]
  pcols <- intersect(c("HR_cptac", "lo", "hi", "p_cptac"), names(rp))
  rp <- rp[, c("module", pcols), with = FALSE]
  setnames(rp, pcols, paste0(pcols, "_parsimonious"))
  hdr <- merge(hdr, rp, by = "module", all.x = TRUE)
}
save_tsv(hdr, "10_lncRNA_headline_table.tsv")
print(hdr)

write_session_info("10_paper_analyses")
banner("10 | done")
