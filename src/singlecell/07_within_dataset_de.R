# 07_within_dataset_de.R: within-dataset differential expression (run after 06).
# edgeR quasi-likelihood on donor x compartment pseudobulk counts, one model per
# dataset x contrast (de$contrasts). Datasets are pooled only in 08.
# Inferential when both groups have >= min_donors_inferential donors: paired
# (~ donor + group) when >= min_paired_donors donors have both groups, else
# unpaired, with TMM, filterByExpr and robust QL dispersion. Otherwise descriptive:
# log2 ratio of mean log-CPM and donor detection, without p-values.
# Genes not mapping uniquely to GENCODE v36 are fitted but not reported.
# Outputs: results/singlecell/07_de_all_genes.tsv.gz, 07_de_lncRNA.tsv.gz,
#          07_contrast_summary.tsv, figures/singlecell/07_de_summary

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
suppressPackageStartupMessages({ library(edgeR); library(ggplot2) })
banner("07  within-dataset differential expression")
start_log("07_within_dataset_de")

DE <- CFG$de
MIN_DONORS <- CFG$pseudobulk$min_donors_inferential
BIOTYPES <- CFG$annotation$lncrna_biotypes
man <- read_manifest()$manifest
pb_files <- list.files(DERIVED_DIR, pattern = "[.]pseudobulk[.]rds$", full.names = TRUE)
if (!length(pb_files)) stop("no pseudobulk objects: run 06 first")

group_stats <- function(cpm, st, levels) {
  one <- function(lvl, tag) {
    x <- cpm[, st$grp == lvl, drop = FALSE]
    out <- data.table(rowMeans(x), apply(x, 1, stats::median), rowMeans(x >= DE$detection_cpm))
    setnames(out, paste0(c("mean_cpm_", "median_cpm_", "donor_detection_"), tag))
    out
  }
  cbind(one(levels[1], "ref"), one(levels[2], "test"))
}

fit_contrast <- function(y_counts, st, levels, paired) {
  X <- if (paired) model.matrix(~ factor(donor_id) + grp, data = st) else model.matrix(~ grp, data = st)
  y <- DGEList(y_counts, group = st$grp)
  keep <- filterByExpr(y, design = X, min.count = DE$filter_min_count)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, X)
  fit <- glmQLFit(y, X, robust = isTRUE(DE$robust))
  res <- glmQLFTest(fit, coef = ncol(X))
  list(y = y, keep = keep, table = res$table, df = res$df.total, design = X)
}

all_rows <- list(); summ <- list()
for (f in pb_files) {
  PB <- readRDS(f); ds <- PB$dataset
  role <- man[dataset_name == ds, analysis_role]
  counts <- PB$counts; st0 <- as.data.table(PB$samples); g <- as.data.table(PB$genes)
  require_raw_counts(counts, paste(ds, "pseudobulk"))
  assert_no_outcome_fields(st0, paste(ds, "pseudobulk samples"))

  # reportable genes: unique mapping to a v36 gene
  dup_key <- g$key %in% g$key[duplicated(g$key) & !is.na(g$key)]
  reportable <- !is.na(g$key) & !dup_key
  msg(ds, ": ", sum(!reportable), " matrix rows not reportable (unmapped or duplicated key)")

  for (cn in names(DE$contrasts)) {
    levels <- unlist(DE$contrasts[[cn]])
    st <- st0[compartment %in% levels]
    n_by <- st[, uniqueN(donor_id), by = compartment]
    n_ref <- n_by[compartment == levels[1], V1]; n_test <- n_by[compartment == levels[2], V1]
    n_ref <- if (length(n_ref)) n_ref else 0L; n_test <- if (length(n_test)) n_test else 0L
    base <- data.table(dataset = ds, analysis_role = role, contrast = cn, n_donors_ref = n_ref, n_donors_test = n_test)
    if (n_ref == 0 || n_test == 0) { summ[[paste(ds, cn)]] <- base[, inference := "not_evaluable"]; next }

    d <- design_for_contrast(st, levels = levels, min_paired = DE$min_paired_donors)
    st <- d$table
    cnt <- counts[, match(paste(st$donor_id, st$compartment, sep = "|"), colnames(counts)), drop = FALSE]
    lib <- Matrix::colSums(cnt)
    cpm_all <- t(t(as.matrix(cnt)) / lib * 1e6)
    gs <- group_stats(cpm_all, st, levels)
    inferential <- min(d$n_group) >= MIN_DONORS

    if (inferential) {
      fc <- fit_contrast(as.matrix(cnt), st, levels, d$paired)
      tab <- as.data.table(fc$table, keep.rownames = "row_id")
      tab[, row := which(fc$keep)]
      tab[, `:=`(log2FC = logFC, SE = se_from_f(logFC, F), PValue = PValue)]
      tab[, `:=`(CI_low = log2FC - qt(0.975, fc$df) * SE, CI_high = log2FC + qt(0.975, fc$df) * SE)]
      # leave-one-donor-out direction stability
      tab[, lodo_direction_flip := NA]
      if (isTRUE(DE$leave_one_donor_out)) {
        flips <- rep(FALSE, nrow(tab))
        for (dn in unique(st$donor_id)) {
          sti <- st[donor_id != dn]
          if (min(table(sti$grp)) < 2) next
          pairedi <- d$paired && sti[, uniqueN(grp), by = donor_id][V1 == 2, .N] >= 2
          if (d$paired && !pairedi) next
          ri <- tryCatch(fit_contrast(as.matrix(cnt[, st$donor_id != dn, drop = FALSE]), sti, levels, d$paired),
                         error = function(e) NULL)
          if (is.null(ri)) next
          li <- ri$table$logFC[match(tab$row, which(ri$keep))]
          flips <- flips | (!is.na(li) & sign(li) != sign(tab$log2FC))
        }
        tab[, lodo_direction_flip := flips]
      }
    } else {
      pc <- DE$descriptive_prior_count
      lc <- log2(t(t(as.matrix(cnt) + pc) / (lib + 2 * pc) * 1e6))
      mref <- rowMeans(lc[, st$grp == levels[1], drop = FALSE]); mtest <- rowMeans(lc[, st$grp == levels[2], drop = FALSE])
      tab <- data.table(row = seq_len(nrow(cnt)), log2FC = mtest - mref, SE = NA_real_, CI_low = NA_real_,
                        CI_high = NA_real_, PValue = NA_real_, lodo_direction_flip = NA)
      tab <- tab[Matrix::rowSums(cnt) >= DE$filter_min_count]
    }

    # paired donors: fraction whose own log2 ratio has the pooled direction
    tab[, paired_direction_consistency := NA_real_]
    if (d$paired) {
      pc <- DE$descriptive_prior_count
      lc <- log2(t(t(as.matrix(cnt) + pc) / (lib + 2 * pc) * 1e6))
      dn <- unique(st$donor_id)
      per <- sapply(dn, function(x) lc[, st$donor_id == x & st$grp == levels[2]] - lc[, st$donor_id == x & st$grp == levels[1]])
      tab[, paired_direction_consistency := rowMeans(sign(per[row, , drop = FALSE]) == sign(log2FC))]
    }

    tab <- cbind(tab, gs[tab$row])
    tab <- tab[reportable[row]]
    tab[, `:=`(gene_key = g$key[row], symbol = g$ref_symbol[row], gene_type = g$gene_type[row],
               is_lncRNA = is_lncrna(g$gene_type[row], BIOTYPES))]
    tab[, FDR_all_genes := if (inferential) p.adjust(PValue, "BH") else NA_real_]
    tab[, FDR_lncRNA := NA_real_]
    if (inferential) tab[is_lncRNA == TRUE, FDR_lncRNA := p.adjust(PValue, "BH")]
    tab[, `:=`(dataset = ds, analysis_role = role, contrast = cn, inference = if (inferential) "inferential" else "descriptive",
               paired = d$paired, n_donors_ref = d$n_group[1], n_donors_test = d$n_group[2])]
    all_rows[[paste(ds, cn)]] <- tab[, .(dataset, analysis_role, contrast, inference, paired, gene_key, symbol, gene_type,
                                         is_lncRNA, log2FC, SE, CI_low, CI_high, PValue, FDR_lncRNA, FDR_all_genes,
                                         n_donors_ref, n_donors_test, donor_detection_ref, donor_detection_test,
                                         mean_cpm_ref, mean_cpm_test, median_cpm_ref, median_cpm_test,
                                         paired_direction_consistency, lodo_direction_flip)]
    sig <- if (inferential) tab[is_lncRNA == TRUE & FDR_lncRNA < DE$fdr] else tab[0]
    summ[[paste(ds, cn)]] <- base[, `:=`(n_donors_ref = d$n_group[1], n_donors_test = d$n_group[2],
      inference = if (inferential) "inferential" else "descriptive", paired = d$paired,
      genes_tested = nrow(tab), lnc_tested = tab[is_lncRNA == TRUE, .N],
      lnc_up_fdr_lfc = sig[log2FC >= DE$min_abs_log2fc, .N], lnc_down_fdr_lfc = sig[log2FC <= -DE$min_abs_log2fc, .N],
      lnc_up_fdr_lfc_lodo_stable = sig[log2FC >= DE$min_abs_log2fc & lodo_direction_flip %in% FALSE, .N])]
    msg("  ", cn, ": ", summ[[paste(ds, cn)]]$inference, if (d$paired) " paired" else "",
        " (", d$n_group[1], " vs ", d$n_group[2], " donors); lncRNA up ", summ[[paste(ds, cn)]]$lnc_up_fdr_lfc %||% NA)
  }
  rm(PB, counts); gc(verbose = FALSE)
}

A <- rbindlist(all_rows, use.names = TRUE)
S <- rbindlist(summ, use.names = TRUE, fill = TRUE)
save_tsv(A, "07_de_all_genes.tsv.gz")
save_tsv(A[is_lncRNA == TRUE], "07_de_lncRNA.tsv.gz")
save_tsv(S, "07_contrast_summary.tsv")
print(S)

pd <- melt(S[inference == "inferential"], id.vars = c("dataset", "contrast"),
           measure.vars = c("lnc_up_fdr_lfc", "lnc_down_fdr_lfc"), variable.name = "direction", value.name = "n")
if (nrow(pd)) {
  p <- ggplot(pd, aes(x = dataset, y = n, fill = direction)) + geom_col(position = "dodge") +
    facet_wrap(~ contrast, ncol = 1, scales = "free_y") + coord_flip() +
    labs(x = NULL, y = sprintf("lncRNAs at FDR < %.2f and |log2FC| >= %.1f", DE$fdr, DE$min_abs_log2fc),
         title = "Within-dataset donor-level DE (edgeR QL on pseudobulk)",
         subtitle = "Inferential contrasts only; each dataset analysed separately") +
    theme_bw(base_size = 8)
  save_fig(p, "07_de_summary", width = 8, height = 9)
}
write_session_info("07_within_dataset_de")
