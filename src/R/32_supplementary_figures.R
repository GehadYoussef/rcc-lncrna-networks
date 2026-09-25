# 32_supplementary_figures.R: ten supplementary composite figures.
# Sensitivity analyses: network diagnostics, module preservation, matched normal libraries,
# published signatures against random nulls, RIN and plate, normalisation, alternative
# endpoints, metric correlates, hub-lncRNA cross-validation, and calibration at the
# secondary horizon. As in 13_figures.R, every number in a title, subtitle or annotation
# is computed from the table the panel plots, and a panel whose source is missing is
# skipped and logged. Writes SupplementaryFigureS1 to S5, S7, S8 and S11 to S13 to
# figures/. Log labels S1 to S10 follow build order, not the file numbers.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
banner("32 | Supplementary figures")

# ---- helpers (same contract as 13_figures.R) --------------------------------
MISSING <- character(0)
SKIPPED <- character(0)
R_  <- function(f) fread(file.path(RESULTS_DIR, f))
Rx  <- function(f) {
  p <- file.path(RESULTS_DIR, f)
  if (file.exists(p)) return(fread(p))
  MISSING <<- unique(c(MISSING, f)); NULL
}
skip <- function(panel, why) {
  SKIPPED <<- c(SKIPPED, sprintf("%s (%s)", panel, why))
  msg("  SKIPPED ", panel, ": ", why); NULL
}
wrap_lab <- function(x, w = 30)
  vapply(x, function(s) paste(strwrap(s, width = w), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
# `word` names what the event column counts: deaths for overall survival, other
# events for the alternative endpoints.
n_ev_text <- function(n, ev, word = "deaths") {
  n <- n[is.finite(n)]; ev <- ev[is.finite(ev)]
  if (!length(n)) return("")
  if (uniqueN(n) == 1L && uniqueN(ev) == 1L)
    sprintf("%d patients, %d %s", as.integer(n[1]), as.integer(ev[1]), word)
  else sprintf("%d-%d patients, %d-%d %s", min(n), max(n), min(ev), max(ev), word)
}
ci_h <- function(mapping, width = 0.2, ...)
  geom_errorbar(mapping, orientation = "y", width = width, ...)
bt <- function(x) ifelse(grepl("^lnc", x, ignore.case = TRUE), "lncRNA", "protein-coding")

W <- 183 / 25.4
.fam <- unique(systemfonts::system_fonts()$family)
FIG_FONT <- if ("Arial" %in% .fam) "Arial" else
            if ("Liberation Sans" %in% .fam) "Liberation Sans" else "sans"
base <- theme_bw(base_size = 8, base_family = FIG_FONT) +
  theme(plot.title    = element_text(face = "bold", size = 8.5,
                                     margin = margin(t = 0, b = 2)),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 6.8, colour = "grey25",
                                     margin = margin(b = 3)),
        axis.title    = element_text(size = 7.5),
        axis.text     = element_text(size = 7, colour = "grey15"),
        strip.text    = element_text(size = 7, margin = margin(2, 2, 2, 2)),
        strip.background = element_rect(fill = "grey93", colour = NA),
        legend.title  = element_text(size = 7),
        legend.text   = element_text(size = 6.8),
        legend.key.size = unit(0.32, "cm"),
        legend.margin = margin(t = 0, b = 0),
        panel.grid.minor = element_blank(),
        plot.margin   = margin(4, 5, 3, 4))
theme_set(base)

OI <- c(orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00",
        purple = "#CC79A7", grey = "#666666")
COH_COL <- c("TCGA-KIRC" = unname(OI["vermillion"]), "CPTAC-3" = unname(OI["blue"]),
             "TCGA-KIRP" = unname(OI["green"]),      "TCGA-KICH" = unname(OI["purple"]))
BIO_COL <- c(lncRNA = unname(OI["vermillion"]), `protein-coding` = "#333333")
# Tissue colours (ColorBrewer BrBG endpoints), distinct from the cohort and
# biotype colours.
TISSUE_COL <- c(tumour = "#8C510A", normal = "#01665E",
                `matched normal` = "#01665E")
# Explicit breaks for log10 hazard-ratio axes. The defaults label too few values here.
HR_BRK <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.25, 1.5, 1.75, 2, 3)
# Six well-separated hues for small points with thin whiskers (no pure yellow).
SIG_PAL <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#666666", "#E69F00")
ANN <- function(lab, x = -Inf, y = Inf, hjust = -0.12, vjust = 1.3, size = 2.2)
  annotate("text", x = x, y = y, hjust = hjust, vjust = vjust, size = size,
           family = FIG_FONT, colour = "grey15", label = lab, lineheight = 0.95)
TG <- function(p, ch) p + labs(tag = ch) +
  theme(plot.tag = element_text(face = "bold", size = 9.5, family = FIG_FONT,
                                hjust = 0, vjust = 1),
        plot.tag.position = "topleft", plot.tag.location = "margin",
        plot.margin = margin(t = 13, r = 6, b = 3, l = 4))
ROW <- function(x) wrap_elements(full = x)
# Compose only the panels that were built. Letters follow the surviving order.
compose <- function(panels) {
  panels <- panels[!vapply(panels, is.null, logical(1))]
  if (!length(panels)) return(NULL)
  for (i in seq_along(panels)) panels[[i]] <- TG(panels[[i]], letters[i])
  panels
}
msg("Figure font: ", FIG_FONT, "; width ", round(W, 3), " in (183 mm)")

# ==== Supplementary Figure S1: network construction diagnostics ==============
banner("S1 | network diagnostics")
stl <- Rx("02_lncRNA_soft_threshold.tsv"); stm <- Rx("02_mRNA_soft_threshold.tsv")
szl <- Rx("02_lncRNA_module_sizes.tsv");   szm <- Rx("02_mRNA_module_sizes.tsv")
nsum <- Rx("02_network_summary.tsv")
P <- list()
if (!is.null(stl) && !is.null(stm)) {
  sft <- rbind(stl[, network := "lncRNA"], stm[, network := "protein-coding"])
  chosen <- if (!is.null(nsum))
    data.table(network = bt(nsum$biotype), power = nsum$soft_power) else NULL
  P$a <- ggplot(sft, aes(Power, SFT_R2, colour = network)) +
    geom_hline(yintercept = RSQ_CUT, linetype = 2, colour = "grey45",
               linewidth = 0.35) +
    geom_line(linewidth = 0.4) + geom_point(size = 1.2) +
    { if (!is.null(chosen)) geom_vline(data = chosen, aes(xintercept = power,
                                                          colour = network),
                                       linetype = 3, linewidth = 0.45) } +
    scale_colour_manual(values = BIO_COL, name = NULL) +
    labs(title = "Scale-free topology fit",
         subtitle = if (!is.null(chosen))
           sprintf("Dotted lines: chosen power (%s); dashed line: target R2 = %.2f",
                   paste(sprintf("%s %d", chosen$network, chosen$power),
                         collapse = ", "), RSQ_CUT) |> wrap_lab(58)
         else sprintf("Dashed line: target R2 = %.2f", RSQ_CUT),
         # Quoted: an unquoted hyphen inside a plotmath expression is parsed as
         # the minus operator and renders as "Scale - free".
         x = "Soft-thresholding power", y = expression("Scale-free fit"~R^2)) +
    theme(legend.position = "bottom")
  # mean connectivity reaches zero at high powers, and log10 of zero is dropped
  P$b <- ggplot(sft[mean_k > 0], aes(Power, mean_k, colour = network)) +
    geom_line(linewidth = 0.4) + geom_point(size = 1.2) +
    { if (!is.null(chosen)) geom_vline(data = chosen, aes(xintercept = power,
                                                          colour = network),
                                       linetype = 3, linewidth = 0.45) } +
    scale_y_log10() +
    scale_colour_manual(values = BIO_COL, guide = "none") +
    labs(title = "Mean connectivity",
         subtitle = "Log scale; same colours and chosen powers as a",
         x = "Soft-thresholding power", y = "Mean connectivity k")
} else skip("S1a/S1b", "soft-threshold tables missing")
if (!is.null(szl) && !is.null(szm)) {
  sz <- rbind(szl[, network := "lncRNA"], szm[, network := "protein-coding"])
  sz[, grey := module == "grey"]
  sz[, key := paste(network, module)]
  sz[, key := factor(key, levels = sz[order(network, n_genes), key])]
  tot <- sz[, .(total = sum(n_genes), grey = sum(n_genes[module == "grey"]),
                k = sum(module != "grey")), by = network]
  P$c <- ggplot(sz, aes(n_genes, key, fill = grey)) +
    geom_col(width = 0.72) +
    facet_wrap(~ network, scales = "free", ncol = 2) +
    scale_y_discrete(labels = function(x) sub("^(lncRNA|protein-coding) ", "", x)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = c(`FALSE` = unname(OI["blue"]), `TRUE` = "grey72"),
                      labels = c(`FALSE` = "assigned module", `TRUE` = "grey (unassigned)"),
                      name = NULL) +
    # "Genes", not "transcripts": the GDC STAR-Counts matrices are gene-level.
    labs(title = "Module sizes",
         subtitle = paste(sprintf("%s: %d modules, %s of %s genes unassigned",
                                  tot$network, tot$k,
                                  format(tot$grey, big.mark = ",", trim = TRUE),
                                  format(tot$total, big.mark = ",", trim = TRUE)),
                          collapse = "; ") |> wrap_lab(108),
         x = "Genes", y = NULL) +
    theme(axis.text.y = element_text(size = 6), legend.position = "bottom")
} else skip("S1c", "module size tables missing")
pp <- compose(P)
if (length(pp)) {
  s1 <- if (length(pp) == 3) ROW(pp[[1]] | pp[[2]]) / ROW(pp[[3]]) +
                             plot_layout(heights = c(1, 1.25))
        else Reduce(`|`, pp)
  save_fig(s1, "SupplementaryFigureS1_network_diagnostics", W, 6.6)
  msg("  S1 written")
}

# ==== Supplementary Figure S8: module preservation in every test cohort ======
banner("S2 | module preservation")
pl <- Rx("24_module_preservation_lncRNA.tsv"); pm <- Rx("24_module_preservation_mRNA.tsv")
ps_ <- Rx("24_module_preservation_summary.tsv")
P <- list()
if (!is.null(pl) || !is.null(pm)) {
  # modulePreservation always emits "grey" (unassigned) and "gold" (random reference)
  # rows, which are not modules.
  pres <- rbindlist(list(pl, pm), fill = TRUE)[!module %in% c("grey", "gold")]
  pres[, network := bt(network)]
  pres[, test_cohort := factor(test_cohort, levels = names(COH_COL))]
  # modulePreservation caps moduleSize at maxModuleSize (1,000), so the true
  # discovery sizes are taken from 02_*_module_sizes.tsv where available.
  szl2 <- Rx("02_lncRNA_module_sizes.tsv"); szm2 <- Rx("02_mRNA_module_sizes.tsv")
  n_capped <- 0L
  if (!is.null(szl2) && !is.null(szm2)) {
    sizes <- rbind(szl2[, .(network = "lncRNA", module, true_size = n_genes)],
                   szm2[, .(network = "protein-coding", module, true_size = n_genes)])
    pres <- merge(pres, sizes, by = c("network", "module"), all.x = TRUE)
    n_capped <- pres[is.finite(true_size) & true_size > moduleSize, uniqueN(module)]
    pres[, size_plot := fifelse(is.finite(true_size), as.numeric(true_size),
                                as.numeric(moduleSize))]
  } else pres[, size_plot := as.numeric(moduleSize)]
  P$a <- ggplot(pres, aes(size_plot, Zsummary, colour = test_cohort)) +
    geom_hline(yintercept = 2, linetype = 3, colour = "grey45", linewidth = 0.35) +
    geom_hline(yintercept = 10, linetype = 2, colour = "grey45", linewidth = 0.35) +
    annotate("text", x = Inf, y = 2, label = "Z = 2 (weak)", size = 1.9,
             hjust = 1.04, vjust = -0.45, family = FIG_FONT, colour = "grey35") +
    annotate("text", x = Inf, y = 10, label = "Z = 10 (strong)", size = 1.9,
             hjust = 1.04, vjust = -0.45, family = FIG_FONT, colour = "grey35") +
    geom_point(size = 1.4, alpha = 0.9) +
    facet_grid(network ~ test_cohort) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = COH_COL, guide = "none") +
    labs(title = sprintf("Preservation of all %d modules in each test cohort",
                         uniqueN(pres[, paste(network, module)])),
         subtitle = sprintf("WGCNA modulePreservation; %d module x cohort pairs (grey and gold excluded); sizes are discovery module sizes%s",
                            nrow(pres),
                            if (n_capped) sprintf(", %d of which modulePreservation reports capped at %s",
                                                  n_capped,
                                                  format(max(pres$moduleSize), big.mark = ",",
                                                         trim = TRUE)) else "") |> wrap_lab(108),
         x = "Module size (genes, log scale)", y = "Zsummary (log scale)")
}
if (!is.null(ps_)) {
  sm <- melt(ps_, id.vars = c("network", "test_cohort"),
             measure.vars = c("n_strong", "n_weak_moderate", "n_none"),
             variable.name = "class", value.name = "n")
  sm[, network := bt(network)]
  sm[, class := factor(c(n_strong = "strong (Z > 10)",
                         n_weak_moderate = "weak to moderate (2-10)",
                         n_none = "none (Z < 2)")[as.character(class)],
                       levels = c("strong (Z > 10)", "weak to moderate (2-10)",
                                  "none (Z < 2)"))]
  sm[, test_cohort := factor(test_cohort, levels = names(COH_COL))]
  P$b <- ggplot(sm, aes(test_cohort, n, fill = class)) +
    geom_col(width = 0.68) +
    facet_wrap(~ network) +
    scale_fill_manual(values = c(unname(OI["green"]), unname(OI["yellow"]),
                                 unname(OI["vermillion"])), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Preservation class counts",
         subtitle = sprintf("%d permutations for lncRNA, %d for protein-coding",
                            max(ps_[grepl("lnc", network), n_permutations]),
                            max(ps_[!grepl("lnc", network), n_permutations])),
         x = NULL, y = "Modules") +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 6.4))
}
pp <- compose(P)
if (length(pp)) {
  s2 <- if (length(pp) == 2) ROW(pp[[1]]) / ROW(pp[[2]]) +
                             plot_layout(heights = c(1.35, 1)) else ROW(pp[[1]])
  save_fig(s2, "SupplementaryFigureS8_module_preservation", W, 6.8)
  msg("  S2 written")
} else skip("S2", "24_module_preservation_* missing")

# ==== Supplementary Figure S3: matched normal libraries ======================
banner("S3 | matched normal libraries")
tvn <- Rx("25_tumour_vs_normal_distribution.tsv")
pss <- Rx("25_paired_scores_per_sample.tsv")
nap <- Rx("25_normal_axis_projection.tsv")
nms <- Rx("25_normal_metric_survival.tsv")
P <- list()
if (!is.null(tvn)) {
  d <- tvn[normal_set == "all normal libraries" &
           variable %in% c("pct_noFeature", "pct_multimapping")]
  if (nrow(d)) {
    dl <- rbind(d[, .(cohort, variable, description, tissue = "tumour",
                      n = n_tumour, med = median_tumour, lo = q25_tumour, hi = q75_tumour)],
                d[, .(cohort, variable, description, tissue = "matched normal",
                      n = n_normal, med = median_normal, lo = q25_normal, hi = q75_normal)])
    dl[, tissue := factor(tissue, levels = c("tumour", "matched normal"))]
    dl[, cohort := factor(cohort, levels = names(COH_COL))]
    # The unpaired contrast mixes tissue with plate: normal libraries occupy few tumour
    # plates, and plate is the largest source of variance in the metric. The
    # plate-restricted and within-patient contrasts are quoted beside it.
    nconf <- d[unpaired_plate_confounded == TRUE, .N]
    nf    <- d[unpaired_plate_confounded == TRUE & variable == "pct_noFeature"]
    # The plate-restricted contrast is reported for both normal sets (all normal
    # libraries, and normals of the analysed tumours) with its library counts.
    nf_pat <- if (nrow(nf)) tvn[cohort == nf$cohort[1] &
                                normal_set == "normals of analysed-tumour patients" &
                                variable == "pct_noFeature"] else NULL
    conf_clause <- if (nconf && nrow(nf)) {
      sprintf(". The unpaired contrast mixes tissue with plate in %d of %d rows: in %s the non-feature difference falls from %.2f percentage points unpaired to %.2f on the plates both tissues occupy (%d tumour and %d normal libraries, rank-sum p %.3f) for the normals of the analysed tumours, the figure the text quotes%s, and is %.2f within patient on %d pairs (signed-rank p %.2f)",
              nconf, nrow(d), nf$cohort[1], nf$median_diff_t_minus_n[1],
              if (!is.null(nf_pat) && nrow(nf_pat))
                nf_pat$plate_matched_median_diff_t_minus_n[1] else
                nf$plate_matched_median_diff_t_minus_n[1],
              if (!is.null(nf_pat) && nrow(nf_pat))
                nf_pat$plate_matched_n_tumour[1] else nf$plate_matched_n_tumour[1],
              if (!is.null(nf_pat) && nrow(nf_pat))
                nf_pat$plate_matched_n_normal[1] else nf$plate_matched_n_normal[1],
              if (!is.null(nf_pat) && nrow(nf_pat))
                nf_pat$plate_matched_rank_sum_p[1] else nf$plate_matched_rank_sum_p[1],
              if (!is.null(nf_pat) && nrow(nf_pat))
                sprintf(" (%.2f, p %.3f, on the %d plate-matched libraries of all normal libraries, which is the set these bars show)",
                        nf$plate_matched_median_diff_t_minus_n[1],
                        nf$plate_matched_rank_sum_p[1],
                        nf$plate_matched_n_normal[1]) else "",
              nf$matched_median_diff_t_minus_n[1], nf$matched_n_pairs[1],
              nf$matched_signed_rank_p[1])
    } else ""
    P$a <- ggplot(dl, aes(cohort, med, fill = tissue)) +
      geom_col(width = 0.66, position = position_dodge(0.72)) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.35,
                    colour = "grey25", position = position_dodge(0.72)) +
      facet_wrap(~ description, scales = "free_y") +
      scale_fill_manual(values = TISSUE_COL, name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      labs(title = "Library metrics in tumour and matched normal libraries",
           subtitle = sprintf("Median (IQR); n = %s tumour and %s normal libraries; %d of %d unpaired rank-sum tests at FDR < %.2f%s",
                              paste(unique(d$n_tumour), collapse = "/"),
                              paste(unique(d$n_normal), collapse = "/"),
                              d[fdr < FDR_ALPHA, .N], nrow(d), FDR_ALPHA,
                              conf_clause) |> wrap_lab(112),
           x = NULL, y = "Per cent of reads") +
      theme(legend.position = "bottom")
  } else skip("S3a", "no rows for the all-normal set")
}
if (!is.null(pss)) {
  d <- pss[paired == TRUE & is.finite(pct_noFeature) & is.finite(lnc_axis_proj)]
  if (nrow(d)) {
    d[, cohort := factor(cohort, levels = names(COH_COL))]
    d[, tissue := factor(tissue, levels = c("tumour", "normal"))]
    rho <- d[, .(rho = cor(pct_noFeature, lnc_axis_proj, method = "spearman"),
                 n = .N), by = .(cohort, tissue)]
    P$b <- ggplot(d, aes(pct_noFeature, lnc_axis_proj, colour = tissue)) +
      geom_point(size = 0.7, alpha = 0.6, stroke = 0) +
      facet_wrap(~ cohort, scales = "free_x") +
      scale_colour_manual(values = TISSUE_COL, name = NULL) +
      labs(title = "The projected lncRNA axis in paired libraries",
           subtitle = paste(sprintf("%s %s: rho %.2f (n = %d)", rho$cohort, rho$tissue,
                                    rho$rho, rho$n), collapse = "; ") |> wrap_lab(96),
           x = "Reads in no annotated feature (%)",
           y = "lncRNA axis projected on\nthe discovery PC1") +
      theme(legend.position = "bottom")
  } else skip("S3b", "no paired per-sample rows")
}
if (!is.null(nap)) {
  d <- nap[axis_type == "own_PC1"]
  d[, matrix := bt(matrix)]
  # One row per sample set: the two matrices are computed on slightly different
  # sample counts where a library lacks one of them.
  d[, lab := sprintf("%s, %s\n(n = %s)", cohort, tissue,
                     if (min(n) == max(n)) as.character(min(n))
                     else sprintf("%d-%d", min(n), max(n))), by = .(cohort, tissue)]
  d[, lab := factor(lab, levels = unique(d[order(cohort, -(tissue == "tumour")), lab]))]
  # Reference: median and range of the PC1 share of random discovery subsets of
  # the same n (n_matched_pc1_* columns). Discovery rows have no reference.
  d[, `:=`(ref_med = 100 * n_matched_pc1_share,
           ref_lo  = 100 * n_matched_pc1_share_min,
           ref_hi  = 100 * n_matched_pc1_share_max)]
  n_ref   <- d[is.finite(ref_med), .N]
  n_above <- d[is.finite(ref_hi) & 100 * var_share > ref_hi, .N]
  # A variance share does not depend on PC1 orientation, so every row is shown. The
  # subtitle counts rows whose sign is a convention (pc1_sign_interpretable).
  n_conv <- if ("pc1_sign_interpretable" %in% names(d))
    d[pc1_sign_interpretable == FALSE, .N] else NA_integer_
  P$c <- ggplot(d, aes(100 * var_share, lab, colour = matrix)) +
    ci_h(aes(xmin = ref_lo, xmax = ref_hi), width = 0.1, linewidth = 0.3,
         colour = "grey55", position = position_dodge(0.5), show.legend = FALSE,
         na.rm = TRUE) +
    geom_point(aes(x = ref_med), shape = 21, fill = "white", colour = "grey35",
               size = 1.8, stroke = 0.5, position = position_dodge(0.5),
               show.legend = FALSE, na.rm = TRUE) +
    geom_point(size = 2.1, position = position_dodge(0.5)) +
    scale_colour_manual(values = BIO_COL, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.1))) +
    labs(title = "Share of variance on each set's own PC1",
         subtitle = sprintf("%d sample sets; open grey points and bars are the median and range of the PC1 share of %s size-matched random subsets of the discovery samples (%d of %d rows); %d observed shares sit above that range%s",
                            uniqueN(d$lab),
                            format(max(d$n_matched_draws, na.rm = TRUE),
                                   big.mark = ",", trim = TRUE),
                            n_ref, nrow(d), n_above,
                            if (is.na(n_conv)) "" else
                              sprintf("; a share carries no sign, and the orientation of %d of these %d components is a convention, so only their magnitudes are quoted elsewhere",
                                      n_conv, nrow(d))) |> wrap_lab(58),
         x = "PC1 variance share (%)", y = NULL) +
    theme(axis.text.y = element_text(size = 6.2, lineheight = 0.9),
          legend.position = "bottom")
}
if (!is.null(nms)) {
  # Linear-scale terms only (nf_n / nf_t). The file also carries a log10 version
  # of every model, which would otherwise plot as a second point on each row.
  d <- nms[exposure == "pct_noFeature" & term %in% c("nf_n", "nf_t")]
  if (nrow(d)) {
    # Short model labels with the library as a facet strip: the full model names
    # are three lines each and collide when stacked on one axis.
    MODS <- c(`exposure alone` = "Exposure alone",
              `+ age, sex, T, N, M1, ordinal grade` = "+ clinical terms",
              `tumour and normal metric jointly + clinical terms` =
                "Both libraries\n+ clinical terms")
    d <- d[model %in% names(MODS)]
    d[, lab := factor(MODS[model], levels = rev(unname(MODS)))]
    d[, strip := factor(paste(exposure_library, "library"),
                        levels = c("tumour library", "normal library"))]
    d[, cohort := factor(cohort, levels = names(COH_COL))]
    P$d <- ggplot(d, aes(HR, lab, colour = cohort, alpha = !underpowered)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.3) +
      ci_h(aes(xmin = lo, xmax = hi), width = 0.16, linewidth = 0.35,
           position = position_dodge(0.62), show.legend = FALSE) +
      geom_point(size = 1.9, position = position_dodge(0.62)) +
      facet_wrap(~ strip, ncol = 1) +
      # Sparser breaks than HR_BRK, whose labels collide at half-page width.
      scale_x_log10(breaks = c(0.6, 0.8, 1, 1.25, 1.5, 2)) +
      scale_colour_manual(values = COH_COL, name = NULL) +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
      labs(title = "The metric as exposure in paired patients",
           subtitle = sprintf("Non-feature fraction, linear scale, per SD; %s; faded where events per parameter < %d",
                              n_ev_text(d$n, d$events), MIN_EPV) |> wrap_lab(58),
           x = "Hazard ratio per 1 SD (95% CI)", y = NULL) +
      theme(axis.text.y = element_text(size = 6.2, lineheight = 0.9),
            legend.position = "bottom")
  } else skip("S3d", "no matching normal-metric survival rows")
}
pp <- compose(P)
if (length(pp) >= 4) {
  s3 <- ROW(pp[[1]]) / ROW(pp[[2]]) / ROW(pp[[3]] | pp[[4]]) +
    plot_layout(heights = c(0.95, 0.95, 1.4))
  save_fig(s3, "SupplementaryFigureS3_matched_normal", W, 9.2)
  msg("  S3 written")
} else if (length(pp)) {
  save_fig(Reduce(function(a, b) a / b, lapply(pp, ROW)),
           "SupplementaryFigureS3_matched_normal", W, 2.6 * length(pp))
  msg("  S3 written (", length(pp), " panels)")
} else skip("S3", "25_* tables missing")

# ==== Supplementary Figure S7: published lncRNA signatures against random nulls
banner("S4 | published signatures")
sm_ <- Rx("26_published_signatures_models.tsv")
rn  <- Rx("26_random_signature_null.tsv")
sc  <- Rx("26_published_signatures_correlations.tsv")
hrg <- Rx("26_headline_ranges_by_provenance.tsv")
# Weight provenance decides what may be pooled. Some signatures publish
# coefficients, others had their weights signed by a univariate Cox fit in this
# cohort, which inflates their hazard ratios. Every quoted range is read from
# 26_headline_ranges_by_provenance.tsv, computed within one provenance class.
PROV <- c(external_coefficients       = "External coefficients",
          in_sample_signed            = "Weights signed in TCGA-KIRC",
          in_sample_refit             = "Coefficients refitted in TCGA-KIRC",
          footnote_incomplete_mapping = "Incomplete mapping (footnote)")
prov_lab <- function(x) fifelse(x %in% names(PROV), unname(PROV[x]), x)
# Subtitles break the signature count down by weight provenance.
PROV_PHRASE <- c(external_coefficients = "with coefficients published outside this cohort",
                 in_sample_signed = "whose per-gene weights had to be signed by a univariate Cox fit in this cohort",
                 in_sample_refit = "with coefficients refitted in this cohort",
                 footnote_incomplete_mapping = "held back for incomplete mapping")
prov_breakdown <- function(d) {
  z <- unique(d[, .(signature, weight_provenance)])[, .N, by = weight_provenance]
  z <- z[order(match(weight_provenance, names(PROV)))]
  paste(sprintf("%d %s", z$N,
                fifelse(z$weight_provenance %in% names(PROV_PHRASE),
                        unname(PROV_PHRASE[z$weight_provenance]),
                        z$weight_provenance)),
        collapse = " and ")
}
prov_fac <- function(x) factor(prov_lab(x),
                               levels = unname(PROV[names(PROV) %in% unique(x)]))
# One per-provenance range as stage 26 wrote it, NA if the quantity is absent.
rng <- function(qpat, spat, prov) {
  if (is.null(hrg)) return(NA_character_)
  v <- hrg[grepl(qpat, quantity) & grepl(spat, stratum, fixed = TRUE) &
           weight_provenance == prov, range_text]
  if (length(v)) v[1] else NA_character_
}
att3_ext <- rng("^Attenuation", "3_score_plus_clinical_plus_STAR | published weights",
                "external_coefficients")
att3_sig <- rng("^Attenuation", "3_score_plus_clinical_plus_STAR | published weights",
                "in_sample_signed")
att4_ext <- rng("^Attenuation", "4_residualised_score_plus_clinical | published weights",
                "external_coefficients")
att4_sig <- rng("^Attenuation", "4_residualised_score_plus_clinical | published weights",
                "in_sample_signed")
psr_ext  <- rng("Percentile of \\|rho\\|.*sign-randomised", "published weights",
                "external_coefficients")
psr_sig  <- rng("Percentile of \\|rho\\|.*sign-randomised", "published weights",
                "in_sample_signed")
pos_ext  <- rng("Percentile of \\|rho\\|.*outcome-signed", "published weights",
                "external_coefficients")
pos_sig  <- rng("Percentile of \\|rho\\|.*outcome-signed", "published weights",
                "in_sample_signed")
phr_ext  <- rng("Percentile of the univariate HR", "published weights",
                "external_coefficients")
phr_sig  <- rng("Percentile of the univariate HR", "published weights",
                "in_sample_signed")
rho_ext  <- rng("^Spearman rho", "observed expression | published weights",
                "external_coefficients")
rho_sig  <- rng("^Spearman rho", "observed expression | published weights",
                "in_sample_signed")
P <- list()
if (!is.null(rn)) {
  # Two nulls over the same gene sets: outcome-signed and sign-randomised.
  # Outcome-signing itself induces correlation with the metric, so both are drawn.
  NULLLAB <- c(outcome_signed = "outcome-signed",
               sign_randomised = "sign-randomised")
  draws   <- rn[row_type == "random_set" & sign_scheme %in% names(NULLLAB)]
  obs_all <- rn[row_type == "signature_percentile"]
  # Footnote-tier signatures mapped too few members to render the published
  # model. They are counted in the subtitle, not drawn.
  obs    <- obs_all[reporting_tier == "primary"]
  n_foot <- uniqueN(obs_all[reporting_tier != "primary", signature])
  if (nrow(draws)) {
    draws[, null_lab := factor(NULLLAB[sign_scheme], levels = unname(NULLLAB))]
    NULL_COL <- setNames(c("grey82", "#CBD9E6"), unname(NULLLAB))
    n_per <- draws[, .N, by = .(set_size, sign_scheme)][, max(N)]
    # The two nulls must share gene sets. The overlap is counted, not assumed.
    n_shared <- length(intersect(
      draws[sign_scheme == "outcome_signed",  unique(gene_set_key)],
      draws[sign_scheme == "sign_randomised", unique(gene_set_key)]))
    med_sr   <- draws[sign_scheme == "sign_randomised", median(HR)]
    pct_os_p <- 100 * draws[sign_scheme == "outcome_signed", mean(HR_p < 0.05)]
    P$a <- ggplot(draws, aes(abs(rho_noFeature), factor(set_size),
                             fill = null_lab)) +
      geom_violin(colour = "grey55", linewidth = 0.25, scale = "width",
                  position = position_dodge(0.8)) +
      { if (nrow(obs)) geom_point(data = obs, aes(abs(rho_noFeature),
                                                  factor(set_size), shape = weights),
                                  colour = unname(OI["vermillion"]), size = 2.1,
                                  inherit.aes = FALSE) } +
      scale_fill_manual(values = NULL_COL, name = "Random-set null") +
      scale_shape_manual(values = c(published = 16, refit_KIRC = 17),
                         labels = c(published = "published weights",
                                    refit_KIRC = "weights refitted in TCGA-KIRC"),
                         name = NULL) +
      guides(fill = guide_legend(order = 1, nrow = 2),
             shape = guide_legend(order = 2, nrow = 2)) +
      labs(title = "Published signatures against two random-signature nulls" |>
             wrap_lab(46),
           # Percentile ranges are quoted within a weight-provenance class.
           subtitle = sprintf("%s random sets per size and null; the two nulls are the same %s gene sets re-signed two ways. Points are %d signatures x %d weightings at the primary tier (%d held back for incomplete mapping). Published weights sit at percentiles %s of the sign-randomised null with external coefficients and %s with weights signed in TCGA-KIRC; outcome-signing itself induces the correlation, so that null sits higher (%s and %s)",
                              format(n_per, big.mark = ",", trim = TRUE),
                              format(n_shared, big.mark = ",", trim = TRUE),
                              uniqueN(obs$signature), uniqueN(obs$weights), n_foot,
                              psr_ext, psr_sig, pos_ext, pos_sig) |>
             wrap_lab(70),
           x = "|Spearman rho| with the non-feature fraction",
           y = "Signature size (genes)") +
      theme(legend.position = "bottom", legend.text = element_text(size = 6))
    P$b <- ggplot(draws, aes(HR, factor(set_size), fill = null_lab)) +
      geom_violin(colour = "grey55", linewidth = 0.25, scale = "width",
                  position = position_dodge(0.8)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey45", linewidth = 0.3) +
      { if (nrow(obs)) geom_point(data = obs, aes(HR, factor(set_size),
                                                  shape = weights),
                                  colour = unname(OI["vermillion"]), size = 2.1,
                                  inherit.aes = FALSE) } +
      scale_x_log10() +
      scale_fill_manual(values = NULL_COL, name = "Random-set null") +
      scale_shape_manual(values = c(published = 16, refit_KIRC = 17),
                         labels = c(published = "published weights",
                                    refit_KIRC = "weights refitted in TCGA-KIRC"),
                         name = NULL) +
      guides(fill = guide_legend(order = 1, nrow = 2),
             shape = guide_legend(order = 2, nrow = 2)) +
      labs(title = "Hazard ratio against two random-signature nulls" |> wrap_lab(46),
           subtitle = sprintf("Random sets, %s; the sign-randomised null has median HR %.2f while the outcome-signed null's sets are prognostic by construction (%.0f%% at p < 0.05). Published weights sit at percentiles %s of the outcome-signed null with external coefficients and %s with weights signed in TCGA-KIRC",
                              n_ev_text(draws$n, draws$events), med_sr, pct_os_p,
                              phr_ext, phr_sig) |> wrap_lab(60),
           x = "Hazard ratio per 1 SD (log scale)",
           y = "Signature size (genes)") +
      theme(legend.position = "bottom", legend.text = element_text(size = 6))
  } else skip("S4a/S4b", "no random-set rows")
}
if (!is.null(sm_)) {
  d <- sm_[reporting_tier == "primary" & weights == "published" &
           specification %in% c("1_score_alone", "2_score_plus_clinical",
                                "3_score_plus_clinical_plus_STAR",
                                "4_residualised_score_plus_clinical")]
  if (nrow(d)) {
    d[, spec := factor(wrap_lab(specification_label, 30),
                       levels = rev(unique(wrap_lab(specification_label, 30))))]
    # Faceted by weight provenance so that external and in-sample weights are not
    # read as one group.
    d[, prov := prov_fac(weight_provenance)]
    P$c <- ggplot(d, aes(HR, spec, colour = label)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.3) +
      ci_h(aes(xmin = lo, xmax = hi), width = 0.14, linewidth = 0.3,
           position = position_dodge(0.7), show.legend = FALSE) +
      geom_point(size = 1.7, position = position_dodge(0.7)) +
      facet_wrap(~ prov, nrow = 1) +
      scale_x_log10(breaks = c(0.75, 1, 1.25, 1.5, 2, 2.5)) +
      # Contrast-checked palette: no pure yellow and no near-neighbour hues.
      scale_colour_manual(values = SIG_PAL[seq_len(uniqueN(d$label))],
                          name = NULL) +
      labs(title = "Published signatures under four specifications, by weight provenance" |>
             wrap_lab(70),
           subtitle = sprintf("%d primary-tier signatures scored with their published gene sets (%s), %s; %d of %d intervals exclude 1. Attenuation of the clinically adjusted log hazard ratio on adding the three STAR metrics: %s per cent with external coefficients and %s with weights signed in this cohort; on residualising the expression instead, %s and %s",
                              uniqueN(d$signature), prov_breakdown(d),
                              n_ev_text(d$n, d$events),
                              d[lo > 1, .N], nrow(d),
                              att3_ext, att3_sig, att4_ext, att4_sig) |> wrap_lab(112),
           x = "Hazard ratio per 1 SD (95% CI)", y = NULL) +
      theme(axis.text.y = element_text(size = 6, lineheight = 0.9),
            legend.position = "bottom", legend.text = element_text(size = 6)) +
      guides(colour = guide_legend(nrow = 2, byrow = TRUE))
  } else skip("S4c", "no primary published-weight rows")
}
if (!is.null(sc)) {
  d <- sc[reporting_tier == "primary" & weights == "published" &
          variable == "pct_noFeature"]
  if (nrow(d)) {
    d[, expr := factor(expression, levels = c("observed", "residualised"))]
    d[, prov := prov_fac(weight_provenance)]
    res_max <- d[expr == "residualised", max(abs(spearman_rho))]
    P$d <- ggplot(d, aes(spearman_rho, reorder(label, spearman_rho), colour = expr)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
      geom_point(size = 2, position = position_dodge(0.5)) +
      facet_grid(prov ~ ., scales = "free_y", space = "free_y") +
      scale_colour_manual(values = c(observed = unname(OI["vermillion"]),
                                     residualised = unname(OI["blue"])),
                          name = "Expression matrix") +
      labs(title = "Signature score against the metric, by weight provenance" |>
             wrap_lab(70),
           subtitle = sprintf("%d primary-tier signatures scored with their published gene sets (%s), Spearman rho, n = %s; observed %s with external coefficients and %s with weights signed in TCGA-KIRC; after residualisation every |rho| is below %.2f, which is the check that the adjustment does what it claims",
                              uniqueN(d$signature), prov_breakdown(d),
                              paste(unique(d$n), collapse = "/"),
                              rho_ext, rho_sig, ceiling(res_max * 100) / 100) |>
             wrap_lab(112),
           x = "Spearman rho with the non-feature fraction", y = NULL) +
      theme(axis.text.y = element_text(size = 6), legend.position = "bottom",
            strip.text.y = element_text(size = 6.2, angle = 0))
  } else skip("S4d", "no primary published-weight correlation rows")
}
pp <- compose(P)
if (length(pp) >= 4) {
  # Panels c and d are faceted by weight provenance and need the full width. The
  # two null panels share the top row.
  s4 <- ROW(pp[[1]] | pp[[2]]) / ROW(pp[[3]]) / ROW(pp[[4]]) +
    plot_layout(heights = c(1.15, 1, 0.8))
  save_fig(s4, "SupplementaryFigureS7_published_signatures", W, 9.4)
  msg("  S4 written")
} else if (length(pp)) {
  save_fig(Reduce(function(a, b) a / b, lapply(pp, ROW)),
           "SupplementaryFigureS7_published_signatures", W, 3.0 * length(pp))
  msg("  S4 written (", length(pp), " panels)")
} else skip("S4", "26_* tables missing")

# ==== Supplementary Figure S4: the non-feature fraction against RIN and plate ==
banner("S5 | RIN and plate")
bio <- Rx("23_biospecimen_kirc.tsv"); axs <- Rx("27_axis_scores_per_sample.tsv")
avb <- Rx("23_axis_vs_biospecimen.tsv")
# Panels are built in source order and lettered in reading order (S5_ORDER).
P <- list()
if (!is.null(axs) && !is.null(bio)) {
  # RIN and plate come from the BCR biospecimen record (23), the axis and the
  # metric from the per-sample axis table (27), joined on the sample barcode.
  axs <- merge(axs[, .(sample_barcode, pct_noFeature, lnc_axis_z)],
               bio[, .(sample_barcode, rin, plate)], by = "sample_barcode")
  # One library has RIN recorded as 0.0, a coded missing value. It is dropped
  # here and the exclusion is stated on the panel, so n differs by one from
  # 23_axis_vs_biospecimen.tsv.
  n_rin0 <- axs[is.finite(rin) & rin <= 0, .N]
  rin_rho <- function(v) {
    d <- axs[is.finite(rin) & rin > 0 & is.finite(get(v))]
    list(d = d, rho = cor(d$rin, d[[v]], method = "spearman"), n = nrow(d))
  }
  rin_note <- if (n_rin0)
    sprintf("; %d library with RIN recorded as 0 excluded as a coded missing value",
            n_rin0) else ""
  # No trend line: the rank correlation is weak.
  a <- rin_rho("pct_noFeature")
  P$a <- ggplot(a$d, aes(rin, pct_noFeature)) +
    geom_point(size = 0.7, alpha = 0.45, colour = unname(OI["grey"]), stroke = 0) +
    ANN(sprintf("Spearman rho = %.3f\nn = %d", a$rho, a$n)) +
    labs(title = "Metric against RIN",
         subtitle = sprintf("TCGA-KIRC, RIN from the BCR record%s", rin_note) |>
           wrap_lab(40),
         x = "RIN", y = "Non-feature reads (%)")
  b <- rin_rho("lnc_axis_z")
  P$b <- ggplot(b$d, aes(rin, lnc_axis_z)) +
    geom_point(size = 0.7, alpha = 0.45, colour = unname(OI["vermillion"]), stroke = 0) +
    ANN(sprintf("Spearman rho = %.3f\nn = %d", b$rho, b$n)) +
    labs(title = "lncRNA axis against RIN",
         subtitle = sprintf("Same axis as Figure 2a%s", rin_note) |> wrap_lab(40),
         x = "RIN", y = "lncRNA PC1 (z)")
  pd <- axs[is.finite(pct_noFeature) & !is.na(plate) & plate != ""]
  ord <- pd[, .(m = median(pct_noFeature), n = .N), by = plate][order(m)]
  pd[, plate := factor(plate, levels = ord$plate)]
  P$c <- ggplot(pd, aes(plate, pct_noFeature)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.3, fill = "grey90",
                 colour = "grey35") +
    labs(title = "The metric by sequencing plate",
         subtitle = sprintf("%d plates, %d libraries, %d to %d per plate",
                            nrow(ord), nrow(pd), min(ord$n), max(ord$n)) |> wrap_lab(70),
         x = "Plate (ordered by median)", y = "Reads in no annotated feature (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8))
  ad <- axs[is.finite(lnc_axis_z) & !is.na(plate) & plate != ""]
  ad[, plate := factor(plate, levels = ord$plate)]
  P$d <- ggplot(ad, aes(plate, lnc_axis_z)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.3,
                 fill = alpha(unname(OI["vermillion"]), 0.25), colour = "grey35") +
    labs(title = "The lncRNA axis by sequencing plate",
         subtitle = sprintf("Plates in the same order; %d libraries with an axis score",
                            nrow(ad)) |> wrap_lab(70),
         x = "Plate (ordered by median non-feature fraction)", y = "lncRNA PC1 (z)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8))
} else skip("S5a-S5d", "27_axis_scores_per_sample.tsv or 23_biospecimen_kirc.tsv missing")
if (!is.null(avb)) {
  MEAS <- c(rin = "RIN", a260_a280 = "A260/A280",
            pct_necrosis_mean = "Necrosis (slide, %)",
            pct_tumor_nuclei_mean = "Tumour nuclei (slide, %)",
            pct_stromal_mean = "Stroma (slide, %)",
            StromalScore = "ESTIMATE stromal", ImmuneScore = "ESTIMATE immune")
  d <- avb[variable %in% c("pct_noFeature", "lnc_axis") & measure %in% names(MEAS)]
  d[, variable := c(pct_noFeature = "Non-feature fraction",
                    lnc_axis = "lncRNA PC1")[variable]]
  d[, mlab := MEAS[measure]]
  # Filled where p < 0.05 (unadjusted), open otherwise.
  d[, sig := p < 0.05]
  P$e <- ggplot(d, aes(spearman_rho, reorder(mlab, spearman_rho),
                       colour = variable, shape = sig)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_point(size = 2, stroke = 0.7, fill = "white",
               position = position_dodge(0.5)) +
    scale_colour_manual(values = c(`Non-feature fraction` = unname(OI["grey"]),
                                   `lncRNA PC1` = unname(OI["vermillion"])),
                        name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), name = NULL,
                       breaks = c("TRUE", "FALSE"),
                       labels = c(`TRUE` = "p < 0.05", `FALSE` = "p >= 0.05")) +
    labs(title = "Biospecimen correlates",
         subtitle = sprintf("Spearman rho, n = %d-%d; %d of %d tests at p < 0.05, unadjusted",
                            min(d$n), max(d$n), d[sig == TRUE, .N], nrow(d)) |>
           wrap_lab(44),
         x = "Spearman rho", y = NULL) +
    theme(axis.text.y = element_text(size = 6.4), legend.position = "bottom",
          legend.box = "vertical", legend.spacing.y = unit(0, "cm")) +
    guides(colour = guide_legend(nrow = 2, order = 1),
           shape = guide_legend(nrow = 1, order = 2))
}
# Reading order: RIN scatters and biospecimen correlates on the top row, then
# the two full-width plate panels. compose() letters by position.
S5_ORDER <- c("a", "b", "e", "c", "d")
P <- P[intersect(S5_ORDER, names(P))]
pp <- compose(P)
if (length(pp) >= 5) {
  s5 <- ROW(pp[[1]] | pp[[2]] | pp[[3]]) / ROW(pp[[4]]) / ROW(pp[[5]]) +
    plot_layout(heights = c(1.05, 1, 1))
  save_fig(s5, "SupplementaryFigureS4_rin_and_plate", W, 8.0)
  msg("  S5 written")
} else if (length(pp)) {
  save_fig(Reduce(function(a, b) a / b, lapply(pp, ROW)),
           "SupplementaryFigureS4_rin_and_plate", W, 2.6 * length(pp))
  msg("  S5 written (", length(pp), " panels)")
} else skip("S5", "source tables missing")

# ==== Supplementary Figure S2: normalisation check ===========================
banner("S6 | normalisation")
nz  <- Rx("29_normalisation_pc1.tsv")
mcn <- Rx("29_normalisation_module_contrast.tsv")
P <- list()
if (!is.null(nz)) {
  nz[, lab := factor(wrap_lab(matrix, 26), levels = rev(wrap_lab(matrix, 26)))]
  # Each PC1 is oriented against the cached FPKM axis, so signed correlations
  # are comparable across matrices. The subtitle counts matrices whose
  # mean-expression anchor is indeterminate.
  n_indet <- if ("mean_expr_anchor_determinate" %in% names(nz))
    nz[mean_expr_anchor_determinate == FALSE, .N] else NA_integer_
  p6a <- ggplot(nz, aes(rho_noFeature, lab)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = rho_noFeature, yend = lab), linewidth = 0.4,
                 colour = "grey60") +
    geom_point(size = 2.2, colour = unname(OI["vermillion"])) +
    geom_text(aes(label = sprintf("%.3f", rho_noFeature)), hjust = -0.35,
              size = 2.1, family = FIG_FONT, colour = "grey25") +
    scale_x_continuous(limits = c(0, 1.06), breaks = seq(0, 1, 0.25)) +
    labs(title = sprintf("The axis tracks the metric under all %d normalisations tested",
                         nrow(nz)) |> wrap_lab(36),
         subtitle = sprintf("Spearman rho of each matrix's PC1 with the non-feature fraction, n = %d; each sign is anchored to the cached FPKM axis%s",
                            max(nz$n),
                            if (is.na(n_indet) || !n_indet) "" else
                              sprintf(", the mean-expression anchor being indeterminate for %d of %d",
                                      n_indet, nrow(nz))) |> wrap_lab(64),
         x = "Spearman rho with the non-feature fraction", y = NULL) +
    theme(axis.text.y = element_text(size = 6.2, lineheight = 0.9))
  p6b <- ggplot(nz, aes(100 * pc1_var_share, lab)) +
    geom_segment(aes(x = 0, xend = 100 * pc1_var_share, yend = lab),
                 linewidth = 0.4, colour = "grey60") +
    geom_point(size = 2.2, colour = unname(OI["blue"])) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * pc1_var_share)), hjust = -0.3,
              size = 2.1, family = FIG_FONT, colour = "grey25") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(title = "PC1 variance share", subtitle = "Same matrices, same samples",
         x = "PC1 variance share (%)", y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p6c <- ggplot(nz, aes(rho_axis_fpkm, lab)) +
    geom_segment(aes(x = 0, xend = rho_axis_fpkm, yend = lab), linewidth = 0.4,
                 colour = "grey60") +
    geom_point(size = 2.2, colour = unname(OI["green"])) +
    geom_text(aes(label = sprintf("%.3f", rho_axis_fpkm)), hjust = -0.3,
              size = 2.1, family = FIG_FONT, colour = "grey25") +
    scale_x_continuous(limits = c(0, 1.14), breaks = seq(0, 1, 0.25)) +
    labs(title = "Agreement with the FPKM axis",
         subtitle = "Against the cached log2(FPKM+1) axis",
         x = "Spearman rho with the FPKM axis", y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  P$a <- p6a; P$b <- p6b; P$c <- p6c
} else skip("S6a-S6c", "29_normalisation_pc1.tsv missing")
# Module hazard ratios under TMM and FPKM: stage 29 bootstraps the difference in
# log hazard ratio, with BH control within each arm. Controlled counts are
# quoted, with the raw counts beside them.
if (!is.null(mcn)) {
  SCL <- c(fpkm = "log2(FPKM+1)", tmm = "log2 CPM TMM")
  w <- mcn[, .(module, n, events, fpkm = HR_fpkm_observed, tmm = HR_tmm_observed,
               sig_obs = as.logical(boot_sig_fdr_observed),
               ci_obs  = as.logical(boot_ci_excludes_0_observed),
               sig_res = as.logical(boot_sig_fdr_residualised),
               ci_res  = as.logical(boot_ci_excludes_0_residualised),
               nboot   = boot_n_observed)]
  w[, lab := factor(module, levels = w[order(fpkm), module])]
  SHIFT <- c(`TRUE`  = sprintf("shift at FDR < %.2f", FDR_ALPHA),
             `FALSE` = sprintf("no shift at FDR < %.2f", FDR_ALPHA))
  w[, shift := factor(SHIFT[as.character(sig_obs)], levels = unname(SHIFT))]
  lg <- melt(w, id.vars = c("module", "lab"), measure.vars = c("fpkm", "tmm"),
             variable.name = "scale", value.name = "HR")
  lg[, scale_lab := factor(SCL[as.character(scale)], levels = unname(SCL))]
  P$d <- ggplot() +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_segment(data = w, aes(x = fpkm, xend = tmm, y = lab, yend = lab,
                               colour = shift), linewidth = 0.9) +
    geom_point(data = lg, aes(HR, lab, shape = scale_lab), size = 1.9,
               colour = "grey15", fill = "white", stroke = 0.55) +
    scale_x_log10(breaks = HR_BRK) +
    scale_colour_manual(values = setNames(c(unname(OI["vermillion"]), "grey75"),
                                          unname(SHIFT)), name = NULL,
                        drop = FALSE) +
    scale_shape_manual(values = setNames(c(21, 16), unname(SCL)),
                       name = "Normalisation") +
    labs(title = sprintf("%d of %d lncRNA module hazard ratios shift detectably between the two scales",
                         w[sig_obs == TRUE, .N], nrow(w)) |> wrap_lab(70),
         subtitle = sprintf("Paired patient bootstrap of the trimmed-mean minus FPKM difference in the log hazard ratio, %s resamples, %s, Benjamini-Hochberg within arm; %d of %d at FDR < %.2f on the observed matrix and %d of %d on the residualised matrix (before multiplicity control, %d and %d intervals excluded zero); same modules, same patients, same covariates on both scales",
                            format(max(w$nboot), big.mark = ",", trim = TRUE),
                            n_ev_text(mcn$n, mcn$events),
                            w[sig_obs == TRUE, .N], nrow(w), FDR_ALPHA,
                            w[sig_res == TRUE, .N], nrow(w),
                            w[ci_obs == TRUE, .N], w[ci_res == TRUE, .N]) |>
           wrap_lab(112),
         x = "Hazard ratio per 1 SD (log scale)", y = NULL) +
    theme(axis.text.y = element_text(size = 6.4), legend.position = "bottom",
          legend.box = "horizontal")
} else skip("S6d", "29_normalisation_module_contrast.tsv missing")
pp <- compose(P)
if (length(pp) >= 4) {
  s6 <- ROW(pp[[1]] | pp[[2]] | pp[[3]]) / ROW(pp[[4]]) +
    plot_layout(heights = c(1, 1.35))
  save_fig(s6, "SupplementaryFigureS2_normalisation", W, 7.2)
  msg("  S6 written")
} else if (length(pp) == 3L) {
  save_fig(ROW(pp[[1]] | pp[[2]] | pp[[3]]) + plot_layout(),
           "SupplementaryFigureS2_normalisation", W, 3.4)
  msg("  S6 written (3 panels)")
} else if (length(pp)) {
  save_fig(Reduce(function(a, b) a / b, lapply(pp, ROW)),
           "SupplementaryFigureS2_normalisation", W, 3.4 * length(pp))
  msg("  S6 written (", length(pp), " panels)")
} else skip("S6", "29_* tables missing")

# ==== Supplementary Figure S13: endpoint sensitivity =========================
banner("S7 | endpoint sensitivity")
em <- Rx("30_endpoint_modules.tsv"); ee <- Rx("30_endpoint_exposure.tsv")
END <- c(OS_pipeline = "Overall survival (pipeline)",
         OS_cdr = "Overall survival (TCGA-CDR)",
         DSS = "Disease-specific survival",
         PFI = "Progression-free interval")
END_COL <- setNames(unname(OI[c("vermillion", "blue", "green", "purple")]), unname(END))
# What the event column counts under each endpoint. Only the overall survival
# endpoints count deaths.
END_EVENT <- c(OS_pipeline = "deaths", OS_cdr = "deaths",
               DSS = "cancer-specific deaths", PFI = "progression or death")
END_SHORT <- c(OS_pipeline = "overall survival", OS_cdr = "overall survival",
               DSS = "disease-specific survival",
               PFI = "the progression-free interval")
ev_gloss <- function(codes) {
  codes <- unique(as.character(codes))
  sprintf("events are %s", paste(unique(sprintf("%s for %s", END_EVENT[codes],
                                                END_SHORT[codes])),
                                 collapse = ", "))
}
P <- list()
if (!is.null(em)) {
  d <- em[specification == "principal"]
  d[, biotype := bt(biotype)]
  d[, lab := paste0(biotype, ": ", module)]
  d[, lab := factor(lab, levels = unique(d[endpoint == "OS_pipeline"][order(HR), lab]))]
  gloss_a <- ev_gloss(d$endpoint)
  d[, endpoint := factor(END[endpoint], levels = unname(END))]
  P$a <- ggplot(d, aes(HR, lab, colour = endpoint)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.3) +
    ci_h(aes(xmin = lo, xmax = hi), width = 0.18, linewidth = 0.32,
         position = position_dodge(0.7), show.legend = FALSE) +
    geom_point(size = 1.8, position = position_dodge(0.7)) +
    scale_x_log10(breaks = HR_BRK) +
    scale_colour_manual(values = END_COL, name = NULL) +
    labs(title = "The five prognostic modules under four endpoints",
         subtitle = sprintf("Principal specification; %s (%s); %d of %d at FDR < %.2f",
                            n_ev_text(d$n, d$events, "events"), gloss_a,
                            d[fdr_5modules < FDR_ALPHA, .N],
                            nrow(d), FDR_ALPHA) |> wrap_lab(92),
         x = "Hazard ratio per 1 SD (95% CI)", y = NULL) +
    theme(axis.text.y = element_text(size = 6.6), legend.position = "bottom")
}
if (!is.null(ee)) {
  d <- copy(ee)
  gloss_b <- ev_gloss(d$endpoint)
  # The subtitle counts intervals excluding 1 in each endpoint group.
  pfi_row <- d$endpoint == "PFI"
  k_surv <- sum(d$lo[!pfi_row] > 1); n_surv <- sum(!pfi_row)
  k_pfi  <- sum(d$lo[pfi_row]  > 1); n_pfi  <- sum(pfi_row)
  d[, endpoint := factor(END[endpoint], levels = unname(END))]
  # Strip the file's (a)/(b) prefixes, which clash with panel letters, and use
  # the plain name of the metric.
  d[, model_lab := sub("^\\([a-z]\\)\\s*", "", model)]
  d[, model_lab := sub("pct_noFeature alone", "Non-feature fraction alone",
                       model_lab, fixed = TRUE)]
  d[, lab := factor(wrap_lab(model_lab, 34),
                    levels = rev(unique(wrap_lab(model_lab, 34))))]
  P$b <- ggplot(d, aes(HR, lab, colour = endpoint)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.3) +
    ci_h(aes(xmin = lo, xmax = hi), width = 0.16, linewidth = 0.32,
         position = position_dodge(0.7), show.legend = FALSE) +
    geom_point(size = 1.8, position = position_dodge(0.7)) +
    facet_wrap(~ exposure_scale) +
    scale_x_log10(breaks = HR_BRK) +
    scale_colour_manual(values = END_COL, name = NULL) +
    labs(title = "The non-feature fraction as exposure under four endpoints",
         subtitle = sprintf("%s (%s); the interval excludes 1 in %d of %d estimates for overall and disease-specific survival and in %d of %d for the progression-free interval",
                            n_ev_text(d$n, d$events, "events"), gloss_b,
                            k_surv, n_surv, k_pfi, n_pfi) |> wrap_lab(92),
         x = "Hazard ratio per 1 SD (95% CI)", y = NULL) +
    theme(axis.text.y = element_text(size = 6.2, lineheight = 0.9),
          legend.position = "bottom")
}
# -- c: the cross-validated increment, endpoint by endpoint -------------------
# Only the bootstrap of the repeat-averaged predictor carries an interval. The
# mean over single repeats is drawn beside it with its range, labelled as such.
ci30 <- Rx("30_endpoint_cv_increment.tsv")
if (!is.null(ci30) && nrow(ci30)) {
  d <- copy(ci30)
  KIND <- c(boot = "bootstrap 95% CI of the predictor averaged over repeats",
            rep  = "mean over repeats (range, not a confidence interval)")
  dd <- rbind(
    d[, .(endpoint, kind = unname(KIND["boot"]), est = boot_delta,
          lo = boot_lo, hi = boot_hi)],
    d[, .(endpoint, kind = unname(KIND["rep"]), est = delta_mean,
          lo = range_lo, hi = range_hi)])
  dd[, ep := factor(END[endpoint], levels = rev(unname(END)))]
  dd[, kind := factor(kind, levels = unname(KIND))]
  # All annotations start at the same x, past the widest interval.
  ann <- d[, .(ep = factor(END[endpoint], levels = rev(unname(END))),
               x = max(dd$hi, na.rm = TRUE),
               txt = sprintf("  p %.3f, BH %.3f", boot_p, boot_p_BH_4endpoints))]
  P$c <- ggplot(dd, aes(est, ep, colour = ep, shape = kind)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    ci_h(aes(xmin = lo, xmax = hi), width = 0.16, linewidth = 0.32,
         position = position_dodge(0.6), show.legend = FALSE) +
    geom_point(size = 1.9, position = position_dodge(0.6)) +
    geom_text(data = ann, aes(x, ep, label = txt), hjust = 0, size = 1.9,
              family = FIG_FONT, colour = "grey25", inherit.aes = FALSE) +
    scale_colour_manual(values = END_COL, guide = "none") +
    scale_shape_manual(values = c(16, 1), name = NULL,
                       labels = function(x) wrap_lab(x, 36)) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.72))) +
    labs(title = "Increment over the clinical model, endpoint by endpoint",
         subtitle = sprintf("%d x %d cross-validation; %d module scores added to age, sex, T, N, M1 and ordinal grade; %s. The interval excludes zero for %d of %d endpoints and %d of %d survive Benjamini-Hochberg over the four (smallest adjusted p %.3f); the mean over repeats is smaller than the bootstrap estimate in %d of %d",
                            ML_N_REPEATS, ML_N_FOLDS, max(d$n_module_scores),
                            n_ev_text(d$n, d$events, "events"),
                            d[boot_lo > 0, .N], nrow(d),
                            d[boot_p_BH_4endpoints < 0.05, .N], nrow(d),
                            min(d$boot_p_BH_4endpoints),
                            d[delta_mean < boot_delta, .N], nrow(d)) |> wrap_lab(62),
         x = "Change in concordance index", y = NULL) +
    theme(axis.text.y = element_text(size = 6.4), legend.position = "bottom",
          legend.direction = "vertical", legend.text = element_text(size = 6),
          legend.key.height = unit(0.26, "cm"))
} else skip("S7c", "30_endpoint_cv_increment.tsv missing")
# -- d: modules rebuilt inside each fold against a matched fixed-module arm ----
# If the increment depended on modules defined on the whole cohort, rebuilding
# them inside each training fold would remove it. The fixed-module arm is
# matched (same patients, folds and inner splits).
fw <- Rx("30_endpoint_foldwise.tsv")
if (!is.null(fw) && nrow(fw)) {
  d <- copy(fw)
  ANA <- c("modules fixed on full cohort (apparent, locked coefficients of 09)" =
             "Apparent, locked coefficients",
           "modules fixed on full cohort (10 x 10 CV, section 6)" =
             "10 x 10 cross-validation, modules fixed",
           "modules fixed on full cohort (3 x 5 assignments, identical folds)" =
             "3 folds x 5 assignments, modules fixed",
           "modules recomputed inside each fold (3 x 5 assignments)" =
             "3 folds x 5 assignments, modules rebuilt in fold")
  IVL <- c("paired patient bootstrap 95% (in-sample)" = "in-sample bootstrap 95% CI",
           "paired patient bootstrap 95% on the repeat-averaged out-of-fold predictor" =
             "bootstrap 95% CI, averaged predictor",
           "range across fold assignments" = "range across assignments, not a CI")
  d[, alab := fifelse(analysis %in% names(ANA), unname(ANA[analysis]), analysis)]
  d[, ilab := fifelse(interval_type %in% names(IVL), unname(IVL[interval_type]),
                      interval_type)]
  # Each row label carries its own n and events, because only the two three-fold
  # rows are matched arms. n and events are collapsed separately.
  rng_txt <- function(x) if (min(x) == max(x)) as.character(as.integer(min(x)))
             else sprintf("%d-%d", as.integer(min(x)), as.integer(max(x)))
  d[, lab := wrap_lab(sprintf("%s (%s; n = %s, %s events)", alab, ilab,
                              rng_txt(n), rng_txt(events)), 46), by = analysis]
  d[, lab := factor(lab, levels = rev(unique(lab)))]
  d[, ep := factor(END[endpoint], levels = unname(END))]
  reb  <- d[grepl("recomputed inside each fold", analysis)]
  fxd  <- d[grepl("identical folds", analysis)]
  cv10 <- d[grepl("10 x 10", analysis)]
  opt  <- merge(fxd[, .(endpoint, fixed = delta_C)],
                reb[, .(endpoint, rebuilt = delta_C)], by = "endpoint")
  opt[, cost := fixed - rebuilt]
  cost_txt <- paste(sprintf("%+.4f for %s", opt$cost, END_SHORT[opt$endpoint]),
                    collapse = " and ")
  dss <- opt[endpoint == "DSS"]
  # Optimism is fixed minus rebuilt (the file's definition), shown signed.
  ttl_d <- if (nrow(dss)) {
    sprintf("Module-definition optimism for %s: %+.4f in the concordance index",
            END_SHORT[["DSS"]], dss$cost)
  } else {
    "Modules rebuilt inside each fold against a matched fixed-module arm"
  }
  # For disease-specific survival, the matched arm separates the effect of
  # coarser cross-validation from that of module definition.
  dss_clause <- if (nrow(dss) && nrow(cv10[endpoint == "DSS"])) {
    sprintf(" For %s both %d-fold rows (%+.4f fixed, %+.4f rebuilt) sit below the %d x %d value (%+.4f), which the matched fixed-module arm attributes to the coarser cross-validation rather than to module definition.",
            END_SHORT[["DSS"]], FOLDWISE_K, dss$fixed, dss$rebuilt,
            ML_N_REPEATS, ML_N_FOLDS, cv10[endpoint == "DSS", delta_C])
  } else ""
  P$d <- ggplot(d, aes(delta_C, lab, colour = ep)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    ci_h(aes(xmin = delta_lo, xmax = delta_hi), width = 0.16, linewidth = 0.32,
         position = position_dodge(0.62), show.legend = FALSE) +
    geom_point(size = 1.9, position = position_dodge(0.62)) +
    scale_colour_manual(values = END_COL, name = NULL, drop = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))) +
    labs(title = ttl_d |> wrap_lab(58),
         subtitle = sprintf("Only the two %d-fold rows are matched arms (the same patients, folds and inner splits, differing only in how the modules are defined), so rebuilding costs %s.%s Ranges span %d fold assignments and are not confidence intervals; every row carries its own n",
                            FOLDWISE_K, cost_txt, dss_clause,
                            max(d[grepl("^range", interval_type),
                                  n_assignments], na.rm = TRUE)) |> wrap_lab(66),
         x = "Change in concordance index", y = NULL) +
    theme(axis.text.y = element_text(size = 5.8, lineheight = 0.9),
          legend.position = "bottom", legend.text = element_text(size = 6))
} else skip("S7d", "30_endpoint_foldwise.tsv missing")
pp <- compose(P)
if (length(pp) >= 4) {
  # Row 3 carries panel d's long labels and subtitle, so it is taller.
  s7 <- ROW(pp[[1]]) / ROW(pp[[2]]) / ROW(pp[[3]] | pp[[4]]) +
    plot_layout(heights = c(0.95, 0.9, 1.35))
  save_fig(s7, "SupplementaryFigureS13_endpoint_sensitivity", W, 10.4)
  msg("  S7 written")
} else if (length(pp)) {
  save_fig(Reduce(function(a, b) a / b, lapply(pp, ROW)),
           "SupplementaryFigureS13_endpoint_sensitivity", W, 3.2 * length(pp))
  msg("  S7 written (", length(pp), " panels)")
} else skip("S7", "30_* tables missing")

# ==== Supplementary Figure S5: what the non-feature fraction correlates with ==
banner("S8 | metric correlates")
mb <- Rx("27_metric_biology_correlations.tsv"); mu <- Rx("28_noFeature_vs_mutations.tsv")
P <- list()
if (!is.null(mb)) {
  d <- mb[exposure == "pct_noFeature" & is.finite(spearman_rho)]
  # Fill encodes the source matrix. Hypoxia metagenes and the lncRNA FPKM share
  # use the observed matrix, ESTIMATE scores and slide measures use none.
  MATLAB <- c(observed = "observed matrix", residualised = "residualised matrix",
              none = "not computed from an expression matrix")
  d[, mat := fifelse(matrix == "", unname(MATLAB["none"]),
                     unname(MATLAB[matrix]))]
  d[, mat := factor(mat, levels = unname(MATLAB))]
  d <- d[order(-abs(spearman_rho))][, head(.SD, 24)]
  # Readable names. An eigengene can appear once per matrix, so the matrix is
  # also written into every eigengene label.
  pretty_cand <- function(x) {
    y <- sub("^mRNA_ME", "protein-coding ", x)
    y <- sub("^lnc_ME", "lncRNA ", y)
    y <- ifelse(y != x, paste0(y, " eigengene"), y)
    DICT <- c(lnc_fpkm_share = "lncRNA share of FPKM",
              StromalScore = "ESTIMATE stromal", ImmuneScore = "ESTIMATE immune",
              ESTIMATEScore = "ESTIMATE combined",
              slide_stromal_pct = "Stroma (slide, %)",
              slide_tumour_nuclei_pct = "Tumour nuclei (slide, %)",
              proliferation_score_10genes = "Proliferation score (10 genes)",
              hypoxia_score_Buffa2010 = "Hypoxia score (Buffa 2010)",
              hypoxia_score_Buffa2010_both_tubulins = "Hypoxia (Buffa, both tubulins)",
              hypoxia_score_Buffa2010_median_FPKM_ge1 = "Hypoxia (Buffa, FPKM >= 1 genes)",
              hypoxia_score_Buffa2010_TUBA1A_mapping = "Hypoxia (Buffa, TUBA1A mapping)")
    ifelse(x %in% names(DICT), DICT[x], y)
  }
  d[, cand_lab := pretty_cand(candidate)]
  d[, cand_lab := fifelse(grepl("eigengene$", cand_lab),
                          sprintf("%s (%s)", cand_lab,
                                  sub(" matrix$", "", as.character(mat))),
                          cand_lab)]
  d[, key := paste(candidate, mat)]
  d[, key := factor(key, levels = d[order(spearman_rho), key])]
  KEY_LAB <- setNames(d$cand_lab, as.character(d$key))
  P$a <- ggplot(d, aes(spearman_rho, key, fill = mat)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_col(width = 0.72) +
    scale_y_discrete(labels = function(x) KEY_LAB[x]) +
    # Breaks every 0.1 across the whole range, negative side included.
    scale_x_continuous(breaks = seq(-1, 1, 0.1),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_manual(values = setNames(unname(OI[c("vermillion", "blue", "grey")]),
                                        unname(MATLAB)),
                      name = "Source of the candidate", drop = FALSE) +
    labs(title = "What the non-feature fraction correlates with",
         subtitle = sprintf("%d strongest of %d candidates, n = %d-%d; observed eigengenes carry the metric by construction",
                            nrow(d), nrow(mb[exposure == "pct_noFeature"]),
                            min(d$n), max(d$n)) |> wrap_lab(92),
         x = "Spearman rho with the non-feature fraction", y = NULL) +
    theme(axis.text.y = element_text(size = 5.8), legend.position = "bottom")
}
if (!is.null(mu)) {
  d <- mu[test == "Wilcoxon rank-sum" & is.finite(difference_of_medians)]
  if (nrow(d)) {
    d[, out := fifelse(grepl("noFeature", outcome), "Non-feature fraction (%)",
                       "lncRNA axis (PC1, z)")]
    # One y entry per gene: n_mutated differs slightly between the two outcomes
    # because the axis is defined on fewer libraries, so the label carries the
    # range and the ordering is taken from the metric panel.
    d[, glab := sprintf("%s\n(%s mutated)", variable,
                        if (min(n_mutated) == max(n_mutated)) as.character(min(n_mutated))
                        else sprintf("%d-%d", min(n_mutated), max(n_mutated))),
      by = variable]
    ordg <- d[grepl("noFeature", outcome)][order(difference_of_medians), glab]
    d[, glab := factor(glab, levels = ordg)]
    P$b <- ggplot(d, aes(difference_of_medians, glab, fill = fdr < FDR_ALPHA)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
      geom_col(width = 0.68) +
      facet_wrap(~ out, scales = "free_x") +
      scale_fill_manual(values = c(`FALSE` = "grey72",
                                   `TRUE` = unname(OI["vermillion"])),
                        labels = c(`FALSE` = sprintf("FDR >= %.2f", FDR_ALPHA),
                                   `TRUE` = sprintf("FDR < %.2f", FDR_ALPHA)),
                        name = NULL) +
      # The subtitle gives the smallest mutated group, which limits power.
      labs(title = sprintf("No driver mutation is associated with the metric at FDR < %.2f",
                           FDR_ALPHA),
           subtitle = sprintf("Difference of medians, mutated minus wild type; n = %s patients with a mutation call; %d of %d at FDR < %.2f, on as few as %d mutated patients",
                              if (min(d$n) == max(d$n)) as.character(max(d$n))
                              else sprintf("%d-%d", min(d$n), max(d$n)),
                              d[fdr < FDR_ALPHA, .N], nrow(d), FDR_ALPHA,
                              min(d$n_mutated, na.rm = TRUE)) |>
             wrap_lab(92),
           x = "Difference of medians (mutated - wild type)", y = NULL) +
      theme(axis.text.y = element_text(size = 6, lineheight = 0.9),
            legend.position = "bottom")
  } else skip("S8b", "no rank-sum rows in 28_noFeature_vs_mutations.tsv")
}
pp <- compose(P)
if (length(pp)) {
  s8 <- if (length(pp) == 2) ROW(pp[[1]]) / ROW(pp[[2]]) +
                             plot_layout(heights = c(1.3, 1)) else ROW(pp[[1]])
  save_fig(s8, "SupplementaryFigureS5_metric_correlates", W, 7.4)
  msg("  S8 written")
} else skip("S8", "27/28 tables missing")

# ==== Supplementary Figure S11: hub-lncRNA cross-validation ==================
banner("S9 | hub cross-validation")
cvr <- Rx("05_cv_cindex_per_repeat.tsv"); hub <- Rx("05_hub_module_selection_frequency.tsv")
P <- list()
if (!is.null(cvr)) {
  MODS <- c(clinical = "Clinical", clinical_eig = "Clinical + modules",
            clinical_hub = "Clinical + hub lncRNAs", augmented = "Augmented",
            augmented_eig = "Augmented + modules",
            augmented_hub = "Augmented + hub lncRNAs")
  d <- melt(cvr, id.vars = "repeat_id", measure.vars = names(MODS),
            variable.name = "model", value.name = "C")
  d[, mlab := factor(MODS[as.character(model)], levels = unname(MODS))]
  P$a <- ggplot(d, aes(mlab, C)) +
    geom_line(aes(group = repeat_id), colour = "grey78", linewidth = 0.25) +
    geom_point(aes(colour = grepl("^Augmented", mlab)), size = 1.3, alpha = 0.85) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.45, linewidth = 0.3,
                 colour = "grey15") +
    scale_colour_manual(values = c(`FALSE` = unname(OI["orange"]),
                                   `TRUE` = unname(OI["skyblue"])),
                        labels = c(`FALSE` = "clinical comparator",
                                   `TRUE` = "augmented comparator"), name = NULL) +
    # The subtitle defines both comparators as built by clinical_design().
    labs(title = sprintf("Cross-validated concordance, %d repeats of %d-fold",
                         nrow(cvr), ML_N_FOLDS),
         subtitle = paste("Lines join the same repeat; bars are the mean over repeats.",
                          "Clinical = age, sex, T, N, M1, ordinal grade;",
                          "Augmented = clinical + ESTIMATE stromal and immune scores",
                          "+ the three STAR library metrics") |> wrap_lab(74),
         x = NULL, y = "Out-of-fold concordance index") +
    theme(axis.text.x = element_text(size = 6, angle = 20, hjust = 1),
          legend.position = "bottom")
  DEL <- c(delta_clinical_eig_vs_clinical = "Modules over clinical",
           delta_clinical_hub_vs_clinical = "Hub lncRNAs over clinical",
           delta_augmented_eig_vs_augmented = "Modules over augmented",
           delta_augmented_hub_vs_augmented = "Hub lncRNAs over augmented")
  dd <- melt(cvr, id.vars = "repeat_id", measure.vars = names(DEL),
             variable.name = "comparison", value.name = "delta")
  dd[, clab := factor(DEL[as.character(comparison)], levels = rev(unname(DEL)))]
  pos <- dd[, .(k = sum(delta > 0), n = .N), by = clab]
  P$b <- ggplot(dd, aes(delta, clab)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_point(size = 1.3, alpha = 0.7, colour = unname(OI["grey"]),
               position = position_jitter(height = 0.12, seed = SEED)) +
    stat_summary(fun = mean, geom = "point", size = 2.4,
                 colour = unname(OI["vermillion"])) +
    geom_text(data = pos, aes(x = Inf, y = clab, label = sprintf("%d/%d positive", k, n)),
              hjust = 1.05, size = 2.0, family = FIG_FONT, colour = "grey25") +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.28))) +
    labs(title = "Increment per cross-validation repeat",
         subtitle = sprintf("Grey points are repeats, red points the mean over %d repeats",
                            nrow(cvr)) |> wrap_lab(60),
         x = "Change in concordance index", y = NULL) +
    theme(axis.text.y = element_text(size = 6.4))
} else skip("S9a/S9b", "05_cv_cindex_per_repeat.tsv missing")
if (!is.null(hub)) {
  P$c <- ggplot(hub, aes(100 * frac_folds_selected, reorder(module, frac_folds_selected))) +
    geom_col(width = 0.68, fill = unname(OI["blue"])) +
    geom_text(aes(label = sprintf("%d", n_folds_selected)), hjust = -0.25,
              size = 2.1, family = FIG_FONT, colour = "grey25") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(title = "How often each module supplies a selected hub",
         subtitle = sprintf("Elastic-net selection inside %d folds (%d repeats x %d folds); labels are fold counts",
                            round(max(hub$n_folds_selected) / max(hub$frac_folds_selected)),
                            ML_N_REPEATS, ML_N_FOLDS) |> wrap_lab(60),
         x = "Folds in which a hub from the module was selected (%)", y = NULL) +
    theme(axis.text.y = element_text(size = 6.4))
} else skip("S9c", "05_hub_module_selection_frequency.tsv missing")
pp <- compose(P)
if (length(pp) >= 3) {
  s9 <- ROW(pp[[1]]) / ROW(pp[[2]] | pp[[3]]) + plot_layout(heights = c(1, 1))
  save_fig(s9, "SupplementaryFigureS11_hub_cross_validation", W, 6.4)
  msg("  S9 written")
} else if (length(pp)) {
  save_fig(Reduce(function(a, b) a / b, lapply(pp, ROW)),
           "SupplementaryFigureS11_hub_cross_validation", W, 3.2 * length(pp))
  msg("  S9 written (", length(pp), " panels)")
} else skip("S9", "05_* tables missing")

# ==== Supplementary Figure S12: calibration and utility sensitivity ==========
# CPTAC-3 at the secondary horizon and under discovery-fixed standardisation.
banner("S10 | calibration and utility sensitivity")
cb <- Rx("12_calibration_bins.tsv"); cs2 <- Rx("12_calibration_summary.tsv")
dc <- Rx("12_validation_decision_curve.tsv")
MOD2 <- c(comparator = "Clinical", `comparator + modules` = "Clinical + modules")
MOD2_COL <- c(Clinical = "#4D4D4D", `Clinical + modules` = unname(OI["orange"]))
P <- list()
cal_panel <- function(yrs, std, ttl) {
  d <- cb[cohort == "validation" & comparator == "clinical" &
          standardisation == std & years == yrs]
  s <- cs2[cohort == "validation" & comparator == "clinical" &
           standardisation == std & years == yrs]
  if (!nrow(d)) return(skip(ttl, "no calibration bins for this arm"))
  d[, mlab := factor(MOD2[model], levels = unname(MOD2))]
  s[, mlab := factor(MOD2[model], levels = unname(MOD2))]
  s[, short := c(Clinical = "Clinical:",
                 `Clinical + modules` = "+ modules:")[as.character(mlab)]]
  L <- c(0, max(c(d$predicted, d$obs_hi), na.rm = TRUE) * 1.05)
  ggplot(d, aes(predicted, observed, colour = mlab)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55",
                linewidth = 0.3) +
    geom_errorbar(aes(ymin = obs_lo, ymax = obs_hi), width = 0.015, linewidth = 0.3,
                  show.legend = FALSE) +
    geom_line(linewidth = 0.35) + geom_point(size = 1.5) +
    # Equal limits and a square panel (aspect.ratio, since coord_fixed clips the title
    # under patchwork), so calibration reads as distance from the diagonal.
    coord_cartesian(xlim = L, ylim = L) +
    scale_colour_manual(values = MOD2_COL, name = NULL) +
    labs(title = ttl,
         subtitle = paste(c(sprintf("CPTAC-3, %d at risk at %d years; bars are Kaplan-Meier 95%% CIs; slope (95%% CI), observed/expected:",
                                    max(s$n_at_risk), yrs) |> wrap_lab(52),
                            sprintf("%s %.2f (%.2f-%.2f), O/E %.2f", s$short, s$slope,
                                    s$slope_lo, s$slope_hi, s$OE)), collapse = "\n"),
         x = "Predicted risk of death", y = "Observed risk") +
    theme(legend.position = "bottom", aspect.ratio = 1)
}
dca_panel <- function(yrs, std, ttl) {
  STRAT <- c(comparator = "Clinical", `comparator + modules` = "Clinical + modules",
             `treat all` = "Treat all", `treat none` = "Treat none")
  d <- dc[baseline == "discovery" & comparator == "clinical" &
          standardisation == std & years == yrs]
  if (!nrow(d)) return(skip(ttl, "no decision-curve rows for this arm"))
  d[, strat := factor(STRAT[strategy], levels = unname(STRAT))]
  nb <- dcast(d[strategy %in% c("comparator", "comparator + modules")],
              threshold ~ strategy, value.var = "net_benefit")
  nb[, dd := `comparator + modules` - comparator]
  win <- nb[threshold >= 0.10 & threshold <= 0.40]
  # The y floor clears the lowest model curve. The treat-all reference may run
  # off the panel.
  mod_nb <- d[strategy %in% c("comparator", "comparator + modules"), net_benefit]
  flo <- min(-0.02, min(mod_nb, na.rm = TRUE) * 1.15)
  ggplot(d, aes(threshold, net_benefit, colour = strat, linetype = strat)) +
    geom_line(linewidth = 0.5) +
    coord_cartesian(ylim = c(flo, max(d$net_benefit, na.rm = TRUE) * 1.08)) +
    scale_colour_manual(values = c(unname(MOD2_COL), "grey40", "grey72"), name = NULL) +
    scale_linetype_manual(values = c(1, 1, 2, 3), name = NULL) +
    labs(title = ttl,
         subtitle = sprintf("Median change in net benefit over thresholds 0.10-0.40: %+.3f (positive at %d of %d); the treat-all reference runs below the panel",
                            median(win$dd), sum(win$dd > 0), nrow(win)) |> wrap_lab(52),
         x = "Threshold probability", y = "Net benefit") +
    theme(legend.position = "bottom", legend.direction = "vertical",
          legend.text = element_text(size = 6),
          legend.key.height = unit(0.26, "cm"))
}
if (!is.null(cb) && !is.null(cs2)) {
  P$a <- cal_panel(SECONDARY_HORIZON_YR, "cohort",
                   sprintf("Calibration at %d years", SECONDARY_HORIZON_YR))
  P$b <- cal_panel(PRIMARY_HORIZON_YR, "discovery",
                   sprintf("Calibration at %d years,\ndiscovery-fixed scores",
                           PRIMARY_HORIZON_YR))
} else skip("S10a/S10b", "12_calibration_* missing")
if (!is.null(dc)) {
  # Both decision curves are scored in CPTAC-3, which the titles name.
  P$c <- dca_panel(SECONDARY_HORIZON_YR, "cohort",
                   sprintf("Decision curve at %d years, CPTAC-3",
                           SECONDARY_HORIZON_YR))
  P$d <- dca_panel(PRIMARY_HORIZON_YR, "discovery",
                   sprintf("Decision curve at %d years, CPTAC-3,\ndiscovery-fixed scores",
                           PRIMARY_HORIZON_YR))
} else skip("S10c/S10d", "12_validation_decision_curve.tsv missing")
pp <- compose(P)
if (length(pp) >= 4) {
  # The calibration row is the taller one: those two panels fix their aspect
  # ratio, so the only way to make them wider is to give the row more height.
  s10 <- ROW(pp[[1]] | pp[[2]]) / ROW(pp[[3]] | pp[[4]]) +
    plot_layout(heights = c(1.25, 1))
  save_fig(s10, "SupplementaryFigureS12_calibration_sensitivity", W, 7.6)
  msg("  S10 written")
} else if (length(pp)) {
  save_fig(Reduce(function(a, b) a | b, lapply(pp, ROW)),
           "SupplementaryFigureS12_calibration_sensitivity", W, 3.6)
  msg("  S10 written (", length(pp), " panels)")
} else skip("S10", "12_* tables missing")

# ---- report missing sources and skipped panels ----
if (length(MISSING))
  msg("NOTE: source files absent: ", paste(MISSING, collapse = ", "))
if (length(SKIPPED))
  msg("NOTE: panels skipped: ", paste(SKIPPED, collapse = "; "))
msg("Supplementary figures written to ", FIG_DIR, " as SVG and ", FIG_DPI, " dpi PNG")
write_session_info("32_supplementary_figures")
banner("32 | done")
