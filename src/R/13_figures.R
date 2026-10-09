# 13_figures.R: the seven composite main figures.
#
# Each panel carries its letter and a short title. Sample sizes, model
# definitions and summary statistics are given in the figure legends, and every
# number drawn inside a panel is computed from the table the panel plots. A
# panel whose source table is absent is skipped with a message.
# Reads result tables from stages 01 to 43. Writes to figures/ (SVG and PNG):
# Figure1_design_cohorts to Figure7_prediction and
# SupplementaryFigureS9_lncRNA_module_annotation. Sections are in build order, and
# object names (fig3, p3a, ...) do not follow the figure numbers.
# Style: lower-case panel letters, 183 mm width, Okabe-Ito palette, hazard
# ratios per 1 SD, and the clinical model as primary comparator.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
banner("13 | Main figures")

# ---- small helpers ----------------------------------------------------------
R_  <- function(f) fread(file.path(RESULTS_DIR, f))
# Optional source: NULL if the file is absent, so the panel can be skipped.
Rx  <- function(f) if (file.exists(file.path(RESULTS_DIR, f))) R_(f) else NULL
MISSING <- character(0)                       # source files that were absent

# Wrap a long axis label onto several lines at `w` characters.
wrap_lab <- function(x, w = 30)
  vapply(x, function(s) paste(strwrap(s, width = w), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
# Horizontal 95% CI (geom_errorbarh is deprecated in ggplot2 >= 3.5).
ci_h <- function(mapping, width = 0.2, ...)
  geom_errorbar(mapping, orientation = "y", width = width, ...)
# One biotype label across source tables.
bt <- function(x) ifelse(grepl("^lnc", x, ignore.case = TRUE), "lncRNA", "protein-coding")

# ---- page geometry, fonts, palette -----------------------------------------
W <- 183 / 25.4                               # 183 mm double-column width
.fam <- unique(systemfonts::system_fonts()$family)
FIG_FONT <- if ("Arial" %in% .fam) "Arial" else
            if ("Liberation Sans" %in% .fam) "Liberation Sans" else "sans"

base <- theme_bw(base_size = 8, base_family = FIG_FONT) +
  theme(plot.title    = element_text(face = "bold", size = 8.5,
                                     margin = margin(t = 0, b = 2)),
        plot.title.position = "plot",
        axis.title   = element_text(size = 7.5),
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

# Okabe-Ito: distinguishable under the common forms of colour vision deficiency.
OI <- c(orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00",
        purple = "#CC79A7", grey = "#666666")
DIR_COL   <- c(adverse = unname(OI["vermillion"]), protective = unname(OI["blue"]),
               neutral = "grey60")
COH_COL   <- c("TCGA-KIRC" = unname(OI["vermillion"]), "CPTAC-3" = unname(OI["blue"]),
               "TCGA-KIRP" = unname(OI["green"]),      "TCGA-KICH" = unname(OI["purple"]))
BIO_COL   <- c(lncRNA = unname(OI["vermillion"]), `protein-coding` = "#333333")
# Fixes the order of the two library protocols.
PROTO_COL <- c(`polyA` = unname(OI["skyblue"]),
               `ribo-depleted total RNA` = unname(OI["orange"]))
STATUS_LEVELS <- c("robust", "LOST on adjustment", "gained on adjustment",
                   "not significant")
STATUS_COL <- c(robust = unname(OI["green"]),
                `LOST on adjustment` = unname(OI["vermillion"]),
                `gained on adjustment` = unname(OI["blue"]),
                `not significant` = "grey72")
MODEL_COL <- c(`Clinical` = "#4D4D4D", `Clinical + modules` = unname(OI["orange"]),
               `Augmented` = unname(OI["skyblue"]),
               `Augmented + modules` = unname(OI["green"]),
               `Null (no covariates)` = "grey78")
# Explicit breaks for log10 hazard-ratio axes. The defaults label too few
# values over this range.
HR_BRK <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 2)
ANN <- function(lab, x = -Inf, y = Inf, hjust = -0.12, vjust = 1.3, size = 2.2)
  annotate("text", x = x, y = y, hjust = hjust, vjust = vjust, size = size,
           family = FIG_FONT, colour = "grey15", label = lab, lineheight = 0.95)
# Lower-case panel letters, set per panel because a wrapped row is a single
# element for automatic tagging. The tag sits in the top margin, clear of the
# panel title.
TG <- function(p, ch) p + labs(tag = ch) +
  theme(plot.tag = element_text(face = "bold", size = 9.5, family = FIG_FONT,
                                hjust = 0, vjust = 1),
        plot.tag.position = "topleft", plot.tag.location = "margin",
        plot.margin = margin(t = 13, r = 6, b = 3, l = 4))
# Wrap a row so that panels align within the row but not across rows, and one
# panel with long labels does not shift every other row.
ROW <- function(x) wrap_elements(full = x)

msg("Figure font: ", FIG_FONT, "; width ", round(W, 3), " in (183 mm)")

# ==== Figure 1: study design and cohorts ====================================
banner("Figure 1 | design and cohorts")
cons <- R_("01_consort_flow.tsv")
cn   <- function(k) as.integer(cons[step == k, n])
t1d  <- R_("01_table1_cohort.tsv")
t1v  <- R_("08_validation_table1.tsv")
sst  <- R_("11_subtype_cohort_summary.tsv")
cq   <- R_("11_cohort_quality_metrics.tsv")

# -- a: CONSORT flow, every count read from 01_consort_flow.tsv ---------------
n_files   <- cn("gdc_files_total")
n_notprim <- cn("excluded_not_primary_tumour")
n_fail    <- cn("excluded_failed_library")
n_dup     <- cn("excluded_duplicate_patient_samples")
n_onepp   <- cn("tumour_samples_one_per_patient")
n_bados   <- cn("excluded_nonpositive_os_time")
n_final   <- cn("cohort_survival_analysis")
n_prim    <- n_files - n_notprim
n_pass    <- n_prim - n_fail
stopifnot(n_pass - n_dup == n_onepp, n_onepp - n_bados == n_final)

main <- data.table(
  y = c(5, 4, 3, 2, 1),
  lab = c(sprintf("GDC STAR-Counts files,\nTCGA-KIRC (n = %d)", n_files),
          sprintf("Primary tumour libraries\n(n = %d)", n_prim),
          sprintf("Libraries passing quality\n(n = %d)", n_pass),
          sprintf("One tumour per patient\n(n = %d)", n_onepp),
          sprintf("Discovery cohort analysed\n(n = %d, %d deaths)",
                  n_final, as.integer(t1d$deaths))))
main[, fill := c(rep("step", 4), "final")]
excl <- data.table(
  y = c(4.5, 3.5, 2.5, 1.5),
  eh = c(0.30, 0.42, 0.30, 0.30),
  lab = c(sprintf("Not primary tumour\n(n = %d)", n_notprim),
          sprintf("Failed library (n = %d):\n< %g M assigned reads or\n> %g%% non-feature",
                  n_fail, MIN_ASSIGNED_READS / 1e6, MAX_NOFEATURE_PCT),
          sprintf("Duplicate patient aliquot\n(n = %d)", n_dup),
          sprintf("Non-positive survival time\n(n = %d)", n_bados)))
# Two columns that do not overlap: the main flow on [0.2, 5.5], the exclusions
# on [5.85, 10.4], joined by an elbow from the spine at XM.
BW <- 2.65; BH <- 0.33; EW <- 2.275
XM <- 2.85; XE <- 8.125
p1a <- ggplot() +
  geom_segment(data = main[y > 1], linewidth = 0.35, colour = "grey35",
               aes(x = XM, xend = XM, y = y - BH, yend = y - 1 + BH),
               arrow = arrow(length = unit(1.3, "mm"), type = "closed")) +
  geom_segment(data = excl, linewidth = 0.3, colour = "grey55",
               aes(x = XM, xend = XE - EW, y = y, yend = y)) +
  geom_rect(data = excl, fill = "grey95", colour = "grey60", linewidth = 0.25,
            aes(xmin = XE - EW, xmax = XE + EW, ymin = y - eh, ymax = y + eh)) +
  geom_text(data = excl, aes(XE, y, label = lab), size = 2.25,
            family = FIG_FONT, lineheight = 0.95, colour = "grey15") +
  geom_rect(data = main, colour = "grey25", linewidth = 0.3,
            aes(xmin = XM - BW, xmax = XM + BW, ymin = y - BH, ymax = y + BH,
                fill = fill)) +
  geom_text(data = main, aes(XM, y, label = lab), size = 2.4,
            family = FIG_FONT, lineheight = 0.95, colour = "grey10") +
  scale_fill_manual(values = c(step = "white", final = "#DDEAF5"), guide = "none") +
  scale_x_continuous(limits = c(0, 10.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 5.5), expand = c(0, 0)) +
  labs(title = "Discovery cohort", x = NULL, y = NULL) +
  theme_void(base_size = 8, base_family = FIG_FONT) +
  theme(plot.title = element_text(face = "bold", size = 8.5,
                                  margin = margin(t = 2, b = 4)),
        plot.title.position = "plot", plot.margin = margin(4, 4, 3, 4))

# -- b: the four cohorts, patients and deaths ---------------------------------
coh <- rbind(
  data.table(cohort = "TCGA-KIRC", role = "discovery (clear cell)",
             n = as.integer(t1d$n), events = as.integer(t1d$deaths)),
  data.table(cohort = "CPTAC-3", role = "external validation (clear cell)",
             n = as.integer(t1v$n), events = as.integer(t1v$deaths)),
  sst[cohort %in% c("TCGA-KIRP", "TCGA-KICH"),
      .(cohort, role = "renal subtype", n = as.integer(n_scored),
        events = as.integer(events))])
coh[, cohort := factor(cohort, levels = names(COH_COL))]
coh[, role := factor(role, levels = unique(role))]
coh[, tick := paste0(cohort, "\n", wrap_lab(as.character(role), 14))]
coh[, tick := factor(tick, levels = coh[order(cohort), tick])]
p1b <- ggplot(coh, aes(tick, n, fill = cohort)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%d\n(%d deaths)", n, events)), vjust = -0.2,
            size = 2.2, family = FIG_FONT, lineheight = 0.9, colour = "grey10") +
  scale_fill_manual(values = COH_COL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.26))) +
  labs(title = "Cohorts", x = NULL, y = "Patients") +
  theme(axis.text.x = element_text(size = 6.6, lineheight = 0.9))

# -- c: library protocol and the non-feature fraction per cohort --------------
# Bars are filled by cohort, as in panel b. The protocol is in the tick label.
cq[, cohort := factor(cohort, levels = names(COH_COL))]
cq[, protocol := factor(protocol, levels = names(PROTO_COL))]
cq[, tick := paste0(cohort, "\n", wrap_lab(as.character(protocol), 13))]
cq[, tick := factor(tick, levels = cq[order(cohort), tick])]
# Each cohort used one protocol, so protocol and cohort are confounded here.
p1c <- ggplot(cq, aes(tick, noFeature_median, fill = cohort)) +
  geom_col(width = 0.62) +
  geom_errorbar(aes(ymin = noFeature_q25, ymax = noFeature_q75), width = 0.16,
                linewidth = 0.35, colour = "grey25") +
  geom_text(aes(y = noFeature_q75, label = sprintf("%.1f%%", noFeature_median)),
            vjust = -0.6, size = 2.1, family = FIG_FONT, colour = "grey10") +
  scale_fill_manual(values = COH_COL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Non-feature read fraction", x = NULL,
       y = "Reads in no annotated\nfeature (%)") +
  theme(axis.text.x = element_text(size = 6.4, lineheight = 0.9))

fig1 <- (TG(p1a, "a") | (TG(p1b, "b") / TG(p1c, "c"))) +
  plot_layout(widths = c(1.1, 1))
save_fig(fig1, "Figure1_design_cohorts", W, 4.7)
msg("  Figure 1 written")

# ==== Figure 2: the library-quality axis ===================================
banner("Figure 2 | lncRNA axis")
ax <- R_("27_axis_scores_per_sample.tsv")
pg <- R_("07_per_gene_quality_correlation.tsv")
ps <- R_("07_per_gene_quality_summary.tsv")
nx <- R_("23_axis_nested_extended.tsv")
bv <- R_("23_batch_variance.tsv")
nd <- R_("23_noFeature_determinants.tsv")

# -- a and b: the two leading axes against the metric, on identical axes ------
alnc <- ax[is.finite(lnc_axis_z) & is.finite(pct_noFeature)]
apc  <- ax[is.finite(pc_axis_z)  & is.finite(pct_noFeature)]
XL <- range(c(alnc$pct_noFeature, apc$pct_noFeature))
YL <- range(c(alnc$lnc_axis_z, apc$pc_axis_z))
rho_l <- cor(alnc$lnc_axis_z, alnc$pct_noFeature, method = "spearman")
rho_p <- cor(apc$pc_axis_z,   apc$pct_noFeature,  method = "spearman")

axis_panel <- function(d, yv, col, ttl, rho, ylab) {
  ggplot(d, aes(pct_noFeature, .data[[yv]])) +
    geom_point(alpha = 0.45, size = 0.6, colour = col, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey15",
                linewidth = 0.4) +
    coord_cartesian(xlim = XL, ylim = YL) +
    ANN(sprintf("Spearman rho = %.3f\nn = %d", rho, nrow(d))) +
    labs(title = ttl, x = "Reads in no annotated feature (%)", y = ylab)
}
p2a <- axis_panel(alnc, "lnc_axis_z", unname(BIO_COL["lncRNA"]),
                  "lncRNA axis", rho_l, "lncRNA PC1 (z)")
p2b <- axis_panel(apc, "pc_axis_z", unname(BIO_COL["protein-coding"]),
                  "Protein-coding axis", rho_p, "Protein-coding PC1 (z)")

# -- c: per-gene |rho| with the metric, lncRNA against protein-coding ---------
pg[, biotype := bt(biotype)]
ps[, biotype := bt(biotype)]
p2c <- ggplot(pg, aes(abs(rho_noFeature), colour = biotype, fill = biotype)) +
  geom_density(alpha = 0.22, linewidth = 0.5, adjust = 1.1) +
  geom_vline(data = ps, aes(xintercept = median_abs_rho, colour = biotype),
             linetype = 2, linewidth = 0.4, show.legend = FALSE) +
  scale_colour_manual(values = BIO_COL, name = NULL) +
  scale_fill_manual(values = BIO_COL, name = NULL) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.36))) +
  # A padded white label keeps the text clear of the dashed median lines.
  annotate("label", x = Inf, y = Inf, hjust = 1.06, vjust = 1.06, size = 2.05,
           family = FIG_FONT, colour = "grey15", lineheight = 0.95,
           fill = "white", linewidth = 0, label.padding = unit(1.1, "mm"),
           label = paste(sprintf("%s (%s genes):\nmedian absolute rho %.2f, %.0f%% above 0.3",
                                 ps$biotype,
                                 format(ps$n_genes, big.mark = ",", trim = TRUE),
                                 ps$median_abs_rho, 100 * ps$frac_abs_rho_gt_0.3),
                         collapse = "\n")) +
  labs(title = "Per-gene correlation",
       x = "Absolute rho with the non-feature fraction", y = "Density") +
  theme(legend.position = "bottom")

# -- d: eight nested adjustments of the axis ---------------------------------
nx[, lab := factor(wrap_lab(model, 30), levels = rev(wrap_lab(model, 30)))]
p2d <- ggplot(nx, aes(HR, lab)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.3) +
  ci_h(aes(xmin = lo, xmax = hi), width = 0.16, linewidth = 0.4, colour = "grey30") +
  geom_point(aes(colour = p < 0.05), size = 1.9) +
  scale_colour_manual(values = c(`TRUE` = unname(DIR_COL["adverse"]),
                                 `FALSE` = "grey55"), name = NULL,
                      labels = c(`TRUE` = "p < 0.05", `FALSE` = "p >= 0.05"),
                      breaks = c("TRUE", "FALSE")) +
  scale_x_log10(breaks = c(0.6, 0.8, 1, 1.25, 1.6)) +
  labs(title = "Nested Cox models", x = "Hazard ratio per 1 SD (95% CI)", y = NULL) +
  theme(axis.text.y = element_text(size = 6.2, lineheight = 0.9),
        legend.position = "bottom")

# -- e: variance explained by batch structure and RIN -------------------------
FAC <- c(plate = "Plate", batch = "Batch", tss = "Tissue\nsource site",
         rin = "RIN")
VAR <- c(pct_noFeature = "Non-feature fraction", lnc_axis = "lncRNA PC1")
ve <- bv[factor %in% names(FAC), .(variable, term = factor, R2, n)]
ve <- rbind(ve, nd[term == "rin", .(variable = "pct_noFeature", term, R2 = R2_alone, n)])
joint <- nd[term == "rin + plate + tss"]
ve <- rbind(ve, data.table(variable = "pct_noFeature", term = "joint",
                           R2 = joint$R2_alone, n = joint$n))
ve[, term_lab := factor(c(FAC, joint = "RIN + plate\n+ site")[term],
                        levels = c(unname(FAC), "RIN + plate\n+ site"))]
ve[, var_lab := factor(VAR[variable], levels = unname(VAR))]
p2e <- ggplot(ve, aes(term_lab, 100 * R2, fill = var_lab)) +
  geom_col(width = 0.7, position = position_dodge2(preserve = "single",
                                                   padding = 0.12)) +
  geom_text(aes(label = sprintf("%.0f", 100 * R2)), vjust = -0.35, size = 2.0,
            family = FIG_FONT, colour = "grey15",
            position = position_dodge2(width = 0.7, preserve = "single",
                                       padding = 0.12)) +
  # Blue for the metric, because dark grey means protein-coding elsewhere in
  # this figure.
  scale_fill_manual(values = c(unname(OI["blue"]), unname(BIO_COL["lncRNA"])),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  # The bars mix quantities: one-way ANOVA R2 for plate, batch and site, and
  # regression R2 for RIN and the joint model (metric only). The joint bar is
  # one model, not a sum of the others.
  labs(title = "Variance explained", x = NULL, y = "Variance explained (%)") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 6.4, lineheight = 0.9))

# -- f: the axis against the metric in every sample set -----------------------
# Two axis definitions are drawn in separate facets:
#   projected: the discovery PC1 direction applied to the set, with its sign
#     fixed on the discovery matrix. This is the axis used in other analyses.
#   own: each set's own leading component, with arbitrary sign. Stages 25 and
#     35 anchor it and set pc1_sign_interpretable to FALSE where the anchors
#     are weak or disagree. Those rows are plotted as magnitudes (open symbols).
# The stage 25 and 35 tables score the same discovery rotation
# (35_axis_recompute_check.tsv), so they are combined into one table.
nrm <- Rx("25_normal_axis_projection.tsv")
sbt <- Rx("35_axis_projection_subtypes.tsv")
if (is.null(nrm)) MISSING <- unique(c(MISSING, "25_normal_axis_projection.tsv"))
if (is.null(sbt)) MISSING <- unique(c(MISSING, "35_axis_projection_subtypes.tsv"))
axp <- rbindlist(list(nrm, sbt), use.names = TRUE, fill = TRUE)
AXT <- c(projected = "Discovery axis projected onto the set",
         own       = "The set's own leading component")
SGN <- c(`TRUE`  = "signed rho (sign anchored)",
         `FALSE` = "absolute rho (the sign is a convention)")
PROTO_SHORT <- c(`polyA` = "poly(A)",
                 `ribo-depleted total RNA` = "ribo-depleted")
if (nrow(axp)) {
  f2 <- axp[, .(cohort, tissue, matrix = bt(matrix), n, rho = rho_noFeature,
                p = p_noFeature,
                interpretable = as.logical(pc1_sign_interpretable),
                axis = unname(AXT[fifelse(axis_type == "projected_discovery_PC1",
                                          "projected", "own")]))]
} else {
  # Without the projection tables only each set's own component can be drawn,
  # and those tables carry no sign flag.
  f2 <- rbind(
    R_("07_axis_vs_purity_qc_correlations.tsv")[variable == "pct_noFeature",
      .(cohort = "TCGA-KIRC", tissue = "tumour", matrix = bt(matrix), n,
        rho = spearman_rho, p, interpretable = TRUE, axis = unname(AXT["own"]))],
    R_("09_axis_replication_cptac.tsv")[
      , .(cohort, tissue = "tumour", matrix = bt(matrix), n = n_samples,
          rho = rho_noFeature, p = p_noFeature, interpretable = TRUE,
          axis = unname(AXT["own"]))],
    R_("11_axis_replication_subtypes.tsv")[
      , .(cohort, tissue = "tumour", matrix = bt(matrix), n,
          rho = rho_noFeature, p = p_noFeature, interpretable = TRUE,
          axis = unname(AXT["own"]))])
}
f2 <- merge(f2, cq[, .(cohort = as.character(cohort),
                       protocol = as.character(protocol))],
            by = "cohort", all.x = TRUE)
f2[, rho_plot := fifelse(interpretable, rho, abs(rho))]
f2[, sgn := factor(SGN[as.character(interpretable)], levels = unname(SGN))]
f2[, tis := fifelse(tissue == "normal", "matched normal", "tumour")]
# Protocol is carried in the row label so that shape can mark magnitudes. The
# label gives an n range where the two matrices differ in sample count.
f2[, lab := sprintf("%s, %s\n%s, n = %s", cohort, tis,
                    PROTO_SHORT[protocol],
                    if (min(n) == max(n)) as.character(min(n))
                    else sprintf("%d-%d", min(n), max(n))),
   by = .(cohort, tis)]
# Rows are ordered by the projected lncRNA correlation.
ordq <- f2[axis == AXT["projected"] & matrix == "lncRNA"]
if (!nrow(ordq)) ordq <- f2[matrix == "lncRNA"]
f2[, lab := factor(lab, levels = unique(ordq[order(rho_plot), lab]))]
f2[, matrix := factor(matrix, levels = names(BIO_COL))]
f2[, axis := factor(axis, levels = unname(AXT))]
p2f <- ggplot(f2, aes(rho_plot, lab, colour = matrix, shape = sgn,
                      group = matrix)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
  # group = matrix stops position_dodge from splitting rows by the sign flag,
  # which would pull points off their row centres.
  geom_point(size = 2.4, stroke = 0.7, fill = "white",
             position = position_dodge(0.55)) +
  facet_wrap(~ axis, nrow = 1) +
  scale_colour_manual(values = BIO_COL, name = "Expression matrix") +
  scale_shape_manual(values = setNames(c(16, 21), unname(SGN)),
                     name = "Quantity plotted", drop = FALSE) +
  scale_x_continuous(limits = c(min(-0.15, min(f2$rho_plot) - 0.05), 1),
                     breaks = seq(-0.25, 1, 0.25)) +
  labs(title = "Axis and metric across sample sets",
       x = "Spearman rho with the non-feature fraction (magnitude where the sign is a convention)",
       y = NULL) +
  theme(axis.text.y = element_text(size = 6.4, lineheight = 0.9),
        legend.position = "bottom", legend.box = "horizontal")

# Panels are lettered in order of first citation in the text: p2f is d, p2d is
# e and p2e is f.
fig2 <- (ROW(TG(p2a, "a") | TG(p2b, "b") | TG(p2c, "c")) /
         ROW(TG(p2f, "d")) /
         ROW(TG(p2d, "e") | TG(p2e, "f"))) +
  plot_layout(heights = c(1, 1.2, 1.05))
save_fig(fig2, "Figure2_quality_axis", W, 9.0)
msg("  Figure 2 written")

# ==== Figure 4: modules (objects fig3, p3*) ================================
banner("Figure 4 | modules")
sens <- R_("12_technical_adjustment_sensitivity.tsv")[covariate_set == "clinical"]
sens[, biotype := bt(biotype)]
sens[, status := factor(status, levels = STATUS_LEVELS)]
p3a <- ggplot(sens, aes(HR_unadjusted, HR_adjusted, colour = status,
                        shape = biotype)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60",
              linewidth = 0.3) +
  geom_hline(yintercept = 1, linewidth = 0.2, colour = "grey88") +
  geom_vline(xintercept = 1, linewidth = 0.2, colour = "grey88") +
  geom_point(size = 1.9) +
  # Explicit breaks so the adverse half of the range is labelled.
  scale_x_log10(breaks = HR_BRK) + scale_y_log10(breaks = HR_BRK) +
  scale_colour_manual(values = STATUS_COL, name = NULL, drop = FALSE) +
  scale_shape_manual(values = c(lncRNA = 16, `protein-coding` = 17), name = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, order = 1),
         shape = guide_legend(nrow = 2, order = 2)) +
  labs(title = "Observed against adjusted",
       x = "HR per 1 SD without adjustment",
       y = "HR per 1 SD with adjustment") +
  theme(legend.position = "bottom", legend.spacing.x = unit(0.1, "cm"))

# -- b: lncRNA module-trait correlations --------------------------------------
mt  <- R_("02_lncRNA_module_trait.tsv")
mtl <- melt(mt, id.vars = "module", measure.vars = patterns("^r_"),
            variable.name = "trait", value.name = "r")
mtl[, trait := sub("^r_", "", trait)]
mtp <- melt(mt, id.vars = "module", measure.vars = patterns("^p_"),
            variable.name = "trait", value.name = "p")
mtp[, trait := sub("^p_", "", trait)]
mtl <- merge(mtl, mtp, by = c("module", "trait"))
mtl <- mtl[module != "MEgrey"]
mtl[, mod := sub("^ME", "", module)]
TR <- c(age = "Age", male = "Male sex", stage = "Stage", grade = "Grade")
mtl[, trait := factor(TR[trait], levels = unname(TR))]
ord <- mtl[trait == "Grade"][order(r), mod]
mtl[, mod := factor(mod, levels = ord)]
p3b <- ggplot(mtl, aes(trait, mod, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(p.adjust(p, "BH") < FDR_ALPHA, "*", "")),
            size = 2.4, vjust = 0.72, family = FIG_FONT, colour = "grey10") +
  scale_fill_gradient2(low = unname(OI["blue"]), mid = "white",
                       high = unname(OI["vermillion"]), midpoint = 0,
                       limits = c(-0.5, 0.5), name = "Pearson r") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Module-trait correlation", x = NULL, y = NULL) +
  theme(axis.text.y = element_text(size = 6.4), legend.position = "right",
        legend.key.width = unit(0.25, "cm"), panel.grid = element_blank())

# -- c: principal-specification forest, both networks --------------------------
sv <- rbind(R_("03_mRNA_module_survival.tsv")[, biotype := "protein-coding"],
            R_("03_lncRNA_module_survival.tsv")[, biotype := "lncRNA"], fill = TRUE)
sv[, sig := fdr_full < FDR_ALPHA]
sv[, cls := fcase(sig & HR_full > 1, "adverse (FDR < 0.05)",
                  sig & HR_full < 1, "protective (FDR < 0.05)",
                  default = "not significant")]
sv[, cls := factor(cls, levels = c("adverse (FDR < 0.05)",
                                   "protective (FDR < 0.05)", "not significant"))]
sv[, biotype := factor(biotype, levels = names(BIO_COL))]
# Module colours are reused between networks, so the y key is biotype-qualified.
# The printed label is the module name alone.
sv[, key := paste(biotype, module)]
sv[, key := factor(key, levels = sv[order(biotype, HR_full), key])]
p3c <- ggplot(sv, aes(HR_full, key, colour = cls)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.3) +
  ci_h(aes(xmin = HR_full_lo, xmax = HR_full_hi), width = 0.2, linewidth = 0.32,
       colour = "grey45") +
  geom_point(size = 1.6) +
  facet_grid(biotype ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^(lncRNA|protein-coding) ", "", x)) +
  scale_x_log10(breaks = c(0.6, 0.8, 1, 1.25, 1.6)) +
  scale_colour_manual(values = c(`adverse (FDR < 0.05)` = unname(DIR_COL["adverse"]),
                                 `protective (FDR < 0.05)` = unname(DIR_COL["protective"]),
                                 `not significant` = "grey65"), name = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "Principal specification", x = "Hazard ratio per 1 SD (95% CI)",
       y = NULL) +
  theme(axis.text.y = element_text(size = 6.2), legend.position = "bottom",
        legend.text = element_text(size = 6.4))

# Lettered in order of citation: p3c is b and p3b is c.
fig3 <- (ROW(TG(p3a, "a") / TG(p3b, "c")) | ROW(TG(p3c, "b"))) +
  plot_layout(widths = c(1, 0.92))
save_fig(fig3, "Figure4_modules", W, 7.2)
msg("  Figure 4 (modules) written")

# ==== Figure 5: replication and subtypes (objects fig4, p4*) ===============
banner("Figure 5 | replication and subtypes")
rep12 <- R_("09_module_replication.tsv")
rep6  <- R_("09_module_replication_parsimonious.tsv")
for (d in list(rep12, rep6)) d[, biotype := bt(biotype)]
disc <- sv[, .(biotype = as.character(biotype), module, HR = HR_full,
               lo = HR_full_lo, hi = HR_full_hi)]
r4 <- rbind(
  merge(rep12[, .(biotype, module)], disc, by = c("biotype", "module"))[
    , est := "TCGA-KIRC discovery (12-parameter)"],
  rep12[, .(biotype, module, HR = HR_cptac, lo, hi,
            est = "CPTAC-3 (12-parameter)")],
  rep6[,  .(biotype, module, HR = HR_cptac, lo, hi,
            est = "CPTAC-3 (6-parameter)")])
EST_LEV <- c("TCGA-KIRC discovery (12-parameter)", "CPTAC-3 (12-parameter)",
             "CPTAC-3 (6-parameter)")
r4[, est := factor(est, levels = EST_LEV)]
r4[, lab := paste0(biotype, ": ", module)]
r4[, lab := factor(lab, levels = r4[est == EST_LEV[1]][order(HR), lab])]
p4a <- ggplot(r4, aes(HR, lab, colour = est)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.3) +
  ci_h(aes(xmin = lo, xmax = hi), width = 0.22, linewidth = 0.35,
       position = position_dodge(0.62), show.legend = FALSE) +
  geom_point(size = 1.9, position = position_dodge(0.62)) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2.2)) +
  scale_colour_manual(values = c(unname(COH_COL["TCGA-KIRC"]),
                                 unname(COH_COL["CPTAC-3"]), unname(OI["skyblue"])),
                      name = NULL) +
  labs(title = "Estimates in CPTAC-3", x = "Hazard ratio per 1 SD (95% CI)",
       y = NULL) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6.6)) +
  guides(colour = guide_legend(nrow = 3))

# -- b: estimates across the renal subtypes -----------------------------------
st  <- R_("11_subtype_module_effects.tsv")
het <- R_("11_subtype_heterogeneity.tsv")
st[, biotype := bt(biotype)]; het[, biotype := bt(biotype)]
st <- merge(st, het[, .(biotype, module, p_heterogeneity, I2, k_het = k,
                        het_cohorts = cohorts)],
            by = c("biotype", "module"), all.x = TRUE)
# Q and I2 come from the cohorts stage 11 could power, so the row label carries
# k and the legend names the cohorts inside and outside the test.
n_coh_plot <- uniqueN(st$cohort)
st[, lab := sprintf("%s: %s\n(Q p = %.2g, I2 = %g%%, %d of %d cohorts)", biotype,
                    module, p_heterogeneity, I2, k_het, n_coh_plot)]
st[, lab := factor(lab, levels = unique(st[cohort == "TCGA-KIRC"][order(HR), lab]))]
st[, cohort := factor(cohort, levels = c("TCGA-KIRC", "TCGA-KIRP", "TCGA-KICH"))]
# The underpowered subtype has very wide intervals, so the axis is capped at
# twice the widest adequately powered bound and truncation is declared.
XCAP <- c(min(st$lo) / 1.1, ceiling(2 * max(st[adequate == TRUE, hi])))
# Raise the cap so the edge marker clears the largest estimate still on the
# panel.
XCAP[2] <- round(max(XCAP[2], 1.35 * max(st[HR <= XCAP[2], HR])), 1)
st[, cut_hi := hi > XCAP[2]]
trunc_hi <- st[cut_hi == TRUE]
p4b <- ggplot(st, aes(HR, lab, colour = cohort, alpha = adequate)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.3) +
  # Truncated intervals keep their true bound and are clipped at the panel edge,
  # where a ">" marker shows that they continue. The marker layer uses the full
  # table with NA elsewhere so that position_dodge sees every group. Bar width
  # is constant so that dodging stays aligned.
  ci_h(aes(xmin = lo, xmax = hi), width = 0.2, linewidth = 0.32,
       position = position_dodge(0.62), show.legend = FALSE) +
  { if (nrow(trunc_hi))
      geom_point(aes(x = fifelse(cut_hi, XCAP[2] * 0.996, NA_real_)), shape = 62,
                 size = 2.6, stroke = 0.7, show.legend = FALSE,
                 position = position_dodge(0.62), na.rm = TRUE) } +
  geom_point(size = 1.8, position = position_dodge(0.62)) +
  coord_cartesian(xlim = XCAP) +
  # No break above 4: the 4 and 6 labels collide at this width.
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  scale_colour_manual(values = COH_COL, name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none") +
  labs(title = "Renal subtypes", x = "Hazard ratio per 1 SD (95% CI)", y = NULL) +
  theme(axis.text.y = element_text(size = 6, lineheight = 0.9),
        legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1))

# -- c: module preservation of the prognostic modules -------------------------
pres <- rbind(R_("24_module_preservation_lncRNA.tsv"),
              R_("24_module_preservation_mRNA.tsv"), fill = TRUE)
pres[, biotype := bt(network)]
keep5 <- unique(rep12[, .(biotype, module)])
pr5 <- merge(pres, keep5, by = c("biotype", "module"))
pr5[, lab := paste0(biotype, ": ", module)]
pr5[, lab := factor(lab, levels = levels(r4$lab))]
pr5[, test_cohort := factor(test_cohort, levels = c("CPTAC-3", "TCGA-KIRP", "TCGA-KICH"))]
p4c <- ggplot(pr5, aes(Zsummary, lab, colour = test_cohort)) +
  geom_vline(xintercept = c(2, 10), linetype = c(3, 2), colour = "grey45",
             linewidth = 0.35) +
  geom_point(size = 2.2, position = position_dodge(0.6)) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
  scale_colour_manual(values = COH_COL, name = NULL) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  # Right of its line: the left panel edge would clip a label placed before it.
  annotate("text", x = 2, y = Inf, label = "Z = 2 (weak)", size = 2.0, vjust = 1.3,
           hjust = -0.06, family = FIG_FONT, colour = "grey35") +
  annotate("text", x = 10, y = Inf, label = "Z = 10 (strong)", size = 2.0,
           vjust = 1.3, hjust = -0.06, family = FIG_FONT, colour = "grey35") +
  labs(title = "Module preservation", x = "Zsummary (log scale)", y = NULL) +
  theme(axis.text.y = element_text(size = 6.6), legend.position = "bottom")

fig4 <- (ROW(TG(p4a, "a") | TG(p4b, "b")) / ROW(TG(p4c, "c"))) +
  plot_layout(heights = c(1.2, 1))
save_fig(fig4, "Figure5_replication_subtypes", W, 6.6)
msg("  Figure 5 (replication) written")

# ==== Figure 7: prediction (objects fig5, p5*) ==============================
banner("Figure 7 | prediction")
d09 <- R_("09_delta_cindex_validation.tsv")[comparator == "clinical"]
d05 <- R_("05_delta_cindex.tsv")
fwa <- R_("12_foldwise_per_assignment.tsv")[comparator == "clinical"]
g05 <- function(cmp) d05[comparison == cmp]
gd9 <- function(std, coh) d09[standardisation == std & grepl(coh, cohort)]

BOOT <- "bootstrap 95% CI"; RANGE <- "range across estimates"
row_boot <- function(lab, x) data.table(lab = lab, delta = x[[1]], lo = x[[2]],
                                        hi = x[[3]], kind = BOOT,
                                        n = x[[4]], events = x[[5]])
a <- gd9("cohort", "KIRC"); b <- gd9("cohort", "CPTAC"); cc <- gd9("discovery", "CPTAC")
e1 <- g05("clinical_eig - clinical"); h1 <- g05("clinical_hub - clinical")
d5 <- rbindlist(list(
  row_boot("Discovery, apparent (locked model)",
           list(a$delta_C, a$lo, a$hi, a$n, a$events)),
  data.table(lab = sprintf("Discovery, %d x %d cross-validation,\nmodules fixed on the full cohort",
                           ML_N_REPEATS, ML_N_FOLDS),
             delta = e1$delta_mean, lo = e1$range_lo, hi = e1$range_hi,
             kind = sprintf("range across %d repeats", e1$n_repeats),
             n = e1$n, events = e1$events),
  row_boot(sprintf("Discovery, %d x %d cross-validation, modules fixed\n(bootstrap of the averaged predictor)",
                   ML_N_REPEATS, ML_N_FOLDS),
           list(e1$boot_delta, e1$boot_lo, e1$boot_hi, e1$n, e1$events)),
  data.table(lab = sprintf("Discovery, %d x %d cross-validation,\nhub lncRNAs selected inside each fold",
                           ML_N_REPEATS, ML_N_FOLDS),
             delta = h1$delta_mean, lo = h1$range_lo, hi = h1$range_hi,
             kind = sprintf("range across %d repeats", h1$n_repeats),
             n = h1$n, events = h1$events),
  row_boot(sprintf("Discovery, %d x %d cross-validation, hub lncRNAs\n(bootstrap of the averaged predictor)",
                   ML_N_REPEATS, ML_N_FOLDS),
           list(h1$boot_delta, h1$boot_lo, h1$boot_hi, h1$n, h1$events)),
  data.table(lab = sprintf("Discovery, both networks rebuilt inside each fold\n(%d folds x %d assignments)",
                           FOLDWISE_K, FOLDWISE_N_ASSIGN),
             delta = mean(fwa$delta), lo = min(fwa$delta), hi = max(fwa$delta),
             kind = sprintf("range across %d assignments", nrow(fwa)),
             n = fwa$n[1], events = fwa$events[1]),
  row_boot("CPTAC-3, locked model,\nscores standardised within cohort",
           list(b$delta_C, b$lo, b$hi, b$n, b$events)),
  row_boot("CPTAC-3, locked model,\ndiscovery-fixed standardisation",
           list(cc$delta_C, cc$lo, cc$hi, cc$n, cc$events))))
d5[, lab := factor(lab, levels = rev(lab))]
d5[, kind2 := ifelse(kind == BOOT, BOOT, "range (not a confidence interval)")]
p5a <- ggplot(d5, aes(delta, lab, colour = kind2, shape = kind2)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
  ci_h(aes(xmin = lo, xmax = hi), width = 0.2, linewidth = 0.4,
       show.legend = FALSE) +
  geom_point(size = 2.1) +
  # Four decimals, the precision of the source tables.
  geom_text(aes(x = hi, label = sprintf("  %+.4f (%+.4f to %+.4f)", delta, lo, hi)),
            hjust = 0, size = 1.95, family = FIG_FONT, colour = "grey25",
            show.legend = FALSE) +
  scale_x_continuous(limits = c(min(d5$lo) - 0.005, max(d5$hi) + 0.075)) +
  # Neutral colours: gold means Clinical + modules in panels c to e.
  scale_colour_manual(values = c("#1A1A1A", "grey58"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  labs(title = "Increment over the clinical model",
       x = "Change in concordance index", y = NULL) +
  theme(axis.text.y = element_text(size = 6.2, lineheight = 0.9),
        legend.position = "bottom")

# -- b: discrimination of the four locked models in both cohorts --------------
ci <- R_("09_validation_cindex.tsv")[standardisation == "cohort" &
                                     model %in% c("comparator only", "+ module eigengenes")]
ci[, mlab := fcase(comparator == "clinical"  & model == "comparator only",     "Clinical",
                   comparator == "clinical"  & model == "+ module eigengenes", "Clinical + modules",
                   comparator == "augmented" & model == "comparator only",     "Augmented",
                   comparator == "augmented" & model == "+ module eigengenes", "Augmented + modules")]
ci[, mlab := factor(mlab, levels = rev(names(MODEL_COL)[1:4]))]
ci[, coh := factor(fifelse(grepl("KIRC", cohort), "TCGA-KIRC", "CPTAC-3"),
                   levels = c("TCGA-KIRC", "CPTAC-3"))]
p5b <- ggplot(ci, aes(C, mlab, colour = coh)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.3) +
  ci_h(aes(xmin = lo, xmax = hi), width = 0.18, linewidth = 0.4,
       position = position_dodge(0.5), show.legend = FALSE) +
  geom_point(size = 2, position = position_dodge(0.5)) +
  coord_cartesian(xlim = c(0.5, 0.9)) +
  scale_colour_manual(values = COH_COL, name = NULL) +
  # Augmented is the covariate set built by clinical_design(set = "augmented").
  labs(title = "Discrimination", x = "Concordance index (95% CI)", y = NULL) +
  theme(axis.text.y = element_text(size = 6.6), legend.position = "bottom")

# -- c: calibration at the primary horizon in CPTAC-3 -------------------------
MOD2 <- c(comparator = "Clinical", `comparator + modules` = "Clinical + modules")
cal <- R_("12_calibration_bins.tsv")[cohort == "validation" & comparator == "clinical" &
                                     standardisation == "cohort" & years == PRIMARY_HORIZON_YR]
cal[, mlab := factor(MOD2[model], levels = unname(MOD2))]
LIMC <- c(0, max(c(cal$predicted, cal$obs_hi), na.rm = TRUE) * 1.05)
p5c <- ggplot(cal, aes(predicted, observed, colour = mlab)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55",
              linewidth = 0.3) +
  geom_errorbar(aes(ymin = obs_lo, ymax = obs_hi), width = 0.012, linewidth = 0.3,
                show.legend = FALSE) +
  geom_line(linewidth = 0.35) + geom_point(size = 1.5) +
  # Equal limits and a square panel (aspect.ratio, since coord_fixed clips the
  # title under patchwork), so calibration reads as distance from the diagonal.
  coord_cartesian(xlim = LIMC, ylim = LIMC) +
  scale_colour_manual(values = MODEL_COL[1:2], name = NULL) +
  labs(title = "Calibration, CPTAC-3", x = "Predicted risk of death",
       y = "Observed risk") +
  theme(legend.position = "bottom", aspect.ratio = 1)

# -- d: decision curve at the primary horizon in CPTAC-3 ----------------------
STRAT <- c(comparator = "Clinical", `comparator + modules` = "Clinical + modules",
           `treat all` = "Treat all", `treat none` = "Treat none")
dca <- R_("12_validation_decision_curve.tsv")[baseline == "discovery" &
                                              comparator == "clinical" &
                                              standardisation == "cohort" &
                                              years == PRIMARY_HORIZON_YR]
dca[, strat := factor(STRAT[strategy], levels = unname(STRAT))]
# The y floor clears the lowest model curve. The treat-all reference may run
# off the panel.
dfloor <- function(d) {
  m <- d[strategy %in% c("comparator", "comparator + modules"), net_benefit]
  min(-0.02, min(m, na.rm = TRUE) * 1.15)
}
p5d <- ggplot(dca, aes(threshold, net_benefit, colour = strat, linetype = strat)) +
  geom_line(linewidth = 0.5) +
  coord_cartesian(ylim = c(dfloor(dca), max(dca$net_benefit, na.rm = TRUE) * 1.08)) +
  scale_colour_manual(values = c(unname(MODEL_COL[1:2]), "grey40", "grey72"),
                      name = NULL) +
  scale_linetype_manual(values = c(1, 1, 2, 3), name = NULL) +
  labs(title = "Decision curve, CPTAC-3", x = "Threshold probability",
       y = "Net benefit") +
  theme(legend.position = "bottom", legend.text = element_text(size = 6),
        legend.direction = "vertical", legend.key.height = unit(0.26, "cm"))

# -- e: Brier scores ----------------------------------------------------------
br <- R_("12_brier_scores.tsv")[comparator == "clinical" & standardisation == "cohort" &
                                years %in% EVAL_TIMES_YRS]
brl <- melt(br, id.vars = c("cohort", "years", "n", "events"),
            measure.vars = c("brier_clinical", "brier_combined", "brier_null"),
            variable.name = "model", value.name = "brier")
brl[, mlab := factor(c(brier_clinical = "Clinical",
                       brier_combined = "Clinical + modules",
                       brier_null = "Null (no covariates)")[as.character(model)],
                     levels = c("Clinical", "Clinical + modules", "Null (no covariates)"))]
brl[, coh := factor(fifelse(cohort == "discovery", "TCGA-KIRC", "CPTAC-3"),
                    levels = c("TCGA-KIRC", "CPTAC-3"))]
p5e <- ggplot(brl, aes(factor(years), brier, fill = mlab)) +
  geom_col(position = position_dodge(0.78), width = 0.68) +
  facet_wrap(~ coh) +
  scale_fill_manual(values = MODEL_COL[c("Clinical", "Clinical + modules",
                                         "Null (no covariates)")], name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Brier score", x = "Years since diagnosis", y = "Brier score") +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.key.height = unit(0.26, "cm"), legend.text = element_text(size = 6))

# Row 2 is taller because the calibration panel has a fixed aspect ratio.
# Lettered in order of citation: p5e is d and p5d is e.
fig5 <- (ROW(TG(p5a, "a")) /
         ROW((TG(p5b, "b") | TG(p5c, "c")) + plot_layout(widths = c(1.1, 1))) /
         ROW(TG(p5e, "d") | TG(p5d, "e"))) +
  # Row 1 needs extra height for eight two-line labels.
  plot_layout(heights = c(1.32, 1.38, 1.05))
save_fig(fig5, "Figure7_prediction", W, 9.1)
msg("  Figure 7 (prediction) written")

# ==== Figure 6: biology of the replicating modules ==========================
banner("Figure 6 | biology")
go  <- R_("06_GO_enrichment_top10_per_module.tsv")
svm <- R_("03_mRNA_module_survival.tsv")
svl <- R_("03_lncRNA_module_survival.tsv")
pc_mod <- catabolic_reference_module()
if (!pc_mod %in% go$module)
  pc_mod <- svm[HR_full < 1 & module %in% go$module][which.min(fdr_full), module]
gob <- go[module == pc_mod][order(p.adjust)][seq_len(min(8, .N))]
p6a <- ggplot(gob, aes(-log10(p.adjust), reorder(wrap_lab(Description, 34), -p.adjust))) +
  # Neutral grey: these are adjusted p values, and blue and orange-red carry the
  # sign of a correlation in panel c.
  geom_col(fill = "grey45", width = 0.68) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = sprintf("Protein-coding %s module", pc_mod),
       x = expression(-log[10]~adjusted~italic(p)), y = NULL) +
  theme(axis.text.y = element_text(size = 6, lineheight = 0.9))

# -- b: sign-split guilt-by-association for the adverse lncRNA module ---------
# One panel per lncRNA module: the top GO terms of its positively and of its
# negatively correlated protein-coding partners. Also used for Supplementary
# Figure 9a and b.
go_signed_panel <- function(gs, mod) {
  g <- gs[module == mod][order(sign, p.adjust)][, head(.SD, 6), by = sign]
  g[, lab := factor(ifelse(sign == "negative", "negatively correlated mRNA partners",
                                               "positively correlated mRNA partners"),
                    levels = c("positively correlated mRNA partners",
                               "negatively correlated mRNA partners"))]
  g[, Description := factor(wrap_lab(Description, 32),
                            levels = rev(unique(wrap_lab(Description, 32))))]
  ggplot(g, aes(-log10(p.adjust), Description, fill = lab)) +
    geom_col(width = 0.68) +
    facet_wrap(~ lab, ncol = 1, scales = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = setNames(c(unname(DIR_COL["adverse"]),
                                          unname(DIR_COL["protective"])),
                                        levels(g$lab)), guide = "none") +
    labs(title = sprintf("lncRNA %s partners", mod),
         x = expression(-log[10]~adjusted~italic(p)), y = NULL) +
    theme(axis.text.y = element_text(size = 5.8, lineheight = 0.9))
}
gs <- Rx("10_lncRNA_guilt_by_association_GO_signed.tsv")
p6b <- NULL
if (!is.null(gs)) {
  # The adverse lncRNA module with the smallest FDR that has a signed enrichment
  # table, located by its statistics because module colours change on rebuild.
  cand <- svl[HR_full > 1 & module %in% unique(gs$module)][order(fdr_full)]
  lnc_mod <- cand$module[1]
  p6b <- go_signed_panel(gs, lnc_mod)
} else MISSING <- unique(c(MISSING, "10_lncRNA_guilt_by_association_GO_signed.tsv"))

# -- c: coupling of the lncRNA eigengenes to the catabolic module -------------
cp <- R_("11_lnc_black_coupling_by_cohort.tsv")
cp[, cohort := factor(cohort, levels = names(COH_COL))]
cp[, module := factor(module, levels = cp[cohort == "TCGA-KIRC"][order(pearson_r), module])]
ref_pc <- unique(cp$reference_module)
p6c <- ggplot(cp, aes(pearson_r, module, colour = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
  ci_h(aes(xmin = CI_lo, xmax = CI_hi), width = 0.18, linewidth = 0.35,
       position = position_dodge(0.68), show.legend = FALSE) +
  geom_point(size = 2.1, position = position_dodge(0.68)) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_colour_manual(values = COH_COL, name = NULL) +
  labs(title = sprintf("Coupling to %s", paste(ref_pc, collapse = "/")),
       x = sprintf("Pearson r with the %s eigengene", paste(ref_pc, collapse = "/")),
       y = NULL) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

# -- gene-class composition, two annotation schemes (supplementary figure) ----
# Counts are genes: the GDC STAR-Counts matrices are gene-level.
comp <- R_("10_lncRNA_module_composition.tsv")
cl1 <- melt(comp, id.vars = "module", variable.name = "class", value.name = "n")
cl1[, scheme := "Gene name class"]
cls <- Rx("28_module_composition_by_class.tsv")
if (!is.null(cls)) {
  ncols <- grep("^n_", names(cls), value = TRUE)
  ncols <- setdiff(ncols, "n_genes")
  cl2 <- melt(cls[module %in% comp$module, c("module", ncols), with = FALSE],
              id.vars = "module", variable.name = "class", value.name = "n")
  cl2[, class := gsub("_", " ", sub("^n_", "", class))]
  # Antisense, divergent and intronic occur in both schemes with different
  # definitions, so positional classes are suffixed to keep legend keys distinct.
  cl2[, class := paste0(class, " (positional)")]
  cl2[, scheme := "Positional class (GENCODE v36)"]
  cl1[, class := as.character(class)]
  cl <- rbind(cl1, cl2)
} else {
  MISSING <- unique(c(MISSING, "28_module_composition_by_class.tsv"))
  cl1[, class := as.character(class)]; cl <- cl1
}
cl <- cl[n > 0]
cl[, class := sub("^clone-based \\(unnamed\\)$", "clone-based", class)]
cl[, class := sub("^lincRNA \\(named\\)$", "lincRNA", class)]
cl[, class := sub("^sense overlapping", "sense-overlap", class)]
cl[, scheme := factor(scheme, levels = c("Gene name class",
                                         "Positional class (GENCODE v36)"))]
# Classes are ordered by scheme so that legend blocks match the facets.
CLS <- c(sort(unique(cl[scheme == "Gene name class", class])),
         sort(unique(cl[scheme != "Gene name class", class])))
# Okabe-Ito without yellow, plus five further hues, so that rare classes stay
# visible against the white panel.
CLS_PAL <- c(unname(OI[c("orange", "skyblue", "green", "blue", "vermillion",
                         "purple", "grey")]), "#000000", "#8C6D31",
             "#7570B3", "#66A61E", "#A0165C")
CLS_COL <- setNames(rep(CLS_PAL, length.out = length(CLS)), CLS)
cl[, class := factor(class, levels = CLS)]
p6d <- ggplot(cl, aes(module, n, fill = class)) +
  geom_col(position = "fill", width = 0.68) +
  coord_flip() +
  facet_wrap(~ scheme, ncol = 1) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = CLS_COL, name = NULL, drop = FALSE) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  labs(title = "Gene classes", x = NULL, y = "Proportion of module genes") +
  theme(legend.position = "bottom", legend.text = element_text(size = 6),
        legend.key.size = unit(0.26, "cm"))

# Supplementary Figure 9: a and b are the partner GO terms of the other
# prognostic lncRNA modules (protective first), c the gene-class composition.
s9_mods <- if (is.null(gs)) character(0) else
  svl[fdr_full < FDR_ALPHA & module %in% unique(gs$module) & module != lnc_mod][order(HR_full), module]
s9_go <- lapply(s9_mods[seq_len(min(2, length(s9_mods)))], function(m) go_signed_panel(gs, m))
fig_s9 <- if (length(s9_go) == 2) {
  (ROW(TG(s9_go[[1]], "a") | TG(s9_go[[2]], "b")) / ROW(TG(p6d, "c"))) +
    plot_layout(heights = c(1, 1.1))
} else p6d
save_fig(fig_s9, "SupplementaryFigureS9_lncRNA_module_annotation", W, 9.2)

# -- d: single-cell localisation of the five prognostic modules (stage 36) -----
# Per cohort: the module's median log2 fold change in malignant cells minus the
# mean of 2,000 abundance-matched random gene sets. The per-cohort shift is
# plotted instead of the Stouffer z, which treats genes as exchangeable and
# overstates precision. Colour is hazard direction and shape is cohort.
sc <- Rx("36_sc_module_localisation_per_dataset.tsv")
p6d_sc <- NULL
if (!is.null(sc)) {
  prog <- data.table(network = c("protein-coding", "lncRNA", "protein-coding", "lncRNA", "lncRNA"),
                     module  = c("green", "blue", "purple", "greenyellow", "turquoise"),
                     dir     = c("protective", "protective", "adverse", "adverse", "adverse"))
  sc <- merge(sc, prog, by = c("network", "module"))
  sc[, shift := median_log2FC - null_mean]
  sc[, lab := factor(paste(ifelse(network == "lncRNA", "lncRNA", "PC"), module),
                     levels = rev(paste(ifelse(prog$network == "lncRNA", "lncRNA", "PC"), prog$module)))]
  n_by <- sc[, .(k = uniqueN(dataset)), by = contrast]
  sc[, facet := factor(ifelse(contrast == "malignant_vs_normal_epithelial",
                              sprintf("vs normal epithelium (%d cohorts)", n_by[contrast == "malignant_vs_normal_epithelial", k]),
                              sprintf("vs tumour microenvironment (%d cohorts)", n_by[contrast == "malignant_vs_pooled_tme", k])))]
  sc[, facet := factor(facet, levels = sort(levels(facet)))]
  sc[, cohort := sub("^KIRC_", "", dataset)]
  sc[, sig := p_empirical < 0.05]
  mn <- sc[, .(shift = mean(shift)), by = .(facet, lab, dir)]
  SC_SHAPE <- setNames(c(21, 22, 24, 23), c("GSE159115", "Li2022", "GSE222703", "GSE207493"))
  p6d_sc <- ggplot(sc, aes(shift, lab)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_point(data = mn, aes(colour = dir), shape = "|", size = 5, show.legend = FALSE) +
    geom_point(aes(shape = cohort, colour = dir, fill = ifelse(sig, dir, "none")),
               size = 1.9, stroke = 0.6, position = position_jitter(height = 0.12, width = 0, seed = SEED)) +
    facet_wrap(~ facet, ncol = 1) +
    scale_colour_manual(values = DIR_COL, guide = "none") +
    scale_fill_manual(values = c(DIR_COL, none = "white"), guide = "none") +
    scale_shape_manual(values = SC_SHAPE, name = NULL) +
    labs(title = "Single-cell localisation",
         x = "Shift above matched null (log2 fold change)", y = NULL) +
    theme(legend.position = "bottom") +
    guides(shape = guide_legend(nrow = 2, override.aes = list(fill = "grey40", colour = "grey20")))
} else MISSING <- unique(c(MISSING, "36_sc_module_localisation_per_dataset.tsv"))

# Lettered in order of citation: p6c is b and p6b is c.
fig6 <- (ROW(TG(p6a, "a") | TG(p6c, "b")) /
         ROW((if (is.null(p6b)) plot_spacer() else TG(p6b, "c")) |
             (if (is.null(p6d_sc)) plot_spacer() else TG(p6d_sc, "d")))) +
  plot_layout(heights = c(1, 1.25))
save_fig(fig6, "Figure6_biology", W, 8.4)
msg("  Figure 6 written")

# ==== Figure 3: the library-quality axis across TCGA (objects fig7, p7*) =====
# Drawn only when the stage 39 to 41 tables exist. Tumour rows of eligible
# projects are the inferential set. Reference lines are the locked rule
# thresholds from 00_config.R.
ax7 <- Rx("39_pancancer_axis.tsv"); gr7 <- Rx("39_pancancer_class_gradient.tsv")
nw7 <- Rx("40_pancancer_network.tsv"); sv7 <- Rx("41_pancancer_metric_survival.tsv")
mt7 <- Rx("41_pancancer_metric_meta.tsv")
if (!is.null(ax7) && !is.null(nw7)) {
  banner("Figure 3 | pan-cancer")
  lab7 <- function(x) sub("^TCGA-", "", x)
  # -- a: lncRNA against protein-coding PC1, rho with the metric ---------------
  a <- ax7[group == "tumour" & eligible == TRUE]
  a[, proj := factor(lab7(project), levels = lab7(project[order(abs(rho_lnc_pc1_noFeature))]))]
  al <- melt(a, id.vars = c("proj", "rule_A_pass"),
             measure.vars = c("rho_lnc_pc1_noFeature", "rho_pc_pc1_noFeature"),
             variable.name = "matrix", value.name = "rho")
  al[, matrix := factor(ifelse(matrix == "rho_lnc_pc1_noFeature", "lncRNA PC1", "protein-coding PC1"),
                        levels = c("lncRNA PC1", "protein-coding PC1"))]
  nrm <- ax7[group == "normal"][, proj := factor(lab7(project), levels = levels(a$proj))][!is.na(proj)]
  p7a <- ggplot(al, aes(abs(rho), proj)) +
    geom_vline(xintercept = PAN_AXIS_RHO_MIN, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_line(aes(group = proj), colour = "grey75", linewidth = 0.4) +
    geom_point(aes(colour = matrix, shape = matrix), size = 1.8) +
    geom_point(data = nrm, aes(abs(rho_lnc_pc1_noFeature), proj), shape = 4, size = 1.5,
               colour = BIO_COL[["lncRNA"]], inherit.aes = FALSE) +
    scale_colour_manual(values = c(`lncRNA PC1` = BIO_COL[["lncRNA"]], `protein-coding PC1` = BIO_COL[["protein-coding"]]), name = NULL) +
    scale_shape_manual(values = c(`lncRNA PC1` = 16, `protein-coding PC1` = 1), name = NULL) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(title = "Leading components",
         x = "Absolute rho with the non-feature fraction", y = NULL) +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 6))
  # -- b: positional-class gradient --------------------------------------------
  p7b <- NULL
  if (!is.null(gr7)) {
    g <- gr7[group == "tumour" & project %in% a$project]
    g[, class := factor(gsub("_", " ", class), levels = gsub("_", " ", c("intronic", "exonic_sense", "antisense", "sense_overlapping", "divergent", "intergenic")))]
    gm <- g[, .(med = median(median_abs_rho)), by = class]
    p7b <- ggplot(g, aes(class, median_abs_rho)) +
      geom_line(aes(group = project), colour = "grey80", linewidth = 0.3) +
      geom_point(size = 0.8, colour = "grey55") +
      geom_point(data = gm, aes(class, med), colour = BIO_COL[["lncRNA"]], size = 2.4, shape = 18) +
      labs(title = "Positional classes",
           x = NULL, y = "Median per-gene absolute rho") +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
  }
  # -- c: network reorganisation -------------------------------------------------
  n <- nw7[eligible == TRUE]
  n[, proj := factor(lab7(project), levels = lab7(project[order(largest_frac_observed)]))]
  nl <- melt(n, id.vars = c("proj", "rule_C_pass", "ari_all_genes"),
             measure.vars = c("largest_frac_observed", "largest_frac_adjusted"),
             variable.name = "network", value.name = "frac")
  nl[, network := factor(ifelse(network == "largest_frac_observed", "observed", "adjusted"),
                         levels = c("observed", "adjusted"))]
  p7c <- ggplot(nl, aes(frac, proj)) +
    geom_vline(xintercept = PAN_NET_LARGEST_MIN, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_line(aes(group = proj), colour = "grey75", linewidth = 0.4) +
    geom_point(aes(colour = network), size = 1.8) +
    geom_text(data = n, aes(x = 1, y = proj, label = sprintf("%.2f", ari_all_genes)),
              size = 1.9, hjust = 1, colour = "grey30", inherit.aes = FALSE) +
    scale_colour_manual(values = c(observed = BIO_COL[["lncRNA"]], `adjusted` = unname(OI["blue"])), name = NULL) +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(title = "Largest module",
         x = "Largest module, share of genes", y = NULL) +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 6))
  # -- d: the metric's hazard ratio -------------------------------------------------
  p7d <- NULL
  if (!is.null(sv7) && !is.null(mt7)) {
    d <- sv7[model == "primary (plate stratum)" & in_meta == TRUE]
    d[, proj := factor(lab7(project), levels = c("Pooled", lab7(project[order(hr)])))]
    pm <- mt7[model == "primary (plate stratum)"]
    pool <- data.table(proj = factor("Pooled", levels = levels(d$proj)), hr = pm$pooled_hr, lo = pm$lo, hi = pm$hi)
    p7d <- ggplot(d, aes(hr, proj)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.3) +
      ci_h(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.35, colour = "grey35") +
      geom_point(size = 1.5, colour = "grey20") +
      geom_point(data = pool, shape = 18, size = 3.4, colour = BIO_COL[["lncRNA"]]) +
      ci_h(data = pool, aes(xmin = lo, xmax = hi), width = 0.3, linewidth = 0.6, colour = BIO_COL[["lncRNA"]]) +
      scale_x_log10(breaks = c(0.5, 0.7, 1, 1.4, 2)) +
      labs(title = "Metric and overall survival",
           x = "Hazard ratio per SD (log scale)", y = NULL) +
      theme(axis.text.y = element_text(size = 6))
  }
  # -- e: dose-response (post hoc, stage 43) -------------------------------------
  dr <- Rx("43_pancancer_dose_response_projects.tsv"); dt <- Rx("43_pancancer_dose_response.tsv")
  p7e <- NULL
  if (!is.null(dr) && !is.null(dt)) {
    dr[, lab := lab7(project)]
    p7e <- ggplot(dr, aes(mad_log_metric, largest_frac_observed)) +
      geom_hline(yintercept = PAN_NET_LARGEST_MIN, linetype = 2, colour = "grey55", linewidth = 0.3) +
      geom_point(aes(colour = abs_rho_lnc, shape = rule_C_pass), size = 2.2) +
      ggrepel::geom_text_repel(aes(label = lab), size = 1.9, colour = "grey25", seed = SEED,
                               max.overlaps = Inf, min.segment.length = 0.1, segment.size = 0.2) +
      scale_colour_gradient(low = "grey75", high = BIO_COL[["lncRNA"]], limits = c(0, 1), breaks = c(0, 0.5, 1), name = "Absolute rho, lncRNA PC1",
                            guide = guide_colourbar(title.position = "top")) +
      scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), labels = c(`TRUE` = "Rule C met", `FALSE` = "Rule C not met"), name = NULL) +
      scale_y_continuous(labels = scales::percent) +
      labs(title = "Metric spread (exploratory)",
           x = "MAD of log non-feature fraction", y = "Largest observed module") +
      theme(legend.position = "bottom", legend.box = "vertical",
            legend.key.width = unit(0.8, "cm"))
  }
  # Lettered in order of citation: p7e (dose-response) is d and p7d (survival
  # forest) is e.
  fig7 <- (ROW(TG(p7a, "a") | (if (is.null(p7b)) plot_spacer() else TG(p7b, "b"))) /
           ROW(TG(p7c, "c") | (if (is.null(p7e)) plot_spacer() else TG(p7e, "d"))) /
           (if (is.null(p7d)) plot_spacer() else ROW(TG(p7d, "e")))) +
    plot_layout(heights = c(1, 1.05, 0.7))
  save_fig(fig7, "Figure3_pancancer", W, 12.4)
  msg("  Figure 3 (pan-cancer) written")
}

# ---- report skipped panels ----
if (length(MISSING))
  msg("NOTE: source files absent, affected panels skipped: ",
      paste(MISSING, collapse = ", "))
msg("Main figures written to ", FIG_DIR, " as SVG and ", FIG_DPI, " dpi PNG")
write_session_info("13_figures")
banner("13 | done")
