# 03_survival.R: Cox models of module scores against overall survival.
#
# Run after 02 and 07. Scores come from discovery_scores() (unit variance, so
# HRs are per SD) and are checked against the WGCNA eigengenes. Every non-grey
# module is fitted under four specifications, each on its own complete cases:
#   uni   score alone
#   clin  score + age, sex, T, N, M1, ordinal grade
#   full  clin + ESTIMATE stromal and immune + STAR metrics (principal)
#   adj   score + age, sex, stage group, grade group
# BH FDR within network, and jointly across networks for clin and full. Also
# writes proportional-hazards tests, the N/M coding sensitivity, the lncRNA
# axis models, the clinical reference models and cache/survival.rds.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(survminer); library(ggplot2)
})
banner("03 | Survival analysis")

ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
cohort <- as.data.table(ds$cohort_full)

need_cols <- c("M1_as_coded", "N_pos_imputed", "M1_imputed", "libsize",
               "pct_noFeature", "pct_multimapping", "T_stage", "N_pos", "M1",
               "grade_num", "stage_group", "grade_group")
.miss <- setdiff(need_cols, names(cohort))
if (length(.miss))
  stop("cohort_full lacks ", paste(.miss, collapse = ", "),
       ": re-run 01_build_data.R before 03")

# ESTIMATE scores and the lncRNA axis are cached by 07.
need_cache <- function(f) {
  p <- file.path(CACHE_DIR, f)
  if (!file.exists(p)) stop(f, " not found in cache: run 07_purity_qc.R before 03")
  as.data.table(readRDS(p))
}
EST  <- need_cache("estimate_scores.rds")
AXIS <- need_cache("lnc_global_axis.rds")

NET_TAG      <- c(mRNA = "mrna",    lncRNA = "lnc")      # tag -> element of nets / SC
SCORE_PREFIX <- c(mRNA = "mRNA_ME", lncRNA = "lnc_ME")   # discovery_loadings() names
NET_LABEL    <- c(mRNA = "Protein-coding", lncRNA = "lncRNA")

# ---- module scores against the WGCNA eigengenes ----
banner("Module scores: discovery_scores() vs WGCNA eigengenes")
SC <- discovery_scores(nets)

score_check <- rbindlist(lapply(names(NET_TAG), function(tag) {
  S  <- SC[[NET_TAG[[tag]]]]
  ME <- nets[[NET_TAG[[tag]]]]$MEs
  stopifnot(all(rownames(S) %in% rownames(ME)))
  ME <- ME[rownames(S), , drop = FALSE]
  mods <- sub(paste0("^", SCORE_PREFIX[[tag]]), "", colnames(S))
  data.table(biotype = tag, module = mods,
             r = vapply(seq_along(mods), function(j)
               cor(S[, j], ME[[paste0("ME", mods[j])]]), numeric(1)))
}))
print(score_check[, .(biotype, module, r = round(r, 5))])
stopifnot(all(abs(score_check$r) > 0.99))
if (any(score_check$r < 0))
  msg("NOTE: ", sum(score_check$r < 0), " module score(s) are sign-flipped ",
      "relative to the WGCNA eigengene: ",
      paste(score_check[r < 0, paste(biotype, module)], collapse = ", "))
save_tsv(score_check, "03_eigengene_definition_check.tsv")

# ---- model frames ----
# One frame per network. STAR and ESTIMATE covariates are standardised over all
# network samples before complete-case restriction, so units match across
# specifications.
FRAMES <- lapply(NET_TAG, function(k) {
  S  <- SC[[k]]
  cl <- cohort[match(rownames(S), sample_barcode)]
  stopifnot(all(cl$sample_barcode == rownames(S)))
  e  <- EST[match(cl$sample_barcode, EST$sample_barcode)]
  base <- data.table(
    sample_barcode = cl$sample_barcode,
    os_time   = cl$os_time,
    os_event  = cl$os_event,
    age       = cl$age,
    sex       = factor(cl$sex, levels = c("female", "male")),
    stage     = cl$stage_group,
    grade     = cl$grade_group,
    T_stage   = as.numeric(cl$T_stage),
    N_pos     = as.numeric(cl$N_pos),
    M1        = as.numeric(cl$M1),
    M1_coded  = as.numeric(cl$M1_as_coded),
    grade_o   = as.numeric(cl$grade_num),
    N_unknown = as.logical(cl$N_pos_imputed),
    M_unknown = as.logical(cl$M1_imputed),
    noFeat    = as.numeric(scale(cl$pct_noFeature)),
    multimap  = as.numeric(scale(cl$pct_multimapping)),
    depth     = as.numeric(scale(log10(cl$libsize))),
    stromal   = as.numeric(scale(e$StromalScore)),
    immune    = as.numeric(scale(e$ImmuneScore)))
  # Three-level N and M factors for the coding sensitivity: positives are
  # N1/M1, unassessed negatives are NX/MX.
  base[, N_cat := factor(fifelse(N_pos == 1, "N1", fifelse(N_unknown, "NX", "N0")),
                         levels = c("N0", "N1", "NX"))]
  base[, M_cat := factor(fifelse(M1 == 1, "M1", fifelse(M_unknown, "MX", "M0")),
                         levels = c("M0", "M1", "MX"))]
  list(scores = S, base = base)
})

module_frame <- function(tag, m) {
  fr <- FRAMES[[tag]]
  d  <- copy(fr$base)
  d[, ME := as.numeric(fr$scores[, paste0(SCORE_PREFIX[[tag]], m)])]
  d
}

SPEC_TERMS <- list(
  uni  = character(0),
  clin = c("age", "sex", "T_stage", "N_pos", "M1", "grade_o"),
  full = c("age", "sex", "T_stage", "N_pos", "M1", "grade_o",
           "stromal", "immune", "noFeat", "multimap", "depth"),
  adj  = c("age", "sex", "stage", "grade"))

# Fit one Cox specification on the complete cases of its own covariates.
# cox.zph is reported for the score term and globally.
fit_spec <- function(d, terms, me = "ME") {
  vars <- c("os_time", "os_event", me, terms)
  dc <- d[complete.cases(d[, ..vars])]
  f  <- coxph(as.formula(paste("Surv(os_time, os_event) ~",
                               paste(c(me, terms), collapse = " + "))), data = dc)
  s  <- summary(f)
  tb <- tryCatch(cox.zph(f)$table, error = function(e) NULL)
  list(fit = f, data = dc, n = s$n, events = s$nevent,
       HR = unname(s$conf.int[me, "exp(coef)"]),
       lo = unname(s$conf.int[me, "lower .95"]),
       hi = unname(s$conf.int[me, "upper .95"]),
       p  = unname(s$coefficients[me, "Pr(>|z|)"]),
       C  = unname(s$concordance["C"]),
       ph_ME     = if (is.null(tb) || !me %in% rownames(tb)) NA_real_ else tb[me, "p"],
       ph_global = if (is.null(tb) || !"GLOBAL" %in% rownames(tb)) NA_real_
                   else tb["GLOBAL", "p"])
}

# ---- per-module Cox models ----
analyse_biotype <- function(tag) {
  banner(paste("Survival:", tag))
  fr <- FRAMES[[tag]]
  modules <- setdiff(sub(paste0("^", SCORE_PREFIX[[tag]]), "", colnames(fr$scores)),
                     "grey")
  msg(tag, ": testing ", length(modules), " modules (grey excluded) on n = ",
      nrow(fr$base), ", events = ", sum(fr$base$os_event, na.rm = TRUE))

  mod_size <- nets[[NET_TAG[[tag]]]]$gene_tbl[, .N, by = module]
  setkey(mod_size, module)

  res <- rbindlist(lapply(modules, function(m) {
    d    <- module_frame(tag, m)
    fits <- lapply(SPEC_TERMS, function(tt) fit_spec(d, tt))
    u <- fits$uni; k <- fits$clin; q <- fits$full; a <- fits$adj

    # median-split Kaplan-Meier on the unadjusted sample set
    du <- copy(u$data)
    du[, grp := factor(ifelse(ME > median(ME), "High", "Low"),
                       levels = c("Low", "High"))]
    lr <- survdiff(Surv(os_time, os_event) ~ grp, data = du)
    p_lr <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)

    data.table(
      module     = m,
      n_genes    = mod_size[m, N],
      n          = u$n,
      events     = u$events,
      HR_uni     = u$HR, HR_uni_lo = u$lo, HR_uni_hi = u$hi, p_uni = u$p, C_uni = u$C,
      ph_p_ME_uni     = u$ph_ME,
      ph_p_global_uni = u$ph_global,
      n_adj      = a$n,
      events_adj = a$events,
      HR_adj     = a$HR, HR_adj_lo = a$lo, HR_adj_hi = a$hi, p_adj = a$p, C_adj = a$C,
      ph_p_ME_adj     = a$ph_ME,
      ph_p_global_adj = a$ph_global,
      n_clin      = k$n,
      events_clin = k$events,
      HR_clin     = k$HR, HR_clin_lo = k$lo, HR_clin_hi = k$hi, p_clin = k$p, C_clin = k$C,
      ph_p_ME_clin     = k$ph_ME,
      ph_p_global_clin = k$ph_global,
      n_full      = q$n,
      events_full = q$events,
      HR_full     = q$HR, HR_full_lo = q$lo, HR_full_hi = q$hi, p_full = q$p, C_full = q$C,
      ph_p_ME_full     = q$ph_ME,
      ph_p_global_full = q$ph_global,
      p_logrank   = p_lr
    )
  }))

  # FDR within network across all non-grey modules, per specification
  res[, fdr_uni     := p.adjust(p_uni,     "BH")]
  res[, fdr_adj     := p.adjust(p_adj,     "BH")]
  res[, fdr_clin    := p.adjust(p_clin,    "BH")]
  res[, fdr_full    := p.adjust(p_full,    "BH")]
  res[, fdr_logrank := p.adjust(p_logrank, "BH")]
  setorder(res, p_full)

  msg(tag, ": ", sum(res$fdr_uni < FDR_ALPHA), " FDR-significant unadjusted; ",
      sum(res$fdr_clin < FDR_ALPHA), " under the clinical comparator; ",
      sum(res$fdr_full < FDR_ALPHA), " under the principal specification ",
      "(clinical + composition + library quality)")
  print(res[, .(module, n_genes,
                HR_clin = round(HR_clin, 3), fdr_clin = signif(fdr_clin, 3),
                HR_full = round(HR_full, 3), fdr_full = signif(fdr_full, 3),
                n_full, events_full, ph_ME = signif(ph_p_ME_full, 2))][1:min(15, .N)])

  # ---- forest plot, principal specification ----
  fp <- copy(res)
  fp[, sig := fdr_full < FDR_ALPHA]
  fp[, module := factor(module, levels = rev(module[order(HR_full)]))]
  g <- ggplot(fp, aes(x = HR_full, y = module, colour = sig)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_errorbarh(aes(xmin = HR_full_lo, xmax = HR_full_hi), height = 0.25) +
    geom_point(size = 1.8) +
    scale_x_log10() +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey55"),
                        name = paste0("FDR < ", FDR_ALPHA)) +
    labs(title = paste0(tag, ": module score vs overall survival"),
         subtitle = paste0("Cox model adjusted for age, sex, T, N, M1, grade, ",
                           "ESTIMATE and STAR metrics; HR per 1 SD of score"),
         x = "Adjusted hazard ratio (95% CI)", y = NULL) +
    theme_bw(base_size = 9)
  save_fig(g, paste0("03_", tag, "_forest_principal"),
           6.5, max(3, 0.22 * nrow(fp)))

  # ---- KM curves for modules significant under the principal specification ----
  sig_mods <- res[fdr_full < FDR_ALPHA, module]
  for (m in sig_mods) {
    d <- module_frame(tag, m)[!is.na(os_time) & !is.na(os_event) & !is.na(ME)]
    d[, Eigengene := factor(ifelse(ME > median(ME), "High", "Low"),
                            levels = c("Low", "High"))]
    fit <- survfit(Surv(os_time, os_event) ~ Eigengene, data = d)
    pl <- ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE,
                     palette = c("#2166AC", "#B2182B"),
                     conf.int = TRUE, xlab = "Days from diagnosis",
                     ylab = "Overall survival", legend.title = "ME",
                     title = paste0(tag, " module ", m,
                                    " (median split, n = ", nrow(d), ")"))
    # A ggsurvplot with a risk table is composite, so it cannot go to ggsave().
    save_fig_base(paste0("03_", tag, "_KM_", m), 6, 6.5, print(pl))
  }
  msg(tag, ": wrote ", length(sig_mods), " Kaplan-Meier figures")
  res[]
}

surv <- list(mrna = analyse_biotype("mRNA"), lnc = analyse_biotype("lncRNA"))

# ---- joint FDR across both networks ----
joint <- rbindlist(lapply(surv, function(r) r[, .(module, p_full, p_clin)]),
                   idcol = "net")
joint[, fdr_full_joint := p.adjust(p_full, "BH")]
joint[, fdr_clin_joint := p.adjust(p_clin, "BH")]
for (k in names(surv)) {
  surv[[k]][joint[net == k], on = "module",
            `:=`(fdr_full_joint = i.fdr_full_joint, fdr_clin_joint = i.fdr_clin_joint)]
}
save_tsv(surv$mrna, "03_mRNA_module_survival.tsv")
save_tsv(surv$lnc,  "03_lncRNA_module_survival.tsv")
msg("Joint FDR over 26 modules: ", sum(joint$fdr_full_joint < FDR_ALPHA),
    " significant under full, ", sum(joint$fdr_clin_joint < FDR_ALPHA),
    " under clin")

# Modules significant under the principal specification (within-network FDR)
SIG <- rbindlist(lapply(names(NET_TAG), function(tag) {
  r <- surv[[NET_TAG[[tag]]]][fdr_full < FDR_ALPHA]
  data.table(biotype = tag, module = r$module,
             label = paste(NET_LABEL[[tag]], r$module))
}))
msg("Modules significant under full (within-network FDR < ", FDR_ALPHA, "): ",
    if (nrow(SIG)) paste(SIG$label, collapse = ", ") else "none")

# ---- all modules under the principal specification ----
all_principal <- rbindlist(lapply(names(NET_TAG), function(tag) {
  r <- surv[[NET_TAG[[tag]]]]
  data.table(biotype = tag, module = r$module, n_genes = r$n_genes,
             n = r$n_full, events = r$events_full,
             HR = r$HR_full, lo = r$HR_full_lo, hi = r$HR_full_hi, p = r$p_full,
             fdr_within = r$fdr_full, fdr_joint = r$fdr_full_joint,
             ph_p_ME = r$ph_p_ME_full, ph_p_global = r$ph_p_global_full)
}))
setorder(all_principal, p)
save_tsv(all_principal, "03_all_modules_principal.tsv")

# ---- proportional hazards ----
banner("Proportional hazards (cox.zph) for every module and specification")
ph_all <- rbindlist(lapply(names(NET_TAG), function(tag) {
  r <- surv[[NET_TAG[[tag]]]]
  data.table(biotype = tag, module = r$module,
             n_uni = r$n, events_uni = r$events,
             ph_p_ME_uni = r$ph_p_ME_uni, ph_p_global_uni = r$ph_p_global_uni,
             n_clin = r$n_clin, events_clin = r$events_clin,
             ph_p_ME_clin = r$ph_p_ME_clin, ph_p_global_clin = r$ph_p_global_clin,
             n_full = r$n_full, events_full = r$events_full,
             ph_p_ME_full = r$ph_p_ME_full, ph_p_global_full = r$ph_p_global_full,
             n_adj = r$n_adj, events_adj = r$events_adj,
             ph_p_ME_adj = r$ph_p_ME_adj, ph_p_global_adj = r$ph_p_global_adj,
             fdr_full = r$fdr_full)
}))
save_tsv(ph_all, "03_proportional_hazards_all_modules.tsv")
print(ph_all[, .(biotype, module, ph_ME_uni = signif(ph_p_ME_uni, 2),
                 ph_ME_clin = signif(ph_p_ME_clin, 2),
                 ph_ME_full = signif(ph_p_ME_full, 2),
                 ph_global_full = signif(ph_p_global_full, 2))])

# Significant modules only. The score-alone test uses the principal
# specification's complete-case set. 13_figures.R reads the 22_* tables.
ph_principal <- rbindlist(lapply(seq_len(nrow(SIG)), function(i) {
  d  <- module_frame(SIG$biotype[i], SIG$module[i])
  ff <- fit_spec(d, SPEC_TERMS$full)
  fu <- fit_spec(ff$data, SPEC_TERMS$uni)
  data.table(module = SIG$label[i], n = ff$n, events = ff$events,
             HR = round(ff$HR, 3),
             ph_ME_principal     = signif(ff$ph_ME, 2),
             ph_global_principal = signif(ff$ph_global, 2),
             ph_ME_unadjusted    = signif(fu$ph_ME, 2))
}))
print(ph_principal, row.names = FALSE)
save_tsv(ph_principal, "22_proportional_hazards_principal.tsv")

# ---- nodal and metastasis coding sensitivity (principal covariate set) ----
banner("Nodal / metastasis coding sensitivity")
msg("N status imputed (NX or missing): ", sum(cohort$N_pos_imputed),
    "; M status imputed: ", sum(cohort$M1_imputed),
    if ("M1_stage_reconciled" %in% names(cohort))
      paste0("; M1 reconciled from stage IV: ",
             sum(cohort$M1_stage_reconciled == 1, na.rm = TRUE)) else "",
    "; of ", nrow(cohort))

swap <- function(x, from, to) { x[x == from] <- to; x }
NM_SPECS <- list(
  A = list(terms = SPEC_TERMS$full, assessed = FALSE,
           label = "primary: NX -> N0, MX -> M0, M1 reconciled from stage IV"),
  B = list(terms = swap(swap(SPEC_TERMS$full, "N_pos", "N_cat"), "M1", "M_cat"),
           assessed = FALSE, label = "NX and MX as own categories"),
  C = list(terms = SPEC_TERMS$full, assessed = TRUE,
           label = "assessed N and M only"),
  D = list(terms = swap(SPEC_TERMS$full, "M1", "M1_coded"), assessed = FALSE,
           label = "M1 as coded in pM (no stage-IV reconciliation)"))

for (sp in names(NM_SPECS)) msg("Specification ", sp, ": ", NM_SPECS[[sp]]$label)

nxmx <- rbindlist(lapply(seq_len(nrow(SIG)), function(i) {
  rbindlist(lapply(names(NM_SPECS), function(sp) {
    s <- NM_SPECS[[sp]]
    d <- module_frame(SIG$biotype[i], SIG$module[i])
    if (s$assessed) d <- d[!N_unknown & !M_unknown]
    f <- fit_spec(d, s$terms)
    data.table(module = SIG$label[i], spec = sp, n = f$n, events = f$events,
               HR = round(f$HR, 3), lo = round(f$lo, 3), hi = round(f$hi, 3),
               p = signif(f$p, 3))
  }))
}))
print(nxmx, row.names = FALSE)
save_tsv(nxmx, "22_nodal_metastasis_sensitivity.tsv")
save_tsv(nxmx, "03_nodal_metastasis_sensitivity.tsv")

# ---- leading lncRNA axis (PC1 of the observed lncRNA matrix, from 07) ----
banner("Leading lncRNA expression axis")
gl <- cohort[match(AXIS$sample_barcode, sample_barcode)]
stopifnot(all(gl$sample_barcode == AXIS$sample_barcode))
gd <- data.table(os_time = gl$os_time, os_event = gl$os_event,
                 axis = as.numeric(scale(AXIS$lnc_axis)),
                 age = gl$age, sex = factor(gl$sex, levels = c("female", "male")),
                 stage = gl$stage_group, grade = gl$grade_group)
axis_fits <- list(
  `unadjusted`                          = fit_spec(gd, character(0), me = "axis"),
  `adjusted for age, sex, stage, grade` = fit_spec(gd, SPEC_TERMS$adj, me = "axis"))
axis_tbl <- rbindlist(lapply(names(axis_fits), function(nm) {
  f <- axis_fits[[nm]]
  data.table(model = nm, n = f$n, events = f$events, HR = f$HR, lo = f$lo, hi = f$hi,
             p = f$p, C = f$C, ph_p = f$ph_ME)
}))
save_tsv(axis_tbl, "03_lncRNA_global_axis.tsv")
print(axis_tbl, row.names = FALSE)

saveRDS(list(mrna = surv$mrna, lnc = surv$lnc, axis = axis_tbl),
        file.path(CACHE_DIR, "survival.rds"))

# ---- clinical reference models on cohort_adj ----
# Apparent concordance of the clinical and augmented comparators, on complete
# cases of clinical_design().
banner("Clinical reference models (cohort_adj)")
cl_adj <- as.data.table(ds$cohort_adj)
e_adj  <- EST[match(cl_adj$sample_barcode, EST$sample_barcode)]
qc_adj <- cl_adj[, .(pct_noFeature, pct_multimapping, assigned_reads = libsize)]

ref <- lapply(c(clinical = "clinical", augmented = "augmented"), function(set) {
  X  <- clinical_design(cl_adj, set = set, est = e_adj, qc = qc_adj)
  ok <- complete.cases(X) & !is.na(cl_adj$os_time) & !is.na(cl_adj$os_event)
  dd <- cbind(data.frame(os_time = cl_adj$os_time[ok], os_event = cl_adj$os_event[ok]),
              as.data.frame(X[ok, , drop = FALSE]))
  f  <- coxph(as.formula(paste("Surv(os_time, os_event) ~",
                               paste(colnames(X), collapse = " + "))), data = dd)
  s  <- summary(f)
  cc <- cindex_ci(Surv(dd$os_time, dd$os_event), predict(f, type = "lp"))
  list(
    terms = data.table(set = set, term = rownames(s$conf.int),
                       HR = unname(s$conf.int[, "exp(coef)"]),
                       lo = unname(s$conf.int[, "lower .95"]),
                       hi = unname(s$conf.int[, "upper .95"]),
                       p  = unname(s$coefficients[, "Pr(>|z|)"]),
                       n = s$n, events = s$nevent),
    cindex = data.table(set = set, n = s$n, events = s$nevent, n_terms = ncol(X),
                        C = cc[["C"]], se = cc[["se"]], lo = cc[["lo"]], hi = cc[["hi"]]))
})
clin_tbl <- rbindlist(lapply(ref, `[[`, "terms"))
clin_c   <- rbindlist(lapply(ref, `[[`, "cindex"))
save_tsv(clin_tbl, "03_clinical_reference_model.tsv")
save_tsv(clin_c,   "03_clinical_reference_cindex.tsv")
print(clin_tbl[, .(set, term, HR = round(HR, 3), lo = round(lo, 3), hi = round(hi, 3),
                   p = signif(p, 3), n, events)], row.names = FALSE)
print(clin_c[, .(set, n, events, n_terms, C = round(C, 3), lo = round(lo, 3),
                 hi = round(hi, 3))], row.names = FALSE)

write_session_info("03_survival")
banner("03 | done")
