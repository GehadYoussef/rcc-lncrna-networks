# 09_validate.R: external validation of the locked models in CPTAC-3 ccRCC.
#
# Module membership (02), eigengene loadings (discovery_loadings) and model
# coefficients are fixed in TCGA-KIRC and applied unchanged to CPTAC-3. Two
# comparators: clinical (primary: age, sex, T, N, M1, grade) and augmented
# (plus ESTIMATE and the STAR metrics). Scores are standardised within cohort,
# with a discovery-fixed standardisation arm for the clinical comparator.
# Technical residualisation is fitted within CPTAC-3.
# Inputs: caches from 01 to 08. Outputs: locked_model.rds (read by 10, 11, 12,
# 20 and 22) and the validation C-index, delta C, module replication and axis tables.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(glmnet)
  library(estimate); library(ggplot2)
})
banner("09 | External validation in CPTAC-3")
set.seed(SEED)

ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
surv <- readRDS(file.path(CACHE_DIR, "survival.rds"))
val  <- readRDS(file.path(CACHE_DIR, "validation_dataset.rds"))
vco  <- as.data.table(val$cohort)
if (is.null(val$qc) || !"M1_as_coded" %in% names(vco))
  stop("validation_dataset.rds was written by an older 08_validation_data.R: re-run it first")
msg("Validation cohort (clear cell): n = ", nrow(vco), ", events = ", sum(vco$os_event))

lpv <- function(X, b) as.numeric(X[, names(b), drop = FALSE] %*% b)

# ---- 1. CPTAC-3 library-quality metrics (from 08) ----
vqc <- as.data.table(val$qc)[match(vco$sample_barcode, sample_barcode)]
stopifnot(identical(vqc$sample_barcode, vco$sample_barcode),
          !anyNA(vqc$pct_noFeature), !anyNA(vqc$assigned_reads))

# ---- 2. CPTAC-3 ESTIMATE scores ----
vest_f <- file.path(CACHE_DIR, "valid_estimate.rds")
vest <- if (file.exists(vest_f)) as.data.table(readRDS(vest_f)) else NULL
if (is.null(vest) || !all(vco$sample_barcode %in% vest$sample_barcode)) {
  msg("Computing ESTIMATE scores for ", nrow(vco), " CPTAC-3 samples ...")
  pc  <- which(val$gene_ann$gene_type == "protein_coding")
  m   <- val$fpkm[pc, vco$sample_barcode, drop = FALSE]
  sym <- val$gene_ann$gene_name[pc]
  k   <- !is.na(sym) & nzchar(sym) & !duplicated(sym)
  m <- log2(m[k, , drop = FALSE] + 1); rownames(m) <- sym[k]
  inf  <- file.path(CACHE_DIR, "v_est_in.txt")
  filf <- file.path(CACHE_DIR, "v_est_f.gct")
  scf  <- file.path(CACHE_DIR, "v_est_s.gct")
  write.table(data.frame(GeneSymbol = rownames(m), m, check.names = FALSE),
              inf, sep = "\t", quote = FALSE, row.names = FALSE)
  filterCommonGenes(input.f = inf, output.f = filf, id = "GeneSymbol")
  estimateScore(filf, scf, platform = "illumina")
  g  <- fread(scf, skip = 2, header = TRUE)
  mm <- as.matrix(g[, -c(1, 2), with = FALSE]); rownames(mm) <- g[[1]]
  vest <- data.table(sample_barcode = colnames(m),
                     StromalScore = as.numeric(mm["StromalScore", ]),
                     ImmuneScore  = as.numeric(mm["ImmuneScore", ]))
  saveRDS(vest, vest_f)
}
vest <- vest[match(vco$sample_barcode, sample_barcode)]
stopifnot(!anyNA(vest$StromalScore))

# ---- 3. harmonised expression on the network genes ----
vfp    <- val$fpkm[, vco$sample_barcode, drop = FALSE]
gm     <- intersect(colnames(nets$mrna$expr), rownames(vfp))
gl     <- intersect(colnames(nets$lnc$expr),  rownames(vfp))
Vm_obs <- t(log2(vfp[gm, , drop = FALSE] + 1))
Vl_obs <- t(log2(vfp[gl, , drop = FALSE] + 1))
# Residualised within CPTAC-3: its non-feature fraction lies far outside the
# TCGA range (different library protocol), so discovery coefficients would
# extrapolate.
Vm <- remove_technical(Vm_obs, tech_covariates(vqc))
Vl <- remove_technical(Vl_obs, tech_covariates(vqc))
msg("Gene coverage -- mRNA ", ncol(Vm), "/", ncol(nets$mrna$expr),
    ", lncRNA ", ncol(Vl), "/", ncol(nets$lnc$expr))

# ---- 4. module scores by projection of the discovery loadings ----
loadings <- discovery_loadings(nets)
Lm <- loadings[grep("^mRNA_", names(loadings))]
Ll <- loadings[grep("^lnc_",  names(loadings))]

# Discovery scores on each network's samples. CPTAC-3 is scored under both
# standardisations and restricted to the evaluated set below.
St_m <- score_modules(nets$mrna$expr, Lm)
St_l <- score_modules(nets$lnc$expr,  Ll)
Sv_m <- score_modules(Vm, Lm); Sv_l <- score_modules(Vl, Ll)
Sd_m <- score_modules(Vm, Lm, gene_standardise = "discovery", score_standardise = "discovery")
Sd_l <- score_modules(Vl, Ll, gene_standardise = "discovery", score_standardise = "discovery")

# Gene coverage and variance explained by the projected score, per cohort.
coverage <- rbindlist(lapply(names(loadings), function(nm) {
  L <- loadings[[nm]]
  V <- if (grepl("^mRNA_", nm)) Vm else Vl
  g <- intersect(L$genes, colnames(V))
  ve_c <- NA_real_
  if (length(g) >= 10) {
    # as in score_modules(): zero-variance genes get scale 1
    Xg   <- V[, g, drop = FALSE]
    sdv  <- apply(Xg, 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
    Z    <- sweep(sweep(Xg, 2, colMeans(Xg), "-"), 2, sdv, "/")
    s    <- as.numeric(Z %*% L$rot[g])
    ve_c <- var(s) / sum(apply(Z, 2, var))
  }
  data.table(module = nm, n_genes_discovery = length(L$genes), n_genes_cptac = length(g),
             var_explained_discovery = round(L$var_explained, 4),
             var_explained_cptac = round(ve_c, 4))
}))
save_tsv(coverage, "09_module_gene_coverage.tsv")
print(coverage)

# ---- 5. analysis sets: complete clinical inputs and survival ----
cf   <- as.data.table(ds$cohort_full)
test <- as.data.table(readRDS(file.path(CACHE_DIR, "estimate_scores.rds")))
tqc  <- as.data.table(readRDS(file.path(CACHE_DIR, "star_qc.rds")))
clin_cols <- c("age", "sex", "T_stage", "N_pos", "M1", "grade_num", "os_time", "os_event")

net_samp <- intersect(rownames(St_m), rownames(St_l))
flow <- list(discovery_in_both_networks = length(net_samp))
tcl0 <- cf[match(net_samp, sample_barcode)]
ok_t <- complete.cases(tcl0[, ..clin_cols])
flow$discovery_complete_clinical_set <- sum(ok_t)
# Requiring the augmented inputs keeps one locked set for both comparators.
ok_t <- ok_t & net_samp %in% test$sample_barcode & net_samp %in% tqc$sample_barcode
flow$discovery_complete_augmented_set <- sum(ok_t)
tsamp <- net_samp[ok_t]
tcl   <- cf[match(tsamp, sample_barcode)]
test  <- test[match(tsamp, sample_barcode)]
tqc   <- tqc[match(tsamp, sample_barcode)]
stopifnot(!anyNA(test$StromalScore), !anyNA(tqc$pct_noFeature))
flow$discovery_events <- sum(tcl$os_event)

flow$cptac_clear_cell <- nrow(vco)
flow$cptac_T_missing <- sum(is.na(vco$T_stage))
flow$cptac_grade_missing <- sum(is.na(vco$grade_num))
ok_v  <- complete.cases(vco[, ..clin_cols])
vsamp <- vco$sample_barcode[ok_v]
vcl   <- vco[match(vsamp, sample_barcode)]
vqm   <- vqc[match(vsamp, sample_barcode)]
vem   <- vest[match(vsamp, sample_barcode)]
flow$cptac_complete_clinical_set <- length(vsamp)
flow$cptac_events <- sum(vcl$os_event)
flow$cptac_M1_reconciled_from_stage_iv <- sum(vcl$M1_stage_reconciled)
save_tsv(data.frame(step = names(flow), n = unlist(flow), row.names = NULL),
         "09_analysis_set_flow.tsv")
msg("Locked set: TCGA-KIRC n = ", length(tsamp), " (events ", sum(tcl$os_event),
    ");  CPTAC-3 evaluated n = ", length(vsamp), " (events ", sum(vcl$os_event),
    "); CPTAC-3 excluded for missing T: ", flow$cptac_T_missing)

yt <- Surv(tcl$os_time, tcl$os_event)
yv <- Surv(vcl$os_time, vcl$os_event)

# ---- module score matrices, identical columns in every arm ----
Et      <- cbind(St_m[tsamp, , drop = FALSE], St_l[tsamp, , drop = FALSE])
Ev      <- cbind(Sv_m[vsamp, , drop = FALSE], Sv_l[vsamp, , drop = FALSE])
Ev_disc <- cbind(Sd_m[vsamp, , drop = FALSE], Sd_l[vsamp, , drop = FALSE])
me_cols <- Reduce(intersect, list(colnames(Et), colnames(Ev), colnames(Ev_disc)))
if (length(me_cols) < ncol(Et))
  msg("Modules without adequate CPTAC-3 coverage dropped from the locked model: ",
      paste(setdiff(colnames(Et), me_cols), collapse = ", "))
Et <- Et[, me_cols, drop = FALSE]; Ev <- Ev[, me_cols, drop = FALSE]
Ev_disc <- Ev_disc[, me_cols, drop = FALSE]
stopifnot(ncol(Et) > 0, !anyNA(Et), !anyNA(Ev), !anyNA(Ev_disc))

# ---- design matrices ----
age_const <- c(mean = mean(tcl$age), sd = sd(tcl$age))
X <- list(
  clinical = list(
    Xt = clinical_design(tcl, "clinical"),
    Xv = clinical_design(vcl, "clinical"),
    Xv_disc = clinical_design(vcl, "clinical", age_center = age_const[["mean"]],
                              age_scale = age_const[["sd"]])),
  augmented = list(
    Xt = clinical_design(tcl, "augmented", est = test, qc = tqc),
    Xv = clinical_design(vcl, "augmented", est = vem, qc = vqm)))
for (cmp in names(X)) for (nm in names(X[[cmp]])) stopifnot(!anyNA(X[[cmp]][[nm]]))

# ---- 6. lock the models in TCGA-KIRC, one pair per comparator ----
models <- lapply(names(X), function(cmp) {
  Xt <- X[[cmp]]$Xt
  b_clin <- coef(coxph(yt ~ ., data = as.data.frame(Xt)))
  XEt <- cbind(Xt, Et)
  pfe <- c(rep(0, ncol(Xt)), rep(1, ncol(Et)))
  set.seed(SEED)
  fit <- cv.glmnet(XEt, yt, family = "cox", alpha = ML_ALPHA, nfolds = 10,
                   nlambda = ML_NLAMBDA, lambda.min.ratio = ML_LAMBDA_MIN_RATIO,
                   thresh = ML_THRESH, maxit = ML_MAXIT, penalty.factor = pfe)
  b_full <- as.matrix(coef(fit, s = "lambda.min"))[, 1]
  b_full <- b_full[b_full != 0]
  sel_me <- names(b_full)[names(b_full) %in% colnames(Et)]
  msg("[", cmp, "] locked model retains ", length(sel_me), " module scores: ",
      paste(sel_me, collapse = ", "))
  list(b_clin = b_clin, b_full = b_full, sel_me = sel_me, lambda = fit$lambda.min)
})
names(models) <- names(X)

# ---- 7. discrimination: apparent in discovery, external in CPTAC-3 ----
crow <- function(cohort, comparator, model, standardisation, y, lp, n_ev)
  data.table(cohort = cohort, comparator = comparator, model = model,
             standardisation = standardisation, n = length(lp), events = n_ev,
             t(cindex_ci(y, lp)))
lps <- list()   # linear predictors, reused by the bootstrap
res <- rbindlist(lapply(names(X), function(cmp) {
  M <- models[[cmp]]; Xt <- X[[cmp]]$Xt; Xv <- X[[cmp]]$Xv
  lps[[cmp]] <<- list(
    t_ref = lpv(Xt, M$b_clin), t_new = lpv(cbind(Xt, Et), M$b_full),
    v_ref = lpv(Xv, M$b_clin), v_new = lpv(cbind(Xv, Ev), M$b_full))
  out <- rbindlist(list(
    crow("TCGA-KIRC (discovery)", cmp, "comparator only", "cohort", yt, lps[[cmp]]$t_ref, sum(tcl$os_event)),
    crow("TCGA-KIRC (discovery)", cmp, "+ module eigengenes", "cohort", yt, lps[[cmp]]$t_new, sum(tcl$os_event)),
    crow("CPTAC-3 (validation)",  cmp, "comparator only", "cohort", yv, lps[[cmp]]$v_ref, sum(vcl$os_event)),
    crow("CPTAC-3 (validation)",  cmp, "+ module eigengenes", "cohort", yv, lps[[cmp]]$v_new, sum(vcl$os_event))))
  if (cmp == "clinical") {
    Xd <- X$clinical$Xv_disc
    lps[[cmp]]$vd_ref <<- lpv(Xd, M$b_clin)
    lps[[cmp]]$vd_new <<- lpv(cbind(Xd, Ev_disc), M$b_full)
    out <- rbind(out,
      crow("CPTAC-3 (validation)", cmp, "comparator only, discovery-fixed standardisation",
           "discovery", yv, lps[[cmp]]$vd_ref, sum(vcl$os_event)),
      crow("CPTAC-3 (validation)", cmp, "+ module eigengenes, discovery-fixed standardisation",
           "discovery", yv, lps[[cmp]]$vd_new, sum(vcl$os_event)))
  }
  out
}))
save_tsv(res, "09_validation_cindex.tsv")
print(res[, .(cohort, comparator, model, n, events, C = round(C, 3),
              lo = round(lo, 3), hi = round(hi, 3))])

# Paired patient bootstrap of delta C. Cohort labels match
# 09_validation_cindex.tsv, which stage 10 selects on. Discovery rows are
# in-sample.
brow <- function(cmp, std, cohort, y, lp_new, lp_ref) {
  b <- paired_boot_delta_c(y, lp_new, lp_ref)
  data.table(comparator = cmp, standardisation = std, cohort = cohort,
             n = length(lp_new), events = sum(y[, 2]),
             delta_C = round(b[["delta"]], 4), lo = round(b[["lo"]], 4),
             hi = round(b[["hi"]], 4), p_boot = signif(b[["p_boot"]], 3),
             prop_positive = round(b[["prop_positive"]], 3), n_boot = b[["n_boot"]])
}
msg("Paired bootstrap of delta C (", BOOT_B, " resamples per arm) ...")
delta <- rbindlist(list(
  brow("clinical",  "cohort",    "TCGA-KIRC (discovery)", yt, lps$clinical$t_new,  lps$clinical$t_ref),
  brow("clinical",  "cohort",    "CPTAC-3 (validation)",  yv, lps$clinical$v_new,  lps$clinical$v_ref),
  brow("clinical",  "discovery", "CPTAC-3 (validation)",  yv, lps$clinical$vd_new, lps$clinical$vd_ref),
  brow("augmented", "cohort",    "TCGA-KIRC (discovery)", yt, lps$augmented$t_new, lps$augmented$t_ref),
  brow("augmented", "cohort",    "CPTAC-3 (validation)",  yv, lps$augmented$v_new, lps$augmented$v_ref)))
save_tsv(delta, "09_delta_cindex_validation.tsv")
print(delta)

# ---- 8. replication of the discovery-significant modules ----
sig <- rbind(
  as.data.table(surv$mrna)[fdr_full < FDR_ALPHA,
    .(biotype = "mRNA", module, HR_tcga = HR_full, fdr_tcga = fdr_full)],
  as.data.table(surv$lnc)[fdr_full < FDR_ALPHA,
    .(biotype = "lncRNA", module, HR_tcga = HR_full, fdr_tcga = fdr_full)])
msg("Discovery-significant modules: ", nrow(sig))

Xv_aug   <- X$augmented$Xv
Xv_clin  <- X$clinical$Xv
# M1 as coded (MX -> M0, no stage-IV reconciliation) for the sensitivity arm.
Xv_aug_coded  <- Xv_aug;  Xv_aug_coded[, "M1"]  <- as.numeric(vcl$M1_as_coded)
Xv_clin_coded <- Xv_clin; Xv_clin_coded[, "M1"] <- as.numeric(vcl$M1_as_coded)
per_sd <- function(x) as.numeric((x - mean(x)) / sd(x))

# Per-module Cox fit returning HR, CI and p for the ME term.
me_fit <- function(me, Xcov, terms) {
  d <- data.frame(ME = me, Xcov[, terms, drop = FALSE])
  f <- coxph(yv ~ ., data = d)
  s <- summary(f)
  c(HR = s$conf.int["ME", "exp(coef)"], lo = s$conf.int["ME", "lower .95"],
    hi = s$conf.int["ME", "upper .95"], p = s$coefficients["ME", "Pr(>|z|)"],
    n = s$n, events = s$nevent)
}
replicate_modules <- function(terms, Xcov, Xcov_coded) {
  rbindlist(lapply(seq_len(nrow(sig)), function(i) {
    cn <- paste0(ifelse(sig$biotype[i] == "mRNA", "mRNA_ME", "lnc_ME"), sig$module[i])
    if (!cn %in% colnames(Ev)) return(NULL)
    a <- me_fit(per_sd(Ev[, cn]),      Xcov,       terms)
    d <- me_fit(per_sd(Ev_disc[, cn]), Xcov,       terms)
    m <- me_fit(per_sd(Ev[, cn]),      Xcov_coded, terms)
    data.table(biotype = sig$biotype[i], module = sig$module[i],
               HR_tcga  = round(sig$HR_tcga[i], 3),
               fdr_tcga = signif(sig$fdr_tcga[i], 3),
               HR_cptac = round(a[["HR"]], 3), lo = round(a[["lo"]], 3),
               hi = round(a[["hi"]], 3), p_cptac = signif(a[["p"]], 3),
               n = as.integer(a[["n"]]), events = as.integer(a[["events"]]),
               epv = round(a[["events"]] / (length(terms) + 1), 1),
               HR_cptac_disc = round(d[["HR"]], 3), lo_disc = round(d[["lo"]], 3),
               hi_disc = round(d[["hi"]], 3), p_disc = signif(d[["p"]], 3),
               HR_cptac_M1coded = round(m[["HR"]], 3), p_M1coded = signif(m[["p"]], 3))
  }))
}
finish_rep <- function(tb, label) {
  if (!nrow(tb)) return(tb)
  tb[, same_direction := (HR_tcga > 1) == (HR_cptac > 1)]
  tb[, fdr_cptac := signif(p.adjust(p_cptac, "BH"), 3)]
  tb[, adequate_epv := epv >= MIN_EPV]
  setcolorder(tb, c("biotype", "module", "HR_tcga", "fdr_tcga", "HR_cptac", "lo", "hi",
                    "p_cptac", "same_direction", "fdr_cptac", "n", "events", "epv",
                    "adequate_epv"))
  bt <- binom.test(sum(tb$same_direction), nrow(tb), 0.5)
  msg("[", label, "] same direction: ", sum(tb$same_direction), "/", nrow(tb),
      " (binomial p = ", signif(bt$p.value, 3), "); same direction and p < 0.05: ",
      sum(tb$same_direction & tb$p_cptac < 0.05), "; EPV = ", tb$epv[1])
  tb
}
aug_terms <- c(CLIN_TERMS, AUG_TERMS)
par_terms <- c("age", "male", "T_stage", "M1", "grade")
rep_tbl <- finish_rep(replicate_modules(aug_terms, Xv_aug, Xv_aug_coded),
                      "12-parameter, augmented covariates")
rep_par <- finish_rep(replicate_modules(par_terms, Xv_clin, Xv_clin_coded),
                      "6-parameter, parsimonious")
if (nrow(rep_tbl)) { save_tsv(rep_tbl, "09_module_replication.tsv"); print(rep_tbl) }
if (nrow(rep_par)) { save_tsv(rep_par, "09_module_replication_parsimonious.tsv"); print(rep_par) }

# ---- 9. axis replication: PC1 of the observed CPTAC-3 matrices ----
# PC1 as in 07, correlated with the STAR metrics and mean expression on the
# full clear-cell cohort.
axis_pc1 <- function(E, label) {
  x  <- scale(E, center = TRUE, scale = FALSE)
  sv <- svd(x, nu = 1, nv = 1)
  pc <- as.numeric(x %*% sv$v[, 1])
  mexp <- rowMeans(E)
  if (cor(pc, mexp) < 0) pc <- -pc
  share <- sv$d^2 / sum(sv$d^2)
  ct <- function(z) { r <- cor.test(pc, z, method = "spearman", exact = FALSE); c(r$estimate, r$p.value) }
  a <- ct(vqc$pct_noFeature); b <- ct(vqc$pct_multimapping)
  d <- ct(log10(vqc$assigned_reads)); e <- ct(mexp)
  list(pc = setNames(pc, rownames(E)),
       row = data.table(cohort = "CPTAC-3", matrix = label, n_samples = nrow(E),
                        n_genes = ncol(E),
                        var_share_PC1 = round(share[1], 4), var_share_PC2 = round(share[2], 4),
                        var_share_PC3 = round(share[3], 4), var_share_PC4 = round(share[4], 4),
                        var_share_PC5 = round(share[5], 4),
                        rho_noFeature = round(a[1], 3), p_noFeature = signif(a[2], 3),
                        rho_multimap = round(b[1], 3), p_multimap = signif(b[2], 3),
                        rho_log_depth = round(d[1], 3), p_log_depth = signif(d[2], 3),
                        rho_mean_expr = round(e[1], 3), p_mean_expr = signif(e[2], 3)))
}
ax_l <- axis_pc1(Vl_obs, "lncRNA")
ax_m <- axis_pc1(Vm_obs, "protein_coding")
save_tsv(rbind(ax_l$row, ax_m$row), "09_axis_replication_cptac.tsv")
print(rbind(ax_l$row, ax_m$row))

# Nested Cox sequence for the CPTAC-3 lncRNA axis on the evaluated set.
axd <- data.frame(axis_z = per_sd(ax_l$pc[vsamp]), X$clinical$Xv,
                  nf  = per_sd(vqm$pct_noFeature),
                  mm  = per_sd(vqm$pct_multimapping),
                  dep = per_sd(log10(vqm$assigned_reads)))
nested <- list(
  `Axis alone`                               = "axis_z",
  `+ age, sex, T, N, M1, ordinal grade`      = "axis_z + age + male + T_stage + N_pos + M1 + grade",
  `+ non-feature, multimapping, log10 depth` = "axis_z + age + male + T_stage + N_pos + M1 + grade + nf + mm + dep")
axis_nested <- rbindlist(lapply(names(nested), function(nm) {
  f <- coxph(as.formula(paste("yv ~", nested[[nm]])), data = axd)
  s <- summary(f)
  ph <- tryCatch(cox.zph(f)$table, error = function(e) NULL)
  data.table(cohort = "CPTAC-3", model = nm, n = s$n, events = s$nevent,
             HR = round(s$conf.int["axis_z", 1], 3), lo = round(s$conf.int["axis_z", 3], 3),
             hi = round(s$conf.int["axis_z", 4], 3), p = signif(s$coefficients["axis_z", 5], 3),
             C = round(s$concordance[1], 3),
             ph_p = if (is.null(ph)) NA_real_ else signif(ph["axis_z", "p"], 3))
}))
save_tsv(axis_nested, "09_axis_nested_cptac.tsv")
print(axis_nested)

# ---- 10. figures ----
res[, model_lab := factor(model, levels = rev(unique(model)))]
g <- ggplot(res, aes(model_lab, C, colour = cohort)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.4)) +
  coord_flip(ylim = c(0.45, 0.9)) +
  facet_grid(comparator ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("#B2182B", "#2166AC")) +
  labs(title = "Locked models: discovery vs external validation",
       subtitle = "Module membership, eigengene loadings and coefficients fixed in TCGA-KIRC",
       x = NULL, y = "Harrell C index (95% CI)", colour = NULL) +
  theme_bw() + theme(legend.position = "bottom")
save_fig(g, "09_validation_cindex", 8, 5.5)
res[, model_lab := NULL]

if (nrow(rep_tbl)) {
  rp <- copy(rep_tbl)[, lab := paste0(biotype, ":", module)]
  g2 <- ggplot(rp, aes(HR_cptac, reorder(lab, HR_cptac), colour = same_direction)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y") +
    geom_point(size = 2) +
    geom_point(aes(x = HR_tcga), shape = 4, size = 2.4, colour = "grey30") +
    scale_x_log10() +
    scale_colour_manual(values = c(`TRUE` = "#1B7837", `FALSE` = "#B2182B"),
                        name = "same direction") +
    labs(title = "Module replication in CPTAC-3",
         subtitle = "Circles = CPTAC-3 (augmented covariates, per within-cohort SD); crosses = TCGA-KIRC estimate",
         x = "Hazard ratio per 1 SD (95% CI)", y = NULL) +
    theme_bw()
  save_fig(g2, "09_module_replication", 7, 4.5)
}

# ---- 11. save the locked objects for 10, 11, 12, 20 and 22 ----
saveRDS(list(version = "v9",
             loadings = loadings,
             tsamp = tsamp, vsamp = vsamp,
             tcl = tcl, vcl = vcl,
             yt = yt, yv = yv,
             Et = Et, Ev = Ev, Ev_disc = Ev_disc,
             X = X,
             age_const = age_const,
             models = models,
             replication = rep_tbl, replication_parsimonious = rep_par),
        file.path(CACHE_DIR, "locked_model.rds"))

write_session_info("09_validate")
banner("09 | done")
