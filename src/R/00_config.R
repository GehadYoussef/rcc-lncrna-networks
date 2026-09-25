# 00_config.R: shared configuration and helpers for the bulk pipeline
# Sourced by every numbered script. Holds the tunable parameters, the output
# paths and the helpers used by more than one stage.

options(stringsAsFactors = FALSE, warn = 1)

# ---- Reproducibility --------------------------------------------------------
# blockwiseModules() pre-clusters genes with a randomised projective k-means,
# so module colours are stable between runs only with a fixed seed.
SEED <- 20260816
set.seed(SEED)

# ---- Paths ------------------------------------------------------------------
# R_DIR is set by the bootstrap at the top of each numbered script. The
# fallbacks allow direct sourcing.
if (!exists("R_DIR") || is.null(R_DIR)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .f    <- sub("^--file=", "", .args[grep("^--file=", .args)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/"))
           else if (dir.exists(file.path(getwd(), "analysis", "R"))) file.path(getwd(), "analysis", "R")
           else file.path(getwd(), "src", "R")
}
# Repository layout (src/R beside src/singlecell): downloads go to data/raw,
# caches to data/derived/cache, outputs to results/ and figures/. A standalone
# <project>/analysis/R layout is also supported. ANALYSIS_DIR is the output
# base in both.
.r_parent    <- normalizePath(dirname(R_DIR), winslash = "/", mustWork = FALSE)
REPO_LAYOUT  <- basename(.r_parent) == "src" && dir.exists(file.path(.r_parent, "singlecell"))
ANALYSIS_DIR <- if (REPO_LAYOUT) normalizePath(dirname(.r_parent), winslash = "/", mustWork = FALSE) else .r_parent
PROJECT_ROOT <- if (REPO_LAYOUT) ANALYSIS_DIR else normalizePath(dirname(ANALYSIS_DIR), winslash = "/", mustWork = FALSE)
DOWNLOAD_DIR <- if (REPO_LAYOUT) file.path(ANALYSIS_DIR, "data", "raw") else ANALYSIS_DIR

GDC_DIR     <- file.path(if (REPO_LAYOUT) DOWNLOAD_DIR else PROJECT_ROOT, "GDCdata", "TCGA-KIRC")
EXPR_DIR    <- file.path(GDC_DIR, "Transcriptome_Profiling", "Gene_Expression_Quantification")
CLIN_XML_DIR<- file.path(GDC_DIR, "Clinical", "Clinical_Supplement")

CACHE_DIR   <- if (REPO_LAYOUT) file.path(ANALYSIS_DIR, "data", "derived", "cache") else file.path(ANALYSIS_DIR, "cache")
RESULTS_DIR <- file.path(ANALYSIS_DIR, "results")
FIG_DIR     <- file.path(ANALYSIS_DIR, "figures")
LOG_DIR     <- file.path(ANALYSIS_DIR, "logs")

for (d in c(CACHE_DIR, RESULTS_DIR, FIG_DIR, LOG_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---- Cohort definition ------------------------------------------------------
PROJECT_ID <- "TCGA-KIRC"

# Only primary tumours enter the analysis, so each patient contributes one row.
KEEP_SAMPLE_TYPES <- "Primary Tumor"

# Overall survival is censored at 10 years, beyond which TCGA follow-up is sparse.
OS_CENSOR_DAYS <- 3650

# ---- Expression filtering ---------------------------------------------------
# WGCNA is run on log2(FPKM + 1). Raw FPKM gives zero-inflated lncRNA modules
# driven by a few high expressors and shared zero patterns.
# Genes must exceed MIN_FPKM in at least MIN_FRAC of samples. lncRNAs have a
# lower threshold because they are intrinsically low-abundance.
MRNA_MIN_FPKM   <- 1.0
MRNA_MIN_FRAC   <- 0.25
LNC_MIN_FPKM    <- 0.30
LNC_MIN_FRAC    <- 0.20

# Keep the most variable genes (by MAD) after filtering so the network fits in
# one block. Inf disables.
MRNA_TOP_N_MAD  <- 12000
LNC_TOP_N_MAD   <- 6000

# ---- WGCNA parameters -------------------------------------------------------
NETWORK_TYPE    <- "signed"        # keeps positively and negatively correlated genes apart
TOM_TYPE        <- "signed"
MIN_MODULE_SIZE <- 30
MERGE_CUT_HEIGHT<- 0.25
RSQ_CUT         <- 0.85            # scale-free topology target
DEFAULT_POWER   <- 12              # fallback for signed networks if fit fails
# One block, so modules may span the whole transcriptome. The default
# maxBlockSize of 5000 splits the genes into blocks and destabilises modules.
MAX_BLOCK_SIZE  <- 25000
DEEP_SPLIT      <- 2               # WGCNA default
# Genes whose correlation with their own module eigengene falls below this
# return to grey. It is the main driver of the grey fraction.
MIN_KME_TO_STAY <- 0.3             # WGCNA default

# lncRNA network parameters, chosen by 02b_lncrna_tuning.R from 96
# configurations on stability and coherence only (never survival, stage or
# grade). NA falls back to the automatic rule.
LNC_POWER       <- 10
LNC_DEEPSPLIT   <- 2
LNC_MERGE       <- 0.15
LNC_MINKME      <- 0.3
# "raw", "resid_PC1" or "resid_PC1_2": project out the leading principal
# component(s) before network construction. Sensitivity option, the primary
# analysis uses "raw".
LNC_VARIANT     <- "raw"

# ---- Technical covariate adjustment before network construction -------------
# The STAR metrics are regressed out of expression before networks are built,
# so they cannot shape module membership. The non-feature fraction is the
# leading axis of the lncRNA matrix (rho = 0.90 with PC1, stage 07) and is
# also a covariate in the principal survival model.
# Tumour composition (ESTIMATE) is kept because immune infiltrate is genuine
# ccRCC biology.
ADJUST_TECHNICAL <- TRUE
TECH_COVARIATES  <- c("pct_noFeature", "pct_multimapping", "log_depth")

# Residualise each gene on a matrix of technical covariates (samples x k).
remove_technical <- function(expr, covars) {
  keep <- apply(covars, 2, function(z) is.finite(sd(z)) && sd(z) > 0)
  covars <- covars[, keep, drop = FALSE]
  if (!ncol(covars)) return(expr)
  X <- cbind(1, scale(covars))
  fit <- qr.lstsq <- qr.solve(crossprod(X), crossprod(X, expr))
  resid <- expr - X %*% fit
  # add gene means back to keep the original scale
  sweep(resid, 2, colMeans(expr), "+")
}

# Project out the first k principal components across samples.
remove_leading_pcs <- function(expr, variant) {
  if (identical(variant, "raw")) return(expr)
  k <- if (identical(variant, "resid_PC1")) 1L else 2L
  x <- scale(expr, center = TRUE, scale = FALSE)
  s <- svd(x, nu = k, nv = k)
  x - s$u[, 1:k, drop = FALSE] %*% diag(s$d[1:k], k, k) %*%
      t(s$v[, 1:k, drop = FALSE])
}

# ---- Statistics -------------------------------------------------------------
FDR_ALPHA       <- 0.05
DCOR_N_PERM     <- 5000            # draws for the global dCor null distribution
DCOR_TOP_HUBS   <- 150             # hub genes per module carried into dCor
DCOR_MAX_PER_SIDE <- 400           # global cap on genes per side (runtime)
DCOR_CONFIRM_TOP  <- 200           # top pairs re-tested with an exact permutation test
ML_MAX_FEATURES <- 200             # candidate lncRNAs entering the penalised model
ML_N_REPEATS    <- 10              # repeats of k-fold CV for the ML model
ML_N_FOLDS      <- 10
# Elastic net with a truncated lambda path. Hub lncRNAs of one module are
# near-collinear, and pure LASSO converges poorly on them.
ML_ALPHA        <- 0.5
ML_NLAMBDA      <- 50
ML_LAMBDA_MIN_RATIO <- 0.05
ML_THRESH       <- 1e-6
ML_MAXIT        <- 50000
EVAL_TIMES_YRS  <- c(1, 3, 5)

# ---- Threads ----------------------------------------------------------------
N_THREADS <- max(1L, min(16L, parallel::detectCores() - 2L))

# ---- Small helpers ----------------------------------------------------------
# Windows MAX_PATH workaround. Deep GDC paths can exceed the Win32 limit, so
# file.exists() and fread() report the files as missing. The \\?\ prefix lifts
# the limit and is honoured by R's file functions and data.table::fread().
winlong <- function(p) {
  if (.Platform$OS.type != "windows") return(p)
  p <- gsub("/", "\\", p, fixed = TRUE)
  ifelse(startsWith(p, "\\\\?\\"), p, paste0("\\\\?\\", p))
}

# Observed (non-residualised) expression for a network. Modules may be defined
# on residualised data, but correlations between transcripts and clinical
# predictors are measured on observed levels.
obs_expr <- function(net) if (!is.null(net$expr_obs)) net$expr_obs else net$expr

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

banner <- function(x) {
  cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
}

# ---- Figure output ----------------------------------------------------------
# Every figure is written as SVG and as a 600 dpi PNG. svglite keeps text
# editable, whereas grDevices::svg converts it to paths.
FIG_DPI     <- 600
FIG_FORMATS <- c("svg", "png")

save_fig <- function(plot, name, width = 7, height = 5) {
  for (fmt in FIG_FORMATS) {
    f <- file.path(FIG_DIR, paste0(name, ".", fmt))
    ggplot2::ggsave(f, plot, width = width, height = height, units = "in",
                    dpi = FIG_DPI,
                    device = if (fmt == "svg") svglite::svglite else ragg::agg_png)
  }
  invisible(name)
}

# For base graphics and composite objects that ggsave cannot handle. The
# expression is re-evaluated once per device.
save_fig_base <- function(name, width, height, expr) {
  e <- substitute(expr); pf <- parent.frame()
  for (fmt in FIG_FORMATS) {
    f <- file.path(FIG_DIR, paste0(name, ".", fmt))
    if (fmt == "svg") svglite::svglite(f, width = width, height = height)
    else ragg::agg_png(f, width = width, height = height, units = "in",
                       res = FIG_DPI)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    eval(e, envir = pf)
    grDevices::dev.off()
  }
  invisible(name)
}

save_tsv <- function(x, file) {
  utils::write.table(x, file.path(RESULTS_DIR, file),
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(file.path(RESULTS_DIR, file))
}

# Writes sessionInfo() for a stage to logs/.
write_session_info <- function(tag) {
  f <- file.path(LOG_DIR, paste0("sessionInfo_", tag, ".txt"))
  con <- file(f, open = "wt"); sink(con); print(sessionInfo()); sink(); close(con)
  invisible(f)
}

# ---- Shared helpers ---------------------------------------------------------
# Used by stages 03 onward so that discovery, validation and subtype cohorts
# are scored identically with the discovery module rotations.

# ---- Analysis parameters ----------------------------------------------------
FOLDWISE_K           <- 3      # folds for the network-rebuild cross-validation (12)
FOLDWISE_N_ASSIGN    <- 5      # random fold assignments, both networks rebuilt per fold
BOOT_B               <- 2000   # bootstrap resamples for concordance increments
PRIMARY_HORIZON_YR   <- 3      # primary horizon for absolute-risk metrics in CPTAC-3
SECONDARY_HORIZON_YR <- 5
MIN_EPV              <- 5      # events per parameter required for inference
# Each cohort is residualised on its own STAR metrics (per-gene OLS). CPTAC-3
# used ribo-depleted total RNA, so its non-feature fraction (about 40%) lies
# outside the polyA TCGA range (about 2-20%).
# Libraries with fewer assigned reads than MIN_ASSIGNED_READS or a non-feature
# fraction above MAX_NOFEATURE_PCT are excluded in 01, before one sample is
# chosen per patient. In TCGA-KIRC this removes one outlying library.
MIN_ASSIGNED_READS <- 1e7
MAX_NOFEATURE_PCT  <- 30
CLIN_TERMS <- c("age", "male", "T_stage", "N_pos", "M1", "grade")
AUG_TERMS  <- c("stromal", "immune", "noFeature", "multimap", "libsize_z")

# The three technical covariates, in the one order used everywhere.
tech_covariates <- function(qc) {
  cbind(pct_noFeature    = as.numeric(qc$pct_noFeature),
        pct_multimapping = as.numeric(qc$pct_multimapping),
        log_depth        = log10(as.numeric(qc$assigned_reads)))
}

# Residualisation split into fit and apply, so a held-out fold uses training
# coefficients. apply_technical(expr, fit_technical(expr, cov), cov) equals
# remove_technical(expr, cov).
fit_technical <- function(expr, covars) {
  mu  <- colMeans(covars); sdv <- apply(covars, 2, sd)
  keep <- is.finite(sdv) & sdv > 0
  X <- cbind(1, scale(covars[, keep, drop = FALSE], center = mu[keep], scale = sdv[keep]))
  B <- qr.solve(crossprod(X), crossprod(X, expr))
  list(cov_center = mu[keep], cov_scale = sdv[keep], keep = keep, B = B,
       gene_mean = colMeans(expr))
}
apply_technical <- function(expr, fit, covars) {
  X <- cbind(1, scale(covars[, fit$keep, drop = FALSE],
                      center = fit$cov_center, scale = fit$cov_scale))
  resid <- expr - X %*% fit$B
  sweep(resid, 2, fit$gene_mean, "+")
}

# STAR summary from one augmented_star_gene_counts.tsv (unstranded column).
# assigned_reads sums counts over ENSG genes. Percentages are of assigned reads
# plus the four STAR summary rows.
read_star_summary <- function(path) {
  d <- data.table::fread(path, skip = 1, showProgress = FALSE,
                         select = c("gene_id", "unstranded"))
  v <- setNames(as.numeric(d$unstranded), d$gene_id)
  assigned <- sum(v[grepl("^ENSG", names(v))], na.rm = TRUE)
  qn  <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
  tot <- assigned + sum(v[qn], na.rm = TRUE)
  data.table::data.table(
    pct_unmapped     = 100 * v[["N_unmapped"]]     / tot,
    pct_multimapping = 100 * v[["N_multimapping"]] / tot,
    pct_noFeature    = 100 * v[["N_noFeature"]]    / tot,
    pct_ambiguous    = 100 * v[["N_ambiguous"]]    / tot,
    assigned_reads   = assigned)
}

# ---- Module loadings and scores ---------------------------------------------
# A module eigengene is PC1 of the module's gene-standardised expression,
# sign-aligned to mean module expression (WGCNA convention).
# fit_module_loadings() fits it once on the residualised discovery matrix and
# score_modules() applies the stored rotation to any cohort. Gene and score
# standardisation come from the target cohort ("cohort", primary) or from
# discovery ("discovery", sensitivity analysis).
fit_module_loadings <- function(expr, gene_tbl, prefix, min_genes = 10) {
  gt   <- data.table::as.data.table(gene_tbl)[module != "grey"]
  sets <- split(gt$gene_id, gt$module)
  out  <- lapply(names(sets), function(m) {
    g <- intersect(sets[[m]], colnames(expr))
    if (length(g) < min_genes) return(NULL)
    X   <- expr[, g, drop = FALSE]
    mu  <- colMeans(X); sdv <- apply(X, 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
    Z   <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
    sv  <- svd(Z, nu = 1, nv = 1)
    rot <- setNames(sv$v[, 1], g)
    s   <- as.numeric(Z %*% rot)
    if (cor(s, rowMeans(Z)) < 0) { rot <- -rot; s <- -s }
    list(module = m, genes = g, gene_center = mu, gene_scale = sdv, rot = rot,
         score_center = mean(s), score_scale = sd(s),
         var_explained = sv$d[1]^2 / sum(Z^2))
  })
  names(out) <- paste0(prefix, names(sets))
  out[!vapply(out, is.null, logical(1))]
}

score_modules <- function(expr, loadings,
                          gene_standardise  = c("cohort", "discovery"),
                          score_standardise = c("cohort", "discovery", "none"),
                          min_genes = 10) {
  gene_standardise  <- match.arg(gene_standardise)
  score_standardise <- match.arg(score_standardise)
  M <- matrix(NA_real_, nrow(expr), length(loadings),
              dimnames = list(rownames(expr), names(loadings)))
  for (nm in names(loadings)) {
    L <- loadings[[nm]]; g <- intersect(L$genes, colnames(expr))
    if (length(g) < min_genes) next
    X <- expr[, g, drop = FALSE]
    if (gene_standardise == "cohort") {
      mu <- colMeans(X); sdv <- apply(X, 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
    } else { mu <- L$gene_center[g]; sdv <- L$gene_scale[g] }
    Z <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
    s <- as.numeric(Z %*% L$rot[g])
    if (score_standardise == "cohort")         s <- (s - mean(s)) / sd(s)
    else if (score_standardise == "discovery") s <- (s - L$score_center) / L$score_scale
    M[, nm] <- s
  }
  M[, colSums(is.na(M)) == 0, drop = FALSE]
}

# Discovery loadings, fitted on the residualised matrices in networks.rds and
# cached in module_loadings.rds.
discovery_loadings <- function(nets, force = FALSE) {
  f <- file.path(CACHE_DIR, "module_loadings.rds")
  if (!force && file.exists(f)) return(readRDS(f))
  L <- c(fit_module_loadings(nets$mrna$expr, nets$mrna$gene_tbl, "mRNA_ME"),
         fit_module_loadings(nets$lnc$expr,  nets$lnc$gene_tbl,  "lnc_ME"))
  saveRDS(L, f); L
}
# Unit-variance discovery module scores on each network's own sample set.
# These equal WGCNA's moduleEigengenes() output up to scale (checked in 03).
discovery_scores <- function(nets, loadings = discovery_loadings(nets)) {
  list(mrna = score_modules(nets$mrna$expr, loadings[grep("^mRNA_", names(loadings))]),
       lnc  = score_modules(nets$lnc$expr,  loadings[grep("^lnc_",  names(loadings))]))
}

# ---- Clinical design matrices -----------------------------------------------
# "clinical": age, sex, T category, nodal status, M1 and ordinal grade.
# "augmented" adds ESTIMATE stromal and immune scores and the three STAR
# metrics. Missing covariates are returned as NA and callers drop those rows.
clinical_design <- function(cl, set = c("clinical", "augmented"), est = NULL, qc = NULL,
                            age_center = NULL, age_scale = NULL) {
  set <- match.arg(set)
  age_z <- if (is.null(age_center)) as.numeric(scale(cl$age))
           else (cl$age - age_center) / age_scale
  X <- cbind(age     = age_z,
             male    = as.numeric(cl$sex == "male"),
             T_stage = as.numeric(cl$T_stage),
             N_pos   = as.numeric(cl$N_pos),
             M1      = as.numeric(cl$M1),
             grade   = as.numeric(cl$grade_num))
  if (set == "augmented") {
    stopifnot(!is.null(est), !is.null(qc))
    X <- cbind(X,
               stromal   = as.numeric(scale(est$StromalScore)),
               immune    = as.numeric(scale(est$ImmuneScore)),
               noFeature = as.numeric(scale(qc$pct_noFeature)),
               multimap  = as.numeric(scale(qc$pct_multimapping)),
               libsize_z = as.numeric(scale(log10(qc$assigned_reads))))
  }
  rownames(X) <- cl$sample_barcode
  X
}

# ---- Discrimination, bootstrap, absolute risk --------------------------------
cindex <- function(y, lp) unname(survival::concordance(y ~ lp, reverse = TRUE)$concordance)
cindex_ci <- function(y, lp) {
  cc <- survival::concordance(y ~ lp, reverse = TRUE); se <- sqrt(cc$var)
  c(C = unname(cc$concordance), se = se,
    lo = unname(cc$concordance) - 1.96 * se, hi = unname(cc$concordance) + 1.96 * se)
}
# Paired patient bootstrap of the difference in concordance between two linear
# predictors evaluated on the same patients (predictors held fixed).
paired_boot_delta_c <- function(y, lp_new, lp_ref, B = BOOT_B, seed = SEED, min_events = 5) {
  set.seed(seed); n <- length(lp_new); ev <- y[, 2]
  b <- vapply(seq_len(B), function(k) {
    i <- sample.int(n, n, replace = TRUE)
    if (sum(ev[i]) < min_events) return(NA_real_)
    cindex(y[i], lp_new[i]) - cindex(y[i], lp_ref[i])
  }, numeric(1))
  b <- b[is.finite(b)]
  c(delta = cindex(y, lp_new) - cindex(y, lp_ref),
    lo = unname(quantile(b, 0.025)), hi = unname(quantile(b, 0.975)),
    p_boot = 2 * min(mean(b <= 0), mean(b >= 0)), prop_positive = mean(b > 0),
    n_boot = length(b))
}
# Breslow baseline cumulative hazard for a locked linear predictor, slope fixed
# at 1: H0(t) = sum over event times t_i <= t of d_i / sum_{j at risk} exp(lp_j).
# basehaz() on an offset-only coxph does not return the curve at offset zero.
baseline_risk_fun <- function(y_train, lp_train) {
  tm <- y_train[, 1]; ev <- y_train[, 2]; w <- exp(lp_train)
  ut <- sort(unique(tm[ev == 1]))
  H0 <- cumsum(vapply(ut, function(ti) sum(ev[tm == ti]) / sum(w[tm >= ti]), numeric(1)))
  function(lp, t) {
    k <- findInterval(t, ut)
    h <- if (k == 0) 0 else H0[k]
    1 - exp(-h * exp(lp))
  }
}
# IPCW Brier score at time t (Graf 1999). The censoring distribution is the
# Kaplan-Meier estimate in the scored cohort, looked up as a step function.
brier_ipcw <- function(y, pred, t) {
  tm <- y[, 1]; ev <- y[, 2]
  cn <- survival::survfit(survival::Surv(tm, 1 - ev) ~ 1)
  G  <- function(u) approx(c(0, cn$time), c(1, cn$surv), xout = pmin(u, max(cn$time)),
                           method = "constant", rule = 2, f = 0)$y
  w <- ifelse(tm <= t & ev == 1, 1 / pmax(G(tm), 1e-6),
       ifelse(tm >  t,           1 / pmax(G(t),  1e-6), 0))
  mean(w * (as.numeric(tm <= t & ev == 1) - pred)^2)
}
# Calibration summary at horizon t: slope of the linear predictor in the
# scored cohort, and observed (Kaplan-Meier) over expected (mean predicted) risk.
calibration_summary <- function(y, lp, pred_t, t) {
  f  <- survival::coxph(y ~ lp)
  km <- summary(survival::survfit(y ~ 1), times = t, extend = TRUE)
  obs <- 1 - km$surv
  c(slope = unname(coef(f)), slope_lo = unname(confint(f)[1]), slope_hi = unname(confint(f)[2]),
    observed = obs, expected = mean(pred_t), OE = obs / mean(pred_t),
    n_at_risk = km$n.risk)
}
# Calibration bins (quintiles of predicted risk) with Kaplan-Meier observed risk.
calibration_bins <- function(y, pred_t, t, n_bins = 5) {
  g <- cut(pred_t, breaks = unique(quantile(pred_t, seq(0, 1, length.out = n_bins + 1))),
           include.lowest = TRUE, labels = FALSE)
  data.table::rbindlist(lapply(sort(unique(g)), function(k) {
    i  <- which(g == k)
    km <- summary(survival::survfit(survival::Surv(y[i, 1], y[i, 2]) ~ 1), times = t, extend = TRUE)
    data.table::data.table(bin = k, n = length(i), events_by_t = sum(y[i, 1] <= t & y[i, 2] == 1),
                           predicted = mean(pred_t[i]), observed = 1 - km$surv,
                           obs_lo = 1 - km$upper, obs_hi = 1 - km$lower)
  }))
}
# Net benefit at horizon t across threshold probabilities (Vickers, censoring by
# Kaplan-Meier within the treated group).
net_benefit <- function(y, risk, t, thr = seq(0.02, 0.60, by = 0.01)) {
  vapply(thr, function(pt) {
    fl <- risk >= pt; if (!any(fl)) return(0)
    s <- summary(survival::survfit(survival::Surv(y[fl, 1], y[fl, 2]) ~ 1), times = t, extend = TRUE)$surv
    p <- mean(fl); (1 - s) * p - s * p * (pt / (1 - pt))
  }, numeric(1))
}
net_benefit_treat_all <- function(y, t, thr = seq(0.02, 0.60, by = 0.01)) {
  s <- summary(survival::survfit(y ~ 1), times = t, extend = TRUE)$surv
  (1 - s) - s * (thr / (1 - thr))
}

# ---- The reference protein-coding programme --------------------------------
# Module colours change when the network is rebuilt, so the protein-coding
# module carrying the fatty-acid and organic-acid catabolic programme ("black")
# is located by its GO enrichment (06). Fallback: the most protective
# significant protein-coding module in the principal survival model (03).
catabolic_reference_module <- function() {
  f_go <- file.path(RESULTS_DIR, "06_GO_enrichment_all_modules.tsv")
  if (file.exists(f_go)) {
    go <- data.table::fread(f_go)
    hit <- go[grepl("organic acid catabolic|fatty acid.*(catabolic|oxidation)",
                    Description, ignore.case = TRUE)]
    if (nrow(hit)) {
      best <- hit[, .(p = min(p.adjust)), by = module][order(p)]
      msg("Reference catabolic protein-coding module (by GO enrichment): ",
          best$module[1], " (adjusted p ", signif(best$p[1], 2), ")")
      return(best$module[1])
    }
  }
  f_sv <- file.path(RESULTS_DIR, "03_mRNA_module_survival.tsv")
  sv <- data.table::fread(f_sv)[HR_full < 1][order(fdr_full)]
  msg("Reference catabolic protein-coding module (fallback, most protective): ", sv$module[1])
  sv$module[1]
}

# ---- Pan-cancer decision rules (stages 38 to 41) -----------------------------
# Fixed before any pan-cancer project was analysed and recorded in
# results/38_decision_rules_lock.json. Stages 39 to 41 stop if these values
# differ from the lock (verify_decision_lock() in pan_helpers.R).
PAN_DATA_DIR          <- file.path(DOWNLOAD_DIR, "pancancer_data")
# Eligibility for inference: QC-passing primary tumours, one per patient.
PAN_MIN_TUMOURS       <- 100
# Solid-tissue normals are analysed as the tumour-free control where at least
# this many QC-passing normal libraries (one per patient) exist.
PAN_MIN_NORMALS       <- 20
# TCGA-LAML (blood, no normals) is reported descriptively only and excluded
# from prevalence counts and the meta-analysis.
PAN_EXCLUDE_INFERENCE <- "TCGA-LAML"
PAN_TUMOUR_TYPES      <- c("Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood")
PAN_NORMAL_TYPE       <- "Solid Tissue Normal"
# Rule A, axis present: the project's lncRNA PC1 has Spearman |rho| >= 0.5 with
# the non-feature fraction, and this exceeds the protein-coding PC1 |rho| by
# at least 0.3.
PAN_AXIS_RHO_MIN      <- 0.5
PAN_AXIS_MARGIN       <- 0.3
# Rule C, network reorganised: in the observed-expression lncRNA network the
# largest non-grey module holds >= 40% of genes, its eigengene has Pearson
# |r| >= 0.7 with the non-feature fraction, and the adjusted Rand index between
# observed and residualised partitions (grey as a class) is < 0.2. Networks use
# the discovery lncRNA parameters in every project.
PAN_NET_LARGEST_MIN   <- 0.40
PAN_NET_EIGEN_R_MIN   <- 0.7
PAN_NET_ARI_MAX       <- 0.2
# Stage 41 survival: TCGA Clinical Data Resource overall survival, censored at
# OS_CENSOR_DAYS. Primary Cox model: metric per SD + age + sex (+ ordinal stage
# where >= 70% of the project has it), stratified by plate (plates with < 10
# patients pooled). Secondary: no plate stratum. Projects with >= 20 events
# enter the REML random-effects meta-analysis.
PAN_SURV_MIN_EVENTS   <- 20
PAN_STAGE_MIN_FRAC    <- 0.70
PAN_PLATE_MIN_N       <- 10
