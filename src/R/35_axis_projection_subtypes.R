# 35_axis_projection_subtypes.R: discovery lncRNA axis projected onto KIRP and KICH
# Projects TCGA-KIRP and TCGA-KICH primary tumours onto the discovery (TCGA-KIRC) lncRNA
# PC1 of 07, centred on discovery gene means, with the discovery protein-coding PC1 as a
# control. Each projected score and each cohort's own PC1 is correlated (Spearman, BH)
# with the non-feature fraction, multimapping fraction, log10 assigned reads and mean
# expression. Gene sets, scaling, sign convention and columns match
# 25_normal_axis_projection.tsv, so the two tables row-bind (checked in section 5 if 25 has run).
# Inputs: cache/dataset.rds, networks.rds, lnc_global_axis.rds,
# subtype_KIRP.rds and subtype_KICH.rds (11), and 07_pc_variance_explained.tsv.
# Outputs: 35_axis_projection_subtypes.tsv, 35_axis_recompute_check.tsv.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages(library(data.table))
banner("35 | The discovery library-quality axis projected onto TCGA-KIRP and TCGA-KICH")
set.seed(SEED)
t_start <- Sys.time()

SUBTYPES <- c("TCGA-KIRP", "TCGA-KICH")

ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
axis <- as.data.table(readRDS(file.path(CACHE_DIR, "lnc_global_axis.rds")))
cohort <- as.data.table(ds$cohort_full)
if (!"assigned_reads" %in% names(cohort)) cohort[, assigned_reads := libsize]

# The matrices 07 decomposed: observed log2(FPKM + 1) on the network samples.
E_lnc_disc <- obs_expr(nets$lnc)
E_pc_disc  <- obs_expr(nets$mrna)
genes_lnc  <- colnames(E_lnc_disc); genes_pc <- colnames(E_pc_disc)
msg("Discovery matrices: lncRNA ", nrow(E_lnc_disc), " x ", ncol(E_lnc_disc),
    "; protein-coding ", nrow(E_pc_disc), " x ", ncol(E_pc_disc))

# ---- helpers, identical to those of 25_matched_normal_control.R ----
spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(c(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}
# Benjamini-Hochberg over finite p-values only (p.adjust() would count NA rows).
bh <- function(p) { out <- rep(NA_real_, length(p)); ok <- is.finite(p)
  if (any(ok)) out[ok] <- p.adjust(p[ok], "BH"); out }

# Leading principal components of a samples x genes matrix, as in 07:
# column-centred, unscaled SVD. u1 is the unit-norm score 07 caches as
# lnc_axis, v1 the rotation, d1 the singular value.
# PC1 sign is arbitrary and is fixed by three anchors, all written to the results table:
#   1. correlation with mean expression over all genes (07's rule)
#   2. correlation with mean expression over the SIGN_TOP_GENES largest-|loading| genes
#   3. correlation with the projection onto the discovery PC1 of the same biotype
#      (own PC1s only)
# Anchor 1 is used when |r| >= SIGN_MIN_R and anchor 2 does not contradict it, anchor 2
# when anchor 1 is weak, and anchor 3 when anchors 1 and 2 are decisive but disagree,
# since its orientation is fixed outside the matrix. Rows with a weak anchor 1 or
# disagreeing anchors have pc1_sign_interpretable = FALSE, and their rho values are
# read as magnitudes.
SIGN_MIN_R     <- 0.30
SIGN_TOP_GENES <- 500
# ref: the libraries' projection onto the discovery PC1 of the same biotype, in
# the row order of E. NULL for the discovery matrices themselves.
pca_lead <- function(E, k = 5, ref = NULL) {
  mu <- colMeans(E)
  Xc <- scale(E, center = TRUE, scale = FALSE)
  s  <- svd(Xc, nu = k, nv = k)
  u1 <- s$u[, 1]; v1 <- s$v[, 1]
  r_all <- stats::cor(u1, rowMeans(E))
  gtop  <- order(abs(v1), decreasing = TRUE)[seq_len(min(SIGN_TOP_GENES, ncol(E)))]
  r_top <- stats::cor(u1, rowMeans(E[, gtop, drop = FALSE]))
  r_ref <- NA_real_
  if (!is.null(ref) && length(ref) == length(u1)) {
    rr <- as.numeric(ref)
    if (sum(is.finite(rr)) >= 5 && is.finite(stats::sd(rr, na.rm = TRUE)) &&
        stats::sd(rr, na.rm = TRUE) > 0)
      r_ref <- stats::cor(u1, rr, use = "complete.obs")
  }
  dec <- function(r) is.finite(r) && abs(r) >= SIGN_MIN_R
  disagree <- dec(r_all) && dec(r_top) && sign(r_all) != sign(r_top)
  if (disagree && dec(r_ref)) {
    rule <- paste0("projection on the discovery PC1 (tie-break: the all-gene and the ",
                   length(gtop), "-largest-|loading| anchors are each decisive and disagree)")
    stat <- r_ref
    resolved <- "third anchor: projection on the discovery PC1"
  } else if (disagree) {
    rule <- "CONTRADICTORY"
    stat <- r_top
    resolved <- "none: anchors 1 and 2 disagree and the projection anchor is not decisive"
  } else if (dec(r_all)) {
    rule <- "mean expression, all genes"
    stat <- r_all
    resolved <- "first anchor: mean expression over all genes"
  } else if (dec(r_top)) {
    rule <- paste0("mean expression, ", length(gtop), " largest |loadings|")
    stat <- r_top
    resolved <- paste0("second anchor: mean expression over the ", length(gtop),
                       " largest |loadings|")
  } else {
    rule <- "UNDETERMINED"; stat <- r_all
    resolved <- "none: no anchor reached the threshold"
  }
  flip <- is.finite(stat) && stat < 0
  if (flip) { u1 <- -u1; v1 <- -v1 }
  orient <- function(x) if (flip) -x else x
  primary_weak <- !dec(r_all)
  list(mu = mu, v1 = setNames(v1, colnames(E)), d1 = s$d[1],
       u1 = setNames(u1, rownames(E)), score = setNames(u1 * s$d[1], rownames(E)),
       var_share = s$d^2 / sum(s$d^2),
       sign_rule = rule, sign_stat = unname(orient(stat)),
       sign_r_mean_all = unname(orient(r_all)),
       sign_r_mean_top = unname(orient(r_top)),
       sign_r_projected = unname(orient(r_ref)),
       sign_resolved_by = resolved,
       # TRUE where 07's all-gene anchor could not fix the sign.
       sign_primary_weak = primary_weak,
       # TRUE where anchors 1 and 2 both reach SIGN_MIN_R with opposite signs.
       sign_anchors_disagree = disagree,
       # FALSE where the sign is a convention: rho values are then magnitudes.
       sign_interpretable = !(primary_weak || disagree))
}
# Prints the anchors of every PC1. Stops if a matrix has no decisive anchor.
# Disagreeing anchors are logged and flagged, not fatal.
assert_sign_convention <- function(objs) {
  objs <- objs[!vapply(objs, is.null, logical(1))]
  gv <- function(f, mode) vapply(objs, function(z) z[[f]], mode)
  tbl <- data.table(matrix = names(objs),
                    sign_rule = gv("sign_rule", character(1)),
                    sign_stat = round(gv("sign_stat", numeric(1)), 3),
                    r_mean_all_genes = round(gv("sign_r_mean_all", numeric(1)), 3),
                    r_mean_top_genes = round(gv("sign_r_mean_top", numeric(1)), 3),
                    r_projected_on_discovery_PC1 = round(gv("sign_r_projected", numeric(1)), 3),
                    primary_anchor_weak = gv("sign_primary_weak", logical(1)),
                    anchors_disagree = gv("sign_anchors_disagree", logical(1)),
                    resolved_by = gv("sign_resolved_by", character(1)),
                    sign_interpretable = gv("sign_interpretable", logical(1)))
  print(tbl, row.names = FALSE)
  bad <- tbl$sign_rule == "UNDETERMINED" | !is.finite(tbl$sign_stat) | tbl$sign_stat < SIGN_MIN_R
  if (any(bad))
    stop("PC1 sign convention undetermined for: ", paste(tbl$matrix[bad], collapse = ", "),
         " (no anchor reached |r| = ", SIGN_MIN_R, ")")
  dis <- tbl[anchors_disagree == TRUE]
  if (nrow(dis)) for (i in seq_len(nrow(dis)))
    msg("SIGN CONVENTION, ANCHORS DISAGREE: ", dis$matrix[i],
        " -- all-gene mean-expression anchor ", dis$r_mean_all_genes[i], ", ",
        SIGN_TOP_GENES, "-largest-|loading| anchor ", dis$r_mean_top_genes[i],
        ", projection on the discovery PC1 ", dis$r_projected_on_discovery_PC1[i],
        "; orientation ", if (grepl("^projection", dis$sign_rule[i]))
          "re-anchored on the projection (which sides with the largest-|loading| anchor)"
        else "left on the largest-|loading| anchor, unresolved",
        ". The SIGN of every correlation reported for this component is a convention: quote |rho| only.")
  wk <- tbl[primary_anchor_weak == TRUE]
  if (nrow(wk)) for (i in seq_len(nrow(wk)))
    msg("SIGN CONVENTION, WEAK PRIMARY ANCHOR: ", wk$matrix[i],
        " -- all-gene mean-expression anchor ", wk$r_mean_all_genes[i],
        " is below |", SIGN_MIN_R, "|; sign fixed by ", wk$resolved_by[i],
        ". Quote |rho| only.")
  msg("PC1 sign convention: an anchor was decisive for all ", nrow(tbl), " matrices; ",
      sum(tbl$primary_anchor_weak), " needed the largest-|loading| anchor (",
      if (any(tbl$primary_anchor_weak)) paste(tbl$matrix[tbl$primary_anchor_weak], collapse = ", ") else "none",
      "); ", nrow(dis), " have contradictory level anchors (",
      if (nrow(dis)) paste(dis$matrix, collapse = ", ") else "none", "); ",
      sum(!tbl$sign_interpretable), " of ", nrow(tbl),
      " matrices therefore carry a sign that is a convention and must be read as a magnitude")
  invisible(tbl)
}
# Projection onto a stored rotation, centred on the discovery gene means and
# divided by d1, so a discovery sample's projection equals its cached u1. All
# rotation genes must be present, since dropping genes would move the direction.
project_on <- function(E, rot) {
  g <- names(rot$v1)
  stopifnot(all(g %in% colnames(E)))
  Xc <- sweep(E[, g, drop = FALSE], 2, rot$mu[g], "-")
  as.numeric(Xc %*% rot$v1) / rot$d1
}
# Share of a set's own (column-centred) variance that lies along a fixed direction.
share_along <- function(E, rot) {
  g <- names(rot$v1)
  Xc <- scale(E[, g, drop = FALSE], center = TRUE, scale = FALSE)
  sum((Xc %*% rot$v1)^2) / sum(Xc^2)
}
# PC1 variance share from the eigenvalues of the n x n Gram matrix. Identical to
# the squared-singular-value share of pca_lead() and far cheaper when genes >>
# samples, which is what the subsampling below needs.
pc1_share <- function(E) {
  Xc <- scale(E, center = TRUE, scale = FALSE)
  ev <- eigen(tcrossprod(Xc), symmetric = TRUE, only.values = TRUE)$values
  ev[1] / sum(ev)
}
# PC1 share rises as n falls, so a set's own share is compared with the
# discovery matrix subsampled to the same n: median and range over N_PC1_DRAWS
# draws, memoised by (matrix, n). The share varies by about 0.02 (SD) between
# subsets at n = 72, so 200 draws keep the median stable to the third decimal.
N_PC1_DRAWS <- 200
.pc1_memo <- new.env(parent = emptyenv())
NM_NA <- c(median = NA_real_, lo = NA_real_, hi = NA_real_, draws = NA_real_)
n_matched_pc1 <- function(E_disc, n_target, key) {
  if (!is.finite(n_target) || n_target < 10 || n_target >= nrow(E_disc)) return(NM_NA)
  k <- paste(key, n_target, sep = "_")
  if (!is.null(.pc1_memo[[k]])) return(.pc1_memo[[k]])
  set.seed(SEED + n_target)
  v <- vapply(seq_len(N_PC1_DRAWS), function(i)
    pc1_share(E_disc[sample.int(nrow(E_disc), n_target), , drop = FALSE]), numeric(1))
  out <- c(median = median(v), lo = min(v), hi = max(v), draws = N_PC1_DRAWS)
  assign(k, out, envir = .pc1_memo); out
}

# =============================================================================
# 1. The discovery axis, recomputed and checked against 07's cache
# =============================================================================
banner("1 | Discovery PC1 recomputed as 07 and checked against the cache")
rot_lnc <- pca_lead(E_lnc_disc)
rot_pc  <- pca_lead(E_pc_disc)
sd_l <- sd(rot_lnc$u1); sd_p <- sd(rot_pc$u1)   # projections are per discovery SD

ax_cached <- axis$lnc_axis[match(names(rot_lnc$u1), axis$sample_barcode)]
r_axis    <- stats::cor(rot_lnc$u1, ax_cached, use = "complete.obs")
maxdiff   <- max(abs(rot_lnc$u1 - ax_cached), na.rm = TRUE)
self_l    <- max(abs(project_on(E_lnc_disc, rot_lnc) - rot_lnc$u1))
self_p    <- max(abs(project_on(E_pc_disc,  rot_pc)  - rot_pc$u1))
msg("Recomputed lncRNA PC1 versus cache/lnc_global_axis.rds: r = ", sprintf("%.6f", r_axis),
    ", max |difference| = ", signif(maxdiff, 3), " (n = ", sum(!is.na(ax_cached)), ")")
stopifnot(r_axis > 0.9999)

pc07   <- fread(file.path(RESULTS_DIR, "07_pc_variance_explained.tsv"))
q_disc <- function(E) cohort[match(rownames(E), sample_barcode)]
rec_row <- function(label, E, rot, self_diff, r_cache, md, n_cache, lab07) {
  q <- q_disc(E); s <- spearman(rot$u1, q$pct_noFeature)
  data.table(matrix = label, n_samples = nrow(E), n_genes = ncol(E),
             pc1_var_share = round(rot$var_share[1], 4),
             pc1_var_share_07 = pc07[matrix == lab07 & pc == 1, var_explained],
             rho_noFeature = round(s[["rho"]], 3),
             rho_noFeature_07 = pc07[matrix == lab07 & pc == 1, rho_noFeature],
             n_cached_axis = n_cache,
             pearson_r_vs_cached_axis = if (is.na(r_cache)) NA_real_ else round(r_cache, 6),
             max_abs_diff_vs_cached_axis = if (is.na(md)) NA_real_ else signif(md, 3),
             self_projection_max_abs_diff = signif(self_diff, 3),
             rotation = "kept: v1 (gene loadings), discovery gene means, d1")
}
recheck <- rbind(
  rec_row("lncRNA", E_lnc_disc, rot_lnc, self_l, r_axis, maxdiff, sum(!is.na(ax_cached)), "lncRNA_observed"),
  rec_row("protein_coding", E_pc_disc, rot_pc, self_p, NA_real_, NA_real_, NA_integer_, "protein_coding_observed"))
recheck[, check := "recomputed discovery PC1 versus cache/lnc_global_axis.rds and 07_pc_variance_explained.tsv"]
recheck[, passed := (matrix != "lncRNA" | (pearson_r_vs_cached_axis >= 0.999999 &
                                           max_abs_diff_vs_cached_axis == 0)) &
                    pc1_var_share == pc1_var_share_07 &
                    rho_noFeature == rho_noFeature_07 &
                    self_projection_max_abs_diff < 1e-10]
recheck[, note := paste("this is the same check 25_axis_recompute_check.tsv performs, on the same",
                        "rotation; the rows of 35_axis_projection_subtypes.tsv are therefore",
                        "directly comparable with those of 25_normal_axis_projection.tsv")]
save_tsv(recheck, "35_axis_recompute_check.tsv"); print(recheck, row.names = FALSE)
if (!all(recheck$passed))
  stop("35: the discovery PC1 recomputed here does not reproduce 07's cached axis or table")

# =============================================================================
# 2. The subtype cohorts
# =============================================================================
banner("2 | TCGA-KIRP and TCGA-KICH: cached expression and STAR metrics")
load_subtype <- function(proj) {
  tag <- sub("TCGA-", "", proj)
  st  <- readRDS(file.path(CACHE_DIR, paste0("subtype_", tag, ".rds")))
  co  <- as.data.table(st$cohort); q <- as.data.table(st$qc)
  # Expression and STAR metrics must come from the same library. Each subtype barcode
  # has one primary-tumour STAR-Counts file, so barcode indexing is unambiguous. In
  # TCGA-KIRC it is not, so 01, 07 and 23 index by file_id.
  stopifnot(!anyDuplicated(colnames(st$fpkm)), !anyDuplicated(co$sample_barcode))
  bc <- co$sample_barcode
  stopifnot(all(bc %in% colnames(st$fpkm)))
  q  <- q[match(bc, sample_barcode)]
  stopifnot(identical(as.character(q$sample_barcode), as.character(bc)))
  need <- c("pct_noFeature", "pct_multimapping", "assigned_reads")
  miss <- need[!need %in% names(q)]
  if (length(miss))
    stop(proj, ": cache/subtype_", tag, ".rds carries no ", paste(miss, collapse = ", "),
         "; recompute the STAR metrics with read_star_summary() over the local ",
         "STAR-Counts files as 11_subtype_specificity.R does, then re-run")
  stopifnot(all(is.finite(q$pct_noFeature)), all(is.finite(q$pct_multimapping)),
            all(is.finite(q$assigned_reads)), all(q$assigned_reads > 0))
  miss_l <- setdiff(genes_lnc, rownames(st$fpkm)); miss_p <- setdiff(genes_pc, rownames(st$fpkm))
  msg(proj, ": ", length(bc), " primary-tumour libraries; production genes present -- lncRNA ",
      length(genes_lnc) - length(miss_l), "/", length(genes_lnc), ", protein-coding ",
      length(genes_pc) - length(miss_p), "/", length(genes_pc))
  if (length(miss_l) || length(miss_p))
    stop(proj, ": ", length(miss_l) + length(miss_p), " production genes are absent from the ",
         "cached FPKM matrix. The discovery rotation cannot be applied to a reduced gene set ",
         "without moving the direction and breaking comparability with 25")
  # The library-failure rule of 01 is reported, not applied, as in 11 and 25.
  fail_reads <- sum(q$assigned_reads < MIN_ASSIGNED_READS)
  fail_nf    <- sum(q$pct_noFeature > MAX_NOFEATURE_PCT)
  msg("  library-failure rule as written (reported, not applied): reads < ",
      MIN_ASSIGNED_READS / 1e6, "M n = ", fail_reads, "; non-feature > ",
      MAX_NOFEATURE_PCT, "% n = ", fail_nf)
  msg("  non-feature fraction: median ", round(median(q$pct_noFeature), 2), "%, range ",
      round(min(q$pct_noFeature), 2), " to ", round(max(q$pct_noFeature), 2), "%")
  list(cohort = proj, n = length(bc),
       lnc = t(log2(st$fpkm[genes_lnc, bc, drop = FALSE] + 1)),
       pc  = t(log2(st$fpkm[genes_pc,  bc, drop = FALSE] + 1)),
       q = q, n_fail_reads = fail_reads, n_fail_noFeature = fail_nf)
}
subs <- lapply(SUBTYPES, load_subtype); names(subs) <- SUBTYPES

# =============================================================================
# 3. Projection onto the discovery rotations, and each cohort's own PC1
# =============================================================================
banner("3 | Projections onto the discovery axes; each cohort's own leading component")
build_set <- function(S) {
  E_l <- S$lnc; E_p <- S$pc
  # Projections first: they are the third sign anchor for the own PC1.
  proj_l <- project_on(E_l, rot_lnc) / sd_l
  proj_p <- project_on(E_p, rot_pc)  / sd_p
  own_l <- if (nrow(E_l) >= 10) pca_lead(E_l, ref = proj_l) else NULL
  own_p <- if (nrow(E_p) >= 10) pca_lead(E_p, ref = proj_p) else NULL
  d <- data.table(
    cohort = S$cohort, tissue = "tumour",
    lnc_axis_proj = proj_l,
    pc_axis_proj  = proj_p,
    mean_lnc_expr = rowMeans(E_l), mean_pc_expr = rowMeans(E_p),
    lnc_own_pc1_z = if (is.null(own_l)) NA_real_ else as.numeric(scale(own_l$score)),
    pc_own_pc1_z  = if (is.null(own_p)) NA_real_ else as.numeric(scale(own_p$score)),
    pct_noFeature    = as.numeric(S$q$pct_noFeature),
    pct_multimapping = as.numeric(S$q$pct_multimapping),
    log_depth        = log10(as.numeric(S$q$assigned_reads)))
  list(tbl = d, own = list(lnc = own_l, pc = own_p),
       share_along = list(lnc = share_along(E_l, rot_lnc), pc = share_along(E_p, rot_pc)),
       n_genes = c(lnc = ncol(E_l), pc = ncol(E_p)))
}
sets <- lapply(subs, build_set)

banner("3b | PC1 sign convention, every matrix")
sign_objs <- c(list(`discovery lncRNA (07 rotation)` = rot_lnc,
                    `discovery protein_coding (07 rotation)` = rot_pc),
               setNames(lapply(SUBTYPES, function(p) sets[[p]]$own$lnc),
                        paste(SUBTYPES, "tumour own lncRNA PC1")),
               setNames(lapply(SUBTYPES, function(p) sets[[p]]$own$pc),
                        paste(SUBTYPES, "tumour own protein_coding PC1")))
sign_tbl <- assert_sign_convention(sign_objs)

# ---- per-cohort axis summaries, in the schema of 25 ----
proj_rows <- rbindlist(lapply(SUBTYPES, function(nm) {
  S <- sets[[nm]]; d <- S$tbl
  covs <- list(noFeature = d$pct_noFeature, multimap = d$pct_multimapping,
               log_depth = d$log_depth)
  mk <- function(matrix, axis_type, x, n_genes, var_share, mean_expr, disc_E,
                 own_share = rep(NA_real_, 5), r_own = NA_real_, n_own = NA_integer_,
                 sgn = NULL) {
    ok <- is.finite(x)
    # Size-matched baseline for the own-PC1 share (see n_matched_pc1()).
    nmp <- if (axis_type == "own_PC1") n_matched_pc1(disc_E, sum(ok), matrix) else NM_NA
    out <- data.table(cohort = d$cohort[1], tissue = d$tissue[1], matrix = matrix,
                      axis_type = axis_type, n = sum(ok), n_genes = n_genes,
                      var_share = round(var_share, 4))
    for (v in names(covs)) { s <- spearman(x, covs[[v]])
      set(out, j = paste0("rho_", v), value = round(s[["rho"]], 3))
      set(out, j = paste0("p_", v),   value = signif(s[["p"]], 3)) }
    s <- spearman(x, mean_expr)
    # All three sign anchors are written out. A projected row inherits the
    # discovery rotation's orientation and carries its anchor values.
    out[, `:=`(pc1_sign_rule = if (is.null(sgn)) NA_character_ else
                 paste0(sgn$sign_rule, if (axis_type != "own_PC1") " (inherited from the discovery rotation)" else ""),
               pc1_sign_stat = if (is.null(sgn)) NA_real_ else round(sgn$sign_stat, 3),
               pc1_sign_r_mean_all_genes = if (is.null(sgn)) NA_real_ else round(sgn$sign_r_mean_all, 3),
               pc1_sign_r_mean_top_genes = if (is.null(sgn)) NA_real_ else round(sgn$sign_r_mean_top, 3),
               pc1_sign_r_projected_on_discovery = if (is.null(sgn)) NA_real_ else round(sgn$sign_r_projected, 3),
               pc1_sign_resolved_by = if (is.null(sgn)) NA_character_ else sgn$sign_resolved_by,
               pc1_sign_primary_anchor_weak = if (is.null(sgn)) NA else sgn$sign_primary_weak,
               pc1_sign_anchors_disagree = if (is.null(sgn)) NA else sgn$sign_anchors_disagree,
               pc1_sign_interpretable = if (is.null(sgn)) NA else sgn$sign_interpretable)]
    out[, `:=`(rho_mean_expr = round(s[["rho"]], 3), p_mean_expr = signif(s[["p"]], 3),
               r_projected_vs_own = round(r_own, 3),
               n_projected_vs_own = if (is.na(r_own)) NA_integer_ else as.integer(n_own),
               n_matched_pc1_share = round(nmp[["median"]], 4),
               n_matched_pc1_share_min = round(nmp[["lo"]], 4),
               n_matched_pc1_share_max = round(nmp[["hi"]], 4),
               n_matched_draws = nmp[["draws"]],
               own_pc2_share = round(own_share[2], 4), own_pc3_share = round(own_share[3], 4),
               own_pc4_share = round(own_share[4], 4), own_pc5_share = round(own_share[5], 4))]
    out
  }
  # r_projected_vs_own is a complete-case Pearson correlation with its own n.
  fin_l <- is.finite(d$lnc_axis_proj) & is.finite(d$lnc_own_pc1_z)
  fin_p <- is.finite(d$pc_axis_proj)  & is.finite(d$pc_own_pc1_z)
  r_l <- if (all(is.na(d$lnc_own_pc1_z))) NA_real_ else stats::cor(d$lnc_axis_proj, d$lnc_own_pc1_z, use = "complete.obs")
  r_p <- if (all(is.na(d$pc_own_pc1_z)))  NA_real_ else stats::cor(d$pc_axis_proj,  d$pc_own_pc1_z,  use = "complete.obs")
  vs_l <- if (is.null(S$own$lnc)) rep(NA_real_, 5) else S$own$lnc$var_share[1:5]
  vs_p <- if (is.null(S$own$pc))  rep(NA_real_, 5) else S$own$pc$var_share[1:5]
  rbind(mk("lncRNA", "projected_discovery_PC1", d$lnc_axis_proj, S$n_genes[["lnc"]], S$share_along$lnc,
           d$mean_lnc_expr, E_lnc_disc, r_own = r_l, n_own = sum(fin_l), sgn = rot_lnc),
        mk("lncRNA", "own_PC1", d$lnc_own_pc1_z, S$n_genes[["lnc"]], vs_l[1],
           d$mean_lnc_expr, E_lnc_disc, vs_l, r_l, sum(fin_l), sgn = S$own$lnc),
        mk("protein_coding", "projected_discovery_PC1", d$pc_axis_proj, S$n_genes[["pc"]], S$share_along$pc,
           d$mean_pc_expr, E_pc_disc, r_own = r_p, n_own = sum(fin_p), sgn = rot_pc),
        mk("protein_coding", "own_PC1", d$pc_own_pc1_z, S$n_genes[["pc"]], vs_p[1],
           d$mean_pc_expr, E_pc_disc, vs_p, r_p, sum(fin_p), sgn = S$own$pc))
}))
set.seed(SEED)   # the subsampling above reseeded per target size

proj_rows[, var_share_definition := fifelse(axis_type == "own_PC1",
  "PC1 share of the set's own centred matrix",
  "share of the set's own centred variance along the discovery PC1 direction")]
proj_rows[, n_matched_pc1_note := fifelse(
  axis_type != "own_PC1", "not applicable: this row is a projection onto a fixed direction, not a PC1 share",
  fifelse(is.na(n_matched_pc1_share),
          "not applicable: n is at or above the discovery matrix, so this row IS (or exceeds) the discovery PCA",
          paste0("median [min-max] PC1 share of ", n_matched_draws,
                 " random subsets of this matrix's discovery samples drawn at this row's n")))]
# Per-row sign note: the anchors consulted, their values in this row's
# orientation and whether they agreed.
anchor_list <- function(all_g, top_g, proj) paste0(
  "anchors in this row's orientation: all-gene mean expression r = ", all_g,
  "; mean expression over the ", SIGN_TOP_GENES, " largest |loadings| r = ", top_g,
  fifelse(is.na(proj), "; projection on the discovery PC1 not applicable to this row",
          paste0("; projection on the discovery PC1 r = ", proj)))
proj_rows[, pc1_sign_note := fifelse(
  is.na(pc1_sign_primary_anchor_weak), NA_character_,
  fifelse(pc1_sign_anchors_disagree & grepl("^projection", pc1_sign_rule),
    paste0("SIGN IS A CONVENTION, NOT A FINDING: the two level anchors DISAGREE. ",
           anchor_list(pc1_sign_r_mean_all_genes, pc1_sign_r_mean_top_genes,
                       pc1_sign_r_projected_on_discovery),
           ". The all-gene and largest-|loading| anchors each clear |", SIGN_MIN_R,
           "| and point opposite ways, which is what a contrast component does when its few hundred ",
           "largest loadings oppose the bulk of the matrix. The printed orientation is the one the ",
           "third anchor gives -- the projection onto the discovery PC1, whose own sign is fixed ",
           "outside this matrix and which here agrees with the largest-|loading| anchor -- so that ",
           "the sign shown is supported by two anchors of three rather than by the weakest one. It ",
           "is still a convention: quote |rho| for this row and do not report the direction"),
  fifelse(pc1_sign_anchors_disagree,
    paste0("SIGN IS NOT DETERMINED: the two level anchors DISAGREE and the projection anchor does ",
           "not settle it. ", anchor_list(pc1_sign_r_mean_all_genes, pc1_sign_r_mean_top_genes,
                                          pc1_sign_r_projected_on_discovery),
           ". The printed orientation follows the largest-|loading| anchor only. Quote |rho| for ",
           "this row and do not report the direction"),
  fifelse(pc1_sign_primary_anchor_weak,
    paste0("SIGN IS A CONVENTION, NOT A FINDING: 07's all-gene mean-expression anchor is below |",
           SIGN_MIN_R, "| (its rank version also disagrees in sign -- see rho_mean_expr on this row) ",
           "because this PC1 is a contrast rather than a level axis. ",
           anchor_list(pc1_sign_r_mean_all_genes, pc1_sign_r_mean_top_genes,
                       pc1_sign_r_projected_on_discovery),
           ". Fixed instead by ", pc1_sign_resolved_by, " (r = ", pc1_sign_stat,
           "). Quote |rho| for this row and do not report the direction"),
    paste0("sign fixed by 07's convention, positive correlation with all-gene mean expression; every ",
           "decisive anchor agrees. ", anchor_list(pc1_sign_r_mean_all_genes, pc1_sign_r_mean_top_genes,
                                                   pc1_sign_r_projected_on_discovery))))))]
# A projected row inherits the discovery rotation's orientation and the discovery
# matrix's anchor values.
proj_rows[axis_type != "own_PC1" & !is.na(pc1_sign_note),
          pc1_sign_note := paste0("orientation inherited from the discovery rotation, which is fixed ",
                                  "once on the discovery matrix and is the same for every sample set; ",
                                  "the anchor values below are the discovery matrix's. ", pc1_sign_note)]
# Multiplicity, as in 25: one BH family per cohort x tissue set, 4 axes x
# 4 covariates = 16 tests.
PCOLS <- c("p_noFeature", "p_multimap", "p_log_depth", "p_mean_expr")
FCOLS <- sub("^p_", "fdr_", PCOLS)
proj_rows[, (FCOLS) := {
  nn <- .N; a <- bh(unlist(.SD, use.names = FALSE))
  lapply(seq_along(.SD), function(j) signif(a[((j - 1L) * nn + 1L):(j * nn)], 3))
}, by = .(cohort, tissue), .SDcols = PCOLS]
proj_rows[, fdr_family := paste0("Benjamini-Hochberg within cohort x tissue over the ",
                                 length(PCOLS) * 4, " reported Spearman tests of that set (4 axes x ",
                                 length(PCOLS), " covariates)")]
# Provenance is logged rather than stored as a column, so the column set
# matches 25.
for (p in SUBTYPES)
  msg("Provenance of the ", p, " rows: TCGA primary tumours, poly(A) libraries, expression and ",
      "STAR metrics from cache/subtype_", sub("TCGA-", "", p), ".rds; computed by the procedure of ",
      "25_matched_normal_control.R on the same rotation, gene sets, scale, centring and sign ",
      "convention, so these rows are directly comparable with, and row-bindable onto, ",
      "25_normal_axis_projection.tsv")

setcolorder(proj_rows, c("cohort", "tissue", "matrix", "axis_type", "n", "n_genes", "var_share",
                         "n_matched_pc1_share", "n_matched_pc1_share_min", "n_matched_pc1_share_max",
                         "n_matched_draws",
                         "rho_noFeature", "p_noFeature", "fdr_noFeature",
                         "rho_multimap", "p_multimap", "fdr_multimap",
                         "rho_log_depth", "p_log_depth", "fdr_log_depth",
                         "rho_mean_expr", "p_mean_expr", "fdr_mean_expr",
                         "r_projected_vs_own", "n_projected_vs_own",
                         "pc1_sign_rule", "pc1_sign_stat", "pc1_sign_r_mean_all_genes",
                         "pc1_sign_r_mean_top_genes", "pc1_sign_r_projected_on_discovery",
                         "pc1_sign_resolved_by", "pc1_sign_primary_anchor_weak",
                         "pc1_sign_anchors_disagree", "pc1_sign_interpretable"))
save_tsv(proj_rows, "35_axis_projection_subtypes.tsv")
print(proj_rows[, .(cohort, matrix, axis_type, n, var_share, n_matched_pc1_share,
                    rho_noFeature, fdr_noFeature, rho_multimap, rho_log_depth, rho_mean_expr,
                    r_projected_vs_own, n_projected_vs_own, pc1_sign_r_mean_all_genes,
                    pc1_sign_r_mean_top_genes, pc1_sign_r_projected_on_discovery,
                    pc1_sign_primary_anchor_weak, pc1_sign_anchors_disagree,
                    pc1_sign_interpretable)],
      row.names = FALSE)

# Log every own-PC1 row whose sign is a convention.
wk <- proj_rows[axis_type == "own_PC1" & pc1_sign_interpretable == FALSE]
if (nrow(wk)) for (i in seq_len(nrow(wk)))
  msg("SIGN IS A CONVENTION: ", wk$cohort[i], " ", wk$tissue[i], " ", wk$matrix[i],
      " own PC1 -- all-gene anchor ", wk$pc1_sign_r_mean_all_genes[i], ", ",
      SIGN_TOP_GENES, "-largest-|loading| anchor ", wk$pc1_sign_r_mean_top_genes[i],
      ", projection anchor ", wk$pc1_sign_r_projected_on_discovery[i],
      " (", if (isTRUE(wk$pc1_sign_anchors_disagree[i])) "the two level anchors DISAGREE" else
        "the primary anchor is below the threshold", "; sign fixed by ", wk$pc1_sign_resolved_by[i],
      "). Its rho with the non-feature fraction (", wk$rho_noFeature[i], ", fdr ", wk$fdr_noFeature[i],
      ") must be quoted as a magnitude: the direction is a convention, not a finding.")
n_conv <- proj_rows[axis_type == "own_PC1" & pc1_sign_interpretable == FALSE, .N]
msg("Own-PC1 rows whose sign is a convention: ", n_conv, " of ",
    proj_rows[axis_type == "own_PC1", .N], "; projected rows: 0 of ",
    proj_rows[axis_type != "own_PC1", .N],
    " (a projected row inherits the discovery rotation's fixed orientation)")

# =============================================================================
# 4. Projected versus own leading component, beside 11 and 25
# =============================================================================
banner("4 | Projected versus own leading component, beside 11 and 25")
cmp <- dcast(proj_rows[, .(cohort, matrix, axis_type, rho_noFeature, n)],
             cohort + matrix + n ~ axis_type, value.var = "rho_noFeature")
setnames(cmp, c("projected_discovery_PC1", "own_PC1"),
         c("rho_noFeature_projected", "rho_noFeature_own"))
cmp[, difference := round(rho_noFeature_projected - rho_noFeature_own, 3)]
print(cmp, row.names = FALSE)

f11 <- file.path(RESULTS_DIR, "11_axis_replication_subtypes.tsv")
if (file.exists(f11)) {
  old <- fread(f11)[, .(cohort, matrix, n_11 = n, rho_noFeature_11 = rho_noFeature)]
  chk <- merge(cmp, old, by = c("cohort", "matrix"), all.x = TRUE)
  chk[, own_reproduces_11 := is.finite(rho_noFeature_11) &
        abs(rho_noFeature_own - rho_noFeature_11) < 0.02 & n == n_11]
  msg("Own-PC1 rows reproduce 11_axis_replication_subtypes.tsv: ",
      sum(chk$own_reproduces_11), "/", nrow(chk))
  print(chk[, .(cohort, matrix, n, n_11, rho_noFeature_own, rho_noFeature_11,
                rho_noFeature_projected, own_reproduces_11)], row.names = FALSE)
} else msg("11_axis_replication_subtypes.tsv absent; the own-PC1 cross-check is skipped")

# =============================================================================
# 5. Row-bind onto 25, checked on the files as written
# =============================================================================
# The two TSVs are read back from disk and bound with fill = FALSE, so any
# difference in column sets is an error.
banner("5 | Row-binding onto 25_normal_axis_projection.tsv")
f25 <- file.path(RESULTS_DIR, "25_normal_axis_projection.tsv")
f35 <- file.path(RESULTS_DIR, "35_axis_projection_subtypes.tsv")
if (file.exists(f25)) {
  p25 <- fread(f25); p35 <- fread(f35)
  only25 <- setdiff(names(p25), names(p35)); only35 <- setdiff(names(p35), names(p25))
  msg("Columns: 25 has ", ncol(p25), ", 35 has ", ncol(p35),
      "; in 25 only: ", if (length(only25)) paste(only25, collapse = ", ") else "none",
      "; in 35 only: ", if (length(only35)) paste(only35, collapse = ", ") else "none",
      "; same order: ", identical(names(p25), names(p35)))
  if (length(only25) || length(only35))
    stop("35: the two files no longer share a column set, so they cannot be row-bound without ",
         "filling a column by default. In 25 only: ", paste(only25, collapse = ", "),
         "; in 35 only: ", paste(only35, collapse = ", "))
  if (!identical(names(p25), names(p35)))
    stop("35: the two files share their columns but not their order; fix setcolorder here")
  both <- rbind(p25, p35, fill = FALSE)
  stopifnot(nrow(both) == nrow(p25) + nrow(p35), ncol(both) == ncol(p25))
  # No column may be empty on one side only. n_matched_pc1_* and own_pc*_share
  # are empty on projected rows by definition, so for those columns only
  # own_PC1 rows are compared.
  side_blank <- function(D, nm) {
    v <- D[[nm]]
    if (is.character(v)) all(is.na(v) | !nzchar(v)) else all(is.na(v))
  }
  own_only <- c("n_matched_pc1_share", "n_matched_pc1_share_min", "n_matched_pc1_share_max",
                "n_matched_draws", "own_pc2_share", "own_pc3_share", "own_pc4_share",
                "own_pc5_share")
  blank <- unlist(lapply(names(both), function(nm) {
    rows25 <- if (nm %in% own_only) p25[axis_type == "own_PC1"] else p25
    rows35 <- if (nm %in% own_only) p35[axis_type == "own_PC1"] else p35
    if (side_blank(rows25, nm) != side_blank(rows35, nm)) nm else NULL
  }))
  msg("Row-bind of the two files with fill = FALSE: ", nrow(p25), " + ", nrow(p35), " = ",
      nrow(both), " rows x ", ncol(both), " columns; columns empty on one side only: ",
      if (length(blank)) paste(blank, collapse = ", ") else "none")
  if (length(blank))
    stop("35: these columns are populated in one file and empty in the other, which is a ",
         "placeholder in all but name: ", paste(blank, collapse = ", "))
  msg("Sample sets in the bound table: ", both[, uniqueN(paste(cohort, tissue))],
      " (", paste(unique(both[, paste(cohort, tissue)]), collapse = "; "), ")")
  msg("Six-set comparison of the PROJECTED lncRNA axis against the non-feature fraction:")
  six <- both[matrix == "lncRNA" & axis_type == "projected_discovery_PC1",
              .(cohort, tissue, n, rho_noFeature, fdr_noFeature,
                sign_interpretable = pc1_sign_interpretable)][order(-abs(rho_noFeature))]
  print(six, row.names = FALSE)
  msg("The same six sets on each set's OWN leading component, for contrast (quote a row whose ",
      "sign_interpretable is FALSE as a magnitude only):")
  print(both[matrix == "lncRNA" & axis_type == "own_PC1",
             .(cohort, tissue, n, var_share, n_matched_pc1_share, rho_noFeature,
               sign_interpretable = pc1_sign_interpretable)], row.names = FALSE)
} else msg("25_normal_axis_projection.tsv absent; the row-bind check is skipped")

msg("Elapsed: ", round(difftime(Sys.time(), t_start, units = "mins"), 2), " min")
write_session_info("35")
banner("35 | done")
