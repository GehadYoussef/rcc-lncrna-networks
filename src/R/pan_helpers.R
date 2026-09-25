# pan_helpers.R: helpers shared by the pan-cancer stages 39 to 41
# Every project is filtered, decomposed and adjusted by the same functions.
# Sourced after 00_config.R.

# Stops if the decision rules differ from those locked before any project was
# analysed (results/38_decision_rules_lock.json).
verify_decision_lock <- function() {
  lk <- jsonlite::fromJSON(file.path(RESULTS_DIR, "38_decision_rules_lock.json"))
  now <- list(PAN_MIN_TUMOURS = PAN_MIN_TUMOURS, PAN_MIN_NORMALS = PAN_MIN_NORMALS,
              PAN_EXCLUDE_INFERENCE = PAN_EXCLUDE_INFERENCE,
              PAN_AXIS_RHO_MIN = PAN_AXIS_RHO_MIN, PAN_AXIS_MARGIN = PAN_AXIS_MARGIN,
              PAN_NET_LARGEST_MIN = PAN_NET_LARGEST_MIN, PAN_NET_EIGEN_R_MIN = PAN_NET_EIGEN_R_MIN,
              PAN_NET_ARI_MAX = PAN_NET_ARI_MAX, PAN_SURV_MIN_EVENTS = PAN_SURV_MIN_EVENTS,
              PAN_STAGE_MIN_FRAC = PAN_STAGE_MIN_FRAC, PAN_PLATE_MIN_N = PAN_PLATE_MIN_N)
  for (k in names(now))
    if (!isTRUE(all.equal(now[[k]], lk$rules[[k]])))
      stop("decision rule ", k, " differs from the locked value; report the change explicitly")
  msg("decision rules match the lock of ", lk$locked_at, " (commit ", substr(lk$git_commit, 1, 7), ")")
  invisible(lk)
}

pan_projects <- function() {
  f <- list.files(CACHE_DIR, pattern = "^pan_TCGA-[A-Z]+\\.rds$")
  sort(sub("^pan_(.*)\\.rds$", "\\1", f))
}

# Retained libraries of one group ("tumour" or "normal") and their log2(FPKM + 1)
# matrices after the discovery gene filters of 01 (split_and_filter), applied
# within the group. Returns samples x genes matrices, as WGCNA expects.
pan_group <- function(obj, grp) {
  # grp, not group: inside [.data.table a bare `group` would be the column
  s <- obj$samples[keep == TRUE & group == grp]
  fp <- obj$fpkm[, s$file_id, drop = FALSE]
  filt <- function(biotype, min_fpkm, min_frac, top_n) {
    ann <- obj$ann[gene_type == biotype]
    m <- fp[ann$gene_id, , drop = FALSE]
    ex <- rowMeans(m >= min_fpkm, na.rm = TRUE) >= min_frac
    m <- m[ex, , drop = FALSE]; ann <- ann[ex]
    lg <- log2(m + 1)
    if (nrow(lg) > top_n) {
      sel <- order(matrixStats::rowMads(lg), decreasing = TRUE)[seq_len(top_n)]
      lg <- lg[sel, , drop = FALSE]; ann <- ann[sel]
    }
    ann <- copy(ann)[, mad_v := matrixStats::rowMads(lg)]
    setorder(ann, gene_name, -mad_v)
    keep <- !duplicated(ann$gene_name) & !is.na(ann$gene_name) & nzchar(ann$gene_name)
    t(lg[ann$gene_id[keep], , drop = FALSE])
  }
  list(samples = s,
       pc  = filt("protein_coding", MRNA_MIN_FPKM, MRNA_MIN_FRAC, MRNA_TOP_N_MAD),
       lnc = filt("lncRNA", LNC_MIN_FPKM, LNC_MIN_FRAC, LNC_TOP_N_MAD))
}

# Leading principal component of a samples x genes matrix, column-centred and
# unscaled, sign-aligned to per-sample mean expression. This is the definition
# used for the discovery axis.
pc1 <- function(X) {
  Xc <- scale(X, center = TRUE, scale = FALSE)
  sv <- svd(Xc, nu = 1, nv = 0)
  sc <- sv$u[, 1] * sv$d[1]
  if (cor(sc, rowMeans(X)) < 0) sc <- -sc
  list(score = sc, var_share = sv$d[1]^2 / sum(Xc^2))
}

# A factor with levels carrying fewer than n samples pooled into "other".
pool_levels <- function(x, n = PAN_PLATE_MIN_N) {
  tb <- table(x); x <- as.character(x)
  x[x %in% names(tb)[tb < n] | is.na(x)] <- "other"
  factor(x)
}
r2_factor <- function(y, f) {
  ok <- is.finite(y) & !is.na(f); f <- droplevels(factor(f[ok]))
  if (nlevels(f) < 2) return(NA_real_)
  summary(lm(y[ok] ~ f))$r.squared
}
