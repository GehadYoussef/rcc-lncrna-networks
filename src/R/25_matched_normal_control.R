# 25_matched_normal_control.R: matched normal kidney libraries as a control
#
# Normal kidney libraries carry no tumour biology but share collection and
# sequencing with the tumours. Tests whether the discovery lncRNA axis tracks
# the non-feature read fraction in normal tissue, and whether a patient's
# tumour and normal libraries share that fraction and axis.
#   1. Read the local Solid Tissue Normal STAR-Counts files (TCGA-KIRC and
#      CPTAC-3). TCGA-KIRP and TCGA-KICH availability is checked on the GDC.
#   2. Recompute the discovery PC1 as in 07, project every library onto it, and
#      run each set's own PCA.
#   3. Paired tumour-normal tests, and variance of the metrics by patient,
#      tissue, tissue source site and plate (tumours, normals and pooled).
#   4. Overall survival against the normal and tumour library metric (paired
#      TCGA-KIRC patients, small sample).
#   5. Unpaired tumour versus normal distributions, with matched-pair and
#      plate-restricted versions of each contrast.
# Inputs: cache/dataset.rds, networks.rds, lnc_global_axis.rds, gdc_file_map.rds,
# valid_file_map.rds, validation_dataset.rds, biospecimen_kirc.rds (23).
# Outputs: results/25_*.tsv and caches 25_normal_libraries.rds,
# 25_normal_aliquots.rds and 25_subtype_normal_files.rds.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(httr); library(jsonlite)
})
banner("25 | Matched normal libraries as a tumour-biology-free control")
set.seed(SEED)
t_start <- Sys.time()

ds       <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets     <- readRDS(file.path(CACHE_DIR, "networks.rds"))
axis     <- as.data.table(readRDS(file.path(CACHE_DIR, "lnc_global_axis.rds")))
cohort   <- as.data.table(ds$cohort_full)
fm_kirc  <- as.data.table(readRDS(file.path(CACHE_DIR, "gdc_file_map.rds")))
fm_cptac <- as.data.table(readRDS(file.path(CACHE_DIR, "valid_file_map.rds")))
vd       <- readRDS(file.path(CACHE_DIR, "validation_dataset.rds"))
vcohort  <- as.data.table(vd$cohort); vqc <- as.data.table(vd$qc)
VALID_DIR   <- file.path(DOWNLOAD_DIR, "validation_data", "CPTAC-3")
NORMAL_TYPE <- "Solid Tissue Normal"
SUBTYPES    <- c("TCGA-KIRP", "TCGA-KICH")
if (!"tss" %in% names(cohort)) cohort[, tss := tstrsplit(sample_barcode, "-", keep = 2)[[1]]]
if (!"plate" %in% names(cohort)) cohort[, plate := NA_character_]
if (!"assigned_reads" %in% names(cohort)) cohort[, assigned_reads := libsize]

# The matrices 07 decomposed: observed log2(FPKM + 1) on the network samples.
E_lnc_disc <- obs_expr(nets$lnc)
E_pc_disc  <- obs_expr(nets$mrna)
genes_lnc  <- colnames(E_lnc_disc); genes_pc <- colnames(E_pc_disc)
stopifnot(all(genes_lnc %in% rownames(ds$lnc$expr)), all(genes_pc %in% rownames(ds$mrna$expr)),
          all(genes_lnc %in% rownames(vd$fpkm)), all(genes_pc %in% rownames(vd$fpkm)))
msg("Discovery matrices: lncRNA ", nrow(E_lnc_disc), " x ", ncol(E_lnc_disc),
    "; protein-coding ", nrow(E_pc_disc), " x ", ncol(E_pc_disc))

# ---- helpers ----------------------------------------------------------------
chr1  <- function(x) if (is.null(x)) NA_character_ else as.character(x[[1]])
modal <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA_character_
                       else names(sort(table(x), decreasing = TRUE))[1] }
bc_field <- function(x, k) vapply(strsplit(as.character(x), "-", fixed = TRUE), function(z)
  if (length(z) >= k) z[k] else NA_character_, character(1))
pool_levels <- function(x, min_n = 10) {
  x <- as.character(x); tb <- table(x)
  x[x %in% names(tb)[tb < min_n]] <- "other"
  factor(x)
}
spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(c(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}
qtl <- function(x) { x <- x[is.finite(x)]
  c(median = median(x), q25 = unname(quantile(x, 0.25)), q75 = unname(quantile(x, 0.75))) }
# Benjamini-Hochberg over finite p-values only, since p.adjust() would count NA rows.
bh <- function(p) { out <- rep(NA_real_, length(p)); ok <- is.finite(p)
  if (any(ok)) out[ok] <- p.adjust(p[ok], "BH"); out }

# Leading principal components of a samples x genes matrix, as in 07:
# column-centred, unscaled SVD. u1 is the unit-norm score 07 caches as
# lnc_axis, v1 the rotation, d1 the singular value.
#
# Sign convention. The sign of PC1 is arbitrary and is fixed by three anchors,
# all written to the results table:
#   1. correlation with mean expression over all genes (07's rule)
#   2. correlation with mean expression over the SIGN_TOP_GENES genes with the
#      largest |loadings|
#   3. correlation with the projection onto the discovery PC1 of the same
#      biotype (own PC1s only, undefined for the discovery rotations)
# Anchor 1 is used when decisive (|r| >= SIGN_MIN_R) and not contradicted by
# anchor 2. Anchor 2 is used when anchor 1 is not decisive. Anchor 3 is used
# when anchors 1 and 2 are both decisive and disagree, since its orientation is
# fixed outside the matrix. A row with a weak anchor 1 or disagreeing anchors
# has pc1_sign_interpretable = FALSE and its rho values are read as magnitudes.
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
# Prints the anchors of every PC1 and stops if a matrix has no decisive anchor.
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
# divided by d1, so a discovery sample's projection equals its cached u1.
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
# PC1 variance share of a samples x genes matrix, from the eigenvalues of its
# n x n Gram matrix. Identical to the squared-singular-value share of pca_lead()
# and far cheaper when genes >> samples, as in the subsampling below.
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

# GDC POST returning the hits. Callers wrap it in tryCatch and cache the result.
gdc_raw <- function(endpoint, body) {
  r <- httr::POST(paste0("https://api.gdc.cancer.gov/", endpoint),
                  body = jsonlite::toJSON(body, auto_unbox = TRUE),
                  httr::content_type_json(), httr::timeout(300))
  httr::stop_for_status(r)
  j <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"),
                          simplifyVector = FALSE)
  if (!is.null(j$data$pagination$total) &&
      j$data$pagination$total > length(j$data$hits))
    warning("GDC returned ", length(j$data$hits), " of ",
            j$data$pagination$total, " records from /", endpoint)
  j$data$hits
}
gdc_notes <- character(0)

# =============================================================================
# 1. Normal-tissue libraries
# =============================================================================
banner("1 | Normal libraries: file lists, STAR metrics, expression, library rule")

# ---- 1a. TCGA-KIRP / TCGA-KICH normals: are any present locally? -----------
fetch_normal_files <- function(proj) {
  filt <- list(op = "and", content = list(
    list(op = "in", content = list(field = "cases.project.project_id", value = list(proj))),
    list(op = "in", content = list(field = "data_category", value = list("Transcriptome Profiling"))),
    list(op = "in", content = list(field = "data_type", value = list("Gene Expression Quantification"))),
    list(op = "in", content = list(field = "analysis.workflow_type", value = list("STAR - Counts"))),
    list(op = "in", content = list(field = "cases.samples.sample_type", value = list(NORMAL_TYPE)))))
  hits <- gdc_raw("files", list(filters = filt, format = "JSON", size = "2000",
                                fields = paste("file_id", "file_name", "cases.submitter_id",
                                               "cases.samples.submitter_id",
                                               "cases.samples.sample_type", sep = ",")))
  empty <- data.table(project = character(0), file_id = character(0), file_name = character(0),
                      patient = character(0), sample_barcode = character(0), sample_type = character(0))
  if (!length(hits)) return(empty)
  rbind(empty, rbindlist(lapply(hits, function(h) {
    cs <- if (length(h$cases)) h$cases[[1]] else NULL
    sm <- NULL
    if (!is.null(cs) && length(cs$samples)) {
      # the case's sample list may carry the tumour as well: take the normal
      types <- vapply(cs$samples, function(s) chr1(s$sample_type), character(1))
      sm <- cs$samples[[if (any(types == NORMAL_TYPE)) which(types == NORMAL_TYPE)[1] else 1L]]
    }
    data.table(project = proj, file_id = chr1(h$file_id), file_name = chr1(h$file_name),
               patient = if (is.null(cs)) NA_character_ else chr1(cs$submitter_id),
               sample_barcode = if (is.null(sm)) NA_character_ else chr1(sm$submitter_id),
               sample_type = if (is.null(sm)) NA_character_ else chr1(sm$sample_type))
  })), fill = TRUE)
}
sub_rds   <- file.path(CACHE_DIR, "25_subtype_normal_files.rds")
sub_files <- if (file.exists(sub_rds)) readRDS(sub_rds) else NULL
if (is.null(sub_files)) {
  msg("Querying GDC /files for ", paste(SUBTYPES, collapse = " and "), " normal-tissue files ...")
  got <- tryCatch(rbindlist(lapply(SUBTYPES, fetch_normal_files)), error = function(e) e)
  if (inherits(got, "error")) {
    gdc_notes <- c(gdc_notes, paste0("subtype normal file list: ", conditionMessage(got)))
    msg("  GDC unreachable (", conditionMessage(got), "); availability recorded as unknown")
  } else { sub_files <- got; saveRDS(sub_files, sub_rds) }
}
availability <- rbindlist(lapply(SUBTYPES, function(p) {
  # A live GDC query. 25_normal_consort.tsv counts the cached file map instead.
  src <- "live GDC /files query for Solid Tissue Normal STAR-Counts files (cached in cache/25_subtype_normal_files.rds)"
  if (is.null(sub_files))
    return(data.table(project = p, n_normal_files_gdc = NA_integer_, file_count_source = src,
                      n_on_disk = NA_integer_, analysed = FALSE,
                      status = "GDC unreachable; 11 downloaded primary-tumour files only, so no normals are expected locally"))
  f <- sub_files[project == p]
  nd <- if (!nrow(f)) 0L else {
    f[, path := winlong(file.path(DOWNLOAD_DIR, "validation_data", p, file_id, file_name))]
    sum(file.exists(f$path))
  }
  data.table(project = p, n_normal_files_gdc = nrow(f), file_count_source = src,
             n_on_disk = nd, analysed = nd > 0,
             status = if (nd == 0) "not present locally (11 downloaded primary tumours only); not analysed"
                      else if (nd < nrow(f)) "partially present locally; analysed" else "present locally; analysed")
}))
save_tsv(availability, "25_subtype_normal_availability.tsv")
print(availability)

# ---- 1b. Normal file tables per cohort --------------------------------------
nk <- fm_kirc[sample_type == NORMAL_TYPE, .(file_id, file_name, patient, sample_barcode, sample_type)]
nk[, path := winlong(file.path(EXPR_DIR, file_id, file_name))]
nc <- fm_cptac[sample_type == NORMAL_TYPE, .(file_id, file_name, patient, sample_barcode, sample_type)]
nc[, path := winlong(file.path(VALID_DIR, file_id, file_name))]
normal_files <- list(`TCGA-KIRC` = nk, `CPTAC-3` = nc)
n_gdc_normals <- c(`TCGA-KIRC` = nrow(nk), `CPTAC-3` = nrow(nc))
for (p in SUBTYPES) {
  if (!is.null(sub_files) && availability[project == p, n_on_disk] > 0) {
    f <- sub_files[project == p]
    f[, path := winlong(file.path(DOWNLOAD_DIR, "validation_data", p, file_id, file_name))]
    normal_files[[p]] <- f[file.exists(path), .(file_id, file_name, patient, sample_barcode, sample_type, path)]
    n_gdc_normals[p] <- nrow(f)
  }
}

# ---- 1c. Read every normal library (cached) --------------------------------
# STAR metrics come from read_star_summary() (the production definition), FPKM
# from the unstranded fpkm column as in 01 and 08, and log2(FPKM + 1) is taken
# on the two production gene sets. Rows are keyed on file_id.
read_library <- function(path) {
  qc <- read_star_summary(path)
  d  <- data.table::fread(path, skip = 1, showProgress = FALSE,
                          select = c("gene_id", "fpkm_unstranded"))
  v  <- setNames(as.numeric(d$fpkm_unstranded), d$gene_id)
  list(qc = qc, lnc = log2(v[genes_lnc] + 1), pc = log2(v[genes_pc] + 1))
}
read_set <- function(tbl, label) {
  ok <- file.exists(tbl$path)
  if (!all(ok)) msg("  ", sum(!ok), " ", label, " normal files listed by GDC are not on disk")
  tbl <- tbl[ok]
  L <- matrix(NA_real_, nrow(tbl), length(genes_lnc), dimnames = list(tbl$file_id, genes_lnc))
  P <- matrix(NA_real_, nrow(tbl), length(genes_pc),  dimnames = list(tbl$file_id, genes_pc))
  qc <- vector("list", nrow(tbl))
  for (i in seq_len(nrow(tbl))) {
    r <- read_library(tbl$path[i])
    qc[[i]] <- r$qc; L[i, ] <- r$lnc; P[i, ] <- r$pc
    if (i %% 25 == 0) msg("    ", i, " / ", nrow(tbl))
  }
  stopifnot(!anyNA(L), !anyNA(P))
  list(files = tbl[, .(file_id, file_name, patient, sample_barcode, sample_type)],
       qc = cbind(tbl[, .(file_id)], rbindlist(qc)), lnc = L, pc = P)
}
lib_rds <- file.path(CACHE_DIR, "25_normal_libraries.rds")
fp <- list(files = lapply(normal_files, function(z) sort(z$file_id)),
           genes_lnc = genes_lnc, genes_pc = genes_pc)
normals <- NULL
if (file.exists(lib_rds)) {
  cached <- readRDS(lib_rds)
  if (identical(cached$fingerprint, fp)) { normals <- cached$sets; msg("Normal libraries from cache") }
  else msg("Normal-library cache is stale (file set or gene set changed); re-reading")
}
if (is.null(normals)) {
  normals <- lapply(names(normal_files), function(coh) {
    msg("Reading ", nrow(normal_files[[coh]]), " ", coh, " normal libraries ...")
    read_set(normal_files[[coh]], coh)
  })
  names(normals) <- names(normal_files)
  saveRDS(list(fingerprint = fp, sets = normals), lib_rds)
}

# ---- 1d. Library-failure rule -----------------------------------------------
for (coh in names(normals)) {
  q <- normals[[coh]]$qc
  q[, failed_reads     := assigned_reads < MIN_ASSIGNED_READS]
  q[, failed_noFeature := pct_noFeature > MAX_NOFEATURE_PCT]
  q[, failed_library   := failed_reads | failed_noFeature]
  # CPTAC-3: only the read-depth criterion excludes. 08 applies no rule to the
  # CPTAC-3 tumours, and the ribo-depleted protocol puts every CPTAC-3 library
  # above the non-feature cap calibrated on poly(A) TCGA libraries.
  q[, excluded := if (coh == "CPTAC-3") failed_reads else failed_library]
  normals[[coh]]$qc <- q
  msg(coh, ": ", nrow(q), " normal libraries read; failing the rule as written: ",
      sum(q$failed_library), " (reads < ", MIN_ASSIGNED_READS / 1e6, "M: ", sum(q$failed_reads),
      "; non-feature > ", MAX_NOFEATURE_PCT, "%: ", sum(q$failed_noFeature), "); excluded: ",
      sum(q$excluded))
}

# ---- 1e. Tumour-side tables, one per cohort with normals --------------------
# TCGA-KIRC: the analysed cohort of 01, one library per patient. CPTAC-3: the
# 200 clear-cell patients of 08.
tumours <- list()
tumours[["TCGA-KIRC"]] <- list(
  meta = cohort[, .(patient, sample_barcode, file_id, pct_unmapped, pct_multimapping, pct_noFeature,
                    pct_ambiguous, assigned_reads, tss, plate,
                    os_time, os_event, age, male = as.numeric(sex == "male"),
                    T_stage = as.numeric(T_stage), N_pos = as.numeric(N_pos),
                    M1 = as.numeric(M1), grade_num = as.numeric(grade_num))],
  lnc = t(ds$lnc$expr[genes_lnc, cohort$sample_barcode, drop = FALSE]),
  pc  = t(ds$mrna$expr[genes_pc, cohort$sample_barcode, drop = FALSE]))
tumours[["CPTAC-3"]] <- list(
  meta = merge(vcohort[, .(patient, sample_barcode, file_id, os_time, os_event, age,
                           male = as.numeric(sex == "male"), T_stage = as.numeric(T_stage),
                           N_pos = as.numeric(N_pos), M1 = as.numeric(M1),
                           grade_num = as.numeric(grade_num))],
               vqc, by = "sample_barcode"),
  lnc = NULL, pc = NULL)
tumours[["CPTAC-3"]]$lnc <- t(log2(vd$fpkm[genes_lnc, tumours[["CPTAC-3"]]$meta$sample_barcode, drop = FALSE] + 1))
tumours[["CPTAC-3"]]$pc  <- t(log2(vd$fpkm[genes_pc,  tumours[["CPTAC-3"]]$meta$sample_barcode, drop = FALSE] + 1))
for (p in SUBTYPES) if (p %in% names(normals)) {
  st <- readRDS(file.path(CACHE_DIR, paste0("subtype_", sub("TCGA-", "", p), ".rds")))
  co <- as.data.table(st$cohort); q <- as.data.table(st$qc)[match(co$sample_barcode, sample_barcode)]
  tumours[[p]] <- list(
    meta = data.table(patient = co$patient, sample_barcode = co$sample_barcode,
                      file_id = NA_character_, pct_multimapping = q$pct_multimapping,
                      pct_noFeature = q$pct_noFeature, assigned_reads = q$assigned_reads,
                      tss = tstrsplit(co$sample_barcode, "-", keep = 2)[[1]], plate = NA_character_,
                      os_time = co$os_time, os_event = co$os_event, age = co$age,
                      male = as.numeric(co$sex == "male"), T_stage = as.numeric(co$T_stage),
                      N_pos = as.numeric(co$N_pos), M1 = as.numeric(co$M1),
                      grade_num = as.numeric(co$grade_num)),
    lnc = t(log2(st$fpkm[genes_lnc, co$sample_barcode, drop = FALSE] + 1)),
    pc  = t(log2(st$fpkm[genes_pc,  co$sample_barcode, drop = FALSE] + 1)))
}
for (coh in names(tumours)) {
  m <- tumours[[coh]]$meta
  if (anyNA(m$file_id)) m[, file_id := paste0(coh, "_tumour_", sample_barcode)]
  rownames(tumours[[coh]]$lnc) <- m$file_id; rownames(tumours[[coh]]$pc) <- m$file_id
  m[, `:=`(excluded = FALSE, failed_library = FALSE, failed_reads = FALSE, failed_noFeature = FALSE)]
  # columns every set must carry so the pooled table is not ragged
  for (cc in c("tss", "plate", "file_name")) if (!cc %in% names(m)) set(m, j = cc, value = NA_character_)
  for (cc in c("pct_unmapped", "pct_ambiguous")) if (!cc %in% names(m)) set(m, j = cc, value = NA_real_)
  tumours[[coh]]$meta <- m
}

# =============================================================================
# 2. The discovery axis recomputed and checked, projections, own PCAs
# =============================================================================
banner("2 | Discovery PC1 recomputed as 07; projections; fresh PCAs on normals")
rot_lnc <- pca_lead(E_lnc_disc)
rot_pc  <- pca_lead(E_pc_disc)
sd_l <- sd(rot_lnc$u1); sd_p <- sd(rot_pc$u1)     # the projected axis is reported per discovery SD

ax_cached <- axis$lnc_axis[match(names(rot_lnc$u1), axis$sample_barcode)]
r_axis    <- stats::cor(rot_lnc$u1, ax_cached, use = "complete.obs")
maxdiff   <- max(abs(rot_lnc$u1 - ax_cached), na.rm = TRUE)
self_l    <- max(abs(project_on(E_lnc_disc, rot_lnc) - rot_lnc$u1))
self_p    <- max(abs(project_on(E_pc_disc,  rot_pc)  - rot_pc$u1))
msg("Recomputed lncRNA PC1 versus cache/lnc_global_axis.rds: r = ", sprintf("%.6f", r_axis),
    ", max |difference| = ", signif(maxdiff, 3), " (n = ", sum(!is.na(ax_cached)), ")")
stopifnot(r_axis > 0.9999)
pc07 <- fread(file.path(RESULTS_DIR, "07_pc_variance_explained.tsv"))
q_disc <- function(E) cohort[match(rownames(E), sample_barcode)]
rec_row <- function(label, E, rot, self_diff, r_cache, md, n_cache, lab07) {
  q <- q_disc(E); s <- spearman(rot$u1, q$pct_noFeature)
  data.table(matrix = label, n_samples = nrow(E), n_genes = ncol(E),
             pc1_var_share = round(rot$var_share[1], 4),
             pc1_var_share_07 = pc07[matrix == lab07 & pc == 1, var_explained],
             rho_noFeature = round(s[["rho"]], 3),
             rho_noFeature_07 = pc07[matrix == lab07 & pc == 1, rho_noFeature],
             n_cached_axis = n_cache, pearson_r_vs_cached_axis = if (is.na(r_cache)) NA_real_ else round(r_cache, 6),
             max_abs_diff_vs_cached_axis = if (is.na(md)) NA_real_ else signif(md, 3),
             self_projection_max_abs_diff = signif(self_diff, 3),
             rotation = "kept: v1 (gene loadings), discovery gene means, d1")
}
recheck <- rbind(
  rec_row("lncRNA", E_lnc_disc, rot_lnc, self_l, r_axis, maxdiff, sum(!is.na(ax_cached)), "lncRNA_observed"),
  rec_row("protein_coding", E_pc_disc, rot_pc, self_p, NA_real_, NA_real_, NA_integer_, "protein_coding_observed"))
save_tsv(recheck, "25_axis_recompute_check.tsv"); print(recheck, row.names = FALSE)

# Score every library of every set: projection onto the discovery rotations
# (per discovery SD, so the discovery tumours have mean 0 and SD 1), and the
# set's own PC1 (z within set). For the TCGA-KIRC tumours the "own" PCA is the
# discovery PCA itself (network samples), so that row reproduces 07.
build_set <- function(coh, tissue, meta, E_l, E_p, own_rows_l = NULL, own_rows_p = NULL) {
  stopifnot(identical(rownames(E_l), meta$file_id), identical(rownames(E_p), meta$file_id))
  keep <- !meta$excluded
  if (is.null(own_rows_l)) own_rows_l <- keep
  if (is.null(own_rows_p)) own_rows_p <- keep
  sc <- data.table(file_id = meta$file_id,
                   lnc_axis_proj = project_on(E_l, rot_lnc) / sd_l,
                   pc_axis_proj  = project_on(E_p, rot_pc) / sd_p,
                   mean_lnc_expr = rowMeans(E_l), mean_pc_expr = rowMeans(E_p),
                   lnc_own_pc1_z = NA_real_, pc_own_pc1_z = NA_real_)
  own_l <- own_p <- NULL; sh_l <- sh_p <- NA_real_
  # The third sign anchor: the same libraries' projection onto the discovery
  # PC1, computed above and independent of the PCA being oriented.
  if (sum(own_rows_l) >= 10) {
    own_l <- pca_lead(E_l[own_rows_l, , drop = FALSE], ref = sc$lnc_axis_proj[own_rows_l])
    sc[own_rows_l, lnc_own_pc1_z := as.numeric(scale(own_l$score))]
  }
  if (sum(own_rows_p) >= 10) {
    own_p <- pca_lead(E_p[own_rows_p, , drop = FALSE], ref = sc$pc_axis_proj[own_rows_p])
    sc[own_rows_p, pc_own_pc1_z := as.numeric(scale(own_p$score))]
  }
  if (sum(keep) >= 10) { sh_l <- share_along(E_l[keep, , drop = FALSE], rot_lnc)
                         sh_p <- share_along(E_p[keep, , drop = FALSE], rot_pc) }
  tbl <- copy(meta)
  set(tbl, j = "cohort", value = rep(coh, nrow(tbl)))
  set(tbl, j = "tissue", value = rep(tissue, nrow(tbl)))
  for (cc in setdiff(names(sc), "file_id")) set(tbl, j = cc, value = sc[[cc]])
  tbl[, log_depth := log10(assigned_reads)]
  list(tbl = tbl, own = list(lnc = own_l, pc = own_p), share_along = list(lnc = sh_l, pc = sh_p),
       n_genes = c(lnc = ncol(E_l), pc = ncol(E_p)))
}
sets <- list()
for (coh in names(tumours)) {
  tm <- tumours[[coh]]
  own_l <- own_p <- NULL
  if (coh == "TCGA-KIRC") {
    own_l <- tm$meta$sample_barcode %in% rownames(E_lnc_disc)
    own_p <- tm$meta$sample_barcode %in% rownames(E_pc_disc)
  }
  sets[[paste(coh, "tumour")]] <- build_set(coh, "tumour", tm$meta, tm$lnc, tm$pc, own_l, own_p)
}
for (coh in names(normals)) {
  nn <- normals[[coh]]
  meta <- merge(nn$files, nn$qc, by = "file_id")
  meta <- meta[match(rownames(nn$lnc), file_id)]
  meta[, tss := if (startsWith(coh, "TCGA")) tstrsplit(sample_barcode, "-", keep = 2)[[1]] else NA_character_]
  meta[, plate := NA_character_]
  sets[[paste(coh, "normal")]] <- build_set(coh, "normal", meta, nn$lnc, nn$pc)
}
# Check: projected scores of the discovery samples equal the cached axis.
chk <- sets[["TCGA-KIRC tumour"]]$tbl
chk_r <- stats::cor(chk$lnc_axis_proj[match(axis$sample_barcode, chk$sample_barcode)] * sd_l, axis$lnc_axis)
msg("Projected TCGA-KIRC tumour axis versus cached axis on the 511 discovery samples: r = ",
    sprintf("%.6f", chk_r))
stopifnot(chk_r > 0.9999)

# Check that every reported PC1 has a sign fixed by a decisive anchor.
banner("2b | PC1 sign convention, every matrix")
sign_objs <- c(list(`discovery lncRNA (07 rotation)` = rot_lnc,
                    `discovery protein_coding (07 rotation)` = rot_pc),
               setNames(lapply(names(sets), function(nm) sets[[nm]]$own$lnc),
                        paste(names(sets), "own lncRNA PC1")),
               setNames(lapply(names(sets), function(nm) sets[[nm]]$own$pc),
                        paste(names(sets), "own protein_coding PC1")))
sign_tbl <- assert_sign_convention(sign_objs)

lib <- rbindlist(lapply(sets, `[[`, "tbl"), fill = TRUE)
lib[, in_discovery_pca := tissue == "tumour" & cohort == "TCGA-KIRC" & sample_barcode %in% rownames(E_lnc_disc)]
tumour_key <- lib[tissue == "tumour", paste(cohort, patient)]
lib[, tumour_in_analysed_cohort := paste(cohort, patient) %in% tumour_key]

# ---- per-cohort axis summaries ---------------------------------------------
proj_rows <- rbindlist(lapply(names(sets), function(nm) {
  S <- sets[[nm]]; d <- S$tbl[excluded == FALSE]
  covs <- list(noFeature = d$pct_noFeature, multimap = d$pct_multimapping, log_depth = d$log_depth)
  mk <- function(matrix, axis_type, x, n_genes, var_share, mean_expr, disc_E,
                 own_share = rep(NA_real_, 5), r_own = NA_real_, n_own = NA_integer_,
                 sgn = NULL) {
    ok <- is.finite(x)
    # Size-matched baseline for the own-PC1 share (see n_matched_pc1()).
    nmp <- if (axis_type == "own_PC1") n_matched_pc1(disc_E, sum(ok), matrix) else NM_NA
    out <- data.table(cohort = d$cohort[1], tissue = d$tissue[1], matrix = matrix, axis_type = axis_type,
                      n = sum(ok), n_genes = n_genes, var_share = round(var_share, 4))
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
  # r_projected_vs_own is a complete-case correlation with its own n (the
  # samples in the set's own PCA).
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
proj_rows[, var_share_definition := fifelse(axis_type == "own_PC1", "PC1 share of the set's own centred matrix",
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
# A projected row inherits the discovery rotation's orientation. Its anchor
# values are the discovery matrix's.
proj_rows[axis_type != "own_PC1" & !is.na(pc1_sign_note),
          pc1_sign_note := paste0("orientation inherited from the discovery rotation, which is fixed ",
                                  "once on the discovery matrix and is the same for every sample set; ",
                                  "the anchor values below are the discovery matrix's. ", pc1_sign_note)]
# Multiplicity: one BH family per cohort x tissue set, 4 axes x 4 covariates
# = 16 tests.
PCOLS <- c("p_noFeature", "p_multimap", "p_log_depth", "p_mean_expr")
FCOLS <- sub("^p_", "fdr_", PCOLS)
proj_rows[, (FCOLS) := {
  nn <- .N; a <- bh(unlist(.SD, use.names = FALSE))
  lapply(seq_along(.SD), function(j) signif(a[((j - 1L) * nn + 1L):(j * nn)], 3))
}, by = .(cohort, tissue), .SDcols = PCOLS]
proj_rows[, fdr_family := paste0("Benjamini-Hochberg within cohort x tissue over the ",
                                 length(PCOLS) * 4, " reported Spearman tests of that set (4 axes x ",
                                 length(PCOLS), " covariates)")]
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
save_tsv(proj_rows, "25_normal_axis_projection.tsv")
print(proj_rows[, .(cohort, tissue, matrix, axis_type, n, var_share, n_matched_pc1_share,
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
# 3. Plate of each normal aliquot (GDC, guarded), pairs, variance partition
# =============================================================================
banner("3 | Normal aliquot plates; paired tumour-normal analysis; variance partition")
fetch_file_aliquots <- function(file_ids) {
  body <- list(filters = list(op = "in", content = list(field = "file_id", value = as.list(unique(file_ids)))),
               fields = paste("file_id", "associated_entities.entity_submitter_id",
                              "associated_entities.entity_type",
                              "cases.samples.submitter_id", "cases.samples.sample_type",
                              "cases.samples.portions.analytes.analyte_type",
                              "cases.samples.portions.analytes.aliquots.submitter_id", sep = ","),
               format = "JSON", size = "600")
  hits <- gdc_raw("files", body)
  rbindlist(lapply(hits, function(h) {
    ae <- h$associated_entities
    a_id <- if (length(ae)) vapply(ae, function(e) chr1(e$entity_submitter_id), character(1)) else character(0)
    a_ty <- if (length(ae)) vapply(ae, function(e) chr1(e$entity_type), character(1)) else character(0)
    nested <- character(0)
    for (cs in h$cases) for (s in cs$samples) for (p in s$portions) for (a in p$analytes)
      for (q in a$aliquots) nested <- c(nested, chr1(q$submitter_id))
    data.table(file_id = chr1(h$file_id),
               aliquot_assoc = if (any(a_ty == "aliquot")) a_id[a_ty == "aliquot"][1] else NA_character_,
               aliquots_nested = paste(unique(nested[!is.na(nested)]), collapse = ";"))
  }))
}
al_rds <- file.path(CACHE_DIR, "25_normal_aliquots.rds")
al_cache <- if (file.exists(al_rds)) readRDS(al_rds) else NULL
tcga_normal_ids <- unlist(lapply(names(normals)[startsWith(names(normals), "TCGA")],
                                 function(coh) normals[[coh]]$files$file_id))
if (is.null(al_cache) || !all(tcga_normal_ids %in% al_cache$file_id)) {
  msg("Querying GDC /files for the aliquot barcode of ", length(tcga_normal_ids), " TCGA normal files ...")
  got <- tryCatch({
    am <- fetch_file_aliquots(tcga_normal_ids)
    if (!nrow(am)) stop("GDC returned no records")
    am
  }, error = function(e) e)
  if (inherits(got, "error")) {
    gdc_notes <- c(gdc_notes, paste0("normal aliquot lookup: ", conditionMessage(got)))
    msg("  GDC aliquot lookup FAILED (", conditionMessage(got), "); falling back to the biospecimen cache of 23")
  } else { al_cache <- got; saveRDS(al_cache, al_rds); msg("  aliquot recovered for ", nrow(got), " files") }
}
# Plate = field 6 of the 28-character aliquot barcode. Sources in order: the
# aliquot the file is associated with, then an RNA aliquot (analyte letter R,
# centre 07) nested under the file's own sample, then the RNA aliquots of that
# normal sample in cache/biospecimen_kirc.rds (modal plate).
bio_rds <- file.path(CACHE_DIR, "biospecimen_kirc.rds")
bio_an  <- if (file.exists(bio_rds)) readRDS(bio_rds)$analytes else NULL
normal_plate <- function(fid, sb) {
  al <- NA_character_; src <- NA_character_
  if (!is.null(al_cache) && fid %in% al_cache$file_id) {
    idx <- match(fid, al_cache$file_id)
    r <- al_cache[idx]
    if (!is.na(r$aliquot_assoc)) { al <- r$aliquot_assoc; src <- "file_aliquot" }
    else if (nzchar(r$aliquots_nested)) {
      cand <- strsplit(r$aliquots_nested, ";", fixed = TRUE)[[1]]
      cand <- cand[startsWith(cand, sb) & sub("^[0-9]+", "", bc_field(cand, 5)) == "R" & bc_field(cand, 7) == "07"]
      if (length(cand)) { al <- modal(cand); src <- "file_nested_rna_aliquot" }
    }
  }
  if (is.na(al) && !is.null(bio_an)) {
    ba <- as.data.table(bio_an)
    b <- ba[ba$sample_barcode == sb & !is.na(ba$analyte_type) & startsWith(ba$analyte_type, "RNA") & !is.na(ba$aliquot_id)]
    b <- b[sub("^[0-9]+", "", bc_field(b$aliquot_id, 5)) == "R" & bc_field(b$aliquot_id, 7) == "07"]
    if (nrow(b)) { al <- modal(b$aliquot_id); src <- "biospecimen_cache_fallback" }
  }
  list(aliquot_id = al, plate = if (is.na(al)) NA_character_ else bc_field(al, 6), plate_source = src)
}
lib[, `:=`(aliquot_id = NA_character_, plate_source = NA_character_)]
if (length(tcga_normal_ids)) {
  pl <- rbindlist(lapply(seq_along(tcga_normal_ids), function(i) {
    fid <- tcga_normal_ids[i]; sb <- lib$sample_barcode[match(fid, lib$file_id)]
    z <- as.data.table(normal_plate(fid, sb)); z[, file_id := fid]; z
  }))
  lib[pl, on = "file_id", `:=`(aliquot_id = i.aliquot_id, plate = i.plate, plate_source = i.plate_source)]
  msg("Plate for TCGA normals: ", sum(!is.na(pl$plate)), " of ", nrow(pl), " (",
      paste(names(table(pl$plate_source)), table(pl$plate_source), sep = " = ", collapse = "; "), ")")
}
lib[tissue == "tumour" & cohort == "TCGA-KIRC", plate_source := fifelse(is.na(plate), NA_character_, "cohort_full (23)")]

# ---- normal library metrics (one row per normal library) --------------------
norm_out <- lib[tissue == "normal", .(cohort, file_id, file_name, sample_barcode, patient, sample_type,
                                      pct_unmapped, pct_multimapping, pct_noFeature, pct_ambiguous,
                                      assigned_reads, log_depth, failed_reads, failed_noFeature,
                                      failed_library, excluded, tumour_in_analysed_cohort,
                                      lnc_axis_proj = round(lnc_axis_proj, 4),
                                      lnc_own_pc1_z = round(lnc_own_pc1_z, 4),
                                      pc_axis_proj = round(pc_axis_proj, 4),
                                      pc_own_pc1_z = round(pc_own_pc1_z, 4),
                                      mean_lnc_expr = round(mean_lnc_expr, 4),
                                      mean_pc_expr = round(mean_pc_expr, 4),
                                      tss, aliquot_id, plate, plate_source)]
setorder(norm_out, cohort, patient)
save_tsv(norm_out, "25_normal_library_metrics.tsv")

# ---- pairs ------------------------------------------------------------------
# pct_unmapped and pct_ambiguous are carried only for the matched columns of
# 25_tumour_vs_normal_distribution.tsv.
LIB_COLS  <- c("sample_barcode", "file_id", "pct_noFeature", "pct_multimapping", "assigned_reads",
               "log_depth", "lnc_axis_proj", "pc_axis_proj", "lnc_own_pc1_z", "pc_own_pc1_z",
               "mean_lnc_expr", "pct_unmapped", "pct_ambiguous", "tss", "plate")
CLIN_COLS <- c("os_time", "os_event", "age", "male", "T_stage", "N_pos", "M1", "grade_num")
pair_cohorts <- intersect(names(tumours), names(normals))
pairs <- rbindlist(lapply(pair_cohorts, function(coh) {
  tu <- lib[cohort == coh & tissue == "tumour" & excluded == FALSE, c("cohort", "patient", LIB_COLS, CLIN_COLS), with = FALSE]
  no <- lib[cohort == coh & tissue == "normal" & excluded == FALSE, c("cohort", "patient", LIB_COLS), with = FALSE]
  stopifnot(!anyDuplicated(tu$patient))
  if (anyDuplicated(no$patient)) {
    msg(coh, ": ", sum(duplicated(no$patient)), " patients with two normal libraries; keeping the deeper one")
    setorder(no, patient, -assigned_reads); no <- no[!duplicated(patient)]
  }
  merge(tu, no, by = c("cohort", "patient"), suffixes = c("_t", "_n"))
}))
stopifnot(nrow(pairs) > 0)
# The paired set is defined by library: where a patient has two normal
# libraries only the deeper one is paired, keeping two libraries per patient.
paired_ids <- c(pairs$file_id_t, pairs$file_id_n)
lib[, paired := file_id %in% paired_ids & excluded == FALSE]
for (coh in pair_cohorts) msg(coh, ": ", sum(pairs$cohort == coh), " patients with an analysed tumour and a normal library")

PAIR_VARS <- c(pct_noFeature = "non-feature read fraction (%)",
               pct_multimapping = "multimapping read fraction (%)",
               log_depth = "log10 assigned reads",
               lnc_axis_proj = "lncRNA axis, projected on discovery PC1 (per discovery SD)",
               pc_axis_proj = "protein-coding PC1, projected (per discovery SD; INTENDED negative control, not a clean one: in the discovery tumours this projection is largely the multimapping axis)")
paired_tbl <- rbindlist(lapply(pair_cohorts, function(coh) rbindlist(lapply(names(PAIR_VARS), function(v) {
  d <- pairs[cohort == coh]; tv <- d[[paste0(v, "_t")]]; nv <- d[[paste0(v, "_n")]]
  ok <- is.finite(tv) & is.finite(nv); tv <- tv[ok]; nv <- nv[ok]
  if (length(tv) < 5) return(NULL)
  s <- spearman(tv, nv)
  w <- suppressWarnings(wilcox.test(tv, nv, paired = TRUE, exact = FALSE))
  qt <- qtl(tv); qn <- qtl(nv); qd <- qtl(tv - nv)
  data.table(cohort = coh, variable = v, description = PAIR_VARS[[v]], n_pairs = length(tv),
             median_tumour = round(qt[["median"]], 4), q25_tumour = round(qt[["q25"]], 4), q75_tumour = round(qt[["q75"]], 4),
             median_normal = round(qn[["median"]], 4), q25_normal = round(qn[["q25"]], 4), q75_normal = round(qn[["q75"]], 4),
             median_paired_diff_t_minus_n = round(qd[["median"]], 4),
             n_tumour_higher = sum(tv > nv), n_normal_higher = sum(nv > tv),
             spearman_rho_tumour_vs_normal = round(s[["rho"]], 3), p_rho = signif(s[["p"]], 3),
             wilcoxon_signed_rank_p = signif(w$p.value, 3))
}))))
stopifnot(nrow(paired_tbl) > 0)
paired_tbl[, fdr_rho := signif(bh(p_rho), 3), by = cohort]
paired_tbl[, fdr_wilcoxon := signif(bh(wilcoxon_signed_rank_p), 3), by = cohort]
# shared_within_patient: positive tumour-versus-normal rank correlation at
# fdr < FDR_ALPHA, i.e. a patient-level component. The pre-specified tests are
# pct_noFeature and lnc_axis_proj. The projected protein-coding PC1 is the
# intended negative control, but in the discovery tumours it tracks the
# multimapping fraction, so it is not a clean biological control.
paired_tbl[, shared_within_patient := is.finite(fdr_rho) & fdr_rho < FDR_ALPHA &
             spearman_rho_tumour_vs_normal > 0]
PAIR_ROLE <- c(pct_noFeature    = "pre-specified sharing test",
               lnc_axis_proj    = "pre-specified sharing test",
               pc_axis_proj     = "intended negative control",
               pct_multimapping = "companion library metric",
               log_depth        = "companion library metric")
paired_tbl[, test_role := unname(PAIR_ROLE[variable])]
paired_tbl[, test_outcome := fifelse(shared_within_patient,
  "POSITIVE: shared within patient",
  fifelse(is.finite(fdr_rho) & spearman_rho_tumour_vs_normal < 0 & fdr_rho < FDR_ALPHA,
          "NEGATIVE: tumour and normal ranks oppose within patient",
          "NULL: no within-patient sharing detected"))]
# Discovery-tumour multimapping correlation of the projected protein-coding
# PC1, quoted in the reading column.
rho_pc_mm_disc <- proj_rows[cohort == "TCGA-KIRC" & tissue == "tumour" &
                            matrix == "protein_coding" & axis_type == "projected_discovery_PC1",
                            rho_multimap][1]
paired_tbl[, reading := fifelse(
  test_role == "pre-specified sharing test" & !shared_within_patient,
  paste0("PRE-SPECIFIED TEST FAILED TO DISCRIMINATE. The within-patient sharing test is null here ",
         "(rho ", spearman_rho_tumour_vs_normal, ", p ", p_rho, ", fdr ", fdr_rho, ", ", n_pairs,
         " pairs). Report it as a test that returned nothing, not as evidence for or against a ",
         "patient-level cause; the positive evidence for the library reading is the normals' own ",
         "correlation with the non-feature fraction (25_normal_axis_projection.tsv) and the plate ",
         "partition (25_paired_variance_partition.tsv), neither of which uses the pairing."),
  fifelse(test_role == "pre-specified sharing test" & shared_within_patient,
  paste0("Pre-specified test positive: this quantity has a patient-level component (rho ",
         spearman_rho_tumour_vs_normal, ", fdr ", fdr_rho, ")."),
  fifelse(test_role == "intended negative control" & shared_within_patient,
  paste0("NOT A CLEAN BIOLOGICAL CONTROL. This intended negative control is itself shared within ",
         "patient (rho ", spearman_rho_tumour_vs_normal, ", fdr ", fdr_rho, "), and in the discovery ",
         "tumours the projected protein-coding PC1 is largely the multimapping axis (Spearman rho ",
         rho_pc_mm_disc, " with the multimapping fraction, n 528). Its sharing is therefore evidence ",
         "of a patient-level component in the LIBRARIES, not evidence that it carries tumour biology, ",
         "and the contrast with the unshared lncRNA axis is not a clean one."),
  fifelse(test_role == "intended negative control",
  "Intended negative control, not shared within patient in this cohort; descriptive only.",
  "Companion library metric, reported for context; not a pre-specified test."))))]
paired_tbl[, note := paste0(
  "Benjamini-Hochberg within cohort over the ", length(PAIR_VARS),
  " variables, separately for the rank correlation and for the signed-rank test; ",
  "shared_within_patient = positive tumour-versus-normal rho at fdr < ", FDR_ALPHA,
  " (a patient-level component); the two pre-specified rows are pct_noFeature and lnc_axis_proj ",
  "and both are null in TCGA-KIRC, so the pairing discriminates nothing and is reported as a ",
  "failed test; see the reading column")]
save_tsv(paired_tbl, "25_paired_tumour_normal.tsv")
print(paired_tbl[, .(cohort, variable, test_role, n_pairs, spearman_rho_tumour_vs_normal,
                     p_rho, fdr_rho, shared_within_patient, test_outcome,
                     wilcoxon_signed_rank_p)], row.names = FALSE)
for (coh in pair_cohorts) {
  ps <- paired_tbl[cohort == coh & test_role == "pre-specified sharing test"]
  sh <- paired_tbl[cohort == coh & shared_within_patient == TRUE]
  msg(coh, ": PRE-SPECIFIED within-patient sharing test -- ",
      paste(sprintf("%s %s (rho = %.3f, p = %s, fdr = %s)", ps$variable,
                    ifelse(ps$shared_within_patient, "SHARED", "NULL"),
                    ps$spearman_rho_tumour_vs_normal, format(ps$p_rho), format(ps$fdr_rho)),
            collapse = "; "))
  if (!any(ps$shared_within_patient))
    msg(coh, ": the pre-specified test discriminates NOTHING (no quantity it was specified on is ",
        "shared within patient). Report it as a test that failed, not as evidence; the positive ",
        "evidence for the library reading is the normals' own axis-versus-metric correlation and ",
        "the plate partition, neither of which uses the pairing.")
  msg(coh, ": quantities shared within patient (positive rho at fdr < ", FDR_ALPHA, "): ",
      if (!nrow(sh)) "none" else paste(sprintf("%s [%s] (rho = %.3f, fdr = %s)", sh$variable,
                                               sh$test_role,
                                               sh$spearman_rho_tumour_vs_normal,
                                               format(sh$fdr_rho, scientific = TRUE)), collapse = "; "))
  if (nrow(sh[test_role == "intended negative control"]))
    msg(coh, ": the ONLY shared quantity of interest is the intended negative control ",
        "(projected protein-coding PC1), which in the discovery tumours is largely the multimapping ",
        "axis (rho = ", rho_pc_mm_disc, "), so its sharing is a patient-level LIBRARY component and ",
        "not tumour biology. Do not report the protein-coding-shared / lncRNA-not-shared contrast ",
        "as a clean biological control.")
}

# ---- variance partition ------------------------------------------------------
# One-way R2 (and adjusted R2) of each metric on a factor. Tissue source site
# and plate levels with fewer than 10 libraries are pooled, as in 23. In the
# paired set the patient factor also gets a one-way intraclass correlation,
# because R2 on a factor with n/2 levels is inflated by its degrees of freedom.
vp_row <- function(coh, set_label, v, f_label, y, g, pooled_min = NULL, want_icc = FALSE) {
  ok <- is.finite(y) & !is.na(g) & nzchar(as.character(g))
  y <- y[ok]; g0 <- factor(as.character(g[ok]))
  if (length(y) < 10 || nlevels(g0) < 2 || nlevels(g0) == length(y)) return(NULL)
  fit0 <- lm(y ~ g0); r2u <- summary(fit0)$r.squared
  g1 <- if (!is.null(pooled_min)) pool_levels(g0, pooled_min) else g0
  base <- data.table(cohort = coh, library_set = set_label, variable = v, factor = f_label,
                     n_libraries = length(y), n_levels = nlevels(g1),
                     n_pooled_other = if (!is.null(pooled_min)) sum(g1 == "other") else 0L,
                     R2 = NA_real_, adj_R2 = NA_real_, anova_p = NA_real_, kruskal_p = NA_real_,
                     R2_unpooled = round(r2u, 4), n_levels_unpooled = nlevels(g0), icc = NA_real_)
  if (nlevels(g1) < 2) return(base)
  fit <- lm(y ~ g1); s <- summary(fit)
  # Computed outside the := expression: inside it a bare `icc` would resolve to
  # the column of that name, not to the function's argument.
  r2v <- round(s$r.squared, 4); adjv <- round(s$adj.r.squared, 4)
  apv <- signif(anova(fit)[["Pr(>F)"]][1], 3)
  kpv <- if (want_icc) NA_real_ else signif(kruskal.test(y ~ g1)$p.value, 3)
  iccv <- NA_real_
  if (want_icc) {
    a <- anova(fit0); msb <- a[["Mean Sq"]][1]; msw <- a[["Mean Sq"]][2]; k <- length(y) / nlevels(g0)
    iccv <- round((msb - msw) / (msb + (k - 1) * msw), 4)
  }
  base[, `:=`(R2 = r2v, adj_R2 = adjv, anova_p = apv, kruskal_p = kpv, icc = iccv)]
  base
}
VP_VARS <- c("pct_noFeature", "lnc_axis_proj", "pct_multimapping")
vp <- rbindlist(lapply(pair_cohorts, function(coh) {
  is_tcga <- startsWith(coh, "TCGA")
  pooled <- lib[cohort == coh & paired == TRUE]
  all_l  <- lib[cohort == coh & excluded == FALSE]
  # Tumours and normals are also partitioned separately, since the plate
  # effect can differ between the two sides.
  tum_l  <- lib[cohort == coh & tissue == "tumour" & excluded == FALSE]
  nor_l  <- lib[cohort == coh & tissue == "normal" & excluded == FALSE]
  rbindlist(lapply(VP_VARS, function(v) rbindlist(c(
    list(vp_row(coh, "paired patients (tumour + normal)", v, "patient", pooled[[v]], pooled$patient, want_icc = TRUE),
         vp_row(coh, "paired patients (tumour + normal)", v, "tissue_type", pooled[[v]], pooled$tissue),
         { # tissue type after patient: partial R2 of tissue given patient (two-way additive)
           ok <- is.finite(pooled[[v]])
           if (sum(ok) >= 10) {
             y <- pooled[[v]][ok]; pf <- factor(pooled$patient[ok]); tf <- factor(pooled$tissue[ok])
             f1 <- lm(y ~ pf); f2 <- lm(y ~ pf + tf)
             # This row reports partial and two-way R2. R2, adj_R2 and
             # R2_unpooled stay empty because they are not comparable.
             data.table(cohort = coh, library_set = "paired patients (tumour + normal)", variable = v,
                        factor = "tissue_type given patient (partial R2)", n_libraries = length(y),
                        n_levels = 2L, n_pooled_other = 0L,
                        R2 = NA_real_, adj_R2 = NA_real_,
                        anova_p = signif(anova(f1, f2)[["Pr(>F)"]][2], 3), kruskal_p = NA_real_,
                        R2_unpooled = NA_real_, n_levels_unpooled = nlevels(pf) + 1L,
                        icc = NA_real_,
                        partial_R2 = round((deviance(f1) - deviance(f2)) / deviance(f1), 4),
                        twoway_R2 = round(summary(f2)$r.squared, 4),
                        twoway_adj_R2 = round(summary(f2)$adj.r.squared, 4))
           } else NULL },
         if (is_tcga) vp_row(coh, "paired patients (tumour + normal)", v, "tss", pooled[[v]], pooled$tss, pooled_min = 10),
         if (is_tcga) vp_row(coh, "paired patients (tumour + normal)", v, "plate", pooled[[v]], pooled$plate, pooled_min = 10),
         vp_row(coh, "all analysed libraries", v, "tissue_type", all_l[[v]], all_l$tissue),
         if (is_tcga) vp_row(coh, "all analysed libraries", v, "tss", all_l[[v]], all_l$tss, pooled_min = 10),
         if (is_tcga) vp_row(coh, "all analysed libraries", v, "plate", all_l[[v]], all_l$plate, pooled_min = 10),
         if (is_tcga) vp_row(coh, "analysed tumour libraries only", v, "tss", tum_l[[v]], tum_l$tss, pooled_min = 10),
         if (is_tcga) vp_row(coh, "analysed tumour libraries only", v, "plate", tum_l[[v]], tum_l$plate, pooled_min = 10),
         if (is_tcga) vp_row(coh, "analysed normal libraries only", v, "tss", nor_l[[v]], nor_l$tss, pooled_min = 10),
         if (is_tcga) vp_row(coh, "analysed normal libraries only", v, "plate", nor_l[[v]], nor_l$plate, pooled_min = 10)))
  , fill = TRUE)), fill = TRUE)
}), fill = TRUE)
vp[, note := fifelse(factor == "patient",
       "R2 on a factor with one level per patient; icc = one-way intraclass correlation of the two libraries",
     fifelse(factor %in% c("tss", "plate"),
       "levels with fewer than 10 libraries pooled as other; R2_unpooled uses every level",
     fifelse(grepl("partial R2", factor, fixed = TRUE),
       "partial_R2 = (deviance(y ~ patient) - deviance(y ~ patient + tissue)) / deviance(y ~ patient); twoway_R2 and twoway_adj_R2 describe the two-way model y ~ patient + tissue; R2, adj_R2 and R2_unpooled are not defined for this row",
       "")))]
# Record that site and plate are not evaluated outside TCGA.
non_tcga <- pair_cohorts[!startsWith(pair_cohorts, "TCGA")]
if (length(non_tcga)) {
  vp <- rbind(vp, rbindlist(lapply(non_tcga, function(coh)
    rbindlist(lapply(VP_VARS, function(v) rbindlist(lapply(c("tss", "plate"), function(f)
      data.table(cohort = coh, library_set = "paired patients (tumour + normal)", variable = v,
                 factor = f, n_libraries = 0L, n_levels = 0L, n_pooled_other = 0L,
                 R2 = NA_real_, adj_R2 = NA_real_, anova_p = NA_real_, kruskal_p = NA_real_,
                 R2_unpooled = NA_real_, n_levels_unpooled = 0L, icc = NA_real_,
                 note = "not evaluated: tissue source site and plate are TCGA barcode fields, and CPTAC-3 libraries carry neither")))))))
    , fill = TRUE)
}
# ---- plate partition by side ----
# Every plate row carries the tumour-only and normal-only plate R2, because
# the pooled figure is inflated when the normals sit on few plates.
PLATE_SIDE_COLS <- c("plate_R2_tumour_only", "plate_p_tumour_only", "plate_n_tumour_only",
                     "plate_R2_normal_only", "plate_p_normal_only", "plate_n_normal_only",
                     "plate_n_levels_normal_only")
vp[, (PLATE_SIDE_COLS) := NA_real_]
plate_side <- vp[factor == "plate" & library_set %in% c("analysed tumour libraries only",
                                                        "analysed normal libraries only")]
if (nrow(plate_side)) {
  pt <- plate_side[library_set == "analysed tumour libraries only",
                   .(cohort, variable, R2_t = R2, p_t = anova_p, n_t = n_libraries, lev_t = n_levels)]
  pn <- plate_side[library_set == "analysed normal libraries only",
                   .(cohort, variable, R2_n = R2, p_n = anova_p, n_n = n_libraries, lev_n = n_levels)]
  ps <- merge(pt, pn, by = c("cohort", "variable"), all = TRUE)
  vp[ps, on = .(cohort, variable), `:=`(
    plate_R2_tumour_only = i.R2_t, plate_p_tumour_only = i.p_t, plate_n_tumour_only = i.n_t,
    plate_R2_normal_only = i.R2_n, plate_p_normal_only = i.p_n, plate_n_normal_only = i.n_n,
    plate_n_levels_normal_only = i.lev_n)]
  vp[factor != "plate", (PLATE_SIDE_COLS) := NA_real_]
  # How much of the plate space the normal libraries actually occupy.
  plate_cov <- rbindlist(lapply(pair_cohorts, function(coh) {
    tl <- lib[cohort == coh & tissue == "tumour" & excluded == FALSE & !is.na(plate), plate]
    nl2 <- lib[cohort == coh & tissue == "normal" & excluded == FALSE & !is.na(plate), plate]
    if (!length(nl2)) return(NULL)
    # Libraries on the busiest k plates, with k capped at the number of plates.
    PLATE_TOP_K <- 3L
    tb <- sort(table(nl2), decreasing = TRUE)
    k  <- min(PLATE_TOP_K, length(tb))
    data.table(cohort = coh, n_plates_tumour = uniqueN(tl), n_plates_normal = uniqueN(nl2),
               txt = paste0(sum(tb[seq_len(k)]), " of ", length(nl2),
                            " normal libraries sit on ", k, " of the ", length(tb),
                            " plates they span, against ", uniqueN(tl),
                            " plates among the ", length(tl), " tumour libraries"))
  }))
  vp[, plate_coverage_note := NA_character_]
  vp[plate_cov, on = "cohort", plate_coverage_note := i.txt]
  vp[factor != "plate", plate_coverage_note := NA_character_]
  vp[factor == "plate" & is.finite(plate_R2_tumour_only), note := paste0(
    note, "; PLATE IS A TUMOUR-SIDE PARTITION: within the ", plate_n_tumour_only,
    " analysed tumour libraries plate explains R2 ", plate_R2_tumour_only, " (p ", plate_p_tumour_only,
    "), within the ", plate_n_normal_only, " analysed normal libraries R2 ", plate_R2_normal_only,
    " on ", plate_n_levels_normal_only, " pooled levels (p ", plate_p_normal_only, "); ",
    plate_coverage_note,
    ", so the pooled and paired-set plate figures must never be quoted without the normal-side figure")]
  for (coh in unique(plate_side$cohort)) {
    r <- vp[cohort == coh & factor == "plate" & variable == "pct_noFeature" &
            library_set == "all analysed libraries"]
    if (nrow(r)) msg(coh, ": plate partition of the non-feature fraction -- pooled R2 ", r$R2,
                     " (", r$n_libraries, " libraries); TUMOURS ONLY R2 ", r$plate_R2_tumour_only,
                     " (p ", r$plate_p_tumour_only, "); NORMALS ONLY R2 ", r$plate_R2_normal_only,
                     " (p ", r$plate_p_normal_only, "). The plate partition is a tumour-side result; ",
                     r$plate_coverage_note, ".")
  }
}

# Multiplicity: BH within each metric over cohorts, library sets and factors,
# with the F and Kruskal-Wallis tests as separate families. Tumour-only and
# normal-only rows form their own family, since they re-partition libraries
# already counted in the pooled rows.
vp[, fdr_subfamily := fifelse(library_set %in% c("analysed tumour libraries only",
                                                 "analysed normal libraries only"),
                              "tumour-only and normal-only diagnostic split", "primary partitions")]
vp[, fdr_anova   := signif(bh(anova_p), 3),   by = .(variable, fdr_subfamily)]
vp[, fdr_kruskal := signif(bh(kruskal_p), 3), by = .(variable, fdr_subfamily)]
# The single-family FDR is also computed. Rows whose significance depends on
# the split are named in fdr_family.
vp[, fdr_anova_pooled   := signif(bh(anova_p), 3),   by = variable]
vp[, fdr_kruskal_pooled := signif(bh(kruskal_p), 3), by = variable]
crossed <- rbindlist(lapply(c("anova", "kruskal"), function(tst) {
  s <- vp[[paste0("fdr_", tst)]]; p <- vp[[paste0("fdr_", tst, "_pooled")]]
  i <- which(is.finite(s) & is.finite(p) & ((s < FDR_ALPHA) != (p < FDR_ALPHA)))
  if (!length(i)) return(NULL)
  data.table(test = if (tst == "anova") "one-way F" else "Kruskal-Wallis",
             cohort = vp$cohort[i], library_set = vp$library_set[i], variable = vp$variable[i],
             factor = vp$factor[i], fdr_split = s[i], fdr_pooled = p[i])
}))
n_diag <- vp[fdr_subfamily != "primary partitions", .N]
crossed_txt <- if (!nrow(crossed)) {
  paste0("moves no row across the ", FDR_ALPHA, " line")
} else {
  paste0("moves ", nrow(crossed), " row(s) across the ", FDR_ALPHA, " line: ",
         paste(sprintf("%s, %s, %s, factor %s -- %s FDR %s under the split against %s pooled",
                       crossed$cohort, crossed$library_set, crossed$variable, crossed$factor,
                       crossed$test, format(crossed$fdr_split), format(crossed$fdr_pooled)),
               collapse = "; "))
}
vp[, fdr_family := paste0(
  "Benjamini-Hochberg within variable and within the family named in fdr_subfamily; the primary family ",
  "is every cohort, library set and factor of the pooled and paired partitions, the diagnostic family ",
  "is the tumour-only and normal-only split of the same libraries; anova_p and kruskal_p adjusted as ",
  "separate families. CONSEQUENCE OF THE SPLIT, disclosed rather than left implicit: folding the ",
  n_diag, " diagnostic rows into the primary family (the counterfactual carried in ",
  "fdr_anova_pooled and fdr_kruskal_pooled) ", crossed_txt,
  ". The split is used because the diagnostic rows re-partition libraries already counted in the ",
  "pooled rows and test nothing new, but a reader who prefers one family should read the pooled columns")]
if (nrow(crossed)) {
  msg("fdr_subfamily split: ", nrow(crossed), " row(s) change significance at FDR ", FDR_ALPHA,
      " when the ", n_diag, " diagnostic rows are folded into the primary family")
  print(crossed, row.names = FALSE)
} else msg("fdr_subfamily split: no row changes significance at FDR ", FDR_ALPHA,
           " when the ", n_diag, " diagnostic rows are folded into the primary family")
save_tsv(vp, "25_paired_variance_partition.tsv")
print(vp[variable == "pct_noFeature", .(cohort, library_set, factor, n_libraries, n_levels, R2, adj_R2,
                                        partial_R2, anova_p, fdr_anova, icc)], row.names = FALSE)

# ---- consort / counts per cohort ------------------------------------------
consort <- rbindlist(lapply(names(normals), function(coh) {
  q <- normals[[coh]]$qc; nl <- lib[cohort == coh & tissue == "normal"]
  pr <- pairs[cohort == coh]
  same_plate <- if (nrow(pr)) sum(!is.na(pr$plate_t) & !is.na(pr$plate_n) & pr$plate_t == pr$plate_n) else NA_integer_
  both_plate <- if (nrow(pr)) sum(!is.na(pr$plate_t) & !is.na(pr$plate_n)) else NA_integer_
  data.table(cohort = coh, n_normal_files_in_file_map = n_gdc_normals[[coh]],
             file_count_source = switch(coh,
               `TCGA-KIRC` = "normal rows of cache/gdc_file_map.rds (local file map, not a live GDC query)",
               `CPTAC-3`   = "normal rows of cache/valid_file_map.rds (local file map, not a live GDC query)",
               "live GDC /files query (cache/25_subtype_normal_files.rds)"),
             n_normal_files_on_disk = nrow(q),
             n_failed_rule_as_written = sum(q$failed_library), n_failed_reads = sum(q$failed_reads),
             n_failed_noFeature = sum(q$failed_noFeature), n_excluded = sum(q$excluded),
             exclusion_rule = if (coh == "CPTAC-3") "assigned reads only (08 applies no rule; protocol above the non-feature cap)" else "as 01 (reads and non-feature cap)",
             n_normal_analysed = sum(!q$excluded),
             n_normal_patient_with_analysed_tumour = sum(nl[excluded == FALSE, tumour_in_analysed_cohort]),
             n_tumour_analysed = lib[cohort == coh & tissue == "tumour", .N],
             n_pairs = nrow(pr), n_pairs_with_survival = if (nrow(pr)) sum(is.finite(pr$os_time) & !is.na(pr$os_event)) else 0L,
             deaths_in_pairs = if (nrow(pr)) sum(pr$os_event, na.rm = TRUE) else 0L,
             n_normal_with_plate = sum(!is.na(nl$plate)),
             n_pairs_both_plates_known = both_plate, n_pairs_same_plate = same_plate,
             frac_pairs_same_plate = if (is.na(both_plate) || both_plate == 0) NA_real_ else round(same_plate / both_plate, 3))
}))
save_tsv(consort, "25_normal_consort.tsv"); print(consort, row.names = FALSE)

# =============================================================================
# 4. Survival: the normal library's non-feature fraction as exposure
# =============================================================================
banner("4 | Overall survival against the NORMAL library's metric (paired TCGA patients)")
# One complete-case set per cohort (survival plus the clinical terms), so
# that the alone and adjusted rows describe the same patients. Hazard ratios
# per SD of the exposure, SD taken over that set. Events per parameter below
# MIN_EPV is flagged as underpowered.
surv_rows <- rbindlist(lapply(pair_cohorts, function(coh) {
  d <- pairs[cohort == coh]
  need <- c("os_time", "os_event", "age", "male", "T_stage", "N_pos", "M1", "grade_num",
            "pct_noFeature_t", "pct_noFeature_n", "lnc_axis_proj_t", "lnc_axis_proj_n")
  d <- d[complete.cases(d[, ..need]) & os_time > 0]
  if (nrow(d) < 10 || sum(d$os_event) < 2) { msg(coh, ": too few paired patients with survival for a Cox model"); return(NULL) }
  d[, `:=`(nf_n = as.numeric(scale(pct_noFeature_n)), nf_t = as.numeric(scale(pct_noFeature_t)),
           nflog_n = as.numeric(scale(log10(pct_noFeature_n))), nflog_t = as.numeric(scale(log10(pct_noFeature_t))),
           ax_n = as.numeric(scale(lnc_axis_proj_n)), ax_t = as.numeric(scale(lnc_axis_proj_t)),
           age_z = as.numeric(scale(age)))]
  clin <- "age_z + male + T_stage + N_pos + M1 + grade_num"
  msg(coh, ": complete cases n = ", nrow(d), ", events = ", sum(d$os_event))
  fit_row <- function(label, rhs, term, lib_label, exposure, scale_label) {
    fit <- tryCatch(coxph(as.formula(paste("Surv(os_time, os_event) ~", rhs)), data = d), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    s <- summary(fit)
    zp <- tryCatch(cox.zph(fit)$table[term, "p"], error = function(e) NA_real_)
    npar <- length(coef(fit))
    data.table(cohort = coh, exposure_library = lib_label, exposure = exposure, exposure_scale = scale_label,
               model = label, term = term, n = s$n, events = s$nevent, n_parameters = npar,
               epv = round(s$nevent / npar, 1), underpowered = s$nevent / npar < MIN_EPV,
               HR = round(s$conf.int[term, 1], 3), lo = round(s$conf.int[term, 3], 3),
               hi = round(s$conf.int[term, 4], 3), p = signif(s$coefficients[term, 5], 3),
               C = round(unname(s$concordance[1]), 3), ph_p = signif(zp, 3))
  }
  spec <- list(
    list("nf_n",    "normal", "pct_noFeature", "linear, per SD"),
    list("nf_t",    "tumour", "pct_noFeature", "linear, per SD"),
    list("nflog_n", "normal", "pct_noFeature", "log10, per SD"),
    list("nflog_t", "tumour", "pct_noFeature", "log10, per SD"),
    list("ax_n",    "normal", "lnc_axis_proj", "per SD"),
    list("ax_t",    "tumour", "lnc_axis_proj", "per SD"))
  out <- rbindlist(lapply(spec, function(z) rbind(
    fit_row("exposure alone", z[[1]], z[[1]], z[[2]], z[[3]], z[[4]]),
    fit_row("+ age, sex, T, N, M1, ordinal grade", paste(z[[1]], "+", clin), z[[1]], z[[2]], z[[3]], z[[4]]))))
  # Both libraries' metrics in one model, each adjusted for the other.
  joint <- rbind(
    fit_row("tumour and normal metric jointly + clinical terms", paste("nf_n + nf_t +", clin), "nf_n", "normal", "pct_noFeature", "linear, per SD"),
    fit_row("tumour and normal metric jointly + clinical terms", paste("nf_n + nf_t +", clin), "nf_t", "tumour", "pct_noFeature", "linear, per SD"),
    fit_row("tumour and normal axis jointly + clinical terms", paste("ax_n + ax_t +", clin), "ax_n", "normal", "lnc_axis_proj", "per SD"),
    fit_row("tumour and normal axis jointly + clinical terms", paste("ax_n + ax_t +", clin), "ax_t", "tumour", "lnc_axis_proj", "per SD"))
  rbind(out, joint)
}), fill = TRUE)
if (nrow(surv_rows)) {
  # Flag models whose cox.zph test rejects proportional hazards for the
  # exposure term.
  surv_rows[, ph_violation := is.finite(ph_p) & ph_p < 0.05]
  surv_rows[, fdr := signif(bh(p), 3), by = cohort]
  n_ph <- sum(surv_rows$ph_violation)
  surv_rows[, note := paste0(
    "paired patients only; hazard ratio per SD of the exposure over the complete-case set; ",
    "underpowered = events per parameter below MIN_EPV; ",
    "fdr = Benjamini-Hochberg within cohort over the Cox p-values of this table; ",
    "ph_violation = cox.zph p < 0.05 for the exposure term (", n_ph, " of ", nrow(surv_rows),
    " rows here)")]
  save_tsv(surv_rows, "25_normal_metric_survival.tsv")
  print(surv_rows[, .(cohort, exposure_library, exposure, exposure_scale, model, n, events, epv,
                      underpowered, HR, lo, hi, p, fdr, ph_p, ph_violation)], row.names = FALSE)
  if (n_ph) print(surv_rows[ph_violation == TRUE, .(cohort, exposure_library, exposure, model, HR, p, ph_p)],
                  row.names = FALSE)
} else {
  save_tsv(data.table(cohort = character(0), note = character(0)), "25_normal_metric_survival.tsv")
  msg("No cohort had enough paired patients with survival for a Cox model")
}

# =============================================================================
# 5. Tumour versus normal distributions per cohort (unpaired)
# =============================================================================
banner("5 | Tumour versus normal distributions per cohort")
DIST_VARS <- c(pct_noFeature = "non-feature read fraction (%)", pct_multimapping = "multimapping read fraction (%)",
               pct_unmapped = "unmapped read fraction (%)", pct_ambiguous = "ambiguous read fraction (%)",
               assigned_reads = "assigned reads", log_depth = "log10 assigned reads",
               lnc_axis_proj = "lncRNA axis, projected on discovery PC1 (per discovery SD)",
               pc_axis_proj = "protein-coding PC1, projected (per discovery SD; INTENDED negative control, not a clean one: in the discovery tumours this projection is largely the multimapping axis)",
               mean_lnc_expr = "mean lncRNA log2(FPKM + 1) over the 3,442 network transcripts")
# The tumour side is the analysed cohort (in CPTAC-3 the clear-cell patients of
# 08). "All normal libraries" includes normals whose tumour 08 excluded, so the
# set restricted to analysed-tumour patients is reported alongside.
DIST_NORMAL_SETS <- c(
  "all normal libraries" = "every analysed normal library of the cohort, whether or not that patient's tumour is in the analysed cohort",
  "normals of analysed-tumour patients" = "normal libraries whose patient's tumour library is in the analysed cohort: the same patient population as the tumour side")
# Tumours and normals do not occupy the same plates, so each unpaired row also
# carries (a) the matched-pair signed-rank contrast and (b) the unpaired
# contrast on plates holding both tumour and normal libraries. The note quotes
# the tumour-only plate R2 from the variance partition.
dist_tbl <- rbindlist(lapply(pair_cohorts, function(coh) rbindlist(lapply(names(DIST_NORMAL_SETS), function(nset) {
  T_all <- lib[cohort == coh & tissue == "tumour" & excluded == FALSE]
  N_all <- lib[cohort == coh & tissue == "normal" & excluded == FALSE]
  if (nset != "all normal libraries") N_all <- N_all[tumour_in_analysed_cohort == TRUE]
  pr <- pairs[cohort == coh]
  shared_plates <- if (startsWith(coh, "TCGA"))
    intersect(unique(na.omit(N_all$plate)), unique(na.omit(T_all$plate))) else character(0)
  n_plate_t <- uniqueN(na.omit(T_all$plate)); n_plate_n <- uniqueN(na.omit(N_all$plate))
  rbindlist(lapply(names(DIST_VARS), function(v) {
    tv <- T_all[[v]]; nv <- N_all[[v]]
    tv <- tv[is.finite(tv)]; nv <- nv[is.finite(nv)]
    if (length(tv) < 5 || length(nv) < 5) return(NULL)
    w <- suppressWarnings(wilcox.test(tv, nv, exact = FALSE))
    qt <- qtl(tv); qn <- qtl(nv)
    # (a) the same contrast within patient
    ct <- paste0(v, "_t"); cn <- paste0(v, "_n")
    mn <- NA_integer_; mdiff <- NA_real_; mp <- NA_real_; mrho <- NA_real_
    if (all(c(ct, cn) %in% names(pr))) {
      a <- pr[[ct]]; b <- pr[[cn]]; ok <- is.finite(a) & is.finite(b)
      mn <- sum(ok)
      if (mn >= 5) {
        mdiff <- median(a[ok] - b[ok])
        mp <- suppressWarnings(wilcox.test(a[ok], b[ok], paired = TRUE, exact = FALSE))$p.value
        mrho <- spearman(a[ok], b[ok])[["rho"]]
      }
    }
    # (b) the same unpaired contrast restricted to plates shared by tumours and normals
    pn_t <- NA_integer_; pn_n <- NA_integer_; pmed_t <- NA_real_; pmed_n <- NA_real_; pp <- NA_real_
    if (length(shared_plates)) {
      a <- T_all[plate %in% shared_plates][[v]]; b <- N_all[plate %in% shared_plates][[v]]
      a <- a[is.finite(a)]; b <- b[is.finite(b)]
      pn_t <- length(a); pn_n <- length(b)
      if (pn_t >= 5 && pn_n >= 5) {
        pmed_t <- median(a); pmed_n <- median(b)
        pp <- suppressWarnings(wilcox.test(a, b, exact = FALSE))$p.value
      }
    }
    data.table(cohort = coh, normal_set = nset, normal_set_definition = DIST_NORMAL_SETS[[nset]],
               tumour_set_definition = "every analysed tumour library of the cohort",
               variable = v, description = DIST_VARS[[v]],
               n_tumour = length(tv), median_tumour = round(qt[["median"]], 4), q25_tumour = round(qt[["q25"]], 4), q75_tumour = round(qt[["q75"]], 4),
               n_normal = length(nv), median_normal = round(qn[["median"]], 4), q25_normal = round(qn[["q25"]], 4), q75_normal = round(qn[["q75"]], 4),
               median_diff_t_minus_n = round(qt[["median"]] - qn[["median"]], 4),
               wilcoxon_rank_sum_p = signif(w$p.value, 3),
               matched_n_pairs = as.integer(mn),
               matched_median_diff_t_minus_n = round(mdiff, 4),
               matched_signed_rank_p = signif(mp, 3),
               matched_rho_tumour_vs_normal = round(mrho, 3),
               n_plates_tumour = as.integer(n_plate_t), n_plates_normal = as.integer(n_plate_n),
               plate_matched_n_tumour = as.integer(pn_t), plate_matched_n_normal = as.integer(pn_n),
               plate_matched_median_tumour = round(pmed_t, 4), plate_matched_median_normal = round(pmed_n, 4),
               plate_matched_median_diff_t_minus_n = round(pmed_t - pmed_n, 4),
               plate_matched_rank_sum_p = signif(pp, 3))
  }))
}))))
stopifnot(nrow(dist_tbl) > 0)
# log10 assigned reads is a monotone transform of assigned reads, so the
# rank-sum test is identical: log_depth is left out of the BH family and
# carries the assigned_reads FDR.
DIST_FDR_FAMILY <- setdiff(names(DIST_VARS), "log_depth")
dist_tbl[, in_fdr_family := variable %in% DIST_FDR_FAMILY]
dist_tbl[, fdr := {
  f <- rep(NA_real_, .N); ok <- in_fdr_family
  if (any(ok)) f[ok] <- p.adjust(wilcoxon_rank_sum_p[ok], "BH")
  ar <- f[variable == "assigned_reads"]
  f[variable == "log_depth"] <- if (length(ar) == 1L) ar else NA_real_
  signif(f, 3)
}, by = .(cohort, normal_set)]
dist_tbl[, fdr_family := fifelse(in_fdr_family,
  paste0("Benjamini-Hochberg within cohort x normal set over the ", length(DIST_FDR_FAMILY),
         " distinct rank-sum tests (log_depth excluded as a monotone duplicate of assigned_reads)"),
  paste0("carried unchanged from assigned_reads: the same rank-sum test on a monotone transform, ",
         "excluded from the BH family. This applies to the UNPAIRED column only; matched_fdr on this ",
         "row is that row's own signed-rank test, which is not transform-invariant"))]
# The matched-pair tests have their own BH family. log_depth enters it
# separately, because the signed-rank test ranks absolute differences, which
# log10 changes.
dist_tbl[, in_matched_fdr_family := is.finite(matched_signed_rank_p)]
dist_tbl[, matched_fdr := {
  f <- rep(NA_real_, .N); ok <- in_matched_fdr_family
  if (any(ok)) f[ok] <- p.adjust(matched_signed_rank_p[ok], "BH")
  signif(f, 3)
}, by = .(cohort, normal_set)]
dist_tbl[, matched_fdr_family := paste0(
  "matched_signed_rank_p adjusted in its OWN Benjamini-Hochberg family, within cohort x normal set, ",
  "over every quantity for which a signed-rank test was performed (", length(DIST_VARS),
  " here, one MORE than the ", length(DIST_FDR_FAMILY), " of the unpaired family). assigned_reads ",
  "and log_depth share one unpaired row because the rank-sum test is invariant under a monotone ",
  "transform, but they are two DIFFERENT signed-rank tests: that test ranks absolute paired ",
  "differences, which log10 does not preserve, so log_depth carries its own matched_fdr and is not ",
  "copied from assigned_reads. Only the unpaired column is transform-invariant. This family is ",
  "therefore not comparable with fdr_wilcoxon of 25_paired_tumour_normal.tsv, which adjusts over the ",
  length(PAIR_VARS), " variables of that table")]
# Whether the two sides of the unpaired contrast occupy the same plates. Where
# they do not, the unpaired difference mixes tissue with plate and the matched
# columns are the interpretable ones.
dist_tbl[, unpaired_plate_confounded := is.finite(n_plates_normal) & is.finite(n_plates_tumour) &
           n_plates_normal > 0 & n_plates_normal < n_plates_tumour]
# The note quotes this variable's tumour-only plate R2, since the pooled R2 is
# inflated by the normals' concentration on few plates. Variables without a
# plate partition are pointed at the variance partition table.
dist_tbl[, plate_R2_tumour_only := NA_real_]
plate_r2_tum <- vp[factor == "plate" & library_set == "analysed tumour libraries only" & is.finite(R2),
                   .(cohort, variable, R2_t = R2)]
if (nrow(plate_r2_tum)) dist_tbl[plate_r2_tum, on = .(cohort, variable), plate_R2_tumour_only := i.R2_t]
dist_tbl[, plate_effect_phrase := fifelse(
  is.finite(plate_R2_tumour_only),
  paste0("plate explains R2 ", plate_R2_tumour_only, " of this variable within the analysed TUMOUR ",
         "libraries (25_paired_variance_partition.tsv, library_set 'analysed tumour libraries only'; ",
         "the pooled figure over tumours and normals together is larger only because the normals are ",
         "concentrated on a few plates)"),
  paste0("plate is not partitioned for this variable (25_paired_variance_partition.tsv gives the ",
         "tumour-side plate partition of the non-feature fraction, the projected lncRNA axis and the ",
         "multimapping fraction)"))]
# Three cases: normals on fewer plates than tumours (confounded), normals on at
# least as many plates, and plate not recorded.
dist_tbl[, note := fifelse(unpaired_plate_confounded,
  paste0("DO NOT QUOTE THE UNPAIRED CONTRAST ALONE. The normal libraries span ", n_plates_normal,
         " sequencing plates against ", n_plates_tumour, " for the tumours, and ", plate_effect_phrase,
         ", so this unpaired difference mixes tissue type with plate. Beside it: within patient the ",
         "median tumour-minus-normal difference is ", matched_median_diff_t_minus_n, " on ",
         matched_n_pairs, " pairs (signed-rank p ", matched_signed_rank_p, ", fdr ", matched_fdr,
         "); restricted to the ", plate_matched_n_tumour, " tumour and ", plate_matched_n_normal,
         " normal libraries on the plates both sides occupy, the difference is ",
         plate_matched_median_diff_t_minus_n, " (rank-sum p ", plate_matched_rank_sum_p, ")."),
  fifelse(is.finite(n_plates_normal) & n_plates_normal > 0,
  paste0("Both sides carry a plate: the normal libraries span ", n_plates_normal,
         " sequencing plates against ", n_plates_tumour, " for the tumours, so the unpaired contrast ",
         "is not restricted to a subset of the tumours' plates in the way it is where the normals sit ",
         "on fewer. The plate-restricted columns give the contrast on the plates both sides occupy (",
         plate_matched_n_tumour, " tumour and ", plate_matched_n_normal, " normal libraries, ",
         "difference ", plate_matched_median_diff_t_minus_n, ", rank-sum p ", plate_matched_rank_sum_p,
         "); within patient the difference is ", matched_median_diff_t_minus_n, " on ", matched_n_pairs,
         " pairs (signed-rank p ", matched_signed_rank_p, ", fdr ", matched_fdr, ")."),
  paste0("Plate is not recorded for any library of this cohort, so plate confounding cannot be ruled ",
         "in or out here and the plate-restricted columns are empty; the matched-pair columns give the ",
         "same contrast within patient on ", matched_n_pairs, " pairs (median difference ",
         matched_median_diff_t_minus_n, ", signed-rank p ", matched_signed_rank_p, ", fdr ",
         matched_fdr, ").")))]
dist_tbl[, plate_effect_phrase := NULL]   # folded into note, value kept in plate_R2_tumour_only
save_tsv(dist_tbl, "25_tumour_vs_normal_distribution.tsv")
print(dist_tbl[, .(cohort, normal_set, variable, n_tumour, median_tumour, n_normal, median_normal,
                   median_diff_t_minus_n, wilcoxon_rank_sum_p, fdr,
                   matched_median_diff_t_minus_n, matched_signed_rank_p, matched_fdr,
                   plate_matched_median_diff_t_minus_n, plate_matched_rank_sum_p)], row.names = FALSE)
for (coh in pair_cohorts) {
  r <- dist_tbl[cohort == coh & variable == "pct_noFeature" &
                normal_set == "normals of analysed-tumour patients"]
  if (nrow(r)) msg(coh, ": non-feature fraction, tumour versus normal -- UNPAIRED ", r$median_tumour,
                   " vs ", r$median_normal, " (difference ", r$median_diff_t_minus_n, ", rank-sum p ",
                   r$wilcoxon_rank_sum_p, ", fdr ", r$fdr, "); WITHIN PATIENT difference ",
                   r$matched_median_diff_t_minus_n, " on ", r$matched_n_pairs, " pairs (signed-rank p ",
                   r$matched_signed_rank_p, ", fdr ", r$matched_fdr, ")",
                   if (isTRUE(r$unpaired_plate_confounded))
                     paste0("; ON SHARED PLATES ONLY difference ", r$plate_matched_median_diff_t_minus_n,
                            " (rank-sum p ", r$plate_matched_rank_sum_p, ", ", r$plate_matched_n_tumour,
                            " tumour and ", r$plate_matched_n_normal, " normal libraries). The unpaired ",
                            "contrast is a plate artefact and is not interpretable on its own.") else ".")
}

# ---- per-library scores for figures ------------------------------------------
per_sample <- lib[, .(cohort, tissue, patient, sample_barcode, file_id, paired, excluded,
                      tumour_in_analysed_cohort, in_discovery_pca,
                      pct_noFeature = round(pct_noFeature, 4), pct_multimapping = round(pct_multimapping, 4),
                      assigned_reads, log_depth = round(log_depth, 4),
                      lnc_axis_proj = round(lnc_axis_proj, 4), pc_axis_proj = round(pc_axis_proj, 4),
                      lnc_own_pc1_z = round(lnc_own_pc1_z, 4), pc_own_pc1_z = round(pc_own_pc1_z, 4),
                      mean_lnc_expr = round(mean_lnc_expr, 4), tss, plate, plate_source,
                      os_time, os_event)]
setorder(per_sample, cohort, patient, tissue)
save_tsv(per_sample, "25_paired_scores_per_sample.tsv")

note_file <- file.path(RESULTS_DIR, "25_NOTE_gdc.txt")
if (length(gdc_notes)) {
  writeLines(c("25_matched_normal_control.R: one or more GDC API calls failed; the affected",
               "quantities fell back to caches or were recorded as unknown (see the tables).",
               paste0("  - ", gdc_notes), paste0("Time: ", format(Sys.time()))), note_file)
} else unlink(note_file)

msg("Elapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
write_session_info("25_matched_normal_control")
banner("25 | done")
