# lib_sc.R: shared functions for the single-cell pipeline.
# Apart from the readers and the downloader, these functions have no side
# effects. They enforce integer counts, unique cell IDs, donor identity, count
# preservation under aggregation and the absence of outcome fields.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- manifest ----------------------------------------------------------------
MANIFEST_REQUIRED <- c(
  "tisch_id", "dataset_name", "geo_accession", "pmid", "species", "platform",
  "reported_patients", "reported_cells", "treatment", "primary_metastatic",
  "download_source", "expression_type", "integer_counts_available",
  "metadata_available", "normal_tissue_available", "analysis_role",
  "local_expression_path", "local_metadata_path", "sha256", "notes")

SAMPLES_REQUIRED <- c("dataset_name", "gsm", "sample_id", "donor_id", "tissue",
                      "sample_type", "include", "file_prefix", "format", "notes")

ANALYSIS_ROLES <- c("discovery", "replication", "replication_candidate",
                    "treated_sensitivity", "excluded")
TISSUES <- c("tumour", "normal_kidney", "blood", "metastasis", "other")

# Returns a character vector of problems (empty means valid).
validate_manifest <- function(man, samples = NULL) {
  p <- character()
  miss <- setdiff(MANIFEST_REQUIRED, names(man))
  if (length(miss)) p <- c(p, paste("manifest missing columns:", paste(miss, collapse = ", ")))
  if (anyDuplicated(man$dataset_name))
    p <- c(p, paste("duplicated dataset_name:", paste(man$dataset_name[duplicated(man$dataset_name)], collapse = ", ")))
  bad_role <- setdiff(man$analysis_role, ANALYSIS_ROLES)
  if (length(bad_role)) p <- c(p, paste("unknown analysis_role:", paste(bad_role, collapse = ", ")))
  for (col in c("dataset_name", "geo_accession", "analysis_role", "expression_type"))
    if (col %in% names(man) && any(is.na(man[[col]]) | man[[col]] == ""))
      p <- c(p, paste("empty required field:", col))
  if (!is.null(samples)) {
    miss <- setdiff(SAMPLES_REQUIRED, names(samples))
    if (length(miss)) p <- c(p, paste("sample sheet missing columns:", paste(miss, collapse = ", ")))
    orphan <- setdiff(unique(samples$dataset_name), man$dataset_name)
    if (length(orphan)) p <- c(p, paste("samples reference unknown datasets:", paste(orphan, collapse = ", ")))
    key <- paste(samples$dataset_name, samples$sample_id)
    if (anyDuplicated(key)) p <- c(p, paste("duplicated sample_id within dataset:", paste(key[duplicated(key)], collapse = "; ")))
    if (any(is.na(samples$donor_id) | samples$donor_id == ""))
      p <- c(p, "sample sheet has samples without donor_id")
    bad_t <- setdiff(samples$tissue, TISSUES)
    if (length(bad_t)) p <- c(p, paste("unknown tissue values:", paste(bad_t, collapse = ", ")))
  }
  p
}

# ---- checksums and downloads -------------------------------------------------
sha256_file <- function(path) {
  stopifnot(file.exists(path))
  digest::digest(file = path, algo = "sha256")
}

geo_sample_url <- function(gsm, filename) {
  stem <- paste0(substr(gsm, 1, nchar(gsm) - 3), "nnn")
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/%s", stem, gsm, filename)
}
geo_series_url <- function(gse, filename) {
  stem <- paste0(substr(gse, 1, nchar(gse) - 3), "nnn")
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s", stem, gse, filename)
}

# File names per sample for each supported raw format.
sample_files <- function(prefix, format) {
  switch(format,
    "10x_h5"     = paste0(prefix, ".h5"),
    "10x_mtx_v3" = paste0(prefix, c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")),
    "10x_mtx_v2" = paste0(prefix, c("barcodes.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz")),
    "dense_txt"  = paste0(prefix, ".txt.gz"),
    stop("unknown sample format: ", format))
}

remote_size <- function(url, ua) {
  h <- curl::new_handle(nobody = TRUE, followlocation = TRUE, useragent = ua)
  r <- tryCatch(curl::curl_fetch_memory(url, handle = h), error = function(e) NULL)
  if (is.null(r) || r$status_code >= 400) return(NA_real_)
  hd <- curl::parse_headers_list(r$headers)
  as.numeric(hd[["content-length"]] %||% NA)
}

# Resumable download with retries. A partial file is kept as <dest>.part and
# renamed into place only once its size matches the server's Content-Length,
# when the server reports one.
download_resumable <- function(url, dest, ua, retries = 4, wait = 20, timeout = 3600) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  expected <- remote_size(url, ua)
  if (file.exists(dest) && (is.na(expected) || file.size(dest) == expected))
    return(list(status = "present", bytes = file.size(dest), expected = expected))
  part <- paste0(dest, ".part")
  if (file.exists(part) && !is.na(expected) && file.size(part) == expected) {
    file.rename(part, dest)
    return(list(status = "downloaded", bytes = file.size(dest), expected = expected))
  }
  for (i in seq_len(retries)) {
    have <- if (file.exists(part)) file.size(part) else 0
    h <- curl::new_handle(useragent = ua, followlocation = TRUE, timeout = timeout,
                          failonerror = TRUE)
    if (have > 0) curl::handle_setopt(h, resume_from = have)
    con <- file(part, open = if (have > 0) "ab" else "wb")
    ok <- tryCatch({
      curl::curl_fetch_stream(url, function(x) writeBin(x, con), handle = h)
      TRUE
    }, error = function(e) { message("  attempt ", i, " failed: ", conditionMessage(e)); FALSE },
    finally = close(con))
    if (ok && (is.na(expected) || file.size(part) == expected)) {
      file.rename(part, dest)
      return(list(status = "downloaded", bytes = file.size(dest), expected = expected))
    }
    if (ok && !is.na(expected) && file.size(part) > expected) unlink(part)
    Sys.sleep(wait)
  }
  list(status = "failed", bytes = if (file.exists(part)) file.size(part) else 0, expected = expected)
}

# ---- identifiers -------------------------------------------------------------
strip_ensembl_version <- function(ids) sub("^(ENS[A-Z]*[GT][0-9]+)\\.[0-9]+(_PAR_Y)?$", "\\1\\2", ids)

# Rows whose unversioned ID collides with another row (e.g. PAR_Y copies, or a
# matrix carrying two versions of one gene). These are reported, never merged.
ensembl_collisions <- function(ids) {
  k <- strip_ensembl_version(ids)
  dup <- k %in% k[duplicated(k)]
  data.table(original_id = ids[dup], unversioned = k[dup])
}

make_cell_ids <- function(dataset, sample, barcode) {
  ids <- paste(dataset, sample, barcode, sep = "|")
  if (anyDuplicated(ids)) stop("duplicated cell identifiers after prefixing: ",
                               paste(head(ids[duplicated(ids)], 3), collapse = ", "))
  ids
}

# ---- count-type guards -------------------------------------------------------
is_integer_counts <- function(m, tol = 1e-8) {
  x <- if (inherits(m, "sparseMatrix")) m@x else as.vector(m)
  x <- x[is.finite(x)]
  length(x) == 0 || (all(x >= 0) && all(abs(x - round(x)) < tol))
}

require_raw_counts <- function(m, what = "matrix") {
  if (!is_integer_counts(m))
    stop(what, " contains non-integer or negative values: count-based methods ",
         "(pseudobulk sums, edgeR) refuse log-normalised or scaled input. ",
         "Use the raw UMI count assay.", call. = FALSE)
  invisible(TRUE)
}

# ---- readers -------------------------------------------------------------------
# Cell Ranger h5: v2 stores one group per genome with gene ids/names, v3 stores
# /matrix with /matrix/features/{id,name,feature_type}.
read_10x_h5 <- function(path) {
  ls <- rhdf5::h5ls(path)
  grp <- if (any(ls$group == "/matrix")) "/matrix" else paste0("/", ls$name[ls$group == "/" & ls$otype == "H5I_GROUP"][1])
  rd <- function(x) as.vector(rhdf5::h5read(path, paste0(grp, "/", x)))
  shape <- rd("shape")
  m <- sparseMatrix(i = rd("indices") + 1L, p = as.integer(rd("indptr")),
                    x = as.numeric(rd("data")), dims = shape)
  if (grp == "/matrix") {
    ids <- rd("features/id"); nm <- rd("features/name")
    ft <- tryCatch(rd("features/feature_type"), error = function(e) rep("Gene Expression", length(ids)))
    keep <- ft == "Gene Expression"
    m <- m[keep, , drop = FALSE]; ids <- ids[keep]; nm <- nm[keep]
  } else { ids <- rd("genes"); nm <- rd("gene_names") }
  rhdf5::h5closeAll()
  dimnames(m) <- list(ids, rd("barcodes"))
  list(counts = m, genes = data.table(gene_id = ids, gene_symbol = nm), h5_group = grp)
}

read_10x_mtx <- function(dir, prefix) {
  f <- function(x) file.path(dir, paste0(prefix, x))
  gene_file <- if (file.exists(f("features.tsv.gz"))) f("features.tsv.gz") else f("genes.tsv.gz")
  m <- as(Matrix::readMM(f("matrix.mtx.gz")), "CsparseMatrix")
  g <- fread(gene_file, header = FALSE)
  b <- fread(f("barcodes.tsv.gz"), header = FALSE)[[1]]
  if (ncol(g) >= 3) g <- g[V3 == "Gene Expression" | is.na(V3)]
  if (nrow(g) != nrow(m) || length(b) != ncol(m))
    stop("mtx dimensions (", nrow(m), "x", ncol(m), ") do not match gene (", nrow(g),
         ") / barcode (", length(b), ") files for ", prefix, call. = FALSE)
  dimnames(m) <- list(g$V1, b)
  list(counts = m, genes = data.table(gene_id = g$V1, gene_symbol = g$V2))
}

# Dense genes-by-cells text (first column gene symbol), read in one pass and
# converted to sparse. The files used are under 15 MB compressed.
read_dense_txt <- function(path) {
  d <- fread(path)
  genes <- d[[1]]; d[[1]] <- NULL
  m <- Matrix(as.matrix(d), sparse = TRUE)
  rownames(m) <- genes
  list(counts = m, genes = data.table(gene_id = NA_character_, gene_symbol = genes))
}

# AnnData .h5ad read without Python. Supports CSR/CSC X or a named layer and
# both the pre-0.8 (__categories) and >= 0.8 (categories/codes) obs encodings.
h5ad_read_obs <- function(path) {
  ls <- rhdf5::h5ls(path, recursive = 2)
  cols <- ls$name[ls$group == "/obs" & !ls$name %in% c("__categories", "_index", "index")]
  idx_name <- tryCatch(rhdf5::h5readAttributes(path, "/obs")[["_index"]], error = function(e) "_index")
  obs <- data.table(obs_names = as.vector(rhdf5::h5read(path, paste0("/obs/", idx_name))))
  has_old_cats <- any(ls$group == "/obs" & ls$name == "__categories")
  for (cn in cols) {
    node <- paste0("/obs/", cn)
    is_grp <- ls$otype[ls$group == "/obs" & ls$name == cn] == "H5I_GROUP"
    if (is_grp) {
      codes <- as.vector(rhdf5::h5read(path, paste0(node, "/codes")))
      cats <- as.vector(rhdf5::h5read(path, paste0(node, "/categories")))
      obs[[cn]] <- ifelse(codes < 0, NA, cats[codes + 1L])
    } else {
      v <- as.vector(rhdf5::h5read(path, node))
      cats <- if (has_old_cats) tryCatch(as.vector(rhdf5::h5read(path, paste0("/obs/__categories/", cn))),
                                         error = function(e) NULL) else NULL
      if (!is.null(cats)) v <- ifelse(v < 0, NA, cats[v + 1L])
      obs[[cn]] <- v
    }
  }
  rhdf5::h5closeAll()
  obs
}

h5ad_read_counts <- function(path, matrix_path = "/X") {
  a <- rhdf5::h5readAttributes(path, matrix_path)
  enc <- a[["encoding-type"]] %||% a[["h5sparse_format"]] %||% "csr_matrix"
  shape <- as.integer(a[["shape"]] %||% a[["h5sparse_shape"]])
  rd <- function(x) as.vector(rhdf5::h5read(path, paste0(matrix_path, "/", x)))
  var_idx <- tryCatch(rhdf5::h5readAttributes(path, "/var")[["_index"]], error = function(e) "_index")
  genes <- as.vector(rhdf5::h5read(path, paste0("/var/", var_idx)))
  obs_names <- h5ad_read_obs(path)$obs_names
  # AnnData is cells x genes, so CSR over cells == CSC over the transposed matrix.
  m <- if (grepl("csr", enc)) {
    sparseMatrix(i = rd("indices") + 1L, p = as.integer(rd("indptr")), x = as.numeric(rd("data")),
                 dims = rev(shape))
  } else {
    t(sparseMatrix(i = rd("indices") + 1L, p = as.integer(rd("indptr")), x = as.numeric(rd("data")),
                   dims = shape))
  }
  rhdf5::h5closeAll()
  dimnames(m) <- list(genes, obs_names)
  m
}

# ---- metadata guards ---------------------------------------------------------
# Columns that carry clinical outcome or disease severity. Candidate discovery
# must never see them. Harmonised discovery metadata is checked against this.
OUTCOME_PATTERN <- "(surviv|os_|pfs|dfs|death|dead|vital|relapse|recur|progress|stage|grade|tnm|response|responder|follow)"

drop_outcome_fields <- function(dt) {
  hit <- grep(OUTCOME_PATTERN, names(dt), ignore.case = TRUE, value = TRUE)
  if (length(hit)) dt <- dt[, setdiff(names(dt), hit), with = FALSE]
  attr(dt, "dropped_outcome_fields") <- hit
  dt
}

assert_no_outcome_fields <- function(dt, context = "candidate discovery") {
  hit <- grep(OUTCOME_PATTERN, names(dt), ignore.case = TRUE, value = TRUE)
  if (length(hit)) stop(context, " received outcome-related fields: ",
                        paste(hit, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

assert_cells_match <- function(expr_cells, meta_cells, min_frac = 0.9, what = "") {
  hit <- mean(meta_cells %in% expr_cells)
  if (!is.finite(hit) || hit < min_frac)
    stop(sprintf("%s: only %.1f%% of annotated cells are present in the expression matrix (threshold %.0f%%). Check barcode suffixes and sample mapping.",
                 what, 100 * hit, 100 * min_frac), call. = FALSE)
  hit
}

assert_donors <- function(donor, what = "") {
  if (!length(donor) || any(is.na(donor) | donor == ""))
    stop(what, ": donor identity is absent for some cells", call. = FALSE)
  if (any(grepl(",", donor)))
    stop(what, ": donor field lists several patients (", unique(donor[grepl(",", donor)])[1],
         "); donor identity is unusable", call. = FALSE)
  invisible(TRUE)
}

# ---- annotation ----------------------------------------------------------------
read_gencode_genes <- function(gtf) {
  x <- fread(cmd = NULL, file = gtf, sep = "\t", header = FALSE, skip = "chr",
             select = c(1, 3, 4, 5, 7, 9), col.names = c("chr", "feature", "start", "end", "strand", "attr"))
  x <- x[feature == "gene"]
  ga <- function(key) sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x$attr)
  data.table(gene_id = ga("gene_id"), gene_name = ga("gene_name"), gene_type = ga("gene_type"),
             chr = x$chr, start = x$start, end = x$end, strand = x$strand)
}

is_lncrna <- function(gene_type, biotypes) !is.na(gene_type) & gene_type %in% biotypes

# Map a matrix's genes onto the reference. Ensembl IDs are matched unversioned.
# Symbols are matched only when unique in the reference. Ambiguous or unknown
# symbols stay unmapped and are counted.
#
# `legacy` (optional) is an older GENCODE gene table, e.g. v19/GRCh37. Symbols
# still unmapped against the reference are looked up there when unique, and the
# Ensembl gene ID is carried to the reference by its unversioned ID. Ensembl
# gene IDs are stable across builds, so no liftover is needed. A symbol whose
# legacy ID was retired stays unmapped.
map_genes <- function(genes, ref, legacy = NULL) {
  ref <- copy(ref)[, key := strip_ensembl_version(gene_id)]
  g <- copy(genes)[, row := .I]
  g[, key := ifelse(!is.na(gene_id) & grepl("^ENSG", gene_id), strip_ensembl_version(gene_id), NA_character_)]
  g <- merge(g, ref[, .(key, ref_gene_id = gene_id, ref_symbol = gene_name, gene_type)],
             by = "key", all.x = TRUE, sort = FALSE)
  setorder(g, row)
  g[, map_status := fifelse(!is.na(ref_gene_id), "ensembl", NA_character_)]
  sym_n <- ref[, .N, by = gene_name]
  uniq <- ref[gene_name %in% sym_n[N == 1, gene_name]]
  need <- is.na(g$ref_gene_id) & !is.na(g$gene_symbol)
  hit <- match(g$gene_symbol[need], uniq$gene_name)
  g$ref_gene_id[need] <- uniq$gene_id[hit]
  g$ref_symbol[need] <- uniq$gene_name[hit]
  g$gene_type[need] <- uniq$gene_type[hit]
  g[need & !is.na(ref_gene_id), map_status := "symbol_unique"]
  if (!is.null(legacy)) {
    lg <- copy(legacy)[, lkey := strip_ensembl_version(gene_id)]
    lg_n <- lg[, .N, by = gene_name]
    lg_u <- lg[gene_name %in% lg_n[N == 1, gene_name]]
    need2 <- is.na(g$ref_gene_id) & !is.na(g$gene_symbol) & !(g$gene_symbol %in% sym_n[N > 1, gene_name])
    lkey <- lg_u$lkey[match(g$gene_symbol[need2], lg_u$gene_name)]
    j <- match(lkey, ref$key)
    g$ref_gene_id[need2] <- ref$gene_id[j]
    g$ref_symbol[need2] <- ref$gene_name[j]
    g$gene_type[need2] <- ref$gene_type[j]
    g[need2 & !is.na(ref_gene_id), map_status := "symbol_legacy_bridge"]
    g[is.na(map_status) & !is.na(gene_symbol) & gene_symbol %in% lg_n[N > 1, gene_name] &
        !gene_symbol %in% ref$gene_name, map_status := "symbol_ambiguous"]
  }
  g[is.na(map_status) & !is.na(gene_symbol) & gene_symbol %in% sym_n[N > 1, gene_name], map_status := "symbol_ambiguous"]
  g[is.na(map_status), map_status := "unmapped"]
  g[, key := strip_ensembl_version(ref_gene_id)]
  g[]
}

# ---- compartments --------------------------------------------------------------
# Harmonise free-text labels into the configured compartments (vocabulary only).
# The caller gives author labels priority.
map_compartment <- function(label) {
  l <- tolower(trimws(as.character(label)))
  out <- rep("other_uncertain", length(l))
  rule <- function(pattern, value) out[grepl(pattern, l) & out == "other_uncertain"] <<- value
  rule("^(tumor|tumour|malignant|cancer|ccrcc|rcc|tumor cells|malignant cells)$|malignan", "malignant")
  rule("^(pt|pt-[abc]|tal|tal[12]?|tal-|tal |dct|cnt|pc|ic-a|ic-b|ic-pc|tal|dl|tal|podo|podocyte|pec|loh|ic|ua_epi)$|^epi[_ -]|^epith|^pt[_-]|proximal|tubul|epithel|podocyt|collecting|intercalated|principal|loop of henle|distal", "normal_epithelial")
  rule("t ?cell|t-cell|b-cell|b_cell|^mac|^t$|cd4|cd8|treg|tprolif|nk|b ?cell|^b$|plasma|macro|mono|mono/macro|dc|cdc|pdc|mast|neutro|immune|lymph|myeloid", "immune")
  rule("^ec$|_ec$|endo|avr|aea|dvr|gc$|^gc|glomerular capillar|vascular|plvap|ackr1", "endothelial")
  rule("fibro|stromal|peri|vsmc|smooth|mesang|myofib|caf", "fibroblast_stromal")
  out[l %in% c("", "na", "ua", "unknown", "unassigned", "doublet", "uc", "erythroblasts")] <- "other_uncertain"
  out
}

# ---- donor-level structure ---------------------------------------------------
donor_compartment_table <- function(cd, min_cells) {
  x <- as.data.table(cd)[, .(n_cells = .N), by = .(dataset, donor_id, tissue, compartment)]
  x[, evaluable := n_cells >= min_cells]
  x[]
}

# Sum raw counts by group. Groups must be donor-level (built from donor_id).
# Refuses non-integer input and preserves totals.
aggregate_pseudobulk <- function(counts, group) {
  require_raw_counts(counts, "pseudobulk input")
  if (length(group) != ncol(counts)) stop("group length must equal number of cells")
  if (anyNA(group)) stop("pseudobulk groups contain NA; every cell needs donor and compartment")
  f <- factor(group)
  mm <- sparseMatrix(i = as.integer(f), j = seq_along(f), x = 1,
                     dims = c(nlevels(f), length(f)), dimnames = list(levels(f), NULL))
  pb <- counts %*% t(mm)
  colnames(pb) <- levels(f)
  pb
}

# Refuses designs in which cells, rather than donors, would be the unit.
assert_donor_level <- function(sample_table) {
  st <- as.data.table(sample_table)
  if (!"donor_id" %in% names(st)) stop("design table has no donor_id column")
  key <- paste(st$donor_id, st$compartment)
  if (anyDuplicated(key))
    stop("more than one pseudobulk per donor x compartment: replicate would not be the donor", call. = FALSE)
  invisible(TRUE)
}

# Donor-level design for a two-group contrast (levels[2] vs levels[1]).
# Paired (~ donor + group, restricted to donors contributing both groups) when
# at least `min_paired` donors contribute both, otherwise unpaired (~ group) on
# all donors. Donors seen in one group carry no information about the group
# effect in a donor-blocked model, so dropping them loses nothing.
design_for_contrast <- function(sample_table, group_col = "compartment", levels, min_paired = 3) {
  st <- as.data.table(sample_table)[get(group_col) %in% levels]
  assert_donor_level(st)
  donors_both <- st[, uniqueN(get(group_col)), by = donor_id][V1 == 2, donor_id]
  paired <- length(donors_both) >= min_paired
  if (paired) st <- st[donor_id %in% donors_both]
  st[, grp := factor(get(group_col), levels = levels)]
  X <- if (paired) model.matrix(~ factor(donor_id) + grp, data = st) else model.matrix(~ grp, data = st)
  list(design = X, paired = paired, n_donors = uniqueN(st$donor_id),
       n_group = as.integer(table(st$grp)), table = st)
}

# Standard error of a log fold change recovered from a 1-df quasi-likelihood
# F statistic (F = (logFC/SE)^2). Undefined when F is 0.
se_from_f <- function(logFC, F) ifelse(is.finite(F) & F > 0, abs(logFC) / sqrt(pmax(F, 1e-300)), NA_real_)

# Random-effects meta-analysis of one gene across datasets. REML by default.
# If REML fails to converge, the DerSimonian-Laird estimate is used and the
# method actually applied is returned. k < 2 returns NA estimates.
meta_one <- function(yi, sei, method = "REML") {
  ok <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[ok]; sei <- sei[ok]; k <- length(yi)
  na <- list(k = k, estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, pval = NA_real_,
             tau2 = NA_real_, Q = NA_real_, Q_pval = NA_real_, I2 = NA_real_, method = NA_character_)
  if (k < 2) return(na)
  fit <- tryCatch(metafor::rma(yi = yi, sei = sei, method = method), error = function(e) NULL)
  used <- method
  if (is.null(fit)) {
    fit <- tryCatch(metafor::rma(yi = yi, sei = sei, method = "DL"), error = function(e) NULL)
    used <- "DL_fallback"
  }
  if (is.null(fit)) return(na)
  list(k = k, estimate = as.numeric(fit$beta), se = fit$se, ci_low = fit$ci.lb, ci_high = fit$ci.ub, pval = fit$pval,
       tau2 = fit$tau2, Q = fit$QE, Q_pval = fit$QEp, I2 = fit$I2 / 100, method = used)
}

# Leave-one-dataset-out: pooled estimate and p-value with each dataset removed.
meta_loo <- function(yi, sei, method = "REML") {
  k <- sum(is.finite(yi) & is.finite(sei) & sei > 0)
  if (k < 3) return(list(min_est = NA_real_, max_est = NA_real_, max_p = NA_real_))
  r <- lapply(seq_along(yi), function(i) meta_one(yi[-i], sei[-i], method))
  est <- vapply(r, `[[`, numeric(1), "estimate"); p <- vapply(r, `[[`, numeric(1), "pval")
  list(min_est = min(est, na.rm = TRUE), max_est = max(est, na.rm = TRUE), max_p = max(p, na.rm = TRUE))
}

# ---- bulk handoff ----------------------------------------------------------------
# R2 of each gene (samples x genes matrix) regressed on technical covariates
# (samples x k), within one cohort.
technical_r2 <- function(expr, covars) {
  X <- cbind(1, scale(covars))
  keep <- apply(X, 2, function(z) all(is.finite(z)))
  X <- X[, keep, drop = FALSE]
  B <- qr.solve(crossprod(X), crossprod(X, expr))
  res <- expr - X %*% B
  ss_tot <- colSums(sweep(expr, 2, colMeans(expr))^2)
  ifelse(ss_tot > 0, 1 - colSums(res^2) / ss_tot, NA_real_)
}

# Unweighted score: z-score each gene with supplied (training) centre and scale,
# then average. fit_zscore() returns the constants. Genes with zero variance in
# the training data are dropped and reported.
fit_zscore <- function(expr) {
  mu <- colMeans(expr); sdv <- apply(expr, 2, sd)
  ok <- is.finite(sdv) & sdv > 0
  list(center = mu[ok], scale = sdv[ok], genes = colnames(expr)[ok], dropped = colnames(expr)[!ok])
}
apply_zscore_mean <- function(expr, zfit) {
  Z <- sweep(sweep(expr[, zfit$genes, drop = FALSE], 2, zfit$center, "-"), 2, zfit$scale, "/")
  rowMeans(Z)
}

# ---- candidate lock ------------------------------------------------------------
# Content digest of a table that ignores row order, column order and float
# noise below `digits` significant digits, so the same scientific content
# always yields the same digest.
table_digest <- function(dt, digits = 10) {
  dt <- copy(as.data.table(dt))
  if (!ncol(dt)) return(digest::digest("empty", algo = "sha256", serialize = FALSE))
  setcolorder(dt, sort(names(dt)))
  for (c in names(dt)) {
    if (is.numeric(dt[[c]])) dt[[c]] <- ifelse(is.na(dt[[c]]), "NA", formatC(signif(dt[[c]], digits), digits = digits, format = "g"))
    else dt[[c]] <- ifelse(is.na(dt[[c]]), "NA", as.character(dt[[c]]))
  }
  rows <- sort(do.call(paste, c(as.list(dt), sep = "\u001f")))
  digest::digest(paste(c(paste(names(dt), collapse = "\u001f"), rows), collapse = "\n"), algo = "sha256", serialize = FALSE)
}

# Lock identifier: digest over everything that defines the selection (ordered
# candidate list, evidence, rules, input digests). Timestamps and environment
# details are excluded so identical inputs reproduce the same id.
lock_digest <- function(candidates, evidence, rules, input_digests) {
  oc <- as.data.table(candidates)[order(lock_rank)]
  cand_order <- paste(oc$lock_rank, oc$gene_key, oc$tier, sep = ":", collapse = ";")
  parts <- c(cand = digest::digest(cand_order, algo = "sha256", serialize = FALSE),
             candtab = table_digest(candidates),
             evidence = table_digest(evidence),
             rules = digest::digest(jsonlite::toJSON(rules, auto_unbox = TRUE, digits = NA), algo = "sha256", serialize = FALSE),
             inputs = digest::digest(paste(sort(paste(names(input_digests), unlist(input_digests))), collapse = ";"),
                                     algo = "sha256", serialize = FALSE))
  digest::digest(paste(names(parts), parts, collapse = "|"), algo = "sha256", serialize = FALSE)
}

# Per-cell rescue of malignant cells hidden under a non-malignant label. A cell
# in a flagged compartment becomes malignant when it detects at least
# `min_markers` of the ccRCC panel and none of the exclusion markers (immune,
# endothelial). Uses protein-coding markers only.
rescue_malignant <- function(counts, symbols, compartment, flagged, panel, min_markers, exclude) {
  cand <- compartment %in% flagged
  out <- rep(FALSE, length(compartment))
  if (!any(cand)) return(out)
  pi <- which(symbols %in% panel); ei <- which(symbols %in% exclude)
  n_panel <- if (length(pi)) Matrix::colSums(counts[pi, cand, drop = FALSE] > 0) else 0
  n_excl <- if (length(ei)) Matrix::colSums(counts[ei, cand, drop = FALSE] > 0) else 0
  out[cand] <- n_panel >= min_markers & n_excl == 0
  out
}

# Resolves a bulk-pipeline path from src/config/singlecell.yml against the
# repository root. If absent, it falls back to a standalone analysis directory
# beside the repository (src/R -> ../analysis/R, data/derived/cache ->
# ../analysis/cache).
bulk_path <- function(x, sub_root = SUB_ROOT, must_exist = TRUE) {
  p1 <- file.path(sub_root, x)
  if (file.exists(p1)) return(normalizePath(p1, winslash = "/"))
  alt <- sub("^src/R/", "../analysis/R/", x)
  alt <- sub("^data/derived/cache/", "../analysis/cache/", alt)
  p2 <- file.path(sub_root, alt)
  if (file.exists(p2)) return(normalizePath(p2, winslash = "/"))
  if (!must_exist) return(normalizePath(p1, winslash = "/", mustWork = FALSE))
  stop("bulk-pipeline file not found: ", p1,
       "; run the bulk pipeline (src/R/run_all.R) through stage 08 first")
}
