# 44_ruv3_prps_comparison.R: per-gene residualisation compared with RUV-III with PRPS.
# Applies RUV-III with pseudo-replicates of pseudo-samples (PRPS, Molania et al. 2023,
# via ruv::RUVIII) to TCGA-KIRC tumours, log2(FPKM + 1).
# PRPS are formed within ccA/ccB x stage groups by plate, by extreme non-feature
# fraction and by extreme depth. Negative controls are biology-poor genes most
# associated with the metric and depth. k = 3 (1 and 5 as sensitivity).
# Observed, residualised and RUV-III data are compared on the leading lncRNA
# component, the lncRNA network, the single-cell tumour-specificity test and module HRs.
# Inputs: caches pan_TCGA-KIRC.rds, 30_TCGA-CDR.xlsx, module_loadings.rds, results of 02, 34 and 42.
# Outputs: results 44_prps_{design,axis,network,sc_truth,module_hr}.tsv.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(WGCNA); library(survival); library(readxl); library(ruv); library(matrixStats) })
enableWGCNAThreads(N_THREADS)
banner("44 | Residualisation against RUV-III with PRPS")
set.seed(SEED)
N_NCG <- 1000; K_MAIN <- 3; N_NULL <- 500
MARKERS <- c("CA9", "NDUFA4L2", "ANGPTL4", "VEGFA", "EGLN3")
SC_DIR <- Filter(dir.exists, c(file.path(PROJECT_ROOT, "submission", "results", "singlecell"),
                               file.path(PROJECT_ROOT, "results", "singlecell")))[1]

# ---- 1. data ------------------------------------------------------------------
obj <- readRDS(file.path(CACHE_DIR, "pan_TCGA-KIRC.rds"))
s <- obj$samples[keep == TRUE & group == "tumour"]
cc <- fread(file.path(RESULTS_DIR, "34_ccAB_classification.tsv"))[cohort == "TCGA-KIRC", .(patient, ccAB_class)]
cdr <- as.data.table(read_excel(file.path(CACHE_DIR, "30_TCGA-CDR.xlsx"), sheet = "TCGA-CDR", guess_max = 20000))[
  , .(patient = bcr_patient_barcode, age = as.numeric(age_at_initial_pathologic_diagnosis), sex = tolower(gender),
      stage_raw = ajcc_pathologic_tumor_stage, os = suppressWarnings(as.numeric(OS)), os_time = suppressWarnings(as.numeric(OS.time)))]
roman <- function(x) { r <- sub("[A-C]$", "", sub("^STAGE\\s*", "", toupper(trimws(x))))
  c(I = 1L, II = 2L, III = 3L, IV = 4L)[r] }
cdr[, stage := unname(roman(stage_raw))]
cdr[os_time > OS_CENSOR_DAYS, `:=`(os = 0, os_time = OS_CENSOR_DAYS)]
s <- merge(s, cc, by = "patient"); s <- merge(s, cdr, by = "patient", all.x = TRUE)
s[, bio := paste(ccAB_class, fifelse(is.na(stage), "stageNA", fifelse(stage <= 2, "I-II", "III-IV")), sep = "_")]
msg(nrow(s), " tumours with a ccA/ccB class")

pc_disc  <- fread(file.path(RESULTS_DIR, "02_mRNA_module_genes.tsv"))$gene_id
lnc_disc <- fread(file.path(RESULTS_DIR, "02_lncRNA_module_genes.tsv"))$gene_id
lnc_all  <- obj$ann[gene_type == "lncRNA", gene_id]
fT <- obj$fpkm[, s$file_id]
lnc_expr <- lnc_all[rowMeans(fT[lnc_all, ] >= LNC_MIN_FPKM) >= LNC_MIN_FRAC]
genes <- unique(c(intersect(pc_disc, rownames(fT)), union(lnc_expr, intersect(lnc_disc, rownames(fT)))))
mk_ids <- obj$ann[gene_name %in% MARKERS & gene_type == "protein_coding", gene_id]
genes <- unique(c(genes, mk_ids))
Y <- t(log2(fT[genes, ] + 1))
cv <- tech_covariates(s); nf <- cv[, "pct_noFeature"]; dp <- cv[, "log_depth"]
msg(ncol(Y), " genes")

# ---- 2. PRPS and negative controls -----------------------------------------
prps <- list(); design <- list()
for (b in unique(s$bio[!grepl("stageNA", s$bio)])) {
  idx <- which(s$bio == b)
  pl <- table(s$plate[idx]); pl <- names(pl)[pl >= 3]
  if (length(pl) >= 2) for (p in pl) {
    ii <- idx[s$plate[idx] == p]
    prps[[length(prps) + 1]] <- list(set = paste0("plate|", b), rows = ii)
  }
  for (v in c("metric", "depth")) {
    x <- if (v == "metric") nf[idx] else dp[idx]
    if (length(idx) >= 30) {
      o <- idx[order(x)]
      prps[[length(prps) + 1]] <- list(set = paste0(v, "|", b), rows = head(o, 10))
      prps[[length(prps) + 1]] <- list(set = paste0(v, "|", b), rows = tail(o, 10))
    }
  }
}
design <- rbindlist(lapply(prps, function(p) data.table(set = p$set, n_samples = length(p$rows))))
sets <- unique(design$set)
design <- design[, .(pseudo_samples = .N, samples_averaged = sum(n_samples)), by = set]
save_tsv(design, "44_prps_design.tsv"); msg(nrow(design), " PRPS sets, ", sum(design$pseudo_samples), " pseudo-samples")
P <- do.call(rbind, lapply(prps, function(p) colMeans(Y[p$rows, , drop = FALSE])))
M <- matrix(0L, nrow(Y) + nrow(P), nrow(Y) + length(sets))
for (i in seq_len(nrow(Y))) M[i, i] <- 1L
for (j in seq_along(prps)) M[nrow(Y) + j, nrow(Y) + match(prps[[j]]$set, sets)] <- 1L

bioF <- apply(Y, 2, function(g) { f <- summary(aov(g ~ factor(s$bio)))[[1]][1, "F value"]; if (is.finite(f)) f else 0 })
cand <- names(bioF)[bioF <= median(bioF)]
un <- rank(-abs(cor(Y[, cand], nf, method = "spearman"))[, 1]) + rank(-abs(cor(Y[, cand], dp, method = "spearman"))[, 1])
ctl <- colnames(Y) %in% cand[order(un)][seq_len(N_NCG)]   # below-median biology F, then top summed |Spearman| ranks

ruv3 <- function(k) {
  Z <- RUVIII(Y = rbind(Y, P), M = M, ctl = ctl, k = k, inputcheck = FALSE)
  Z <- Z[seq_len(nrow(Y)), , drop = FALSE]; dimnames(Z) <- dimnames(Y); Z
}
mats <- list(observed = Y, residualised = remove_technical(Y, cv), ruv3_prps = ruv3(K_MAIN))

# ---- 3a. the leading component -----------------------------------------------
lnc_cols <- intersect(lnc_disc, colnames(Y)); pc_cols <- intersect(pc_disc, colnames(Y))
pc1 <- function(X) { Xc <- scale(X, TRUE, FALSE); sc <- svd(Xc, nu = 1, nv = 0)$u[, 1]; if (cor(sc, rowMeans(X)) < 0) -sc else sc }
plate_f <- { tb <- table(s$plate); factor(ifelse(s$plate %in% names(tb)[tb >= 10], s$plate, "other")) }
axis_row <- function(nm, X) data.table(version = nm,
  rho_lnc_pc1_metric = cor(pc1(X[, lnc_cols]), nf, method = "spearman"),
  rho_pc_pc1_metric = cor(pc1(X[, pc_cols]), nf, method = "spearman"),
  r2_plate_lnc_pc1 = summary(lm(pc1(X[, lnc_cols]) ~ plate_f))$r.squared,
  median_abs_gene_rho_lnc = median(abs(cor(X[, lnc_cols], nf, method = "spearman"))))
ax <- rbindlist(Map(axis_row, names(mats), mats))
ksens <- rbindlist(lapply(c(1, 5), function(k) axis_row(paste0("ruv3_prps_k", k), ruv3(k))))
ax <- rbind(ax, ksens)
save_tsv(ax, "44_prps_axis.tsv"); print(ax)

# ---- 3b. networks on the discovery lncRNAs ---------------------------------------
E0 <- Y[, lnc_cols]
gsg <- goodSamplesGenes(E0, verbose = 0); keep_s <- gsg$goodSamples; keep_g <- gsg$goodGenes
A <- adjacency(t(E0[keep_s, keep_g]), type = "distance"); kk <- colSums(A) - 1; rm(A)
keep_s[which(keep_s)[scale(kk)[, 1] < -2.5]] <- FALSE
build <- function(E) { set.seed(SEED)
  bw <- blockwiseModules(E, power = LNC_POWER, networkType = NETWORK_TYPE, TOMType = TOM_TYPE,
                         minModuleSize = MIN_MODULE_SIZE, mergeCutHeight = LNC_MERGE, deepSplit = LNC_DEEPSPLIT,
                         minKMEtoStay = LNC_MINKME, numericLabels = FALSE, pamRespectsDendro = FALSE,
                         maxBlockSize = MAX_BLOCK_SIZE, saveTOMs = FALSE, verbose = 0)
  setNames(bw$colors, colnames(E)) }
cols <- lapply(mats, function(X) build(X[keep_s, lnc_cols[keep_g]]))
ari <- function(x, y) { tab <- table(x, y); n <- sum(tab); sij <- sum(choose(tab, 2)); si <- sum(choose(rowSums(tab), 2)); sj <- sum(choose(colSums(tab), 2))
  e <- si * sj / choose(n, 2); (sij - e) / ((si + sj) / 2 - e) }
nw <- rbindlist(lapply(names(cols), function(nm) { cl <- cols[[nm]]; tb <- table(cl[cl != "grey"])
  data.table(version = nm, n_samples = sum(keep_s), n_lncRNA = length(cl), n_modules = length(tb),
             largest_frac = max(tb) / length(cl), grey_frac = mean(cl == "grey"),
             ari_vs_observed = ari(cl, cols$observed), ari_vs_residualised = ari(cl, cols$residualised),
             ari_vs_ruv3 = ari(cl, cols$ruv3_prps)) }))
save_tsv(nw, "44_prps_network.tsv"); print(nw)

# ---- 3c. single-cell tumour-specificity test (from 42) ------------------------------
lock <- fread(file.path(SC_DIR, "locked_lncRNA_candidates.tsv"))[, gk := sub("[.].*$", "", gene_key)]
meta <- fread(file.path(SC_DIR, "08_meta_lncRNA.tsv.gz"))[contrast == "malignant_vs_normal_epithelial"]
key <- sub("[.].*$", "", colnames(Y))
t1_ids <- colnames(Y)[key %in% lock[tier == "tier1_tumour_specific", gk]]
pool <- colnames(Y)[key %in% meta[abs(pooled_log2FC) < 0.25 & pooled_FDR > 0.5, gene_key] & colnames(Y) %in% lnc_expr]
ab <- colMeans(Y); bins <- setNames(cut(ab, quantile(ab, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE), colnames(Y))
need <- table(bins[t1_ids])
null_sets <- lapply(seq_len(N_NULL), function(k) unlist(lapply(names(need), function(b) {
  cand <- pool[bins[pool] == as.integer(b)]; if (length(cand)) sample(cand, min(length(cand), need[[b]])) })))
marker <- rowMeans(scale(Y[, mk_ids]))
score <- function(X, ids) rowMeans(scale(X[, ids, drop = FALSE]))
sct <- rbindlist(lapply(names(mats), function(nm) { X <- mats[[nm]]
  r1 <- cor(score(X, t1_ids), marker, method = "spearman")
  rn <- vapply(null_sets, function(ids) cor(score(X, ids), marker, method = "spearman"), 1)
  data.table(version = nm, n_tier1 = length(t1_ids), rho_marker_tier1 = r1, rho_marker_null_median = median(rn),
             excess_over_null = r1 - median(rn), rho_metric_tier1 = cor(score(X, t1_ids), nf, method = "spearman")) }))
save_tsv(sct, "44_prps_sc_truth.tsv"); print(sct)

# ---- 3d. prognostic modules with fixed discovery loadings --------------------------
L <- readRDS(file.path(CACHE_DIR, "module_loadings.rds"))
mods <- c("mRNA_MEgreen", "mRNA_MEpurple", "lnc_MEblue", "lnc_MEgreenyellow", "lnc_MEturquoise")
ok <- is.finite(s$os_time) & s$os_time > 0 & !is.na(s$os) & is.finite(s$age) & !is.na(s$stage)
hr <- rbindlist(lapply(names(mats), function(nm) {
  sc_ <- score_modules(mats[[nm]], L[mods])
  rbindlist(lapply(colnames(sc_), function(m) {
    d <- data.table(t = s$os_time, e = s$os, x = sc_[, m], age = scale(s$age)[, 1], male = as.integer(s$sex == "male"), stage = s$stage)[ok]
    f <- summary(coxph(Surv(t, e) ~ x + age + male + stage, data = d))$coefficients["x", ]
    data.table(version = nm, module = m, n = nrow(d), events = sum(d$e), hr = exp(f[["coef"]]),
               lo = exp(f[["coef"]] - 1.96 * f[["se(coef)"]]), hi = exp(f[["coef"]] + 1.96 * f[["se(coef)"]]), p = f[["Pr(>|z|)"]])
  }))
}))
save_tsv(hr, "44_prps_module_hr.tsv"); print(hr)
write_session_info("44_ruv3_prps_comparison")
msg("44 done")
