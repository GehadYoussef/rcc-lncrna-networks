# 39_pancancer_axis.R: is the leading lncRNA axis a library-quality axis across the 33 TCGA projects? (rules A and B)
# For each project from 38, in retained primary tumours and (with at least
# PAN_MIN_NORMALS) solid-tissue normals: PC1 of the lncRNA and protein-coding matrices
# and its Spearman correlation with the STAR metrics (rule A), per-gene |rho| with the
# non-feature fraction by positional lncRNA class (rule B), and variance explained by
# plate and tissue source site. Prevalence counts eligible tumour projects only.
# verify_decision_lock() stops if the rules differ from results/38_decision_rules_lock.json.
# Positional classes apply the stage 28 rules to every GENCODE v36 lncRNA and must
# reproduce stage 28 for all 3,442 discovery lncRNAs. All retained libraries enter PC1
# (discovery used the 511 network samples). The TCGA-KIRC row is the like-for-like check.
# Outputs: results/39_positional_classes_all_lncRNA.tsv,
#   39_positional_class_check.tsv, 39_pancancer_axis.tsv,
#   39_pancancer_class_gradient.tsv, 39_pancancer_prevalence.tsv.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
source(file.path(R_DIR, "pan_helpers.R"))
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(matrixStats)
})
banner("39 | Pan-cancer axis prevalence (rules A and B)")
set.seed(SEED)
verify_decision_lock()

# ---- 1. positional class of every GENCODE v36 lncRNA ----
classify_lncRNA <- function(gtf) {
  genes <- gtf$genes; tx <- gtf$tx; exons <- gtf$exons
  g_l <- genes[gene_type == "lncRNA"]; g_p <- genes[gene_type == "protein_coding"]
  tx_l <- tx[gene_id %in% g_l$gene_id]; tx_p <- tx[gene_id %in% g_p$gene_id]
  ex_l <- exons[gene_id %in% g_l$gene_id]; ex_p <- exons[gene_id %in% g_p$gene_id]
  mkgr <- function(d) GRanges(d$chr, IRanges(d$start, d$end), strand = d$strand)
  lg <- mkgr(g_l); lg$gene_id <- g_l$gene_id
  pg <- mkgr(g_p); pg$gene_id <- g_p$gene_id
  le <- mkgr(ex_l); le$gene_id <- ex_l$gene_id
  pe <- mkgr(ex_p); pe$gene_id <- ex_p$gene_id; pe$transcript_id <- ex_p$transcript_id
  sc <- function(x) as.character(strand(x))
  exl <- split(pe, pe$transcript_id); rng <- range(exl)
  k <- lengths(rng) == 1L; exl <- exl[k]; rng <- rng[k]
  intr_l <- psetdiff(unlist(rng, use.names = TRUE), exl)
  intr <- unlist(intr_l, use.names = FALSE)
  ov <- findOverlaps(le, pe, ignore.strand = TRUE)
  exf <- data.table(gene_id = le$gene_id[queryHits(ov)],
                    same = sc(le)[queryHits(ov)] == sc(pe)[subjectHits(ov)])[
    , .(exon_overlap_sense = any(same), exon_overlap_antisense = any(!same)), by = gene_id]
  ov <- findOverlaps(lg, intr, type = "within", ignore.strand = TRUE)
  inf <- data.table(gene_id = lg$gene_id[queryHits(ov)],
                    same = sc(lg)[queryHits(ov)] == sc(intr)[subjectHits(ov)])[
    , .(within_intron = TRUE, intron_host_same_strand = any(same)), by = gene_id]
  ov <- findOverlaps(lg, pg, ignore.strand = TRUE)
  spf <- data.table(gene_id = lg$gene_id[queryHits(ov)],
                    same = sc(lg)[queryHits(ov)] == sc(pg)[subjectHits(ov)])[
    , .(span_overlap_sense = any(same), span_overlap_antisense = any(!same)), by = gene_id]
  tss_of <- function(d) ifelse(d$strand == "+", d$start, d$end)
  lt <- tss_of(tx_l); pt <- tss_of(tx_p)
  lts <- GRanges(tx_l$chr, IRanges(lt, lt), strand = tx_l$strand); lts$gene_id <- tx_l$gene_id
  pts <- GRanges(tx_p$chr, IRanges(pt, pt), strand = tx_p$strand)
  win <- GRanges(tx_l$chr, IRanges(pmax(1L, lt - 1000L), lt + 1000L))
  ov <- findOverlaps(win, pts, ignore.strand = TRUE)
  q <- queryHits(ov); s <- subjectHits(ov)
  sl <- sc(lts)[q]; sp <- sc(pts)[s]; tl <- start(lts)[q]; tp <- start(pts)[s]
  h2h <- (sl == "+" & sp == "-" & tp <= tl) | (sl == "-" & sp == "+" & tp >= tl)
  dvf <- unique(data.table(gene_id = lts$gene_id[q][h2h], div_tss = TRUE))
  cls <- data.table(gene_id = g_l$gene_id, gene_name = g_l$gene_name)
  for (f in list(exf, inf, spf, dvf)) cls <- merge(cls, f, by = "gene_id", all.x = TRUE)
  for (v in c("exon_overlap_sense", "exon_overlap_antisense", "within_intron",
              "intron_host_same_strand", "span_overlap_sense", "span_overlap_antisense", "div_tss"))
    set(cls, which(is.na(cls[[v]])), v, FALSE)
  cls[, divergent := div_tss & !span_overlap_sense & !span_overlap_antisense]
  cls[, class_detail := fifelse(exon_overlap_sense, "exonic_sense",
                        fifelse(exon_overlap_antisense, "exonic_antisense",
                        fifelse(within_intron & intron_host_same_strand, "intronic_sense",
                        fifelse(within_intron, "intronic_antisense",
                        fifelse(span_overlap_antisense, "antisense_overlapping",
                        fifelse(span_overlap_sense, "sense_overlapping",
                        fifelse(divergent, "divergent", "intergenic")))))))]
  cls[, class := fifelse(class_detail == "exonic_sense", "exonic_sense",
                 fifelse(class_detail %in% c("exonic_antisense", "antisense_overlapping"), "antisense",
                 fifelse(class_detail %in% c("intronic_sense", "intronic_antisense"), "intronic",
                         class_detail)))]
  cls[, .(gene_id, gene_name, class)]
}
cls_rds <- file.path(CACHE_DIR, "39_lncRNA_positional_all.rds")
if (file.exists(cls_rds)) cls <- readRDS(cls_rds) else {
  msg("classifying every GENCODE v36 lncRNA ...")
  cls <- classify_lncRNA(readRDS(file.path(CACHE_DIR, "28_gencode_v36_parsed.rds")))
  saveRDS(cls, cls_rds)
}
cls[, key := sub("\\..*$", "", gene_id)]
ref <- fread(file.path(RESULTS_DIR, "28_lncRNA_positional_classes.tsv"))[, .(key = sub("\\..*$", "", gene_id), class28 = class)]
chk <- merge(ref, cls[, .(key, class)], by = "key", all.x = TRUE)
check <- data.table(n_discovery_lncRNA = nrow(chk), n_found = sum(!is.na(chk$class)),
                    n_same_class = sum(chk$class == chk$class28, na.rm = TRUE))
save_tsv(check, "39_positional_class_check.tsv"); print(check)
stopifnot(check$n_same_class == check$n_discovery_lncRNA)
save_tsv(cls[, .(gene_id, gene_name, class)], "39_positional_classes_all_lncRNA.tsv")
CLASSES <- c("exonic_sense", "antisense", "intronic", "sense_overlapping", "divergent", "intergenic")

# ---- 2. per project x group ----
sp <- function(x, y) { ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
                       c(unname(ct$estimate), ct$p.value) }
rows <- list(); grad <- list()
for (proj in pan_projects()) {
  obj <- readRDS(file.path(CACHE_DIR, paste0("pan_", proj, ".rds")))
  for (grp in c("tumour", "normal")) {
    n_ret <- obj$samples[keep == TRUE & group == grp, .N]
    if (grp == "normal" && n_ret < PAN_MIN_NORMALS) next
    if (n_ret < 10) next
    G <- pan_group(obj, grp); s <- G$samples
    nf <- s$pct_noFeature; mm <- s$pct_multimapping; dp <- log10(s$assigned_reads)
    L <- pc1(G$lnc); P <- pc1(G$pc)
    rl <- sp(L$score, nf); rp <- sp(P$score, nf)
    gl <- abs(cor(G$lnc, nf, method = "spearman"))[, 1]
    gp <- abs(cor(G$pc,  nf, method = "spearman"))[, 1]
    eligible <- if (grp == "tumour") n_ret >= PAN_MIN_TUMOURS && !proj %in% PAN_EXCLUDE_INFERENCE else NA
    rows[[length(rows) + 1]] <- data.table(
      project = proj, group = grp, n = nrow(s), eligible = eligible,
      n_lncRNA = ncol(G$lnc), n_pc = ncol(G$pc),
      median_noFeature = median(nf),
      lnc_pc1_var_share = L$var_share, pc_pc1_var_share = P$var_share,
      rho_lnc_pc1_noFeature = rl[1], p_lnc_pc1_noFeature = rl[2],
      rho_pc_pc1_noFeature = rp[1], p_pc_pc1_noFeature = rp[2],
      rho_lnc_pc1_multimap = sp(L$score, mm)[1], rho_lnc_pc1_depth = sp(L$score, dp)[1],
      median_abs_rho_lnc_genes = median(gl), median_abs_rho_pc_genes = median(gp),
      frac_lnc_genes_gt0.3 = mean(gl > 0.3), frac_pc_genes_gt0.3 = mean(gp > 0.3),
      r2_plate_metric = r2_factor(nf, pool_levels(s$plate)),
      r2_tss_metric = r2_factor(nf, pool_levels(s$tss)),
      r2_plate_lnc_pc1 = r2_factor(L$score, pool_levels(s$plate)),
      n_plates = uniqueN(s$plate))
    # gkey, not key: `key` is an argument of data.table() itself
    gk <- data.table(gkey = sub("\\..*$", "", colnames(G$lnc)), abs_rho = gl)
    gk <- merge(gk, cls[, .(gkey = key, class)], by = "gkey")
    kw <- kruskal.test(abs_rho ~ factor(class), data = gk)
    grad[[length(grad) + 1]] <- gk[, .(n_genes = .N, median_abs_rho = median(abs_rho)), by = class][
      , `:=`(project = proj, group = grp, kw_p = kw$p.value)]
    msg(sprintf("%-10s %-6s n=%4d  lnc rho %.2f  pc rho %.2f", proj, grp, nrow(s), rl[1], rp[1]))
  }
  rm(obj); gc(verbose = FALSE)
}
ax <- rbindlist(rows)
ax[, rule_A_pass := abs(rho_lnc_pc1_noFeature) >= PAN_AXIS_RHO_MIN &
                    abs(rho_lnc_pc1_noFeature) - abs(rho_pc_pc1_noFeature) >= PAN_AXIS_MARGIN]
gr <- rbindlist(grad)
gw <- dcast(gr, project + group ~ class, value.var = "median_abs_rho")
ax <- merge(ax, gw[, .(project, group, intronic_minus_intergenic = intronic - intergenic,
                       exonic_sense_minus_intergenic = exonic_sense - intergenic)],
            by = c("project", "group"), all.x = TRUE)
setorder(ax, group, -rho_lnc_pc1_noFeature)
save_tsv(ax, "39_pancancer_axis.tsv")
gr[, class := factor(class, levels = CLASSES)]
save_tsv(gr[order(project, group, class), .(project, group, class, n_genes, median_abs_rho, kw_p)],
         "39_pancancer_class_gradient.tsv")

el <- ax[group == "tumour" & eligible == TRUE]
prev <- data.table(
  quantity = c("eligible tumour projects", "rule A met",
               "lncRNA PC1 |rho| >= 0.5 (first half of rule A)",
               "intronic median |rho| above intergenic",
               "normal sets analysed", "normal sets meeting rule A thresholds (descriptive)"),
  n = c(nrow(el), sum(el$rule_A_pass), sum(abs(el$rho_lnc_pc1_noFeature) >= PAN_AXIS_RHO_MIN),
        sum(el$intronic_minus_intergenic > 0, na.rm = TRUE),
        ax[group == "normal", .N], ax[group == "normal", sum(rule_A_pass)]),
  of = c(NA, nrow(el), nrow(el), nrow(el), NA, ax[group == "normal", .N]))
save_tsv(prev, "39_pancancer_prevalence.tsv"); print(prev)
write_session_info("39_pancancer_axis")
msg("39 done")
