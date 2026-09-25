# 34_proximal_tubule_reading.R: is the catabolic protein-coding module a proximal-tubule signal?
# The hubs of the catabolic module (catabolic_reference_module()) are proximal-tubule
# (PT) identity genes. Tests whether the module indexes PT differentiation or normal-kidney
# admixture, and whether the coupled lncRNA modules (blue, greenyellow, turquoise) add
# signal beyond it. Sections: 1 hub annotation (PT markers, HPA, ClearCode34), 2 PT marker
# scores and ccA/ccB axis, 3 analysis frames, 4 correlations, 5 Cox models, 6 ESTIMATE
# overlap, 7 projection of adjacent-normal libraries (cached by 25) onto the module.
# Inputs: cache dataset.rds, networks.rds, expr_raw.rds, estimate_scores.rds,
#   biospecimen_kirc.rds, validation_dataset.rds, locked_model.rds,
#   valid_estimate.rds, 25_normal_libraries.rds.
# Outputs: 34_*.tsv (marker set, hub PT flags, ClearCode34 membership and classes,
#   PT score correlations, adjusted models, ESTIMATE overlap, normal admixture).

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(survival) })
banner("34 | Proximal-tubule reading of the catabolic protein-coding module")
set.seed(SEED)
t_start <- Sys.time()

ds   <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets <- readRDS(file.path(CACHE_DIR, "networks.rds"))
est  <- as.data.table(readRDS(file.path(CACHE_DIR, "estimate_scores.rds")))
bio  <- readRDS(file.path(CACHE_DIR, "biospecimen_kirc.rds"))
vd   <- readRDS(file.path(CACHE_DIR, "validation_dataset.rds"))
L    <- readRDS(file.path(CACHE_DIR, "locked_model.rds"))
stopifnot(identical(L$version, "v9"))
cohort <- as.data.table(ds$cohort_full)
if (!"assigned_reads" %in% names(cohort)) cohort[, assigned_reads := libsize]
msg("TCGA-KIRC: ", nrow(cohort), " patients, ", sum(cohort$os_event), " deaths; CPTAC-3: ",
    nrow(vd$cohort), " patients, ", sum(vd$cohort$os_event), " deaths")

REF_MOD  <- catabolic_reference_module()
REF_COL  <- paste0("mRNA_ME", REF_MOD)
LNC_MODS <- c("blue", "greenyellow", "turquoise")
gt_m     <- as.data.table(nets$mrna$gene_tbl)

# ---- 1. Marker set and the hubs of the module ----
banner("1 | Proximal-tubule marker set and the top hubs of the module")
HPA_SRC  <- "HPA single cell type atlas (Karlsson 2021 Sci Adv 7:eabh2169), gene page accessed 2026-09-10"
KPMP_SRC <- "KPMP atlas (Lake 2023 Nature 619:585)"
# core = canonical PT markers, extended = further markers whose HPA
# single-cell specificity names proximal tubule cells.
PT_MARKERS <- rbindlist(list(
  list("LRP2",    TRUE,  "Cell type enriched (Proximal tubule cells)", "pan-PT (megalin)"),
  list("CUBN",    TRUE,  "Group enriched (Proximal tubule cells, Enterocytes)", "pan-PT (cubilin)"),
  list("SLC34A1", TRUE,  "Cell type enriched (Proximal tubule cells)", "PT S1-S2 (NaPi-IIa)"),
  list("SLC22A6", TRUE,  "Cell type enriched (Proximal tubule cells)", "PT S2 (OAT1)"),
  list("SLC22A8", TRUE,  "Group enriched (Retinal pigment epithelial cells, Proximal tubule cells, Choroid plexus epithelial cells)", "PT S2-S3 (OAT3)"),
  list("SLC5A2",  TRUE,  "Cell type enhanced (Late spermatids, Proximal tubule cells, Late primary spermatocytes, Early spermatids, Cardiomyocytes, Epicardial cells)", "PT S1 (SGLT2)"),
  list("SLC13A3", TRUE,  "Cell type enriched (Proximal tubule cells)", "PT S1-S2 (NaDC3)"),
  list("ALDOB",   TRUE,  "Cell type enriched (Enterocytes)", "PT gluconeogenic enzyme; PT not named in HPA single-cell specificity"),
  list("PDZK1",   TRUE,  "Cell type enhanced (Proximal tubule cells, Enterocytes, Epididymal efferent duct absorptive cells, Cholangiocytes, Hepatocytes, Loop of Henle epithelial cells)", "PT brush border scaffold"),
  list("GLYAT",   TRUE,  "Cell type enhanced (Hepatocytes, Proximal tubule cells, Adipocytes, Epididymal efferent duct absorptive cells)", "PT (glycine N-acyltransferase)"),
  list("HNF4A",   TRUE,  "Cell type enhanced (Colonocytes, Enteric stem cells, Paneth cells, Enteric transient amplifying cells, Enterocytes, Hepatocytes, Goblet cells, Proximal tubule cells, Foveolar cells, Tuft cells, Gastric progenitor cells, Neuroendocrine cells)", "PT lineage transcription factor"),
  list("SLC3A1",  TRUE,  "Cell type enhanced (Proximal tubule cells, Enterocytes, Distal convoluted tubule cells, Pancreatic duct cells, Oocytes, Endometrial luminal cells, Loop of Henle epithelial cells, Renal connecting tubule cells, Epididymal efferent duct absorptive cells, Renal collecting duct intercalated cells, Endometrial glandular cells)", "PT (rBAT)"),
  list("SLC17A1", TRUE,  "Cell type enriched (Proximal tubule cells)", "PT (NPT1)"),
  list("FBP1",    TRUE,  "Cell type enhanced (Enterocytes, Hepatocytes, Urothelial cells, Prostatic club cells, Breast lactating cells, Prostatic hillock cells)", "PT gluconeogenic enzyme; PT not named in HPA single-cell specificity"),
  list("PCK1",    TRUE,  "Group enriched (Enterocytes, Hepatocytes, Colonocytes)", "PT gluconeogenic enzyme; PT not named in HPA single-cell specificity"),
  list("SLC7A13", FALSE, "Cell type enriched (Proximal tubule cells)", "PT S3 (AGT1)"),
  list("SLC34A3", FALSE, "Cell type enhanced (Enterocytes, Proximal tubule cells, Esophageal apical cells, Esophageal suprabasal cells)", "PT (NaPi-IIc)"),
  list("SLC36A2", FALSE, "Group enriched (Myonuclei, Proximal tubule cells)", "PT (PAT2)"),
  list("SLC5A12", FALSE, "Cell type enhanced (Proximal tubule cells, Enterocytes, Epididymal principal cells)", "PT S1 (SMCT2)"),
  list("SLC22A2", FALSE, "Cell type enhanced (Proximal tubule cells, Loop of Henle epithelial cells, Renal connecting tubule cells, Early spermatids, Distal convoluted tubule cells, Tuft cells, Podocytes)", "PT S3 (OCT2)"),
  list("MIOX",    FALSE, "Cell type enriched (Proximal tubule cells)", "PT (myo-inositol oxygenase)")))
setnames(PT_MARKERS, c("gene_name", "core", "hpa_single_cell_specificity", "segment_note"))
PT_MARKERS[, hpa_pt_named := grepl("Proximal tubule", hpa_single_cell_specificity)]
PT_MARKERS[, source := paste0(HPA_SRC, "; ", KPMP_SRC)]

# ---- how exclusive is each HPA specificity string? ----
# An HPA string lists every cell type the gene is enhanced in, so "names
# proximal tubule" is permissive. These helpers record the list length, the
# position of proximal tubule and an exclusivity grade. The strings are
# hand-entered for these genes only, so no background rate is computed.
hpa_cell_types <- function(s) {
  if (is.na(s)) return(character(0))
  inner <- sub("\\)\\s*$", "", sub("^[^(]*\\(", "", s))
  trimws(strsplit(inner, ",", fixed = TRUE)[[1]])
}
hpa_n_types <- function(s) vapply(s, function(x) length(hpa_cell_types(x)), integer(1), USE.NAMES = FALSE)
hpa_pt_pos  <- function(s) vapply(s, function(x) {
  ct <- hpa_cell_types(x); i <- grep("Proximal tubule", ct)
  if (length(i)) as.integer(i[1]) else NA_integer_ }, integer(1), USE.NAMES = FALSE)
hpa_grade <- function(s) {
  n <- hpa_n_types(s); p <- hpa_pt_pos(s)
  ifelse(is.na(s), NA_character_,
  ifelse(is.na(p), "PT not named",
  ifelse(n == 1L, "PT-exclusive",
  ifelse(n == 2L, "PT with 1 other cell type",
         paste0("PT within a list of ", n, " cell types")))))
}
PT_MARKERS[, hpa_n_cell_types := hpa_n_types(hpa_single_cell_specificity)]
PT_MARKERS[, hpa_pt_rank := hpa_pt_pos(hpa_single_cell_specificity)]
PT_MARKERS[, hpa_pt_exclusivity := hpa_grade(hpa_single_cell_specificity)]

# HPA single-cell specificity of the module's top 20 hubs (same access date).
HUB_HPA <- rbindlist(list(
  list("GLYATL1",  "Group enriched (Hepatocytes, Proximal tubule cells)"),
  list("PDZK1",    PT_MARKERS[gene_name == "PDZK1", hpa_single_cell_specificity]),
  list("DDC",      "Cell type enhanced (Rod photoreceptor cells, Adrenal medulla cells, Neuroendocrine cells, Enterocytes, Proximal tubule cells, Retinal horizontal cells, Retinal ganglion cells, Retinal bipolar cells)"),
  list("ACAA2",    "Cell type enhanced (Hepatocytes, Colonocytes, Late spermatids, Cholangiocytes, Enterocytes, Enteric transient amplifying cells)"),
  list("C11orf54", "Cell type enhanced (Hepatocytes, Choroid plexus epithelial cells, Oligodendrocyte progenitor cells, Astrocytes, Bergmann glia, Oligodendrocytes, Ependymal cells, Brain inhibitory neurons)"),
  list("FMO1",     "Cell type enhanced (Salivary ionocytes, Enterocytes, Proximal tubule cells, Ocular epithelial cells, Mesothelial cells, Pituitary stem cells, Granulosa cells)"),
  list("AGMAT",    "Cell type enhanced (Hepatocytes, Gastric progenitor cells, Proximal tubule cells, Enteric transient amplifying cells, Enterocytes, Enteric stem cells, Epididymal efferent duct absorptive cells, Paneth cells)"),
  list("EHHADH",   "Cell type enhanced (Hepatocytes, Proximal tubule cells, Retinal pigment epithelial cells, Adipocytes)"),
  list("AGXT2",    "Group enriched (Proximal tubule cells, Hepatocytes)"),
  list("SLC47A1",  "Cell type enhanced (Proximal tubule cells, Early spermatids, Late primary spermatocytes, Epididymal efferent duct absorptive cells, Myosatellite cells, Adrenal cortex cells, Myonuclei)"),
  list("CLCN5",    "Cell type enhanced (Myosatellite cells, Distal convoluted tubule cells, Loop of Henle epithelial cells, Renal connecting tubule cells, Proximal tubule cells)"),
  list("NAT8",     "Cell type enhanced (Epididymal efferent duct absorptive cells, Hepatocytes, Enterocytes, Proximal tubule cells, Epididymal efferent duct ciliated cells, Adrenal medulla cells)"),
  list("BBOX1",    "Cell type enhanced (Proximal tubule cells, Ependymal cells, Loop of Henle epithelial cells, Epicardial cells, Astrocytes, Esophageal suprabasal cells, Breast secretory cells, Fallopian tube ciliated cells, Suprabasal keratinocytes)"),
  list("SLC22A11", "Group enriched (Syncytiotrophoblasts, Cytotrophoblasts)"),
  list("FUT6",     "Cell type enhanced (Esophageal apical cells, Enterocytes, Proximal tubule cells, Goblet cells, Colonocytes, Paneth cells, Respiratory secretory cells, Esophageal suprabasal cells, Enteric transient amplifying cells)"),
  list("LRP2",     PT_MARKERS[gene_name == "LRP2", hpa_single_cell_specificity]),
  list("GIPC2",    "Cell type enhanced (Adrenal cortex cells, Enterocytes, Proximal tubule cells, Epididymal efferent duct absorptive cells, Pancreatic acinar cells)"),
  list("CRYL1",    "Cell type enhanced (Enterocytes, Proximal tubule cells, Hepatocytes, Renal collecting duct intercalated cells)"),
  list("SLC27A2",  "Cell type enhanced (Cytotrophoblasts, Hepatocytes, Proximal tubule cells, Respiratory ciliated cells, Conjunctival goblet cells, Respiratory secretory cells, Migrating cytotrophoblasts, Lactotrophs, Epididymal principal cells, Retinal pigment epithelial cells, Enterocytes, Epididymal efferent duct absorptive cells)"),
  list("GLYAT",    PT_MARKERS[gene_name == "GLYAT", hpa_single_cell_specificity])))
setnames(HUB_HPA, c("gene_name", "hpa_single_cell_specificity"))

# ClearCode34 (Brooks et al. 2014, Table 1). Current symbols: C13orf1 is
# SPRYD7, FLJ23867 is QSOX1 and UNG2 is the nuclear isoform of UNG. The alias
# UNG2 is also carried by CCNO, which is run as a sensitivity mapping.
CC34 <- data.table(
  symbol_paper = c("MAPT", "STK32B", "FZD1", "RGS5", "GIPC2", "PDGFD", "EPAS1", "MAOB", "CDH5",
                   "TCEA3", "LEPROTL1", "BNIP3L", "EHBP1", "VCAM1", "PHYH", "PRKAA2", "SLC4A4",
                   "ESD", "TLR3", "NRP1", "C11orf1", "ST13", "ARNT", "C13orf1",
                   "SERPINA3", "SLC4A3", "MOXD1", "KCNN4", "ROR2", "FLJ23867", "FOXM1", "UNG2",
                   "GALNT10", "GALNT4"),
  class = c(rep("ccA", 24), rep("ccB", 10)))
CC34[, gene_name := symbol_paper]
CC34[symbol_paper == "C13orf1",  gene_name := "SPRYD7"]
CC34[symbol_paper == "FLJ23867", gene_name := "QSOX1"]
CC34[symbol_paper == "UNG2",     gene_name := "UNG"]
CC34_ALT <- copy(CC34)[symbol_paper == "UNG2", gene_name := "CCNO"]

# ---- observed expression of the marker and classifier genes ----
# Columns are selected by the file_id kept in 01 (four barcodes carry two
# files). Symbols are resolved among protein-coding genes only (RGS5 also
# names a lncRNA).
raw <- readRDS(file.path(CACHE_DIR, "expr_raw.rds"))
ann_pc <- as.data.table(raw$gene_ann)[gene_type == "protein_coding"]
ann_pc <- ann_pc[!duplicated(gene_name)]
sym2id <- function(s) ann_pc$gene_id[match(s, ann_pc$gene_name)]
PT_MARKERS[, gene_id := sym2id(gene_name)]
CC34[, gene_id := sym2id(gene_name)]; CC34_ALT[, gene_id := sym2id(gene_name)]
HUB_HPA[, gene_id := sym2id(gene_name)]
if (anyNA(PT_MARKERS$gene_id)) stop("PT markers missing from the annotation: ",
                                    paste(PT_MARKERS[is.na(gene_id), gene_name], collapse = ", "))
msg("ClearCode34 genes resolved in the annotation: ", sum(!is.na(CC34$gene_id)), "/34",
    if (anyNA(CC34$gene_id)) paste0(" (missing: ", paste(CC34[is.na(gene_id), symbol_paper], collapse = ", "), ")") else "")
need_ids <- unique(na.omit(c(PT_MARKERS$gene_id, CC34$gene_id, CC34_ALT$gene_id)))
stopifnot(all(need_ids %in% rownames(raw$fpkm)))
# fpkm columns follow sample_info row order and are labelled by barcode, so
# the file_id row index is the column index (asserted below).
stopifnot(ncol(raw$fpkm) == nrow(raw$sample_info),
          identical(colnames(raw$fpkm), raw$sample_info$sample_barcode))
j <- match(cohort$file_id, raw$sample_info$file_id)
stopifnot(!anyNA(j), identical(colnames(raw$fpkm)[j], cohort$sample_barcode))
X_kirc <- t(log2(raw$fpkm[need_ids, j, drop = FALSE] + 1))
rownames(X_kirc) <- cohort$sample_barcode
rm(raw); invisible(gc(verbose = FALSE))
msg("Observed marker/classifier expression: ", nrow(X_kirc), " TCGA-KIRC samples x ", ncol(X_kirc), " genes")

# CPTAC-3 (same GENCODE v36 identifiers)
vann <- as.data.table(vd$gene_ann)
vmap <- data.table(gene_id = need_ids, v_id = NA_character_)
vmap[gene_id %in% rownames(vd$fpkm), v_id := gene_id]
if (anyNA(vmap$v_id)) {                       # symbol fallback for absent identifiers
  miss <- vmap[is.na(v_id), gene_id]
  sy   <- ann_pc$gene_name[match(miss, ann_pc$gene_id)]
  alt  <- vann[gene_type == "protein_coding" & gene_id %in% rownames(vd$fpkm)][!duplicated(gene_name)]
  vmap[is.na(v_id), v_id := alt$gene_id[match(sy, alt$gene_name)]]
  msg("CPTAC-3: ", length(miss), " identifiers absent; ",
      sum(!is.na(vmap[gene_id %in% miss, v_id])), " recovered by symbol")
}
vmap <- vmap[!is.na(v_id)]
stopifnot(identical(colnames(vd$fpkm), vd$cohort$sample_barcode))
X_cptac <- t(log2(vd$fpkm[vmap$v_id, , drop = FALSE] + 1))
# columns are relabelled with the discovery identifier so the marker and
# classifier sets are addressed identically in both cohorts
colnames(X_cptac) <- vmap$gene_id
rownames(X_cptac) <- vd$cohort$sample_barcode
msg("CPTAC-3: ", nrow(X_cptac), " samples x ", ncol(X_cptac), " genes")

# ---- module membership, expression level and coupling of every marker -------
SC <- discovery_scores(nets)
green_me <- SC$mrna[, REF_COL]
PT_MARKERS[, in_network := gene_id %in% gt_m$gene_id]
PT_MARKERS[, module := gt_m$module[match(gene_id, gt_m$gene_id)]]
PT_MARKERS[, kME := gt_m$kME[match(gene_id, gt_m$gene_id)]]
PT_MARKERS[, green_member := !is.na(module) & module == REF_MOD]
PT_MARKERS[, median_log2fpkm_kirc := apply(X_kirc[, gene_id, drop = FALSE], 2, median)]
PT_MARKERS[, frac_fpkm_gt1_kirc := colMeans(X_kirc[, gene_id, drop = FALSE] > 1)]
PT_MARKERS[, r_green_kirc := vapply(gene_id, function(g) {
  x <- X_kirc[names(green_me), g]; cor(x, green_me) }, numeric(1))]
PT_MARKERS[, in_cptac := gene_id %in% colnames(X_cptac)]
PT_MARKERS[, used_full := TRUE]
PT_MARKERS[, used_hpa := hpa_pt_named]
PT_MARKERS[, used_leaveout := !green_member]
msg("Marker set: ", nrow(PT_MARKERS), " genes; in network ", sum(PT_MARKERS$in_network),
    " (all in module ", REF_MOD, ": ", all(PT_MARKERS[in_network == TRUE, green_member]), "); ",
    "HPA names proximal tubule for ", sum(PT_MARKERS$hpa_pt_named), "; leave-out (non-",
    REF_MOD, ") set ", sum(PT_MARKERS$used_leaveout), " genes: ",
    paste(PT_MARKERS[used_leaveout == TRUE, gene_name], collapse = ", "))
setcolorder(PT_MARKERS, c("gene_name", "gene_id", "core", "hpa_pt_named", "hpa_pt_exclusivity",
                          "hpa_n_cell_types", "hpa_pt_rank", "hpa_single_cell_specificity",
                          "segment_note", "in_network", "module", "kME", "green_member",
                          "median_log2fpkm_kirc", "frac_fpkm_gt1_kirc", "r_green_kirc", "in_cptac",
                          "used_full", "used_hpa", "used_leaveout", "source"))
out1 <- copy(PT_MARKERS)
out1[, `:=`(kME = round(kME, 3), median_log2fpkm_kirc = round(median_log2fpkm_kirc, 3),
            frac_fpkm_gt1_kirc = round(frac_fpkm_gt1_kirc, 3), r_green_kirc = round(r_green_kirc, 3))]
save_tsv(out1, "34_pt_marker_set.tsv")
print(out1[, .(gene_name, core, hpa_pt_named, module, kME, green_member, r_green_kirc, used_leaveout)])
# The leave-out markers failed the network's expression and variance filters,
# so the leave-out score is measured near the detection floor.
msg("Leave-out markers, observed expression in TCGA-KIRC: median log2(FPKM+1) ",
    paste(sprintf("%s %.2f (FPKM>1 in %.0f%%)", PT_MARKERS[used_leaveout == TRUE, gene_name],
                  PT_MARKERS[used_leaveout == TRUE, median_log2fpkm_kirc],
                  100 * PT_MARKERS[used_leaveout == TRUE, frac_fpkm_gt1_kirc]), collapse = "; "),
    "; module markers median log2(FPKM+1) ",
    round(median(PT_MARKERS[used_leaveout == FALSE, median_log2fpkm_kirc]), 2))

# ---- the 20 top hubs ------------------------------------------------------------
hubs <- gt_m[module == REF_MOD][order(-kME)][1:20]
hubs[, rank := .I]
hubs[, in_pt_marker_set  := gene_id %in% PT_MARKERS$gene_id]
hubs[, in_pt_marker_core := gene_id %in% PT_MARKERS[core == TRUE, gene_id]]
hubs[, hpa_single_cell_specificity := HUB_HPA$hpa_single_cell_specificity[match(gene_name, HUB_HPA$gene_name)]]
hubs[, hpa_names_proximal_tubule := ifelse(is.na(hpa_single_cell_specificity), NA,
                                           grepl("Proximal tubule", hpa_single_cell_specificity))]
hubs[, hpa_n_cell_types := hpa_n_types(hpa_single_cell_specificity)]
hubs[, hpa_pt_rank := hpa_pt_pos(hpa_single_cell_specificity)]
hubs[, hpa_pt_exclusivity := hpa_grade(hpa_single_cell_specificity)]
hubs[, hpa_source := ifelse(is.na(hpa_single_cell_specificity), NA_character_, HPA_SRC)]
hubs[, clearcode34_class := CC34$class[match(gene_name, CC34$gene_name)]]
hubs[, in_clearcode34 := !is.na(clearcode34_class)]
n_hub_all <- gt_m[module == REF_MOD, .N]
hub_out <- hubs[, .(rank, gene_name, gene_id, module, kME = round(kME, 3), in_pt_marker_set,
                    in_pt_marker_core, hpa_names_proximal_tubule, hpa_pt_exclusivity,
                    hpa_n_cell_types, hpa_pt_rank, in_clearcode34, clearcode34_class,
                    hpa_single_cell_specificity, hpa_source)]
save_tsv(hub_out, "34_green_hubs_pt_flag.tsv")
msg("Top 20 hubs of ", REF_MOD, " (", n_hub_all, " genes): ", sum(hubs$in_pt_marker_set),
    " in the marker set; HPA single-cell specificity names proximal tubule cells for ",
    sum(hubs$hpa_names_proximal_tubule, na.rm = TRUE), " of ", sum(!is.na(hubs$hpa_names_proximal_tubule)), " looked up; ",
    sum(hubs$in_clearcode34), " also a ClearCode34 gene (",
    paste(hubs[in_clearcode34 == TRUE, paste0(gene_name, " ", clearcode34_class)], collapse = ", "), ")")
# Report the permissive HPA count with its exclusivity breakdown.
excl <- hubs[, .N, by = hpa_pt_exclusivity][order(-N)]
msg("  HPA exclusivity of the ", sum(hubs$hpa_names_proximal_tubule, na.rm = TRUE),
    " hubs whose string names proximal tubule: ",
    paste(excl[hpa_pt_exclusivity != "PT not named", sprintf("%s: %d", hpa_pt_exclusivity, N)], collapse = "; "),
    "; cell types named per hub: median ", median(hubs$hpa_n_cell_types, na.rm = TRUE),
    ", range ", min(hubs$hpa_n_cell_types, na.rm = TRUE), "-", max(hubs$hpa_n_cell_types, na.rm = TRUE),
    " (", hubs[which.max(hpa_n_cell_types), gene_name], " names the most). PT-exclusive: ",
    paste(hubs[hpa_pt_exclusivity == "PT-exclusive", gene_name], collapse = ", "),
    ". No background rate is computed (the strings are hand-entered for these genes only), so this is",
    " composition evidence, not a test of enrichment.")
# ---- where the ClearCode34 genes sit in the protein-coding network ----
# If many ccA genes are inside the module, the ccA/ccB class and the module
# eigengene are not independent readings.
cc_out <- copy(CC34)[, .(symbol_paper, gene_name, class, gene_id)]
cc_out[, in_network := gene_id %in% gt_m$gene_id]
cc_out[, module := gt_m$module[match(gene_id, gt_m$gene_id)]]
cc_out[, kME := round(gt_m$kME[match(gene_id, gt_m$gene_id)], 3)]
cc_out[, in_reference_module := !is.na(module) & module == REF_MOD]
cc_out[, r_green_kirc := round(vapply(gene_id, function(g)
  if (g %in% colnames(X_kirc)) cor(X_kirc[names(green_me), g], green_me) else NA_real_, numeric(1)), 3)]
cc_out[, source := "Brooks et al. 2014 Eur Urol 66:77-84, Table 1"]
save_tsv(cc_out, "34_clearcode34_gene_membership.tsv")
n_cc_in_mod <- cc_out[in_reference_module == TRUE, .N]
msg("ClearCode34 genes inside the network: ", cc_out[in_network == TRUE, .N], "/34; in module ",
    REF_MOD, ": ", n_cc_in_mod, " (", cc_out[in_reference_module == TRUE, sum(class == "ccA")],
    " of the 24 ccA genes, ", cc_out[in_reference_module == TRUE, sum(class == "ccB")],
    " of the 10 ccB genes): ", paste(cc_out[in_reference_module == TRUE, gene_name], collapse = ", "))
print(hub_out[, .(rank, gene_name, kME, in_pt_marker_set, hpa_names_proximal_tubule)])

# ---- 2. Marker scores and the ClearCode34 axis, both cohorts ----
banner("2 | Proximal-tubule marker scores and ClearCode34 ccA/ccB axis")
id2sym <- setNames(ann_pc$gene_name[match(need_ids, ann_pc$gene_id)], need_ids)
# A gene with zero variance in a cohort is dropped rather than set to z = 0,
# so the reported gene counts reflect the genes actually used.
zscore <- function(X, label = "") {
  sdv <- apply(X, 2, sd)
  bad <- !is.finite(sdv) | sdv == 0
  if (any(bad)) msg(label, ": dropping ", sum(bad), " constant / non-finite gene(s) before z scoring: ",
                    paste(ifelse(is.na(id2sym[colnames(X)[bad]]), colnames(X)[bad],
                                 id2sym[colnames(X)[bad]]), collapse = ", "))
  else msg(label, ": no constant or non-finite genes among the ", ncol(X),
           " marker/classifier genes (minimum SD ", signif(min(sdv), 3), ")")
  Z <- scale(X[, !bad, drop = FALSE])
  attr(Z, "scaled:center") <- NULL; attr(Z, "scaled:scale") <- NULL; Z
}
pt_scores <- function(Z) {
  g_full <- intersect(PT_MARKERS[used_full == TRUE, gene_id], colnames(Z))
  g_hpa  <- intersect(PT_MARKERS[used_hpa == TRUE, gene_id], colnames(Z))
  g_lo   <- intersect(PT_MARKERS[used_leaveout == TRUE, gene_id], colnames(Z))
  data.table(sample_barcode = rownames(Z),
             pt_full = rowMeans(Z[, g_full, drop = FALSE]), n_pt_full = length(g_full),
             pt_hpa = rowMeans(Z[, g_hpa, drop = FALSE]),  n_pt_hpa = length(g_hpa),
             pt_leaveout = rowMeans(Z[, g_lo, drop = FALSE]), n_pt_leaveout = length(g_lo))
}
# The published ClearCode34 PAM centroids are not distributed, so ccA/ccB is
# assigned by template correlation: +1 for ccA and -1 for ccB genes on
# within-cohort z scores. With a two-valued template, cor > 0 is equivalent to
# mean z(ccA) > mean z(ccB), a nearest-centroid rule with cohort class means.
# ccB_score = mean z(ccB) - mean z(ccA), higher is more ccB-like (adverse).
ccab_axis <- function(Z, cc) {
  g <- cc[!is.na(gene_id) & gene_id %in% colnames(Z)]
  Zg <- Z[, g$gene_id, drop = FALSE]
  tmpl <- ifelse(g$class == "ccA", 1, -1)
  r  <- apply(Zg, 1, function(z) cor(z, tmpl))
  sA <- rowMeans(Zg[, g$class == "ccA", drop = FALSE]); sB <- rowMeans(Zg[, g$class == "ccB", drop = FALSE])
  list(n_genes = nrow(g), n_ccA = sum(g$class == "ccA"), n_ccB = sum(g$class == "ccB"),
       r = r, score = sB - sA, class = ifelse(r > 0, "ccA", "ccB"))
}
score_cohort <- function(X, cohort_name) {
  Z <- zscore(X, cohort_name)
  pt <- pt_scores(Z)
  a  <- ccab_axis(Z, CC34); b <- ccab_axis(Z, CC34_ALT)
  msg(cohort_name, ": PT score on ", pt$n_pt_full[1], " (full), ", pt$n_pt_hpa[1], " (HPA-named), ",
      pt$n_pt_leaveout[1], " (leave-out) markers; ClearCode34 on ", a$n_genes, " genes (",
      a$n_ccA, " ccA, ", a$n_ccB, " ccB); ccA ", sum(a$class == "ccA"), ", ccB ", sum(a$class == "ccB"),
      " (", round(100 * mean(a$class == "ccB"), 1), "% ccB); UNG->CCNO sensitivity agrees for ",
      round(100 * mean(a$class == b$class), 1), "%")
  pt[, `:=`(cohort = cohort_name, cc34_n_genes = a$n_genes, cc34_r_template = a$r,
            ccB_score = a$score, ccAB_class = a$class, ccAB_class_CCNO = b$class)]
  pt
}
PT_K <- score_cohort(X_kirc,  "TCGA-KIRC")
PT_V <- score_cohort(X_cptac, "CPTAC-3")
per_patient <- rbind(PT_K, PT_V)
# computed outside the data.table [ so that `cohort` is not shadowed by the
# table's own cohort column
pat <- c(cohort$patient[match(PT_K$sample_barcode, cohort$sample_barcode)],
         vd$cohort$patient[match(PT_V$sample_barcode, vd$cohort$sample_barcode)])
per_patient[, patient := pat]
setcolorder(per_patient, c("cohort", "patient", "sample_barcode"))
out2 <- copy(per_patient)
for (v in c("pt_full", "pt_hpa", "pt_leaveout", "cc34_r_template", "ccB_score")) set(out2, j = v, value = round(out2[[v]], 4))
save_tsv(out2, "34_ccAB_classification.tsv")

# ---- 3. Analysis frames ----
# Observed-space eigengenes: the discovery eigengenes are built on matrices
# residualised on the STAR metrics, so they cannot show a quality or
# normal-admixture imprint. The observed-space score keeps the membership and
# refits the rotation on the observed matrix (as in 12 and 27).
banner("3 | Observed-space eigengenes and the analysis frames")
LOBS_M <- fit_module_loadings(obs_expr(nets$mrna), nets$mrna$gene_tbl, "mRNA_ME")
LOBS_L <- fit_module_loadings(obs_expr(nets$lnc),  nets$lnc$gene_tbl,  "lnc_ME")
SO_K   <- list(mrna = score_modules(obs_expr(nets$mrna), LOBS_M),
               lnc  = score_modules(obs_expr(nets$lnc),  LOBS_L))
stopifnot(REF_COL %in% colnames(SO_K$mrna),
          all(paste0("lnc_ME", LNC_MODS) %in% colnames(SO_K$lnc)))
msg("Observed-space rotation refitted: protein-coding r(observed, residualised) for ", REF_MOD,
    " = ", round(cor(SO_K$mrna[, REF_COL], SC$mrna[rownames(SO_K$mrna), REF_COL]), 3),
    "; lncRNA ", paste(sprintf("%s %.3f", LNC_MODS,
      vapply(LNC_MODS, function(m) cor(SO_K$lnc[, paste0("lnc_ME", m)],
             SC$lnc[rownames(SO_K$lnc), paste0("lnc_ME", m)]), numeric(1))), collapse = ", "))
# CPTAC-3 observed matrices, scored through the same rotations with
# cohort gene standardisation (the transport rule used in 09).
XV_M <- t(log2(vd$fpkm[intersect(colnames(obs_expr(nets$mrna)), rownames(vd$fpkm)), , drop = FALSE] + 1))
XV_L <- t(log2(vd$fpkm[intersect(colnames(obs_expr(nets$lnc)),  rownames(vd$fpkm)), , drop = FALSE] + 1))
SO_V <- list(mrna = score_modules(XV_M, LOBS_M), lnc = score_modules(XV_L, LOBS_L))
msg("CPTAC-3 observed-space scoring on ", ncol(XV_M), " protein-coding and ", ncol(XV_L),
    " lncRNA network genes present of ", ncol(obs_expr(nets$mrna)), " and ", ncol(obs_expr(nets$lnc)))
rm(XV_M, XV_L); invisible(gc(verbose = FALSE))

# TCGA-KIRC frame: clinical data, ESTIMATE, slide tumour-nuclei percentage
# (pooled per sample in 23), discovery and observed-space eigengenes and the
# scores above. Covariates are scaled over all rows before complete-case
# restriction, as in 03.
sl <- as.data.table(bio$slides)
i_sl <- match(cohort$sample_barcode, sl$sample_barcode)
if (all(is.na(i_sl))) i_sl <- match(substr(cohort$sample_barcode, 1, 15), substr(sl$sample_barcode, 1, 15))
D <- cohort[, .(sample_barcode, os_time, os_event, age,
                sex = factor(as.character(sex), levels = c("female", "male")),
                T_stage = as.numeric(T_stage), N_pos = as.numeric(N_pos), M1 = as.numeric(M1),
                grade_o = as.numeric(grade_num), pct_noFeature, pct_multimapping, assigned_reads)]
D[, pct_tumor_nuclei := sl$pct_tumor_nuclei_mean[i_sl]]
D <- merge(D, est[, .(sample_barcode, StromalScore, ImmuneScore, ESTIMATEScore)], by = "sample_barcode", all.x = TRUE)
D[, `:=`(stromal = as.numeric(scale(StromalScore)), immune = as.numeric(scale(ImmuneScore)),
         noFeat = as.numeric(scale(pct_noFeature)), multimap = as.numeric(scale(pct_multimapping)),
         depth = as.numeric(scale(log10(assigned_reads))))]
D[, green := SC$mrna[match(sample_barcode, rownames(SC$mrna)), REF_COL]]
D[, green_obs := SO_K$mrna[match(sample_barcode, rownames(SO_K$mrna)), REF_COL]]
for (m in LNC_MODS) {
  set(D, j = m, value = SC$lnc[match(D$sample_barcode, rownames(SC$lnc)), paste0("lnc_ME", m)])
  set(D, j = paste0(m, "_obs"),
      value = SO_K$lnc[match(D$sample_barcode, rownames(SO_K$lnc)), paste0("lnc_ME", m)])
}
D <- merge(D, PT_K[, .(sample_barcode, pt_full, pt_hpa, pt_leaveout, ccB_score, ccAB_class)],
           by = "sample_barcode", all.x = TRUE)
D[, ccB_class := as.numeric(ccAB_class == "ccB")]
msg("TCGA-KIRC frame: ", nrow(D), " rows; green eigengene for ", sum(!is.na(D$green)),
    ", lncRNA eigengenes for ", sum(!is.na(D$blue)), ", tumour-nuclei % for ", sum(!is.na(D$pct_tumor_nuclei)))

# CPTAC-3: locked-model Ev (cohort-standardised, 199 patients), ESTIMATE from
# 08 (ESTIMATE score = stromal + immune), STAR metrics from validation_dataset.rds.
vest <- as.data.table(readRDS(file.path(CACHE_DIR, "valid_estimate.rds")))
vcl  <- as.data.table(L$vcl); vq <- as.data.table(vd$qc)
stopifnot(all(vcl$sample_barcode %in% rownames(L$Ev)))
V <- vcl[, .(sample_barcode, os_time, os_event, age,
             sex = factor(as.character(sex), levels = c("female", "male")),
             T_stage = as.numeric(T_stage), N_pos = as.numeric(N_pos), M1 = as.numeric(M1),
             grade_o = as.numeric(grade_num))]
V[, pct_noFeature := vq$pct_noFeature[match(sample_barcode, vq$sample_barcode)]]
V[, StromalScore := vest$StromalScore[match(sample_barcode, vest$sample_barcode)]]
V[, ImmuneScore  := vest$ImmuneScore[match(sample_barcode, vest$sample_barcode)]]
V[, ESTIMATEScore := StromalScore + ImmuneScore]
V[, pct_tumor_nuclei := NA_real_]
V[, green := L$Ev[match(sample_barcode, rownames(L$Ev)), REF_COL]]
V[, green_obs := SO_V$mrna[match(sample_barcode, rownames(SO_V$mrna)), REF_COL]]
for (m in LNC_MODS) {
  set(V, j = m, value = L$Ev[match(V$sample_barcode, rownames(L$Ev)), paste0("lnc_ME", m)])
  set(V, j = paste0(m, "_obs"),
      value = SO_V$lnc[match(V$sample_barcode, rownames(SO_V$lnc)), paste0("lnc_ME", m)])
}
V <- merge(V, PT_V[, .(sample_barcode, pt_full, pt_hpa, pt_leaveout, ccB_score, ccAB_class)],
           by = "sample_barcode", all.x = TRUE)
V[, ccB_class := as.numeric(ccAB_class == "ccB")]
msg("CPTAC-3 frame: ", nrow(V), " rows, ", sum(V$os_event), " events; ESTIMATE for ", sum(!is.na(V$ESTIMATEScore)))

# ---- class summary ---------------------------------------------------------
class_summary <- function(d, pp, cohort_name) {
  w <- function(x, g) tryCatch(t.test(x ~ g)$p.value, error = function(e) NA_real_)
  lr <- survdiff(Surv(os_time, os_event) ~ ccAB_class, data = d[!is.na(ccAB_class)])
  cu <- summary(coxph(Surv(os_time, os_event) ~ ccB_class, data = d))
  data.table(cohort = cohort_name, n_scored = nrow(pp), n_in_frame = nrow(d),
             n_ccA = sum(pp$ccAB_class == "ccA"), n_ccB = sum(pp$ccAB_class == "ccB"),
             pct_ccB = round(100 * mean(pp$ccAB_class == "ccB"), 1),
             pct_agreement_UNG_vs_CCNO = round(100 * mean(pp$ccAB_class == pp$ccAB_class_CCNO), 1),
             n_with_green = sum(!is.na(d$green)),
             mean_green_ccA = round(mean(d[ccAB_class == "ccA", green], na.rm = TRUE), 3),
             mean_green_ccB = round(mean(d[ccAB_class == "ccB", green], na.rm = TRUE), 3),
             p_green_by_class = signif(w(d$green, d$ccAB_class), 3),
             mean_pt_full_ccA = round(mean(d[ccAB_class == "ccA", pt_full], na.rm = TRUE), 3),
             mean_pt_full_ccB = round(mean(d[ccAB_class == "ccB", pt_full], na.rm = TRUE), 3),
             p_pt_full_by_class = signif(w(d$pt_full, d$ccAB_class), 3),
             n_surv = cu$n, events_surv = cu$nevent,
             HR_ccB_unadjusted = round(cu$conf.int[1, 1], 3), lo = round(cu$conf.int[1, 3], 3),
             hi = round(cu$conf.int[1, 4], 3), p_cox = signif(cu$coefficients[1, 5], 3),
             p_logrank = signif(pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE), 3))
}
cls <- rbind(class_summary(D, PT_K, "TCGA-KIRC"), class_summary(V, PT_V, "CPTAC-3"))
save_tsv(cls, "34_ccAB_class_summary.tsv"); print(cls, row.names = FALSE)

# ---- 4. Correlations ----
banner("4 | Correlations of the eigengenes with the marker scores and composition")
cor_pair <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < 10 || sd(x[ok]) == 0 || sd(y[ok]) == 0)
    return(list(n = n, r = NA_real_, r_lo = NA_real_, r_hi = NA_real_, p_r = NA_real_, rho = NA_real_, p_rho = NA_real_))
  ct <- cor.test(x[ok], y[ok])
  cs <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(n = n, r = unname(ct$estimate), r_lo = ct$conf.int[1], r_hi = ct$conf.int[2], p_r = ct$p.value,
       rho = unname(cs$estimate), p_rho = cs$p.value)
}
FEATURES  <- c(green = paste("Protein-coding", REF_MOD), blue = "lncRNA blue",
               greenyellow = "lncRNA greenyellow", turquoise = "lncRNA turquoise",
               green_obs = paste("Protein-coding", REF_MOD), blue_obs = "lncRNA blue",
               greenyellow_obs = "lncRNA greenyellow", turquoise_obs = "lncRNA turquoise",
               pt_full = "PT marker score (full)")
# residualised = discovery eigengene, observed = same membership with the rotation
# refitted on the observed matrix. The marker score is always observed.
FEAT_SPACE <- ifelse(names(FEATURES) == "pt_full", "observed",
                     ifelse(grepl("_obs$", names(FEATURES)), "observed", "residualised"))
names(FEAT_SPACE) <- names(FEATURES)
VARIABLES <- c(pt_full = "PT marker score (full)", pt_hpa = "PT marker score (HPA-named)",
               pt_leaveout = paste0("PT marker score (leave-out, non-", REF_MOD, ")"),
               ccB_score = "ClearCode34 ccB minus ccA score", ccB_class = "ClearCode34 class (ccB = 1)",
               ESTIMATEScore = "ESTIMATE score", StromalScore = "ESTIMATE stromal score",
               ImmuneScore = "ESTIMATE immune score", pct_tumor_nuclei = "slide tumour nuclei %",
               pct_noFeature = "non-feature read fraction")
cor_tbl <- rbindlist(lapply(list(list(D, "TCGA-KIRC"), list(V, "CPTAC-3")), function(z) {
  d <- z[[1]]
  rbindlist(lapply(names(FEATURES), function(f) rbindlist(lapply(names(VARIABLES), function(v) {
    if (f == v) return(NULL)
    s <- cor_pair(d[[f]], d[[v]])
    data.table(cohort = z[[2]], feature = FEATURES[[f]], space = FEAT_SPACE[[f]],
               variable = VARIABLES[[v]], n = s$n,
               pearson_r = round(s$r, 3), pearson_lo = round(s$r_lo, 3), pearson_hi = round(s$r_hi, 3),
               pearson_p = signif(s$p_r, 3), spearman_rho = round(s$rho, 3), spearman_p = signif(s$p_rho, 3))
  }))))
}))
cor_tbl <- cor_tbl[!is.na(pearson_r) | cohort == "TCGA-KIRC"]
cor_tbl[, fdr_pearson  := signif(p.adjust(pearson_p,  "BH"), 3), by = .(cohort, space)]
cor_tbl[, fdr_spearman := signif(p.adjust(spearman_p, "BH"), 3), by = .(cohort, space)]
save_tsv(cor_tbl, "34_pt_score_correlations.tsv")
print(dcast(cor_tbl, cohort + space + feature ~ variable, value.var = "pearson_r"), row.names = FALSE)
# Each coefficient is printed with the p-value of its own test.
msg("Non-feature fraction (TCGA-KIRC), Spearman rho (spearman_p): ",
    paste(cor_tbl[cohort == "TCGA-KIRC" & variable == VARIABLES[["pct_noFeature"]],
                  sprintf("%s [%s] rho %.3f (p %s)", feature, space, spearman_rho,
                          format(spearman_p, digits = 3, trim = TRUE))], collapse = "; "))
msg("Non-feature fraction (TCGA-KIRC), Pearson r (pearson_p): ",
    paste(cor_tbl[cohort == "TCGA-KIRC" & variable == VARIABLES[["pct_noFeature"]],
                  sprintf("%s [%s] r %.3f (p %s)", feature, space, pearson_r,
                          format(pearson_p, digits = 3, trim = TRUE))], collapse = "; "))
msg("Slide tumour-nuclei % (TCGA-KIRC), Pearson r (pearson_p): ",
    paste(cor_tbl[cohort == "TCGA-KIRC" & variable == VARIABLES[["pct_tumor_nuclei"]],
                  sprintf("%s [%s] r %.3f (p %s, n=%d)", feature, space, pearson_r,
                          format(pearson_p, digits = 3, trim = TRUE), n)], collapse = "; "))
msg("PT marker score (full) vs the module and the lncRNA eigengenes (TCGA-KIRC): ",
    paste(cor_tbl[cohort == "TCGA-KIRC" & variable == VARIABLES[["pt_full"]],
                  sprintf("%s [%s] r %.3f (p %s), rho %.3f (p %s), n=%d", feature, space,
                          pearson_r, format(pearson_p, digits = 3, trim = TRUE),
                          spearman_rho, format(spearman_p, digits = 3, trim = TRUE), n)], collapse = "; "))

# ---- 5. Cox models ----
banner("5 | Cox models: does the lncRNA signal survive the proximal-tubule programme?")
CLIN <- c("age", "sex", "T_stage", "N_pos", "M1", "grade_o")
PRIN <- c(CLIN, "stromal", "immune", "noFeat", "multimap", "depth")
# Model labels include the specification.
SPEC_LABEL <- c(clinical = "clinical", principal = "principal covariates")
# One model on the complete cases of its own covariates. Continuous score terms are
# scaled to unit SD within that set, class indicators stay 0/1.
# Row types (column row_type):
#   term      HR per SD (or ccB vs ccA), with the variance inflation factor
#   contrast  Wald contrast beta_a - beta_b for each pair of score terms (HR
#             column = ratio of the two HRs), with their correlation
#   LRT       nested likelihood-ratio test for adding each term to the others
# Near-collinear score terms can split one shared effect between their Wald p-values,
# so the contrast and the LRT are the interpretable tests.
vif_terms <- function(f) {
  b <- coef(f)
  v <- tryCatch(diag(solve(cov2cor(vcov(f)))), error = function(e) NULL)
  if (is.null(v) || length(v) != length(b)) return(setNames(rep(NA_real_, length(b)), names(b)))
  setNames(as.numeric(v), names(b))
}
fit_cox <- function(d, label_base, score_terms, covs, spec, cohort_name, class_terms = character(0),
                    restrict = character(0)) {
  label <- paste0(label_base, " + ", SPEC_LABEL[[spec]])
  # `restrict` names variables that must be observed but are not fitted, so a
  # one-score model can use the complete-case set of a two-score model (the
  # module eigengene is missing outside the network).
  vars  <- c("os_time", "os_event", score_terms, class_terms, covs, restrict)
  dc <- copy(d[complete.cases(d[, ..vars])])
  for (v in score_terms) set(dc, j = v, value = as.numeric(scale(dc[[v]])))
  terms_all <- c(score_terms, class_terms)
  rhs <- paste(c(terms_all, covs), collapse = " + ")
  f   <- coxph(as.formula(paste("Surv(os_time, os_event) ~", rhs)), data = dc)
  s   <- summary(f); V <- vcov(f); b <- coef(f); vf <- vif_terms(f)
  zp  <- tryCatch(cox.zph(f)$table, error = function(e) NULL)
  # Events per parameter and a MIN_EPV flag on every row. Models below MIN_EPV are
  # labelled descriptive.
  epv_v <- s$nevent / length(b)
  base <- function(...) data.table(cohort = cohort_name, model = label, spec = spec, ...,
                                   n = s$n, events = s$nevent,
                                   epv = round(epv_v, 1), n_parameters = length(b),
                                   min_epv = MIN_EPV, below_min_epv = epv_v < MIN_EPV,
                                   inference_status = if (epv_v < MIN_EPV)
                                     paste0("descriptive: ", round(epv_v, 1),
                                            " events per parameter, below MIN_EPV = ", MIN_EPV)
                                     else "inferential",
                                   C = round(s$concordance[1], 3),
                                   covariates = paste(covs, collapse = "+"))
  row_term <- rbindlist(lapply(terms_all, function(tm) base(
    row_type = "term", term = tm,
    unit = if (tm %in% class_terms) "ccB vs ccA" else "per 1 SD",
    HR = round(s$conf.int[tm, 1], 3), lo = round(s$conf.int[tm, 3], 3), hi = round(s$conf.int[tm, 4], 3),
    estimate = round(unname(b[tm]), 4), se = round(sqrt(V[tm, tm]), 4),
    statistic = round(unname(b[tm]) / sqrt(V[tm, tm]), 3), df = NA_integer_,
    p = signif(s$coefficients[tm, 5], 3), vif = round(unname(vf[tm]), 2), pair_r = NA_real_,
    ph_p = if (!is.null(zp) && tm %in% rownames(zp)) signif(zp[tm, "p"], 3) else NA_real_)))
  row_con <- NULL
  if (length(score_terms) >= 2) row_con <- rbindlist(lapply(
    combn(score_terms, 2, simplify = FALSE), function(p2) {
      e  <- unname(b[p2[1]] - b[p2[2]])
      se <- sqrt(V[p2[1], p2[1]] + V[p2[2], p2[2]] - 2 * V[p2[1], p2[2]])
      z  <- e / se
      base(row_type = "contrast", term = paste(p2[1], "-", p2[2]),
           unit = "difference in log HR per 1 SD (HR column = ratio of the two hazard ratios)",
           HR = round(exp(e), 3), lo = round(exp(e - 1.96 * se), 3), hi = round(exp(e + 1.96 * se), 3),
           estimate = round(e, 4), se = round(se, 4), statistic = round(z, 3), df = NA_integer_,
           p = signif(2 * pnorm(-abs(z)), 3), vif = NA_real_,
           pair_r = round(cor(dc[[p2[1]]], dc[[p2[2]]]), 3), ph_p = NA_real_)
    }))
  row_lrt <- rbindlist(lapply(terms_all, function(tm) {
    rest <- setdiff(terms_all, tm)
    f0 <- coxph(as.formula(paste("Surv(os_time, os_event) ~",
                                 paste(c(rest, covs), collapse = " + "))), data = dc)
    ch <- 2 * as.numeric(logLik(f) - logLik(f0))
    base(row_type = "LRT",
         term = paste0("add ", tm, " to ",
                       if (length(rest)) paste(rest, collapse = "+") else "covariates only"),
         unit = "nested likelihood-ratio test, 1 df",
         HR = NA_real_, lo = NA_real_, hi = NA_real_, estimate = NA_real_, se = NA_real_,
         statistic = round(ch, 3), df = 1L, p = signif(pchisq(ch, 1, lower.tail = FALSE), 3),
         vif = NA_real_,
         pair_r = if (length(rest) == 1L && all(c(tm, rest) %in% score_terms))
                    round(cor(dc[[tm]], dc[[rest]]), 3) else NA_real_,
         ph_p = NA_real_)
  }))
  rbind(row_term, row_con, row_lrt)
}
# Each entry: label, per-SD score terms, class term (or NULL), and optionally
# variables that must be observed but are not fitted (see fit_cox `restrict`).
MODELS <- list(
  list("PT score (full)",                       "pt_full",                      NULL),
  list("PT score (leave-out)",                  "pt_leaveout",                  NULL),
  list("PT score (HPA-named)",                  "pt_hpa",                       NULL),
  list("PT score (full), green analysis set",   "pt_full",                      NULL, "green"),
  list("blue, green analysis set",              "blue",                         NULL, "green"),
  list("green",                                 "green",                        NULL),
  list("green + PT score (full)",               c("green", "pt_full"),          NULL),
  list("green + PT score (leave-out)",          c("green", "pt_leaveout"),      NULL),
  list("blue",                                  "blue",                         NULL),
  list("blue + PT score (full)",                c("blue", "pt_full"),           NULL),
  list("blue + PT score (leave-out)",           c("blue", "pt_leaveout"),       NULL),
  list("blue + green",                          c("blue", "green"),             NULL),
  list("blue + green + PT score (full)",        c("blue", "green", "pt_full"),  NULL),
  list("greenyellow",                           "greenyellow",                  NULL),
  list("greenyellow + PT score (full)",         c("greenyellow", "pt_full"),    NULL),
  list("greenyellow + PT score (leave-out)",    c("greenyellow", "pt_leaveout"), NULL),
  list("greenyellow + green",                   c("greenyellow", "green"),      NULL),
  list("greenyellow + green + PT score (full)", c("greenyellow", "green", "pt_full"), NULL),
  list("turquoise",                             "turquoise",                    NULL),
  list("turquoise + PT score (full)",           c("turquoise", "pt_full"),      NULL),
  list("turquoise + green",                     c("turquoise", "green"),        NULL),
  list("turquoise + green + PT score (full)",   c("turquoise", "green", "pt_full"), NULL),
  list("ccA/ccB class",                         character(0),                   "ccB_class"),
  list("ccB score",                             "ccB_score",                    NULL),
  list("green + ccA/ccB class",                 "green",                        "ccB_class"),
  list("PT score (full) + ccA/ccB class",       "pt_full",                      "ccB_class"),
  list("green + PT score (full) + ccA/ccB class", c("green", "pt_full"),        "ccB_class"))
run_models <- function(d, cohort_name, specs) rbindlist(lapply(names(specs), function(sp)
  rbindlist(lapply(MODELS, function(mm)
    fit_cox(d, mm[[1]], mm[[2]], specs[[sp]], sp, cohort_name,
            class_terms = if (is.null(mm[[3]])) character(0) else mm[[3]],
            restrict    = if (length(mm) < 4L) character(0) else mm[[4]])))))
cox_tbl <- rbind(run_models(D, "TCGA-KIRC", list(clinical = CLIN, principal = PRIN)),
                 run_models(V, "CPTAC-3",   list(clinical = CLIN)))
# "green" is the working label of the reference module. Replace it with the located
# colour, leaving "greenyellow" untouched.
gsub_ref <- function(x) gsub("(?<![[:alnum:]])green(?![[:alnum:]])", REF_MOD, x, perl = TRUE)
cox_tbl[, model := gsub_ref(model)]
cox_tbl[, term  := gsub_ref(term)]
# p-values are nominal. BH within cohort x specification x row_type.
cox_tbl[, fdr := signif(p.adjust(p, "BH"), 3), by = .(cohort, spec, row_type)]
cox_tbl[, fdr_family_n := .N, by = .(cohort, spec, row_type)]
setcolorder(cox_tbl, c("cohort", "model", "spec", "row_type", "term", "unit", "n", "events", "epv",
                       "n_parameters", "min_epv", "below_min_epv", "inference_status",
                       "HR", "lo", "hi", "estimate", "se", "statistic", "df", "p", "fdr",
                       "fdr_family_n", "vif", "pair_r", "C", "ph_p", "covariates"))
save_tsv(cox_tbl, "34_pt_adjusted_models.tsv")
# Events-per-parameter summary per cohort.
for (ch in unique(cox_tbl$cohort)) {
  z <- cox_tbl[cohort == ch]
  msg(ch, ": ", z[below_min_epv == TRUE, .N], " of ", nrow(z), " rows (",
      uniqueN(z[below_min_epv == TRUE, paste(model, spec)]), " of ",
      uniqueN(z[, paste(model, spec)]), " models) sit below MIN_EPV = ", MIN_EPV,
      "; events per parameter ", min(z$epv), " to ", max(z$epv),
      if (z[below_min_epv == TRUE, .N] == nrow(z))
        paste0(". EVERY model in this cohort is descriptive: ", z$events[1],
               " events cannot support ", min(z$n_parameters), " to ", max(z$n_parameters),
               " parameters.") else ".")
}
print(cox_tbl[cohort == "TCGA-KIRC" & spec == "clinical" & row_type == "term",
              .(model, term, n, events, HR, lo, hi, p, fdr, vif, C)], row.names = FALSE, nrows = 200)
print(cox_tbl[row_type == "contrast",
              .(cohort, spec, model, term, pair_r, estimate, se, statistic, p, fdr)],
      row.names = FALSE, nrows = 200)
print(cox_tbl[cohort == "TCGA-KIRC" & row_type == "LRT" & !grepl("covariates only", term),
              .(spec, model, term, statistic, p, fdr)], row.names = FALSE, nrows = 200)
print(cox_tbl[cohort == "CPTAC-3" & row_type == "term",
              .(model, term, n, events, epv, below_min_epv, HR, lo, hi, p, fdr)],
      row.names = FALSE, nrows = 200)

# ---- reading the two-score models ----
report_pair <- function(cohort_name, sp, mod, x, y) {
  cn <- cox_tbl[cohort == cohort_name & spec == sp & model == mod]
  ct <- cn[row_type == "contrast" & (term == paste(x, "-", y) | term == paste(y, "-", x))]
  if (!nrow(ct)) { msg("  [", mod, "] no contrast row (model not fitted)"); return(invisible(NULL)) }
  ab <- trimws(strsplit(ct$term, "-", fixed = TRUE)[[1]]); a <- ab[1]; b <- ab[2]
  la <- cn[row_type == "LRT" & term == paste0("add ", a, " to ", b)]
  lb <- cn[row_type == "LRT" & term == paste0("add ", b, " to ", a)]
  ta <- cn[row_type == "term" & term == a]; tb <- cn[row_type == "term" & term == b]
  if (!nrow(la) || !nrow(lb)) { msg("  [", mod, "] no LRT rows"); return(invisible(NULL)) }
  msg("  [", cohort_name, " / ", sp, "] ", mod, " (n ", ta$n, ", ", ta$events, " events): ",
      a, " HR ", ta$HR, " (", ta$lo, "-", ta$hi, ", p ", ta$p, ", VIF ", ta$vif, "); ",
      b, " HR ", tb$HR, " (", tb$lo, "-", tb$hi, ", p ", tb$p, ", VIF ", tb$vif, "). ",
      "The two terms correlate r = ", ct$pair_r, ". Wald contrast ", a, " - ", b, " = ",
      ct$estimate, " (SE ", ct$se, ", z ", ct$statistic, ", p ", ct$p, "). Nested LRT: adding ",
      a, " chisq ", la$statistic, " p ", la$p, "; adding ", b, " chisq ", lb$statistic,
      " p ", lb$p, ". => ",
      if (ct$p >= 0.05)
        "the model CANNOT separate the two predictors; the split of the Wald p-values between them carries no information."
      else "the two terms are formally distinguishable.")
}
banner("5b | The near-collinear pairs, read correctly")
# Marginal effects on one patient set, so the one-score models are comparable.
same_set <- cox_tbl[cohort == "TCGA-KIRC" & spec == "clinical" & row_type == "term" &
                    model %in% c(paste0(REF_MOD, " + clinical"),
                                 paste0("PT score (full), ", REF_MOD, " analysis set + clinical"),
                                 paste0("blue, ", REF_MOD, " analysis set + clinical"))]
msg("Marginal effects on the module's own analysis set, each adjusted for the clinical comparator ",
    "and for nothing else: ",
    paste(same_set[, sprintf("%s HR %.3f (%.3f-%.3f), p %s, n %d, %d events", term, HR, lo, hi,
                             format(p, digits = 3, trim = TRUE), n, events)], collapse = "; "),
    ". None of these is an unadjusted hazard ratio.")
M_GP_C <- paste0(REF_MOD, " + PT score (full) + clinical")
M_GP_P <- paste0(REF_MOD, " + PT score (full) + principal covariates")
M_BG_C <- paste0("blue + ", REF_MOD, " + clinical")
M_BG_P <- paste0("blue + ", REF_MOD, " + principal covariates")
for (z in list(list("TCGA-KIRC", "clinical",  M_GP_C, REF_MOD, "pt_full"),
               list("TCGA-KIRC", "principal", M_GP_P, REF_MOD, "pt_full"),
               list("TCGA-KIRC", "clinical",  M_BG_C, "blue", REF_MOD),
               list("TCGA-KIRC", "principal", M_BG_P, "blue", REF_MOD)))
  report_pair(z[[1]], z[[2]], z[[3]], z[[4]], z[[5]])
msg("Reading. (a) Where the Wald contrast does not exclude zero, the two terms are statistically ",
    "INDISTINGUISHABLE as predictors on these patients: their point estimates may look very ",
    "different, but the difference is within noise, and the split of the Wald p-values between two ",
    "terms correlated at r > 0.8 is arbitrary. Nothing may be said about one term absorbing the ",
    "other from a Wald p that straddles 0.05. (b) The only asymmetric evidence is the nested LRT, ",
    "and it must be quoted with the specification it came from, because these two specifications ",
    "do not agree. (c) A single nominal p near 0.05 is one of ",
    cox_tbl[cohort == "TCGA-KIRC" & spec == "clinical" & row_type == "term", .N],
    " terms in its cohort x specification family; the fdr column, not p, is what to read.")
msg("Terms with nominal p < 0.05 in TCGA-KIRC / clinical: ",
    cox_tbl[cohort == "TCGA-KIRC" & spec == "clinical" & row_type == "term" & p < 0.05, .N],
    " of ", cox_tbl[cohort == "TCGA-KIRC" & spec == "clinical" & row_type == "term", .N],
    "; surviving BH within that family: ",
    cox_tbl[cohort == "TCGA-KIRC" & spec == "clinical" & row_type == "term" & fdr < FDR_ALPHA, .N])

# ---- 6. Overlap of the protein-coding modules with the ESTIMATE signatures ----
# The ESTIMATE signatures index fibroblast, endothelial and leucocyte content, so this
# section tests stromal/immune contamination only. Normal-kidney admixture is section 7.
banner("6 | ESTIMATE stromal / immune signature genes in each protein-coding module (a stromal/immune contamination test, NOT an admixture test)")
data(SI_geneset, package = "estimate", envir = environment())
sg <- as.matrix(SI_geneset)
sig_genes <- function(row) { g <- as.character(sg[row, -1]); unique(g[!is.na(g) & nzchar(g)]) }
SIGS <- list(stromal = sig_genes("StromalSignature"), immune = sig_genes("ImmuneSignature"))
universe <- unique(gt_m$gene_name)
mods <- c(REF_MOD, setdiff(sort(unique(gt_m$module)), c(REF_MOD, "grey")), "grey")
overlap <- rbindlist(lapply(mods, function(m) {
  gm <- unique(gt_m[module == m, gene_name])
  row <- data.table(module = m, n_module_genes = length(gm), universe_genes = length(universe))
  for (s in names(SIGS)) {
    present <- intersect(SIGS[[s]], universe); ov <- intersect(gm, present)
    p <- phyper(length(ov) - 1, length(present), length(universe) - length(present), length(gm), lower.tail = FALSE)
    row[, paste0("n_", s, "_signature") := length(SIGS[[s]])]
    row[, paste0("n_", s, "_in_universe") := length(present)]
    row[, paste0("n_", s, "_in_module") := length(ov)]
    row[, paste0("pct_module_", s) := round(100 * length(ov) / length(gm), 2)]
    row[, paste0("expected_", s) := round(length(gm) * length(present) / length(universe), 2)]
    row[, paste0("p_hyper_", s) := signif(p, 3)]
    row[, paste0(s, "_genes") := paste(sort(ov), collapse = ";")]
  }
  row
}))
overlap[, `:=`(fdr_hyper_stromal = signif(p.adjust(p_hyper_stromal, "BH"), 3),
               fdr_hyper_immune  = signif(p.adjust(p_hyper_immune,  "BH"), 3))]
save_tsv(overlap, "34_green_estimate_overlap.tsv")
print(overlap[, .(module, n_module_genes, n_stromal_in_module, expected_stromal, p_hyper_stromal,
                  n_immune_in_module, expected_immune, p_hyper_immune)], row.names = FALSE)
msg("Signature genes absent from the annotation-matched universe: stromal ",
    length(setdiff(SIGS$stromal, universe)), "/", length(SIGS$stromal), ", immune ",
    length(setdiff(SIGS$immune, universe)), "/", length(SIGS$immune),
    " (the universe is the 12,000 network genes)")
msg(REF_MOD, " module: stromal genes ", overlap[module == REF_MOD, stromal_genes],
    "; immune genes ", overlap[module == REF_MOD, immune_genes])
msg("This bears on stromal/immune contamination only. Whether the module indexes admixture of ",
    "NORMAL KIDNEY is a different question, tested in section 7 against adjacent-normal libraries.")

# ---- 7. Normal-kidney admixture: where do adjacent-normal libraries sit? ----
# If the module indexed normal PT content, a fully normal library would sit at or above
# the top of the tumour distribution. Adjacent-normal libraries cached by 25 (72 TCGA-KIRC,
# 168 CPTAC-3) are standardised on their cohort's tumour gene means and SDs and scored
# through the tumour rotation, giving tumour-score SD units. The PT marker score uses only
# network markers. The library-failure rule of 01 and 25 applies. The non-feature cap is
# not applied in CPTAC-3, whose ribo-depleted protocol places every library above it.
# Normals and tumours also differ in library quality (non-feature fraction lower in
# TCGA-KIRC normals, higher in CPTAC-3 normals), and observed lncRNA scores largely track
# that axis, so every feature is scored in two spaces (column `space`):
#   observed      observed-space rotation, no technical adjustment
#   residualised  STAR metrics removed from both matrices with tumour-fitted coefficients
#                 applied to the normals at their own metric values. Tumour gene means
#                 anchor both sets, so the tumour-normal contrast is kept.
# The residualised row is conservative: tissue differences collinear with library
# quality are removed too.
banner("7 | Adjacent-normal libraries projected onto the module (the normal-kidney admixture test)")
nl <- readRDS(file.path(CACHE_DIR, "25_normal_libraries.rds"))
MARK_NET <- PT_MARKERS[in_network == TRUE, gene_id]
msg("Marker genes inside the protein-coding network: ", length(MARK_NET), " of ", nrow(PT_MARKERS),
    " (", paste(PT_MARKERS[in_network == TRUE, gene_name], collapse = ", "), ")")

# Tumours and normals are standardised identically, anchored on the tumours
# of the cohort. retarget() keeps membership and rotation but re-anchors gene
# and score centres/scales on the tumour matrix supplied.
retarget <- function(loadings, keys, Tm) {
  out <- lapply(keys, function(k) {
    L  <- loadings[[k]]
    g  <- intersect(L$genes, colnames(Tm))
    X  <- Tm[, g, drop = FALSE]
    mu <- colMeans(X); sdv <- apply(X, 2, sd)
    if (any(!is.finite(sdv) | sdv == 0)) {
      msg("  ", k, ": ", sum(!is.finite(sdv) | sdv == 0),
          " constant gene(s) in the tumour matrix contribute zero deviation")
      sdv[!is.finite(sdv) | sdv == 0] <- 1
    }
    Z <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
    s <- as.numeric(Z %*% L$rot[g])
    list(module = L$module, genes = g, gene_center = mu, gene_scale = sdv, rot = L$rot[g],
         score_center = mean(s), score_scale = sd(s), var_explained = NA_real_)
  })
  names(out) <- keys
  out
}
tumour_units <- function(Tm, Nm, loadings, keys) {
  Lr <- retarget(loadings, keys, Tm)
  list(t = score_modules(Tm, Lr, gene_standardise = "discovery", score_standardise = "discovery"),
       n = score_modules(Nm, Lr, gene_standardise = "discovery", score_standardise = "discovery"))
}
marker_units <- function(Tm, Nm, g) {
  g  <- intersect(g, intersect(colnames(Tm), colnames(Nm)))
  mu <- colMeans(Tm[, g, drop = FALSE]); sdv <- apply(Tm[, g, drop = FALSE], 2, sd)
  keep <- is.finite(sdv) & sdv > 0
  if (any(!keep)) msg("  marker score: dropping ", sum(!keep), " constant gene(s) in the tumour matrix")
  g <- g[keep]; mu <- mu[keep]; sdv <- sdv[keep]
  zt <- rowMeans(sweep(sweep(Tm[, g, drop = FALSE], 2, mu, "-"), 2, sdv, "/"))
  zn <- rowMeans(sweep(sweep(Nm[, g, drop = FALSE], 2, mu, "-"), 2, sdv, "/"))
  m <- mean(zt); s <- sd(zt)
  list(t = (zt - m) / s, n = (zn - m) / s, n_genes = length(g))
}
# Residualise tumours and normals on the STAR metrics with one set of
# tumour-fitted coefficients. Tumour gene means are added back to both.
resid_pair <- function(Tm, Nm, Ct, Cn) {
  stopifnot(nrow(Ct) == nrow(Tm), nrow(Cn) == nrow(Nm), !anyNA(Ct), !anyNA(Cn))
  ft <- fit_technical(Tm, Ct)
  list(t = apply_technical(Tm, ft, Ct), n = apply_technical(Nm, ft, Cn))
}
# Part of the tumour-normal gap explained by library quality: the tumour OLS
# slope on the non-feature fraction times the normal-minus-tumour metric
# difference. A decomposition, not a test. Correlations are Spearman.
quality_diag <- function(st, sn, nft, nfn) {
  ct <- suppressWarnings(cor.test(st, nft, method = "spearman", exact = FALSE))
  cn <- suppressWarnings(cor.test(sn, nfn, method = "spearman", exact = FALSE))
  b  <- unname(coef(stats::lm(st ~ nft))[2])
  dn <- mean(nfn) - mean(nft)
  list(rho_t = unname(ct$estimate), p_t = ct$p.value,
       rho_n = unname(cn$estimate), p_n = cn$p.value,
       slope = b, d_nf = dn, pred = b * dn, beyond = (mean(sn) - mean(st)) - b * dn)
}
# One rule for flagging quality confounding, applied to every row.
QUALITY_RULE <- paste0("flagged when the score tracks the non-feature read fraction in the tumours ",
                       "at |Spearman rho| >= 0.3, or when the tumour-fitted slope on that metric ",
                       "alone accounts for >= 0.1 tumour SD of the normal-minus-tumour gap")
MIXING_NOTE <- paste0("pct_normal_content_for_1_tumour_sd assumes the score is linear in the normal ",
                      "fraction; mixing is linear in FPKM, not in log2(FPKM+1), so the figure is ",
                      "anti-conservative and the conclusion rests on the observed ceiling instead")
admix_row <- function(cohort_name, feature, space, st, sn, n_genes, nft, nfn, note = "") {
  ok_t <- is.finite(st) & is.finite(nft); ok_n <- is.finite(sn) & is.finite(nfn)
  st <- st[ok_t]; nft <- nft[ok_t]; sn <- sn[ok_n]; nfn <- nfn[ok_n]
  w  <- suppressWarnings(wilcox.test(sn, st))
  auc <- unname(w$statistic) / (length(sn) * length(st))   # P(normal > tumour)
  d  <- mean(sn) - mean(st)
  q  <- quality_diag(st, sn, nft, nfn)
  flag <- abs(q$rho_t) >= 0.3 || abs(q$pred) >= 0.1
  why <- sprintf(paste0("%stumour rho with the non-feature fraction %.3f; the normals differ from ",
                        "the tumours by %+.2f percentage points on that metric, which alone predicts ",
                        "%+.3f SD of the %+.3f SD gap, leaving %+.3f SD"),
                 if (space == "residualised")
                   "STAR metrics removed from both matrices with tumour-fitted coefficients; " else "",
                 q$rho_t, q$d_nf, q$pred, d, q$beyond)
  data.table(cohort = cohort_name, feature = feature, space = space, n_genes = n_genes,
             n_tumour = length(st), n_normal = length(sn),
             tumour_mean = round(mean(st), 3), tumour_sd = round(sd(st), 3),
             tumour_min = round(min(st), 3), tumour_max = round(max(st), 3),
             tumour_q95 = round(unname(quantile(st, 0.95)), 3),
             normal_mean = round(mean(sn), 3), normal_sd = round(sd(sn), 3),
             normal_min = round(min(sn), 3), normal_max = round(max(sn), 3),
             delta_tumour_sd = round(d, 3),
             auc_normal_gt_tumour = round(auc, 3),
             p_wilcoxon = signif(w$p.value, 3),
             n_normal_above_tumour_max = sum(sn > max(st)),
             n_normal_above_tumour_q95 = sum(sn > quantile(st, 0.95)),
             pct_normal_content_for_1_tumour_sd = if (d > 0) round(100 / d, 1) else NA_real_,
             note = note,
             rho_noFeature_tumour = round(q$rho_t, 3), p_rho_noFeature_tumour = signif(q$p_t, 3),
             rho_noFeature_normal = round(q$rho_n, 3), p_rho_noFeature_normal = signif(q$p_n, 3),
             delta_noFeature_pp = round(q$d_nf, 2),
             delta_predicted_from_noFeature = round(q$pred, 3),
             delta_beyond_noFeature = round(q$beyond, 3),
             quality_confounded = flag, quality_confound_rule = QUALITY_RULE,
             quality_confound_note = why, mixing_assumption = MIXING_NOTE)
}
TM_PC <- obs_expr(nets$mrna); TM_LN <- obs_expr(nets$lnc)
# Residualised rows use the cached discovery loadings (fitted on residualised
# matrices), observed rows use LOBS_M / LOBS_L from section 3.
LDIS   <- discovery_loadings(nets)
LRES_M <- LDIS[grep("^mRNA_", names(LDIS))]
LRES_L <- LDIS[grep("^lnc_",  names(LDIS))]
adm <- list()
# ---- TCGA-KIRC --------------------------------------------------------------
qk <- as.data.table(nl$sets$`TCGA-KIRC`$qc)
keep_k <- qk$assigned_reads >= MIN_ASSIGNED_READS & qk$pct_noFeature <= MAX_NOFEATURE_PCT
msg("TCGA-KIRC adjacent normals: ", nrow(qk), " cached, ", sum(!keep_k),
    " excluded by the library-failure rule, ", sum(keep_k), " scored")
NK_PC <- nl$sets$`TCGA-KIRC`$pc[keep_k, colnames(TM_PC), drop = FALSE]
NK_LN <- nl$sets$`TCGA-KIRC`$lnc[keep_k, colnames(TM_LN), drop = FALSE]
# the normal matrices are keyed on file_id in the cached row order, so the QC
# rows selected by keep_k are the rows of the matrices, in order
stopifnot(identical(rownames(NK_PC), qk$file_id[keep_k]),
          identical(rownames(NK_LN), qk$file_id[keep_k]))
qkk <- qk[keep_k]
QT_PC <- cohort[match(rownames(TM_PC), sample_barcode)]
QT_LN <- cohort[match(rownames(TM_LN), sample_barcode)]
stopifnot(!anyNA(QT_PC$pct_noFeature), !anyNA(QT_LN$pct_noFeature))
NF_T_PC <- QT_PC$pct_noFeature; NF_T_LN <- QT_LN$pct_noFeature; NF_N_K <- qkk$pct_noFeature
# Fraction of normal libraries whose STAR metrics lie inside the tumour range:
# a tumour-fitted residualisation slope is interpolation only where this is 1.
inside_range <- function(Cn, Ct) mean(apply(
  sapply(seq_len(ncol(Ct)), function(j) Cn[, j] >= min(Ct[, j]) & Cn[, j] <= max(Ct[, j])), 1, all))
msg("TCGA-KIRC non-feature read fraction: tumours mean ", round(mean(NF_T_PC), 2), "% (range ",
    round(min(NF_T_PC), 2), "-", round(max(NF_T_PC), 2), "), scored normals mean ",
    round(mean(NF_N_K), 2), "% (range ", round(min(NF_N_K), 2), "-", round(max(NF_N_K), 2),
    "); normals inside the tumour range on all three STAR metrics: ",
    round(100 * inside_range(tech_covariates(qkk), tech_covariates(QT_PC)), 1), "%")
uk <- tumour_units(TM_PC, NK_PC, LOBS_M, REF_COL)
stopifnot(REF_COL %in% colnames(uk$t), REF_COL %in% colnames(uk$n))
adm[[length(adm) + 1]] <- admix_row("TCGA-KIRC", paste("Protein-coding", REF_MOD), "observed",
                                    uk$t[, REF_COL], uk$n[, REF_COL],
                                    length(LOBS_M[[REF_COL]]$genes), NF_T_PC, NF_N_K,
                                    "observed-space rotation, tumour standardisation")
ul <- tumour_units(TM_LN, NK_LN, LOBS_L, paste0("lnc_ME", LNC_MODS))
for (m in LNC_MODS) adm[[length(adm) + 1]] <- admix_row(
  "TCGA-KIRC", paste("lncRNA", m), "observed",
  ul$t[, paste0("lnc_ME", m)], ul$n[, paste0("lnc_ME", m)],
  length(LOBS_L[[paste0("lnc_ME", m)]]$genes), NF_T_LN, NF_N_K,
  "observed-space rotation, tumour standardisation")
mk <- marker_units(TM_PC, NK_PC, MARK_NET)
adm[[length(adm) + 1]] <- admix_row("TCGA-KIRC", "PT marker score (network markers)", "observed",
                                    mk$t, mk$n, mk$n_genes, NF_T_PC, NF_N_K,
                                    "mean tumour-standardised z of the markers inside the network")
# ---- TCGA-KIRC, residualised on the three STAR metrics -----------------------
RK_PC <- resid_pair(TM_PC, NK_PC, tech_covariates(QT_PC), tech_covariates(qkk))
RK_LN <- resid_pair(TM_LN, NK_LN, tech_covariates(QT_LN), tech_covariates(qkk))
# the residualised tumour matrix is the production one built by 02
stopifnot(max(abs(RK_PC$t - nets$mrna$expr[rownames(RK_PC$t), colnames(RK_PC$t)])) < 1e-8,
          max(abs(RK_LN$t - nets$lnc$expr[rownames(RK_LN$t),  colnames(RK_LN$t)])) < 1e-8)
ukr <- tumour_units(RK_PC$t, RK_PC$n, LRES_M, REF_COL)
adm[[length(adm) + 1]] <- admix_row("TCGA-KIRC", paste("Protein-coding", REF_MOD), "residualised",
                                    ukr$t[, REF_COL], ukr$n[, REF_COL],
                                    length(LRES_M[[REF_COL]]$genes), NF_T_PC, NF_N_K,
                                    "residualised rotation, tumour-fitted STAR coefficients")
ulr <- tumour_units(RK_LN$t, RK_LN$n, LRES_L, paste0("lnc_ME", LNC_MODS))
for (m in LNC_MODS) adm[[length(adm) + 1]] <- admix_row(
  "TCGA-KIRC", paste("lncRNA", m), "residualised",
  ulr$t[, paste0("lnc_ME", m)], ulr$n[, paste0("lnc_ME", m)],
  length(LRES_L[[paste0("lnc_ME", m)]]$genes), NF_T_LN, NF_N_K,
  "residualised rotation, tumour-fitted STAR coefficients")
mkr <- marker_units(RK_PC$t, RK_PC$n, MARK_NET)
adm[[length(adm) + 1]] <- admix_row("TCGA-KIRC", "PT marker score (network markers)", "residualised",
                                    mkr$t, mkr$n, mkr$n_genes, NF_T_PC, NF_N_K,
                                    "mean tumour-standardised z of the markers, STAR metrics removed")
# agreement between the network-marker score and the 21-marker score of section 2
i_pt <- match(rownames(TM_PC), PT_K$sample_barcode)
msg("PT marker score, network markers (", mk$n_genes, ") versus the full ",
    PT_K$n_pt_full[1], "-marker score of section 2 in the tumours: Pearson r ",
    round(cor(mk$t, PT_K$pt_full[i_pt], use = "complete.obs"), 3), ", Spearman rho ",
    round(cor(mk$t, PT_K$pt_full[i_pt], method = "spearman", use = "complete.obs"), 3))
# ---- CPTAC-3 ----------------------------------------------------------------
qv <- as.data.table(nl$sets$`CPTAC-3`$qc)
keep_v <- qv$assigned_reads >= MIN_ASSIGNED_READS
msg("CPTAC-3 adjacent normals: ", nrow(qv), " cached, ", sum(!keep_v),
    " excluded on assigned reads (the non-feature cap is reported but not applied in CPTAC-3: ",
    sum(qv$pct_noFeature > MAX_NOFEATURE_PCT), " libraries exceed it), ", sum(keep_v), " scored")
qvv <- qv[keep_v]
gv_pc <- intersect(colnames(TM_PC), rownames(vd$fpkm)); gv_ln <- intersect(colnames(TM_LN), rownames(vd$fpkm))
TV_PC <- t(log2(vd$fpkm[gv_pc, , drop = FALSE] + 1)); TV_LN <- t(log2(vd$fpkm[gv_ln, , drop = FALSE] + 1))
NV_PC <- nl$sets$`CPTAC-3`$pc[keep_v, gv_pc, drop = FALSE]
NV_LN <- nl$sets$`CPTAC-3`$lnc[keep_v, gv_ln, drop = FALSE]
stopifnot(identical(rownames(NV_PC), qv$file_id[keep_v]), identical(rownames(NV_LN), qv$file_id[keep_v]))
QT_V <- vq[match(rownames(TV_PC), sample_barcode)]
stopifnot(identical(rownames(TV_PC), rownames(TV_LN)), !anyNA(QT_V$pct_noFeature))
NF_T_V <- QT_V$pct_noFeature; NF_N_V <- qvv$pct_noFeature
msg("CPTAC-3 non-feature read fraction: tumours mean ", round(mean(NF_T_V), 2), "% (range ",
    round(min(NF_T_V), 2), "-", round(max(NF_T_V), 2), "), scored normals mean ",
    round(mean(NF_N_V), 2), "% (range ", round(min(NF_N_V), 2), "-", round(max(NF_N_V), 2),
    "); normals inside the tumour range on all three STAR metrics: ",
    round(100 * inside_range(tech_covariates(qvv), tech_covariates(QT_V)), 1),
    "%. The sign of this difference is OPPOSITE to TCGA-KIRC, which is why the observed lncRNA ",
    "rows of the two cohorts point in opposite directions.")
uv <- tumour_units(TV_PC, NV_PC, LOBS_M, REF_COL)
adm[[length(adm) + 1]] <- admix_row("CPTAC-3", paste("Protein-coding", REF_MOD), "observed",
                                    uv$t[, REF_COL], uv$n[, REF_COL],
                                    length(intersect(LOBS_M[[REF_COL]]$genes, gv_pc)), NF_T_V, NF_N_V,
                                    "tumour standardisation within CPTAC-3; module genes present there only")
uvl <- tumour_units(TV_LN, NV_LN, LOBS_L, paste0("lnc_ME", LNC_MODS))
for (m in LNC_MODS) adm[[length(adm) + 1]] <- admix_row(
  "CPTAC-3", paste("lncRNA", m), "observed",
  uvl$t[, paste0("lnc_ME", m)], uvl$n[, paste0("lnc_ME", m)],
  length(intersect(LOBS_L[[paste0("lnc_ME", m)]]$genes, gv_ln)), NF_T_V, NF_N_V,
  "tumour standardisation within CPTAC-3; module genes present there only")
mv <- marker_units(TV_PC, NV_PC, MARK_NET)
adm[[length(adm) + 1]] <- admix_row("CPTAC-3", "PT marker score (network markers)", "observed",
                                    mv$t, mv$n, mv$n_genes, NF_T_V, NF_N_V,
                                    "mean tumour-standardised z of the markers inside the network")
# ---- CPTAC-3, residualised on the three STAR metrics -------------------------
RV_PC <- resid_pair(TV_PC, NV_PC, tech_covariates(QT_V), tech_covariates(qvv))
RV_LN <- resid_pair(TV_LN, NV_LN, tech_covariates(QT_V), tech_covariates(qvv))
uvr <- tumour_units(RV_PC$t, RV_PC$n, LRES_M, REF_COL)
adm[[length(adm) + 1]] <- admix_row("CPTAC-3", paste("Protein-coding", REF_MOD), "residualised",
                                    uvr$t[, REF_COL], uvr$n[, REF_COL],
                                    length(intersect(LRES_M[[REF_COL]]$genes, gv_pc)), NF_T_V, NF_N_V,
                                    "residualised rotation, CPTAC-3 tumour-fitted STAR coefficients")
uvlr <- tumour_units(RV_LN$t, RV_LN$n, LRES_L, paste0("lnc_ME", LNC_MODS))
for (m in LNC_MODS) adm[[length(adm) + 1]] <- admix_row(
  "CPTAC-3", paste("lncRNA", m), "residualised",
  uvlr$t[, paste0("lnc_ME", m)], uvlr$n[, paste0("lnc_ME", m)],
  length(intersect(LRES_L[[paste0("lnc_ME", m)]]$genes, gv_ln)), NF_T_V, NF_N_V,
  "residualised rotation, CPTAC-3 tumour-fitted STAR coefficients")
mvr <- marker_units(RV_PC$t, RV_PC$n, MARK_NET)
adm[[length(adm) + 1]] <- admix_row("CPTAC-3", "PT marker score (network markers)", "residualised",
                                    mvr$t, mvr$n, mvr$n_genes, NF_T_V, NF_N_V,
                                    "mean tumour-standardised z of the markers, STAR metrics removed")
admix <- rbindlist(adm)
# BH over the normal-versus-tumour Wilcoxon tests of this table.
admix[, fdr_wilcoxon := signif(p.adjust(p_wilcoxon, "BH"), 3)]
admix[, fdr_family := paste0("Benjamini-Hochberg within the ", .N, " normal-versus-tumour ",
                             "Wilcoxon tests of this table (", uniqueN(cohort), " cohorts x ",
                             uniqueN(feature), " features x ", uniqueN(space), " spaces)")]
# Direction survives quality adjustment: same sign in both spaces and
# residualised fdr below FDR_ALPHA. Written on both rows of each pair.
admix[, direction_survives_quality_adjustment := {
  o <- delta_tumour_sd[space == "observed"]; r <- delta_tumour_sd[space == "residualised"]
  pr <- fdr_wilcoxon[space == "residualised"]
  rep(length(o) == 1L && length(r) == 1L && sign(o) == sign(r) && pr < FDR_ALPHA, .N)
}, by = .(cohort, feature)]
setcolorder(admix, c("cohort", "feature", "space", "n_genes", "n_tumour", "n_normal",
                     "tumour_mean", "tumour_sd", "tumour_min", "tumour_max", "tumour_q95",
                     "normal_mean", "normal_sd", "normal_min", "normal_max", "delta_tumour_sd",
                     "auc_normal_gt_tumour", "p_wilcoxon", "fdr_wilcoxon",
                     "n_normal_above_tumour_max",
                     "n_normal_above_tumour_q95", "pct_normal_content_for_1_tumour_sd", "note",
                     "quality_confounded", "direction_survives_quality_adjustment",
                     "rho_noFeature_tumour", "p_rho_noFeature_tumour",
                     "rho_noFeature_normal", "p_rho_noFeature_normal", "delta_noFeature_pp",
                     "delta_predicted_from_noFeature", "delta_beyond_noFeature",
                     "quality_confound_rule", "quality_confound_note", "mixing_assumption",
                     "fdr_family"))
# Discovery cohort first, each observed row above its residualised row.
FEAT_ORDER <- c(paste("Protein-coding", REF_MOD), paste("lncRNA", LNC_MODS),
                "PT marker score (network markers)")
admix <- admix[order(match(cohort, c("TCGA-KIRC", "CPTAC-3")), match(feature, FEAT_ORDER),
                     match(space, c("observed", "residualised")))]
save_tsv(admix, "34_normal_admixture.tsv")
print(admix[, .(cohort, feature, space, n_tumour, n_normal, tumour_max, normal_mean,
                delta_tumour_sd, auc_normal_gt_tumour, p_wilcoxon, n_normal_above_tumour_max,
                quality_confounded, direction_survives_quality_adjustment)], row.names = FALSE)
print(admix[, .(cohort, feature, space, rho_noFeature_tumour, delta_noFeature_pp,
                delta_predicted_from_noFeature, delta_tumour_sd, delta_beyond_noFeature)],
      row.names = FALSE)
for (r in seq_len(nrow(admix))) with(admix[r], msg(
  "  ", cohort, " ", feature, " [", space, "]: normals sit ", delta_tumour_sd,
  " tumour SD from the tumour mean (",
  n_normal, " normals, AUC ", auc_normal_gt_tumour, ", Wilcoxon p ", p_wilcoxon, ", BH ", fdr_wilcoxon, "); ",
  n_normal_above_tumour_max, " of ", n_normal, " above the tumour maximum, ",
  n_normal_above_tumour_q95, " above the tumour 95th centile",
  if (is.na(pct_normal_content_for_1_tumour_sd)) "" else
    paste0("; under linear mixing ", pct_normal_content_for_1_tumour_sd,
           "% normal content would move the score by one tumour SD (anti-conservative)"),
  if (quality_confounded) paste0("; CONFOUNDED BY LIBRARY QUALITY: ", quality_confound_note) else ""))
pick <- function(coh, feat, sp = "observed")
  admix[cohort == coh & feature == feat & space == sp]
akg <- pick("TCGA-KIRC", paste("Protein-coding", REF_MOD))
akm <- pick("TCGA-KIRC", "PT marker score (network markers)")
avg <- pick("CPTAC-3",   paste("Protein-coding", REF_MOD))
avm <- pick("CPTAC-3",   "PT marker score (network markers)")
akg_r <- pick("TCGA-KIRC", paste("Protein-coding", REF_MOD), "residualised")
akm_r <- pick("TCGA-KIRC", "PT marker score (network markers)", "residualised")
msg("Admixture reading. Normal kidney is HIGHER on the module and on the PT markers than the ",
    "average tumour in both cohorts, which is what a proximal-tubule programme must do. The ",
    "question is whether that gradient is large enough for normal content to generate the ",
    "tumour-to-tumour spread of the module. In TCGA-KIRC, the discovery cohort in which the module ",
    "was defined and every survival model fitted, it is not: an ENTIRELY normal library sits only ",
    akg$delta_tumour_sd, " tumour SD above the mean tumour on the module (AUC ",
    akg$auc_normal_gt_tumour, "), ", akg$n_normal_above_tumour_max, " of ", akg$n_normal,
    " normals reach the tumour maximum, and the tumours span ",
    round(akg$tumour_max - akg$tumour_min, 1), " SD, so under linear mixing a varying normal ",
    "fraction cannot generate that range. The PT marker score separates further (", akm$delta_tumour_sd,
    " tumour SD, AUC ", akm$auc_normal_gt_tumour, ", ", akm$n_normal_above_tumour_max, " of ",
    akm$n_normal, " above the tumour maximum) but still not enough. In CPTAC-3 the same projection ",
    "separates much more strongly (module ", avg$delta_tumour_sd, " SD, AUC ",
    avg$auc_normal_gt_tumour, "; markers ", avm$delta_tumour_sd, " SD, AUC ",
    avm$auc_normal_gt_tumour, ", with ", avm$n_normal_above_tumour_max, " of ", avm$n_normal,
    " normals above the tumour maximum), so in that cohort normal-kidney content is clearly ",
    "resolvable and a normal-admixture contribution to the module score there cannot be excluded ",
    "on these data. The two cohorts use different library protocols and different normal sampling, ",
    "so the discovery-cohort figure is the one that bears on the discovery models. This is the ",
    "test the ESTIMATE overlap of section 6 cannot do.")
msg("The protein-coding reading does not depend on library quality: the ", REF_MOD,
    " module score barely tracks the non-feature fraction in the tumours (rho ",
    akg$rho_noFeature_tumour, "), the ", akg$delta_noFeature_pp,
    "-percentage-point difference between the normal and tumour libraries predicts only ",
    akg$delta_predicted_from_noFeature, " SD of the ", akg$delta_tumour_sd,
    " SD gap, and residualising both matrices on the three STAR metrics leaves ",
    akg_r$delta_tumour_sd, " SD (AUC ", akg_r$auc_normal_gt_tumour, ", p ", akg_r$p_wilcoxon,
    "); the PT marker score moves from ", akm$delta_tumour_sd, " to ", akm_r$delta_tumour_sd, " SD.")
banner("7b | The lncRNA rows are a library-quality comparison, not an admixture test")
for (m in LNC_MODS) {
  o <- pick("TCGA-KIRC", paste("lncRNA", m)); r <- pick("TCGA-KIRC", paste("lncRNA", m), "residualised")
  ov <- pick("CPTAC-3", paste("lncRNA", m));  rv <- pick("CPTAC-3", paste("lncRNA", m), "residualised")
  msg("  lncRNA ", m, ": TCGA-KIRC observed ", o$delta_tumour_sd, " SD (p ", o$p_wilcoxon,
      ", tumour rho with the metric ", o$rho_noFeature_tumour, ", quality alone predicts ",
      o$delta_predicted_from_noFeature, ") -> residualised ", r$delta_tumour_sd, " SD (p ",
      r$p_wilcoxon, "); CPTAC-3 observed ", ov$delta_tumour_sd, " -> residualised ",
      rv$delta_tumour_sd, " SD (p ", rv$p_wilcoxon, "). Direction survives adjustment: ",
      o$direction_survives_quality_adjustment, " (TCGA-KIRC), ",
      ov$direction_survives_quality_adjustment, " (CPTAC-3).")
}
msg("Reading the lncRNA rows. The normal and tumour libraries differ in library quality as well as ",
    "in tissue (TCGA-KIRC normals ", round(mean(NF_N_K), 2), "% non-feature against ",
    round(mean(NF_T_LN), 2), "% in the tumours; CPTAC-3 normals ", round(mean(NF_N_V), 2),
    "% against ", round(mean(NF_T_V), 2), "%, the opposite sign), and the observed lncRNA scores ",
    "are largely that axis. On the observed matrix the quality difference alone accounts for ",
    "essentially all of the blue gap and most of the turquoise gap in TCGA-KIRC, and both collapse ",
    "when the three STAR metrics are removed from both matrices, so NO CLAIM MAY BE MADE that ",
    "blue or turquoise runs opposite to the protein-coding module in normal kidney. Only ",
    "greenyellow survives: ", pick("TCGA-KIRC", "lncRNA greenyellow")$delta_tumour_sd,
    " SD observed and ", pick("TCGA-KIRC", "lncRNA greenyellow", "residualised")$delta_tumour_sd,
    " SD residualised (p ", pick("TCGA-KIRC", "lncRNA greenyellow", "residualised")$p_wilcoxon,
    "), and it is the one lncRNA module that is lower in normal kidney than ",
    "in tumour. Rows carrying the confound are flagged in quality_confounded with the rule in ",
    "quality_confound_rule and the arithmetic in quality_confound_note.")
msg("Rows flagged as confounded by library quality: ",
    paste(admix[quality_confounded == TRUE, paste0(cohort, " ", feature, " [", space, "]")],
          collapse = "; "), ". Rows in which the direction survives the adjustment: ",
    paste(unique(admix[direction_survives_quality_adjustment == TRUE, paste0(cohort, " ", feature)]),
          collapse = "; "), ".")

msg("Elapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
write_session_info("34_proximal_tubule_reading")
banner("34 | done")
