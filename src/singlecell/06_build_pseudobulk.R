# 06_build_pseudobulk.R: donor x compartment pseudobulk from raw counts.
# Run after 12: it stops unless every data-quality criterion in
# 12_go_no_go_criteria.tsv passes (--force overrides, with a logged warning).
# Raw UMI counts of cells retained by 05 are summed per donor x compartment,
# each compartment taken from its configured tissue. All libraries of a donor
# are summed, so each column is one donor. Donor-compartments with fewer than
# min_cells_per_donor_compartment cells go to the exclusion table.
# Outputs: data/derived/singlecell/<dataset>.pseudobulk.rds, results/singlecell/
#          06_pseudobulk_samples.tsv, 06_pseudobulk_exclusions.tsv,
#          06_count_preservation.tsv, figures/singlecell/06_pseudobulk_qc

SC_DIR <- NULL
source(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "00_config_sc.R") }))
suppressPackageStartupMessages({ library(SingleCellExperiment); library(ggplot2) })
banner("06  donor-level pseudobulk")
start_log("06_build_pseudobulk")

gate <- read_tsv("12_go_no_go_criteria.tsv")
FORCE <- "--force" %in% commandArgs(TRUE)
if (!all(gate$pass) && !FORCE)
  stop("the data-quality criteria in 12_go_no_go_criteria.tsv are not all met; stage 06 will not run (use --force to override)")
if (!all(gate$pass)) msg("WARNING: data-quality criteria not all met; running because --force was given")

PBC <- CFG$pseudobulk_build
MIN_CELLS <- CFG$pseudobulk$min_cells_per_donor_compartment
man <- read_manifest()$manifest
audit <- read_tsv("02_schema_audit.tsv")
datasets <- intersect(audit[status == "imported", dataset],
                      man[analysis_role %in% unlist(PBC$datasets_roles), dataset_name])
tme <- c("immune", "endothelial", "fibroblast_stromal")

samples <- list(); excluded <- list(); preservation <- list()
for (ds in datasets) {
  msg("pseudobulk ", ds)
  sce <- readRDS(file.path(DERIVED_DIR, paste0(ds, ".sce.rds")))
  cd <- readRDS(file.path(DERIVED_DIR, paste0(ds, ".coldata_qc.rds")))
  g <- readRDS(file.path(DERIVED_DIR, paste0(ds, ".genes.rds")))
  stopifnot(identical(cd$cell_id, colnames(sce)), nrow(g) == nrow(sce))
  cnt <- counts(sce); rm(sce)
  require_raw_counts(cnt, paste(ds, "counts"))
  assert_no_outcome_fields(cd, paste(ds, "cell metadata"))

  keep <- cd$qc_reason == ""
  groups <- list()
  for (comp in names(PBC$compartments)) {
    tis <- PBC$compartments[[comp]]
    sel <- keep & cd$tissue == tis &
      (if (comp == "pooled_nonmalignant_tme") cd$compartment %in% tme else cd$compartment == comp)
    if (any(sel)) groups[[comp]] <- which(sel)
  }
  if (!length(groups)) { msg("  no eligible cells"); next }

  mats <- list(); tabs <- list()
  for (comp in names(groups)) {
    idx <- groups[[comp]]
    donors <- cd$donor_id[idx]
    n_by <- table(donors)
    ok_d <- names(n_by)[n_by >= MIN_CELLS]
    bad_d <- names(n_by)[n_by < MIN_CELLS]
    if (length(bad_d)) excluded[[paste(ds, comp)]] <- data.table(dataset = ds, compartment = comp, donor_id = bad_d,
                                                                 n_cells = as.integer(n_by[bad_d]),
                                                                 reason = sprintf("fewer than %d cells", MIN_CELLS))
    idx <- idx[donors %in% ok_d]
    if (!length(idx)) next
    pb <- aggregate_pseudobulk(cnt[, idx, drop = FALSE], cd$donor_id[idx])
    colnames(pb) <- paste(colnames(pb), comp, sep = "|")
    mats[[comp]] <- pb
    tabs[[comp]] <- cd[idx, .(n_cells = .N, n_samples = uniqueN(sample_id),
                              median_umi = as.numeric(median(n_umi)), median_genes = as.numeric(median(n_genes)),
                              median_pct_mito = as.numeric(median(pct_mito)), n_rescued = sum(compartment_before_rescue != compartment)),
                       by = donor_id][, `:=`(dataset = ds, compartment = comp, tissue = PBC$compartments[[comp]])]
    preservation[[paste(ds, comp)]] <- data.table(dataset = ds, compartment = comp,
                                                  cell_total = sum(cnt[, idx, drop = FALSE]), pseudobulk_total = sum(pb))
  }
  pb_all <- do.call(cbind, mats)
  st <- rbindlist(tabs, use.names = TRUE)
  st[, pb_id := paste(donor_id, compartment, sep = "|")]
  st <- st[match(colnames(pb_all), pb_id)]
  st[, library_size := as.numeric(Matrix::colSums(pb_all))]
  assert_donor_level(st)
  assert_no_outcome_fields(st, paste(ds, "pseudobulk samples"))
  require_raw_counts(pb_all, paste(ds, "pseudobulk"))
  saveRDS(list(counts = pb_all, samples = st, genes = g, dataset = ds,
               config_sha256 = CFG_HASH, created = format(Sys.time(), tz = "UTC", usetz = TRUE)),
          file.path(DERIVED_DIR, paste0(ds, ".pseudobulk.rds")))
  samples[[ds]] <- st
  msg("  ", ncol(pb_all), " donor-compartment pseudobulks")
  rm(cnt, pb_all, mats); gc(verbose = FALSE)
}

S <- rbindlist(samples, use.names = TRUE)
P <- rbindlist(preservation)
P[, preserved := cell_total == pseudobulk_total]
if (!all(P$preserved)) stop("pseudobulk totals differ from cell totals: ", paste(P[preserved == FALSE, paste(dataset, compartment)], collapse = ", "))
save_tsv(S[, .(dataset, donor_id, compartment, tissue, n_cells, n_samples, n_rescued, library_size,
               median_umi, median_genes, median_pct_mito)], "06_pseudobulk_samples.tsv")
save_tsv(if (length(excluded)) rbindlist(excluded) else data.table(dataset = character()), "06_pseudobulk_exclusions.tsv")
save_tsv(P, "06_count_preservation.tsv")

p <- ggplot(S, aes(x = donor_id, y = compartment, fill = log10(library_size))) +
  geom_tile(colour = "white") + geom_text(aes(label = n_cells), size = 1.8) +
  facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(name = "log10 library size") +
  labs(x = "Donor", y = NULL, title = "Donor x compartment pseudobulks",
       subtitle = "Numbers are cells summed; each tile is one independent donor-level replicate") +
  theme_bw(base_size = 7) + theme(axis.text.x = element_text(angle = 70, hjust = 1))
save_fig(p, "06_pseudobulk_qc", width = 14, height = 4)
write_session_info("06_build_pseudobulk")
