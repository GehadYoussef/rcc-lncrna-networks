# 26_published_signatures.R: published TCGA-KIRC lncRNA signatures with and without quality adjustment
#
# Re-scores six published TCGA-KIRC lncRNA prognostic signatures on the 528
# discovery tumours: name mapping to GENCODE v36, scores on observed and
# residualised expression, correlation with the STAR metrics and lncRNA axis,
# Cox models with clinical and STAR covariates, paired bootstrap C-index
# increments, and placement in two matched-size random-set nulls.
# Every row carries weight_provenance (see PROV_CLASS_NOTE). Ranges are never
# pooled across provenance classes.
# Inputs: cache dataset.rds, networks.rds, star_qc.rds, lnc_global_axis.rds, expr_raw.rds.
# Outputs: 26_published_signatures_{genes,models,correlations}.tsv,
# 26_random_signature_null.tsv, 26_signature_scores_per_sample.tsv,
# 26_signature_provenance_key.tsv, 26_headline_ranges_by_provenance.tsv, and
# 26_NOTE_ensembl_unreachable.txt when an Ensembl lookup fails.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(httr); library(jsonlite)
})
banner("26 | Published TCGA-KIRC lncRNA signatures with and without quality adjustment")
set.seed(SEED)
t_start <- Sys.time()

N_RANDOM_SETS <- 1000L   # random lncRNA sets per signature size (Venet null)

# ---- 0. Signature definitions, as published ----
# Gene names as printed by the source, published Cox coefficients where
# available (NA otherwise) and the source URL.
# name_convention: "GENCODE19_clone" are Havana clone names of the TCGA legacy
# release, resolved via Ensembl GRCh37. "GDC_current" names are matched to GENCODE v36.

SIGNATURES <- list(

  Jiang_2019_9lnc = list(
    label    = "Jiang 2019, 9-lncRNA nomogram",
    citation = paste("Jiang W, Guo Q, Wang C, Zhu Y. A nomogram based on 9-lncRNAs signature",
                     "for improving prognostic prediction of clear cell renal cell carcinoma.",
                     "Cancer Cell Int. 2019;19:208. PMID 31404170; PMC6683339."),
    doi        = "10.1186/s12935-019-0928-5",
    source_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6683339/fullTextXML",
    coef_source = paste("published risk-score formula, Results:",
                        "'Risk score = (0.77321 x relative expression of RP13-463N16.6)",
                        "- (0.36556 x ... CTD-2201E18.5) + ... - (0.8543 x ... RP11-348J24.2)'"),
    name_convention = "GENCODE19_clone",
    genes = data.table(
      published_name = c("RP13-463N16.6", "CTD-2201E18.5", "RP11-430G17.3", "AC005785.2",
                         "RP11-2E11.9", "TFAP2A-AS1", "RP11-133F8.2", "RP11-297L17.2",
                         "RP11-348J24.2"),
      published_coef = c(0.77321, -0.36556, 0.24349, 0.37839,
                         0.74425, 0.03509, -0.01524, -0.51770,
                         -0.85430))),

  Gui_2021_8autophagy = list(
    label    = "Gui 2021, 8 autophagy-related lncRNAs",
    citation = paste("Gui CP, Cao JZ, Tan L, et al. A panel of eight autophagy-related long",
                     "non-coding RNAs is a good predictive parameter for clear cell renal cell",
                     "carcinoma. Genomics. 2021;113(2):740-754. PMID 33516849."),
    doi        = "10.1016/j.ygeno.2021.01.016",
    # Paywalled source. The panel is quoted from an open-access review.
    source_url = "https://pmc.ncbi.nlm.nih.gov/articles/PMC12403932/",
    coef_source = paste("coefficients NOT retrievable: the source is paywalled (ScienceDirect",
                        "HTTP 403, no PMCID). Panel quoted verbatim from Alimohammadi M et al.,",
                        "Cancer Cell Int 2025 (PMC12403932): 'A risk-score model comprising an",
                        "8-lncRNA signature (AC156455.1, AC107021.2, AC073611.1, AC105446.1,",
                        "SPINT1-AS1, WDFY3-AS2, FOXD2-AS1, and MELTF-AS1)'. Hazard direction",
                        "per lncRNA is not stated there either."),
    name_convention = "GDC_current",
    genes = data.table(
      published_name = c("AC156455.1", "AC107021.2", "AC073611.1", "AC105446.1",
                         "SPINT1-AS1", "WDFY3-AS2", "FOXD2-AS1", "MELTF-AS1"),
      published_coef = rep(NA_real_, 8))),

  Yu_2021_m6A = list(
    label    = "Yu 2021, m6A-related lncRNA signature",
    citation = paste("Yu J, Mao W, Sun S, et al. Identification of an m6A-related lncRNA",
                     "signature for predicting the prognosis in patients with kidney renal",
                     "clear cell carcinoma. Front Oncol. 2021;11:663263.",
                     "PMID 34123820; PMC8187870."),
    doi        = "10.3389/fonc.2021.663263",
    source_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8187870/fullTextXML",
    coef_source = paste("published risk-score formula: 'risk score = 0.935053 * AC012170.2 +",
                        "(-1.93775) * AC025580.3 + 0.416438 * AL157394.1 + 0.291862 *",
                        "AP006621.2 + (-0.35955) * AC124312.5'"),
    name_convention = "GDC_current",
    genes = data.table(
      published_name = c("AC012170.2", "AC025580.3", "AL157394.1", "AP006621.2", "AC124312.5"),
      published_coef = c(0.935053, -1.93775, 0.416438, 0.291862, -0.35955))),

  Hong_2022_cuproptosis = list(
    label    = "Hong 2022, cuproptosis-related lncRNA signature",
    citation = paste("Hong P, Huang W, Du H, et al. Prognostic value and immunological",
                     "characteristics of a novel cuproptosis-related long noncoding RNAs risk",
                     "signature in kidney renal clear cell carcinoma. Front Genet.",
                     "2022;13:1009555. PMID 36406128; PMC9669974."),
    doi        = "10.3389/fgene.2022.1009555",
    source_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9669974/fullTextXML",
    coef_source = paste("published risk-score formula: 'Risk score = (-1.15122394874834 * Exp",
                        "AL161782.1) + (0.4711103719724987 * Exp AC026401.3) +",
                        "(0.678892201655986 * Exp APCDD1L-DT) + (0.468667562066302 * Exp MINCR)'"),
    name_convention = "GDC_current",
    genes = data.table(
      published_name = c("AL161782.1", "AC026401.3", "APCDD1L-DT", "MINCR"),
      published_coef = c(-1.15122394874834, 0.4711103719724987,
                          0.678892201655986, 0.468667562066302))),

  Liu_2023_pyroptosis = list(
    label    = "Liu 2023, pyroptosis-related lncRNA model",
    citation = paste("Liu C, Dai S, Geng H, et al. Development and validation of a kidney renal",
                     "clear cell carcinoma prognostic model relying on pyroptosis-related",
                     "lncRNAs. Eur J Med Res. 2023;28(1):341. PMID 37700389; PMC10498568."),
    doi        = "10.1186/s40001-023-01277-2",
    source_url = paste0("https://static-content.springer.com/esm/art%3A10.1186%2F",
                        "s40001-023-01277-2/MediaObjects/40001_2023_1277_MOESM2_ESM.docx"),
    coef_source = paste("gene list from Additional file 2 (Table S2, 'Names of the six lncRNAs",
                        "in the model that are connected to pyroptosis'). Coefficients are NOT",
                        "published: the main text gives only the generic formula",
                        "'risk score = coef(lncRNA1) x expr(lncRNA1) + ...' and Additional",
                        "file 3 is a table of clinical characteristics. Hazard direction per",
                        "lncRNA is not stated."),
    name_convention = "GDC_current",
    genes = data.table(
      published_name = c("LINC02747", "LUCAT1", "LINC00896", "KLHDC7B-DT",
                         "LINC01138", "LINC01671"),
      published_coef = rep(NA_real_, 6))),

  Feng_2024_disulfidptosis = list(
    label    = "Feng 2024, disulfidptosis-related lncRNAs",
    citation = paste("Feng K, Zhou S, Sheng Y, et al. Disulfidptosis-related lncRNA signatures",
                     "for prognostic prediction in kidney renal clear cell carcinoma.",
                     "Clin Genitourin Cancer. 2024;22(4):102095. PMID 38833825."),
    doi        = "10.1016/j.clgc.2024.102095",
    source_url = paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?",
                        "query=DOI:%2210.1016/j.clgc.2024.102095%22&resultType=core&format=json"),
    coef_source = paste("gene list from the published abstract: 'Six signatures, namely",
                        "FAM83C.AS1, AC136475.2, AC121338.2, AC026401.3, AC254562.3, and",
                        "AC000050.2, were established'. Full text is paywalled (Elsevier, no",
                        "PMCID), so coefficients and per-lncRNA hazard direction were not",
                        "retrievable."),
    name_convention = "GDC_current",
    genes = data.table(
      published_name = c("FAM83C-AS1", "AC136475.2", "AC121338.2", "AC026401.3",
                         "AC254562.3", "AC000050.2"),
      published_coef = rep(NA_real_, 6)))
)

# The Feng source prints names with dots (FAM83C.AS1). The GDC symbol is FAM83C-AS1.
PRINTED_AS <- c("FAM83C-AS1" = "FAM83C.AS1")

# Result of the Ensembl GRCh37 lookup below, recorded on 11 September 2026.
# Used only when the API is unreachable and no cache exists.
LEGACY_FALLBACK <- list(
  "RP13-463N16.6" = "ENSG00000242147",
  "CTD-2201E18.5" = "ENSG00000271788",
  "RP11-430G17.3" = "ENSG00000271200",
  "AC005785.2"    = "ENSG00000268189",
  "RP11-2E11.9"   = c("ENSG00000270953", "ENSG00000270869"),
  "TFAP2A-AS1"    = "ENSG00000229950",
  "RP11-133F8.2"  = "ENSG00000249776",
  "RP11-297L17.2" = "ENSG00000260963",
  "RP11-348J24.2" = "ENSG00000250049")

# ---- 1. Cohort, quality metrics and observed expression of the 528 tumours ----
banner("1 | Cohort, STAR metrics and observed expression")
ds     <- readRDS(file.path(CACHE_DIR, "dataset.rds"))
nets   <- readRDS(file.path(CACHE_DIR, "networks.rds"))
qc     <- as.data.table(readRDS(file.path(CACHE_DIR, "star_qc.rds")))
axis   <- as.data.table(readRDS(file.path(CACHE_DIR, "lnc_global_axis.rds")))
cohort <- as.data.table(ds$cohort_full)
stopifnot(all(c("sample_barcode", "file_id", "os_time", "os_event", "age", "sex",
                "T_stage", "N_pos", "M1", "grade_num") %in% names(cohort)))
msg("Discovery cohort: ", nrow(cohort), " patients, ", sum(cohort$os_event), " deaths")

qc <- qc[match(cohort$sample_barcode, sample_barcode)]
stopifnot(identical(qc$sample_barcode, cohort$sample_barcode), !anyNA(qc$pct_noFeature))
COV <- tech_covariates(qc)                       # noFeature, multimapping, log depth
rownames(COV) <- cohort$sample_barcode

PROD_GENES <- colnames(nets$lnc$expr)            # the 3,442 production lncRNAs
msg("Production lncRNA set: ", length(PROD_GENES), " transcripts (networks.rds)")

# ---- legacy clone-based names -> Ensembl stable IDs ----
# The Ensembl GRCh37 gene set is GENCODE 19, the annotation of the TCGA legacy
# matrices, so the retired clone-based names still resolve there.
legacy_rds  <- file.path(CACHE_DIR, "26_legacy_name_map.rds")
legacy_need <- unique(unlist(lapply(SIGNATURES, function(s)
  if (identical(s$name_convention, "GENCODE19_clone")) s$genes$published_name else NULL)))

# Ensembl xrefs/symbol query for legacy names and for current-style names that
# fail the v36 match. Returns the query outcome:
#   "found"      HTTP 200 with at least one gene cross-reference
#   "not_found"  HTTP 200 with no gene cross-reference
#   "error"      transport failure, non-200 or unparseable body after
#                ENS_ATTEMPTS tries with backoff, never cached
# The REST server fails intermittently, so an error is not evidence of absence.
ENS_ATTEMPTS <- 5L
ens_symbol_probe <- function(symbol, host, attempts = ENS_ATTEMPTS, pause = 0.25) {
  u <- sprintf("https://%s/xrefs/symbol/homo_sapiens/%s?content-type=application/json",
               host, utils::URLencode(symbol, reserved = TRUE))
  mk <- function(status, ids, http, detail, a)
    list(status = status, ids = as.character(ids), http_status = http, detail = detail,
         symbol = symbol, host = host, url = u, queried_on = format(Sys.Date()), attempts = a)
  last <- "no attempt was made"
  for (a in seq_len(attempts)) {
    r <- tryCatch(httr::GET(u, httr::timeout(60)), error = function(e) e)
    if (inherits(r, "error")) {
      last <- paste0("transport error: ", conditionMessage(r))
    } else {
      code <- httr::status_code(r)
      if (code == 200L) {
        txt <- httr::content(r, "text", encoding = "UTF-8")
        j <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = TRUE), error = function(e) e)
        if (inherits(j, "error")) {
          last <- paste0("HTTP 200 but the body did not parse as JSON: ", conditionMessage(j))
        } else if (is.data.frame(j) && nrow(j) > 0) {
          ids <- unique(as.character(j$id[j$type == "gene"]))
          Sys.sleep(pause)
          return(if (length(ids))
            mk("found", ids, 200L, sprintf(
              "HTTP 200, %d cross-reference records, %d of type gene", nrow(j), length(ids)), a)
          else mk("not_found", character(0), 200L, sprintf(
              "HTTP 200, %d cross-reference records, none of type gene", nrow(j)), a))
        } else if (length(j) == 0L) {
          Sys.sleep(pause)
          return(mk("not_found", character(0), 200L, paste(
            "HTTP 200 with an empty result: the server replied and holds no",
            "cross-reference for this symbol"), a))
        } else {
          last <- "HTTP 200 with an unexpected JSON shape"
        }
      } else {
        last <- paste0("HTTP ", code)
      }
    }
    if (a < attempts) Sys.sleep(2^a)
  }
  mk("error", character(0), NA_integer_,
     paste0("lookup FAILED after ", attempts, " attempts; last outcome: ", last), attempts)
}
# Legacy-path wrapper: returns IDs and raises on error, so the fallback map is used.
ens_symbol_lookup <- function(symbol, host) {
  o <- ens_symbol_probe(symbol, host)
  if (identical(o$status, "error")) stop(symbol, " @ ", host, ": ", o$detail, call. = FALSE)
  o$ids
}

# The note records lookups that failed in this run, so a stale copy is removed.
ENS_NOTE <- file.path(RESULTS_DIR, "26_NOTE_ensembl_unreachable.txt")
if (file.exists(ENS_NOTE)) invisible(file.remove(ENS_NOTE))

legacy_map <- if (file.exists(legacy_rds)) readRDS(legacy_rds) else NULL
ens_error  <- NULL
# Accept the cache only if every needed name has a non-empty entry.
if (!is.null(legacy_map) &&
    !(all(legacy_need %in% names(legacy_map)) &&
      all(lengths(legacy_map[legacy_need]) > 0))) legacy_map <- NULL
legacy_from_cache <- !is.null(legacy_map)
if (is.null(legacy_map) && length(legacy_need)) {
  msg("Querying Ensembl GRCh37 REST for ", length(legacy_need), " legacy clone-based names")
  legacy_map <- tryCatch({
    out <- lapply(legacy_need, function(s) ens_symbol_lookup(s, "grch37.rest.ensembl.org"))
    names(out) <- legacy_need
    if (all(lengths(out) > 0)) saveRDS(out, legacy_rds) else
      warning("26: Ensembl GRCh37 returned no gene for ",
              paste(legacy_need[lengths(out) == 0], collapse = ", "),
              "; result NOT cached so the next run re-queries")
    out
  }, error = function(e) { ens_error <<- conditionMessage(e); NULL })
}
LEGACY_SOURCE <- if (!is.null(legacy_map)) {
  if (legacy_from_cache) "Ensembl GRCh37 REST (xrefs/symbol), from cache/26_legacy_name_map.rds"
  else "Ensembl GRCh37 REST (xrefs/symbol), queried in this run"
} else {
  legacy_map <- LEGACY_FALLBACK
  note <- ENS_NOTE
  writeLines(c("Ensembl GRCh37 REST was unreachable when 26 was run.",
               paste("Error:", ens_error),
               "The legacy clone-based names were mapped with the map recorded in the script",
               "(the result of the same xrefs/symbol call on 11 September 2026)."), note)
  warning("26: Ensembl unreachable, using the recorded fallback map (see ", note, ")")
  "recorded fallback (Ensembl GRCh37 REST, 11 September 2026)"
}
msg("Legacy name map source: ", LEGACY_SOURCE)
for (s in legacy_need) msg("  ", s, " -> ", paste(legacy_map[[s]], collapse = " / "))

msg("Reading expr_raw.rds ...")
raw <- readRDS(file.path(CACHE_DIR, "expr_raw.rds"))
ann <- as.data.table(raw$gene_ann)
ann[, ens := sub("\\..*$", "", gene_id)]
jcol <- match(cohort$file_id, raw$sample_info$file_id)
stopifnot(!anyNA(jcol), !anyDuplicated(jcol))

# Candidate rows: possible signature members plus the production set (for the null).
sig_names <- unique(unlist(lapply(SIGNATURES, function(s) s$genes$published_name)))
cand_ens  <- unique(unlist(legacy_map))
keep_rows <- sort(unique(c(which(ann$gene_name %in% sig_names),
                           which(ann$ens %in% cand_ens),
                           which(ann$gene_id %in% PROD_GENES))))
E_OBS <- log2(t(raw$fpkm[keep_rows, jcol, drop = FALSE]) + 1)
rownames(E_OBS) <- cohort$sample_barcode
colnames(E_OBS) <- ann$gene_id[keep_rows]
MED_FPKM <- setNames(apply(raw$fpkm[keep_rows, jcol, drop = FALSE], 1, median),
                     ann$gene_id[keep_rows])
rm(raw); invisible(gc(verbose = FALSE))
msg("Observed matrix: ", nrow(E_OBS), " samples x ", ncol(E_OBS), " genes, log2(FPKM+1)")

# Per-gene OLS, so residualising the whole matrix once is equivalent to
# residualising each signature's genes separately.
TFIT  <- fit_technical(E_OBS, COV)
E_RES <- apply_technical(E_OBS, TFIT, COV)
msg("Residualised matrix built on the three STAR metrics (per-gene OLS, ",
    nrow(E_OBS), " samples)")

# ---- helpers ----------------------------------------------------------------
zc <- function(x) as.numeric(scale(x))
spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(c(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}

# ---- 2. Map published names onto the GDC (GENCODE v36) annotation ----
banner("2 | Mapping published lncRNA names to GENCODE v36")

# ---- probe current-style names that fail the exact v36 match ----
# Both Ensembl hosts are queried and the outcome is recorded per gene. A
# returned stable ID present in the GDC annotation is used for mapping.
gdc_names     <- unique(unlist(lapply(SIGNATURES, function(s)
  if (identical(s$name_convention, "GDC_current")) s$genes$published_name else NULL)))
unmatched_gdc <- sort(setdiff(gdc_names, ann$gene_name))
probe_rds     <- file.path(CACHE_DIR, "26_unmatched_name_probe.rds")
PROBE_HOSTS   <- c(ensembl_current = "rest.ensembl.org",
                   ensembl_grch37  = "grch37.rest.ensembl.org")
# Cache format version 2 stores the outcome of each query per name and host.
# Files in another format are discarded. Only "found" and "not_found" entries
# are reused. "error" entries are never written.
PROBE_CACHE_VERSION <- 2L
probe_entry_ok <- function(e)
  is.list(e) && !is.null(e$status) && e$status %in% c("found", "not_found")

probe_cached <- if (file.exists(probe_rds))
  tryCatch(readRDS(probe_rds), error = function(e) NULL) else NULL
if (!is.null(probe_cached) && !identical(probe_cached$version, PROBE_CACHE_VERSION)) {
  msg("Unmatched-name probe cache is in the pre-v2 format (an empty answer with no record of ",
      "whether the query succeeded); discarding it and re-querying")
  probe_cached <- NULL
}
PROBE <- list()
n_probe_queried <- 0L
if (length(unmatched_gdc))
  msg("Unmatched current-style names to probe against Ensembl (current and GRCh37): ",
      paste(unmatched_gdc, collapse = ", "))
for (nm in unmatched_gdc) {
  PROBE[[nm]] <- list()
  for (hk in names(PROBE_HOSTS)) {
    e <- probe_cached$result[[nm]][[hk]]
    if (probe_entry_ok(e)) {
      e$from_cache <- TRUE
      PROBE[[nm]][[hk]] <- e
      msg("  ", nm, " @ ", hk, ": ", e$status, " (cached outcome of ", e$queried_on, ")")
      next
    }
    e <- ens_symbol_probe(nm, PROBE_HOSTS[[hk]])
    e$from_cache <- FALSE
    n_probe_queried <- n_probe_queried + 1L
    PROBE[[nm]][[hk]] <- e
    msg("  ", nm, " @ ", hk, ": ", e$status, " -- ", e$detail)
  }
}
# Write back ONLY definitive outcomes.
probe_definitive_entries <- lapply(PROBE, function(x) x[vapply(x, probe_entry_ok, logical(1))])
probe_definitive_entries <- probe_definitive_entries[lengths(probe_definitive_entries) > 0]
if (length(probe_definitive_entries))
  saveRDS(list(version = PROBE_CACHE_VERSION, hosts = PROBE_HOSTS,
               written_on = format(Sys.Date()), result = probe_definitive_entries), probe_rds)
probe_failed <- unlist(lapply(names(PROBE), function(nm)
  vapply(names(PROBE[[nm]]), function(hk)
    if (identical(PROBE[[nm]][[hk]]$status, "error")) paste0(nm, " @ ", hk) else NA_character_,
    character(1))), use.names = FALSE)
probe_failed <- probe_failed[!is.na(probe_failed)]
if (length(probe_failed)) {
  note <- ENS_NOTE
  cat(paste0("Ensembl REST did not answer for ", length(probe_failed),
             " name x host combination(s) when 26 was run (", format(Sys.Date()), "): ",
             paste(probe_failed, collapse = ", "), "\n",
             "Those outcomes were NOT cached and are re-queried on the next run. ",
             "For those combinations the absence of a mapping is NOT evidence that no such ",
             "gene exists.\n"),
      file = note, append = file.exists(note))
  warning("26: Ensembl did not answer for ", paste(probe_failed, collapse = ", "),
          "; the failed lookups were not cached (see ", note, ")")
}
PROBE_SOURCE <- paste0("Ensembl REST xrefs/symbol, hosts ",
                       paste(PROBE_HOSTS, collapse = " and "),
                       "; each outcome dated and marked cached or queried in this run")
msg("Unmatched-name probe: ", n_probe_queried, " of ", 2L * length(unmatched_gdc),
    " name x host lookups queried in this run, the rest reused from a recorded outcome; ",
    length(probe_failed), " failed")

probe_status <- function(nm, hk) {
  e <- PROBE[[nm]][[hk]]; if (is.null(e)) "not_probed" else e$status
}
probe_ids <- function(nm, which) {
  e <- PROBE[[nm]][[which]]; if (is.null(e)) character(0) else as.character(e$ids)
}
# TRUE only if every host returned a definitive answer for this name.
probe_definitive <- function(nm) {
  if (is.null(PROBE[[nm]])) return(FALSE)
  all(vapply(names(PROBE_HOSTS), function(hk) probe_status(nm, hk) != "error", logical(1)))
}
probe_outcome_string <- function(nm) {
  if (is.null(PROBE[[nm]])) return(NA_character_)
  paste(vapply(names(PROBE_HOSTS), function(hk)
    paste0(hk, "=", probe_status(nm, hk)), character(1)), collapse = "; ")
}
probe_text <- function(nm) {
  if (is.null(PROBE[[nm]])) return("no Ensembl probe was performed for this name")
  f <- function(hk) {
    e <- PROBE[[nm]][[hk]]
    if (is.null(e)) return("not queried")
    body <- switch(e$status,
      found     = paste0(paste(e$ids, collapse = "/"), " (HTTP 200)"),
      not_found = "no gene of this name (the server answered, HTTP 200, with no gene-level cross-reference)",
      error     = paste0("LOOKUP FAILED, so whether a gene of this name exists is UNKNOWN (",
                         e$detail, ")"))
    paste0(body, " [", if (isTRUE(e$from_cache)) paste0("cached outcome of ", e$queried_on)
                       else paste0("queried in this run, ", e$queried_on), "]")
  }
  paste0("Ensembl xrefs/symbol probe [", PROBE_SOURCE, "]: current release -> ",
         f("ensembl_current"), "; GRCh37/GENCODE 19 -> ", f("ensembl_grch37"))
}
if (length(unmatched_gdc)) for (nm in unmatched_gdc) msg("  ", nm, ": ", probe_text(nm))

# One row per published gene name.
map_one <- function(sig_key, published_name) {
  S <- SIGNATURES[[sig_key]]
  if (identical(S$name_convention, "GENCODE19_clone")) {
    ids <- legacy_map[[published_name]]
    hit <- ann[ens %in% ids]
    method <- paste0(LEGACY_SOURCE, ": ", published_name, " -> ",
                     paste(ids, collapse = "/"), "; stable ID matched in GENCODE v36")
    if (nrow(hit) > 1) {
      hit <- hit[order(-MED_FPKM[gene_id])][1]
      method <- paste0(method, " (several IDs returned; the one present in the GDC",
                       " annotation with the higher median FPKM was taken)")
    }
    if (!nrow(hit)) return(data.table(mapped_gene_id = NA_character_, mapped_gene_name = NA_character_,
                                      gene_type = NA_character_, map_method = paste0(method, " -- NOT FOUND in GENCODE v36")))
    return(data.table(mapped_gene_id = hit$gene_id, mapped_gene_name = hit$gene_name,
                      gene_type = hit$gene_type, map_method = method))
  }
  hit <- ann[gene_name == published_name]
  if (!nrow(hit)) {
    # No exact v36 name: fall back to IDs returned by the Ensembl probe.
    pids <- unique(c(probe_ids(published_name, "ensembl_current"),
                     probe_ids(published_name, "ensembl_grch37")))
    ph   <- ann[ens %in% pids]
    if (nrow(ph)) {
      if (nrow(ph) > 1) ph <- ph[order(-MED_FPKM[gene_id])][1]
      return(data.table(mapped_gene_id = ph$gene_id, mapped_gene_name = ph$gene_name,
                        gene_type = ph$gene_type,
                        map_method = paste0("no exact gene_name match in the GDC GENCODE v36 annotation; ",
                                            probe_text(published_name),
                                            "; a returned stable ID is present in v36 and was used")))
    }
    tail <- if (probe_definitive(published_name))
      paste("; both hosts answered and neither holds a gene of this name, so no stable ID could be",
            "carried into v36 and the gene is not scored")
    else
      paste("; at least one host did NOT answer, so the gene is not scored but its absence from the",
            "annotation is NOT established -- the failed lookup was not cached and is re-queried on",
            "the next run")
    return(data.table(mapped_gene_id = NA_character_, mapped_gene_name = NA_character_,
                      gene_type = NA_character_,
                      map_method = paste0("exact gene_name match against the GDC GENCODE v36 annotation",
                                          " -- NOT FOUND; ", probe_text(published_name), tail)))
  }
  if (nrow(hit) > 1) hit <- hit[order(-MED_FPKM[gene_id])][1]
  data.table(mapped_gene_id = hit$gene_id, mapped_gene_name = hit$gene_name,
             gene_type = hit$gene_type,
             map_method = "exact gene_name match against the GDC GENCODE v36 annotation")
}

GENES <- rbindlist(lapply(names(SIGNATURES), function(k) {
  S <- SIGNATURES[[k]]
  g <- copy(S$genes)
  m <- rbindlist(lapply(g$published_name, function(nm) map_one(k, nm)))
  cbind(data.table(signature = k, label = S$label, citation = S$citation, doi = S$doi,
                   source_url = S$source_url, coefficient_source = S$coef_source,
                   name_convention = S$name_convention,
                   published_name_as_printed = ifelse(g$published_name %in% names(PRINTED_AS),
                                                      PRINTED_AS[g$published_name], g$published_name),
                   # Probe outcome, keeping "no such gene" and "lookup failed"
                   # distinct. NA where no probe was needed.
                   name_probe_outcome = unname(vapply(g$published_name, probe_outcome_string,
                                                      character(1)))),
        g, m)
}))
GENES[, in_expr_raw     := !is.na(mapped_gene_id)]
GENES[, in_production_3442 := !is.na(mapped_gene_id) & mapped_gene_id %in% PROD_GENES]
GENES[, median_fpkm     := ifelse(in_expr_raw, round(MED_FPKM[mapped_gene_id], 3), NA_real_)]

cover <- GENES[, .(n_published = .N, n_mapped = sum(in_expr_raw),
                   n_in_production = sum(in_production_3442),
                   coefficients_published = sum(!is.na(published_coef))), by = .(signature, label)]
print(cover[, .(signature, n_published, n_mapped, n_in_production, coefficients_published)],
      row.names = FALSE)

# ---- 3. Per-gene univariate Cox and correlation with quality ----
banner("3 | Per-gene univariate Cox and quality correlation")
Y528 <- Surv(cohort$os_time, cohort$os_event)
nf   <- qc$pct_noFeature
mm   <- qc$pct_multimapping
ax   <- axis$lnc_axis[match(cohort$sample_barcode, axis$sample_barcode)]
msg("lncRNA axis available for ", sum(!is.na(ax)), " of ", nrow(cohort), " patients")

uni_gene <- function(gid) {
  v <- E_OBS[, gid]
  if (!is.finite(sd(v)) || sd(v) == 0) return(c(HR = NA_real_, p = NA_real_, coef = NA_real_))
  f <- coxph(Y528 ~ zc(v)); s <- summary(f)
  c(HR = unname(s$conf.int[1, 1]), p = unname(s$coefficients[1, 5]), coef = unname(coef(f)[1]))
}
gid_all <- GENES[in_expr_raw == TRUE, unique(mapped_gene_id)]
ug <- rbindlist(lapply(gid_all, function(g) {
  u <- uni_gene(g); s <- spearman(E_OBS[, g], nf); v <- E_OBS[, g]
  data.table(mapped_gene_id = g, uni_HR_per_SD = round(u[["HR"]], 3),
             uni_p = signif(u[["p"]], 3), uni_coef = round(u[["coef"]], 4),
             rho_noFeature = round(s[["rho"]], 3), p_noFeature = signif(s[["p"]], 3),
             # Rounded for display. Tests below use sd_log2fpkm_exact.
             sd_log2fpkm = round(sd(v), 3), sd_log2fpkm_exact = signif(sd(v), 8),
             n_nonzero_fpkm = sum(v > 0))
}))
GENES <- merge(GENES, ug, by = "mapped_gene_id", all.x = TRUE, sort = FALSE)

# ---- 4. Weights: as published, and re-estimated in TCGA-KIRC ----
banner("4 | Signature weights")
# Published arm: published coefficients where available. Otherwise the weight
# is +/-1 from the sign of the gene's univariate Cox coefficient in TCGA-KIRC
# (Venet construction), recorded per gene in weight_source.
GENES[, weight_published_arm := published_coef]
GENES[, weight_source := ifelse(!is.na(published_coef), "published coefficient", NA_character_)]
need_sign <- GENES[is.na(published_coef) & in_expr_raw == TRUE, which = TRUE]
GENES[need_sign, weight_published_arm := ifelse(uni_coef >= 0, 1, -1)]
GENES[need_sign, weight_source := paste("unit weight signed by the univariate Cox coefficient in",
                                        "TCGA-KIRC (neither coefficients nor hazard direction are",
                                        "published/retrievable for this signature)")]
GENES[is.na(mapped_gene_id), weight_source := "not scored (name not mappable to GENCODE v36)"]
# A gene is scored if it maps, has a weight and varies (unrounded SD above
# SCORE_MIN_SD).
SCORE_MIN_SD <- 1e-8
GENES[, used_in_score := in_expr_raw & is.finite(weight_published_arm) &
        !is.na(sd_log2fpkm_exact) & sd_log2fpkm_exact > SCORE_MIN_SD]
GENES[in_expr_raw == TRUE & used_in_score == FALSE & n_nonzero_fpkm == 0,
      weight_source := paste("not scored (annotated in GENCODE v36 but FPKM is identically",
                             "zero across all 528 discovery tumours, so the transcript",
                             "carries no signal)")]
GENES[in_expr_raw == TRUE & used_in_score == FALSE & n_nonzero_fpkm > 0 &
        sd_log2fpkm_exact == 0,
      weight_source := paste("not scored (expressed but constant: identical log2(FPKM+1) in",
                             "all 528 discovery tumours)")]
GENES[in_expr_raw == TRUE & used_in_score == FALSE & n_nonzero_fpkm > 0 &
        sd_log2fpkm_exact > 0,
      weight_source := paste0("not scored (non-constant but below the scoring tolerance: SD = ",
                              signif(sd_log2fpkm_exact, 3), " <= ", SCORE_MIN_SD, ")")]
GENES[used_in_score == FALSE, weight_published_arm := NA_real_]

# ---- provenance and shared membership ----
# Published arms without coefficients or hazard direction are signed in-sample,
# so their p-values are not valid p-values.
PROV <- GENES[, .(n_published_coefs = sum(!is.na(published_coef)),
                  n_scored = sum(used_in_score)), by = signature]
PROV[, weight_provenance_published_arm :=
       ifelse(n_published_coefs > 0, "published coefficients (external to TCGA-KIRC)",
              "unit weights signed by univariate Cox in the same 528 patients (IN-SAMPLE)")]
PROV[, published_arm_is_in_sample := n_published_coefs == 0]
# Reporting tier: fewer than 4 scored lncRNAs makes a signature a footnote.
PROV[, reporting_tier := ifelse(n_scored >= 4, "primary",
                                "footnote (fewer than 4 published lncRNAs scored)")]
# Provenance class attached to every table row. The footnote tier takes
# precedence: Yu 2021 publishes external coefficients, but only 2 of its 5
# lncRNAs are scored and the remnant runs opposite to the published model.
PROV[, weight_provenance_published_arm_class :=
       ifelse(reporting_tier != "primary", "footnote_incomplete_mapping",
       ifelse(n_published_coefs > 0, "external_coefficients", "in_sample_signed"))]
GENES <- merge(GENES, PROV[, .(signature, weight_provenance_published_arm,
                               weight_provenance = weight_provenance_published_arm_class,
                               published_arm_is_in_sample, reporting_tier)],
               by = "signature", all.x = TRUE, sort = FALSE)
GENES[, published_weight_is_external := !is.na(published_coef)]

# A gene may belong to several panels, so the six signatures are correlated
# tests, not independent replications.
shared <- GENES[used_in_score == TRUE, .(sigs = paste(sort(unique(signature)), collapse = "; "),
                                         k = uniqueN(signature)), by = mapped_gene_id][k > 1]
GENES[, shared_with_signatures := NA_character_]
if (nrow(shared)) {
  GENES[shared, on = "mapped_gene_id", shared_with_signatures := i.sigs]
  for (i in seq_len(nrow(shared)))
    msg("SHARED MEMBER: ", GENES[mapped_gene_id == shared$mapped_gene_id[i], mapped_gene_name][1],
        " (", shared$mapped_gene_id[i], ") belongs to ", shared$sigs[i],
        " -- these signatures are not independent replications")
} else msg("No scored gene belongs to more than one signature")

USABLE <- function(k) GENES[signature == k & used_in_score == TRUE]

# Provenance carried onto every model and correlation row.
REFIT_PROV <- "coefficients re-estimated by multivariable Cox in the same 528 patients (IN-SAMPLE)"
prov_detail_of <- function(k, w) if (identical(w, "refit_KIRC")) REFIT_PROV else
  PROV[signature == k, weight_provenance_published_arm][1]
in_sample_of <- function(k, w) identical(w, "refit_KIRC") ||
  isTRUE(PROV[signature == k, published_arm_is_in_sample][1])
tier_of <- function(k) PROV[signature == k, reporting_tier][1]
# Class per signature x weighting: footnote tier first, and refit arms are in-sample.
PROV_CLASSES <- c("external_coefficients", "in_sample_signed", "in_sample_refit",
                  "footnote_incomplete_mapping")
prov_class_of <- function(k, w) {
  if (!identical(tier_of(k), "primary")) return("footnote_incomplete_mapping")
  if (identical(w, "refit_KIRC")) return("in_sample_refit")
  PROV[signature == k, weight_provenance_published_arm_class][1]
}
PROV_CLASS_NOTE <- paste(
  "weight_provenance is one of external_coefficients (coefficients published outside these 528",
  "patients), in_sample_signed (+/-1 signed by the univariate Cox coefficient in these same 528",
  "patients because the source publishes neither coefficients nor a hazard direction),",
  "in_sample_refit (coefficients re-estimated here), footnote_incomplete_mapping (too few published",
  "members scored for the score to render the published model). Ranges must never span two classes.")

# Refit arm: multivariable Cox on the signature's genes, using the expression
# version the score is built from. Fitted and evaluated on the same 528
# patients, so its statistics are apparent (optimistic).
refit_coefs <- function(gids, E) {
  X <- E[, gids, drop = FALSE]
  colnames(X) <- paste0("g", seq_along(gids))
  d <- as.data.frame(X); d$os_time <- cohort$os_time; d$os_event <- cohort$os_event
  f <- coxph(as.formula(paste("Surv(os_time, os_event) ~", paste(colnames(X), collapse = " + "))),
             data = d)
  setNames(unname(coef(f)), gids)
}
score_of <- function(E, gids, w) {
  s <- as.numeric(E[, gids, drop = FALSE] %*% w[gids])
  setNames(s, rownames(E))
}

SIG_KEYS <- names(SIGNATURES)
SC <- list()   # SC[[sig]][[weights]][[expression]] = raw (unstandardised) score
for (k in SIG_KEYS) {
  gg <- USABLE(k); gids <- gg$mapped_gene_id
  w_pub <- setNames(gg$weight_published_arm, gids)
  w_obs <- refit_coefs(gids, E_OBS)
  w_res <- refit_coefs(gids, E_RES)
  SC[[k]] <- list(
    published  = list(observed = score_of(E_OBS, gids, w_pub),
                      residualised = score_of(E_RES, gids, w_pub)),
    refit_KIRC = list(observed = score_of(E_OBS, gids, w_obs),
                      residualised = score_of(E_RES, gids, w_res)))
  GENES[signature == k & mapped_gene_id %in% gids,
        refit_coef_observed := round(w_obs[mapped_gene_id], 4)]
  GENES[signature == k & mapped_gene_id %in% gids,
        refit_coef_residualised := round(w_res[mapped_gene_id], 4)]
  msg(SIGNATURES[[k]]$label, ": scored on ", length(gids), " of ",
      nrow(SIGNATURES[[k]]$genes), " published lncRNAs (",
      sum(gg$in_production_3442), " in the production set)")
}

GENES[, weight_provenance_note := PROV_CLASS_NOTE]
setcolorder(GENES, c("signature", "label", "weight_provenance", "reporting_tier",
                     "published_name", "published_name_as_printed",
                     "name_convention", "mapped_gene_id", "mapped_gene_name", "gene_type",
                     "map_method", "name_probe_outcome",
                     "in_expr_raw", "in_production_3442", "used_in_score",
                     "shared_with_signatures", "published_coef", "published_weight_is_external",
                     "weight_published_arm", "weight_source",
                     "weight_provenance_published_arm", "published_arm_is_in_sample",
                     "refit_coef_observed",
                     "refit_coef_residualised", "median_fpkm", "sd_log2fpkm",
                     "sd_log2fpkm_exact", "n_nonzero_fpkm",
                     "uni_HR_per_SD", "uni_p", "uni_coef", "rho_noFeature", "p_noFeature",
                     "doi", "source_url", "coefficient_source", "citation"))
save_tsv(GENES, "26_published_signatures_genes.tsv")
print(GENES[, .(signature, published_name, mapped_gene_name, in_expr_raw, in_production_3442,
                published_coef, weight_published_arm, uni_HR_per_SD, rho_noFeature)],
      nrows = 60)

# ---- 5. Correlations with the quality metrics and the lncRNA axis ----
banner("5 | Score correlations with library quality and the lncRNA axis")
VARS <- list(pct_noFeature = nf, lnc_axis = ax, pct_multimapping = mm)
# Residualisation removes the STAR metrics, so on residualised expression the
# Pearson correlation of a score with them is zero by construction. Those rows
# are flagged. The lncRNA axis is not a residualisation covariate.
ZERO_BY_CONSTRUCTION <- c("pct_noFeature", "pct_multimapping")
pearson <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(NA_real_)
  unname(cor(x[ok], y[ok]))
}
COR_GRID <- CJ(signature = SIG_KEYS, weights = c("published", "refit_KIRC"),
               expression = c("observed", "residualised"), variable = names(VARS),
               sorted = FALSE, unique = TRUE)
COR <- rbindlist(lapply(seq_len(nrow(COR_GRID)), function(i) {
  g <- COR_GRID[i]
  sc <- SC[[g$signature]][[g$weights]][[g$expression]]
  s  <- spearman(sc, VARS[[g$variable]])
  data.table(signature = g$signature, label = SIGNATURES[[g$signature]]$label,
             weight_provenance = prov_class_of(g$signature, g$weights),
             reporting_tier = tier_of(g$signature),
             weights = g$weights,
             weight_provenance_detail = prov_detail_of(g$signature, g$weights),
             in_sample_weights = in_sample_of(g$signature, g$weights),
             expression = g$expression, variable = g$variable,
             n = s[["n"]], spearman_rho = round(s[["rho"]], 3), p = signif(s[["p"]], 3),
             pearson_r = signif(pearson(sc, VARS[[g$variable]]), 3),
             pearson_zero_by_construction =
               g$expression == "residualised" && g$variable %in% ZERO_BY_CONSTRUCTION)
}))
COR[, fdr := signif(p.adjust(p, "BH"), 3), by = .(variable, expression)]
COR[, n_tests_in_family := .N, by = .(variable, expression)]
COR[, fdr_family_note := "BH over the 12 signature x weighting combinations within this variable x expression; the six signatures share cohort, expression matrix and (for at least one pair) a gene, so these are correlated tests, not independent replications"]
# The pooled family mixes provenances and in-sample p-values are not valid, so
# a within-provenance BH adjustment is given alongside.
COR[, fdr_within_provenance := signif(p.adjust(p, "BH"), 3),
    by = .(variable, expression, weight_provenance)]
COR[, n_tests_in_provenance_family := .N, by = .(variable, expression, weight_provenance)]
COR[, weight_provenance_note := PROV_CLASS_NOTE]
save_tsv(COR, "26_published_signatures_correlations.tsv")
msg("Residualised x STAR-metric rows carry pearson_zero_by_construction = TRUE; ",
    "max |pearson_r| over all such rows = ",
    signif(max(abs(COR[pearson_zero_by_construction == TRUE, pearson_r])), 3),
    " (zero to numerical precision, as the residualisation guarantees). This maximum deliberately ",
    "pools weight provenances because it is a floating-point check on an identity, not an effect ",
    "estimate; every range that IS an effect estimate is split by provenance in section 6b.")
print(COR[expression == "observed", .(signature, weights, variable, n, spearman_rho, p, fdr)],
      nrows = 60)

# ---- 6. Cox models: four specifications x two weightings x six signatures ----
banner("6 | Cox models with and without library-quality adjustment")
X_CLIN <- clinical_design(cohort, set = "clinical")
stopifnot(identical(rownames(X_CLIN), cohort$sample_barcode))
D <- data.table(sample_barcode = cohort$sample_barcode, os_time = cohort$os_time,
                os_event = cohort$os_event)
D <- cbind(D, as.data.table(X_CLIN))
D[, `:=`(noFeature = zc(qc$pct_noFeature), multimap = zc(qc$pct_multimapping),
         libsize_z = zc(log10(qc$assigned_reads)))]
STAR_TERMS <- c("noFeature", "multimap", "libsize_z")
msg("Clinical design complete cases: ", sum(complete.cases(D[, ..CLIN_TERMS])), " of ", nrow(D))

SPECS <- list(
  list(id = "1_score_alone",        expression = "observed",     terms = character(0)),
  list(id = "2_score_plus_clinical", expression = "observed",    terms = CLIN_TERMS),
  list(id = "3_score_plus_clinical_plus_STAR", expression = "observed",
       terms = c(CLIN_TERMS, STAR_TERMS)),
  list(id = "4_residualised_score_plus_clinical", expression = "residualised",
       terms = CLIN_TERMS))
SPEC_LABEL <- c(
  "1_score_alone"                      = "score alone",
  "2_score_plus_clinical"              = "+ age, sex, T, N, M1, ordinal grade",
  "3_score_plus_clinical_plus_STAR"    = "+ clinical + non-feature, multimapping, log10 depth",
  "4_residualised_score_plus_clinical" = "score on technically residualised expression + clinical")

fit_one <- function(k, w, spec) {
  raw_score <- SC[[k]][[w]][[spec$expression]]
  d <- copy(D); d[, score_raw := raw_score[sample_barcode]]
  vars <- c("os_time", "os_event", "score_raw", spec$terms)
  dc <- d[complete.cases(d[, ..vars])]
  sd_raw <- sd(dc$score_raw)
  dc[, score := zc(score_raw)]                      # per SD of the analysed set
  y  <- Surv(dc$os_time, dc$os_event)
  f  <- coxph(as.formula(paste("Surv(os_time, os_event) ~",
                               paste(c("score", spec$terms), collapse = " + "))), data = dc)
  s  <- summary(f)
  ph <- tryCatch(cox.zph(f)$table["score", "p"], error = function(e) NA_real_)

  ref <- NULL
  if (length(spec$terms)) {
    f0 <- coxph(as.formula(paste("Surv(os_time, os_event) ~", paste(spec$terms, collapse = " + "))),
                data = dc)
    lp0 <- predict(f0, type = "lp"); lp1 <- predict(f, type = "lp")
    c0 <- cindex_ci(y, lp0); c1 <- cindex_ci(y, lp1)
    bt <- paired_boot_delta_c(y, lp1, lp0, B = BOOT_B, seed = SEED)
    lrt <- 2 * (f$loglik[2] - f0$loglik[2])
    ref <- list(c0 = c0, c1 = c1, bt = bt,
                lrt_p = pchisq(lrt, df = 1, lower.tail = FALSE))
  }
  data.table(
    signature = k, label = SIGNATURES[[k]]$label, weight_provenance = prov_class_of(k, w),
    reporting_tier = tier_of(k), weights = w,
    weight_provenance_detail = prov_detail_of(k, w), in_sample_weights = in_sample_of(k, w),
    expression = spec$expression, specification = spec$id,
    specification_label = SPEC_LABEL[[spec$id]],
    terms = paste(c("score", spec$terms), collapse = " + "),
    n_genes_scored = nrow(USABLE(k)),
    n_genes_published = nrow(SIGNATURES[[k]]$genes),
    n = s$n, events = s$nevent, score_sd_raw = round(sd_raw, 4),
    HR = round(unname(s$conf.int["score", 1]), 3),
    lo = round(unname(s$conf.int["score", 3]), 3),
    hi = round(unname(s$conf.int["score", 4]), 3),
    p  = signif(unname(s$coefficients["score", 5]), 3),
    C_model = round(unname(s$concordance["C"]), 3),
    ph_p = signif(ph, 3),
    lrt_p_vs_reference = if (is.null(ref)) NA_real_ else signif(ref$lrt_p, 3),
    C_reference        = if (is.null(ref)) NA_real_ else round(ref$c0[["C"]], 3),
    C_reference_lo     = if (is.null(ref)) NA_real_ else round(ref$c0[["lo"]], 3),
    C_reference_hi     = if (is.null(ref)) NA_real_ else round(ref$c0[["hi"]], 3),
    C_with_score       = if (is.null(ref)) NA_real_ else round(ref$c1[["C"]], 3),
    C_with_score_lo    = if (is.null(ref)) NA_real_ else round(ref$c1[["lo"]], 3),
    C_with_score_hi    = if (is.null(ref)) NA_real_ else round(ref$c1[["hi"]], 3),
    delta_C            = if (is.null(ref)) NA_real_ else round(ref$bt[["delta"]], 4),
    delta_C_lo         = if (is.null(ref)) NA_real_ else round(ref$bt[["lo"]], 4),
    delta_C_hi         = if (is.null(ref)) NA_real_ else round(ref$bt[["hi"]], 4),
    delta_C_p_boot     = if (is.null(ref)) NA_real_ else signif(ref$bt[["p_boot"]], 3),
    delta_C_prop_positive = if (is.null(ref)) NA_real_ else round(ref$bt[["prop_positive"]], 3),
    boot_B             = if (is.null(ref)) NA_integer_ else as.integer(ref$bt[["n_boot"]]))
}

MODELS <- rbindlist(lapply(SIG_KEYS, function(k)
  rbindlist(lapply(c("published", "refit_KIRC"), function(w)
    rbindlist(lapply(SPECS, function(sp) fit_one(k, w, sp)))))))
MODELS[, fdr := signif(p.adjust(p, "BH"), 3), by = .(specification, weights)]
MODELS[, n_tests_in_fdr_family := .N, by = .(specification, weights)]
# BH family: the six (correlated) signatures within a specification x weighting.
MODELS[, fdr_family_note := "BH over the six signatures within this specification x weighting; the six signatures share cohort, expression matrix and (for at least one pair) a gene, so they are correlated tests"]
# Within-provenance adjustment, as in the correlation table.
MODELS[, fdr_within_provenance := signif(p.adjust(p, "BH"), 3),
       by = .(specification, weights, weight_provenance)]
MODELS[, n_tests_in_provenance_family := .N,
       by = .(specification, weights, weight_provenance)]
MODELS[, weight_provenance_note := PROV_CLASS_NOTE]
# Attenuation: percentage change in the log HR relative to spec 2 (clinical
# only) when quality is a covariate (spec 3) or removed from the expression
# (spec 4).
MODELS[, HR_clinical_only_reference := HR[specification == "2_score_plus_clinical"],
       by = .(signature, weights)]
MODELS[, pct_change_logHR_vs_clinical := round(
  100 * (log(HR) - log(HR_clinical_only_reference)) / log(HR_clinical_only_reference), 1)]
MODELS[, HR_clinical_only_reference := round(HR_clinical_only_reference, 3)]
save_tsv(MODELS, "26_published_signatures_models.tsv")
print(MODELS[, .(signature, weights, specification, n, events, HR, lo, hi, p,
                 C_reference, C_with_score, delta_C, delta_C_lo, delta_C_hi,
                 pct_change_logHR_vs_clinical)], nrows = 100)

banner("6b | Headline ranges, one weight provenance class at a time")
# Every headline range is built within one weight_provenance class. A class
# with one arm gets a single value. The table is written after section 7 adds
# the null percentiles.
RANGE_ROWS <- list()
add_range <- function(quantity, stratum, d, col, digits = 3) {
  if (!nrow(d)) return(invisible(NULL))
  for (cl in intersect(PROV_CLASSES, unique(d$weight_provenance))) {
    x <- d[weight_provenance == cl]
    v <- suppressWarnings(as.numeric(x[[col]])); ok <- is.finite(v)
    if (!any(ok)) next
    v <- v[ok]
    fmt <- paste0("%.", digits, "f")
    RANGE_ROWS[[length(RANGE_ROWS) + 1L]] <<- data.table(
      quantity = quantity, stratum = stratum, weight_provenance = cl,
      n_arms = sum(ok),
      signatures = paste(sort(unique(sub("_.*", "", x$signature[ok]))), collapse = ", "),
      minimum = min(v), maximum = max(v),
      range_text = if (length(v) == 1L) sprintf(fmt, v)
                   else paste0(sprintf(fmt, min(v)), " to ", sprintf(fmt, max(v))),
      note = PROV_CLASS_NOTE)
  }
  invisible(NULL)
}
rng_txt <- function(d, col, digits = 3) {
  v <- suppressWarnings(as.numeric(d[[col]])); v <- v[is.finite(v)]
  if (!length(v)) return("-")
  fmt <- paste0("%.", digits, "f")
  if (length(v) == 1L) sprintf(fmt, v) else
    paste0(sprintf(fmt, min(v)), " to ", sprintf(fmt, max(v)))
}
SPEC_IDS <- vapply(SPECS, function(s) s$id, character(1))
for (sp in SPEC_IDS) for (w in c("published", "refit_KIRC")) {
  d <- MODELS[specification == sp & weights == w]
  st <- paste0(sp, " | ", w, " weights")
  add_range("HR per SD",                         st, d, "HR", 3)
  add_range("C-index of the model with the score", st, d, "C_with_score", 3)
  # pct_change_logHR_vs_clinical is an attenuation only for specs 3 and 4
  # (spec 2 is the reference), so only those are reported as ranges.
  attn <- sp %in% c("3_score_plus_clinical_plus_STAR", "4_residualised_score_plus_clinical")
  if (attn)
    add_range("Attenuation of the clinically adjusted log HR (%), against specification 2",
              st, d, "pct_change_logHR_vs_clinical", 1)
  add_range("Paired bootstrap increment in C over the reference model", st, d, "delta_C", 4)
  for (cl in intersect(PROV_CLASSES, unique(d$weight_provenance))) {
    x <- d[weight_provenance == cl]
    msg(sp, " | ", w, " | ", cl, " (", paste(sort(sub("_.*", "", x$signature)), collapse = ", "),
        "): HR per SD ", rng_txt(x, "HR"),
        if (attn) paste0("; attenuation of the clinically adjusted log HR ",
                         rng_txt(x, "pct_change_logHR_vs_clinical", 1), " %") else "",
        "; increment in C ", rng_txt(x, "delta_C", 4))
  }
}
# Correlation with the non-feature fraction (observed expression), per class.
for (w in c("published", "refit_KIRC")) {
  d <- COR[expression == "observed" & variable == "pct_noFeature" & weights == w]
  add_range("Spearman rho with the non-feature fraction (observed expression)",
            paste0("observed expression | ", w, " weights"), d, "spearman_rho", 3)
  for (cl in intersect(PROV_CLASSES, unique(d$weight_provenance)))
    msg("rho with the non-feature fraction | observed | ", w, " | ", cl, ": ",
        rng_txt(d[weight_provenance == cl], "spearman_rho"))
}
msg("Footnote tier (fewer than 4 published lncRNAs scored): ",
    paste(PROV[reporting_tier != "primary", signature], collapse = ", "),
    " -- its arms carry weight_provenance = footnote_incomplete_mapping and appear in no",
    " external_coefficients or in_sample_signed range")
# Does refitting raise the HR? Counted within provenance class.
chk <- merge(MODELS[weights == "published",
                    .(signature, specification, weight_provenance, HR_pub = HR, C_pub = C_model)],
             MODELS[weights == "refit_KIRC", .(signature, specification, HR_ref = HR, C_ref = C_model)],
             by = c("signature", "specification"))
chk[, `:=`(HR_up = HR_ref > HR_pub, C_up = C_ref >= C_pub)]
for (cl in intersect(PROV_CLASSES, unique(chk$weight_provenance))) {
  x <- chk[weight_provenance == cl]
  ex <- if (any(!x$HR_up))
    paste(x[HR_up == FALSE, paste0(sub("_.*", "", signature), "/", specification, " ",
                                   sprintf("%.3f", HR_pub), "->", sprintf("%.3f", HR_ref))],
          collapse = ", ") else "none"
  msg("Refit vs published | ", cl, ": HR higher in ", sum(x$HR_up), " of ", nrow(x),
      " signature x specification pairs, model C-index at least as high in ", sum(x$C_up),
      " of ", nrow(x), "; HR exceptions: ", ex)
}
# Proportional hazards for the score term, pooled over all models as a diagnostic.
phb <- MODELS[order(ph_p)]
msg("cox.zph on the score term: ", sum(MODELS$ph_p < 0.05, na.rm = TRUE), " of ", nrow(MODELS),
    " models have p < 0.05; smallest: ",
    paste(phb[1:3, paste0(sub("_.*", "", signature), "/", weights, "/", specification, " p=", ph_p)],
          collapse = ", "))

# ---- 7. Matched-size random lncRNA sets: two nulls ----
banner("7 | Random lncRNA signatures of matched size (outcome-signed and sign-randomised nulls)")
# Sign of every production lncRNA's univariate Cox coefficient in TCGA-KIRC.
msg("Univariate Cox for ", length(PROD_GENES), " production lncRNAs ...")
Z_PROD <- scale(E_OBS[, PROD_GENES, drop = FALSE])
sgn <- vapply(PROD_GENES, function(g) {
  v <- Z_PROD[, g]
  if (!is.finite(sd(v)) || sd(v) == 0) return(NA_real_)
  sign(unname(coef(coxph(Y528 ~ v))[1]))
}, numeric(1))
ok_prod <- names(sgn)[is.finite(sgn) & sgn != 0]
N_RISK <- sum(sgn[ok_prod] > 0); N_PROT <- sum(sgn[ok_prod] < 0)
msg("Signed production lncRNAs: ", length(ok_prod), " of ", length(PROD_GENES), " (",
    N_RISK, " risk-directed, ", N_PROT, " protective; ",
    round(100 * N_RISK / length(ok_prod), 1), "% risk-directed)")
msg("That asymmetry is the mechanism behind the outcome-signed null: with ",
    round(100 * N_RISK / length(ok_prod)), "% of production lncRNAs risk-directed, an ",
    "outcome-signed random sum is close to an unsigned sum of an axis-loaded matrix, ",
    "so the null is itself oriented along the non-feature axis.")
# Sign-flipped observed expression, so a random set's score is a row sum.
E_SGN <- sweep(E_OBS[, ok_prod, drop = FALSE], 2, sgn[ok_prod], "*")

SET_SIZES <- sort(unique(vapply(SIG_KEYS, function(k) nrow(USABLE(k)), integer(1))))
msg("Signature sizes to match: ", paste(SET_SIZES, collapse = ", "))

# Two sign schemes applied to the same draws:
#   outcome_signed   +/-1 from the univariate Cox sign (Venet et al. 2011):
#                    reference for the hazard ratio.
#   sign_randomised  +/-1 at random: reference for the correlation with
#                    library quality, since outcome signing orients every draw
#                    along the cohort's dominant prognostic axis.
SIGN_SCHEMES <- c("outcome_signed", "sign_randomised")

# Gene sets are drawn first (all sizes), then the random signs, so draw b of a
# given size is the same gene set in both nulls.
set.seed(SEED)
DRAW_GENES <- lapply(SET_SIZES, function(sz)
  lapply(seq_len(N_RANDOM_SETS), function(b) sample(ok_prod, sz)))
names(DRAW_GENES) <- as.character(SET_SIZES)
DRAW_SIGNS <- lapply(SET_SIZES, function(sz)
  lapply(seq_len(N_RANDOM_SETS), function(b) sample(c(1, -1), sz, replace = TRUE)))
names(DRAW_SIGNS) <- as.character(SET_SIZES)
stopifnot(all(vapply(as.character(SET_SIZES),
                     function(s) length(DRAW_GENES[[s]]) == N_RANDOM_SETS, logical(1))),
          all(vapply(as.character(SET_SIZES), function(s)
            all(lengths(DRAW_GENES[[s]]) == as.integer(s)), logical(1))))
# gene_set_key: sorted pool indices of the drawn genes, shared by both schemes.
POOL_INDEX  <- setNames(seq_along(ok_prod), ok_prod)
gene_set_key <- function(g) paste(sort(unname(POOL_INDEX[g])), collapse = "|")
DRAW_KEYS <- lapply(as.character(SET_SIZES), function(s) vapply(DRAW_GENES[[s]], gene_set_key, character(1)))
names(DRAW_KEYS) <- as.character(SET_SIZES)
msg("Random draws: ", N_RANDOM_SETS, " gene sets per size (", paste(SET_SIZES, collapse = ", "),
    "), each scored under BOTH weightings; distinct gene sets per size: ",
    paste(vapply(as.character(SET_SIZES), function(s) length(unique(DRAW_KEYS[[s]])), integer(1)),
          collapse = ", "))

draw_score <- function(g, scheme, sgn_rand) {
  if (identical(scheme, "outcome_signed")) rowSums(E_SGN[, g, drop = FALSE])
  else as.numeric(E_OBS[, g, drop = FALSE] %*% sgn_rand)
}

null_stat <- function(s) {
  z  <- zc(s)
  f  <- coxph(Y528 ~ z); su <- summary(f)
  sp <- spearman(s, nf); sa <- spearman(s, ax)
  c(rho_noFeature = sp[["rho"]], p_noFeature = sp[["p"]],
    rho_lnc_axis = sa[["rho"]], n_lnc_axis = sa[["n"]],
    rho_multimapping = spearman(s, mm)[["rho"]],
    HR = unname(su$conf.int[1, 1]), HR_lo = unname(su$conf.int[1, 3]),
    HR_hi = unname(su$conf.int[1, 4]), HR_p = unname(su$coefficients[1, 5]),
    C = unname(su$concordance["C"]))
}
NULL_SETS <- rbindlist(lapply(SIGN_SCHEMES, function(scheme)
  rbindlist(lapply(SET_SIZES, function(sz) {
    msg("  ", scheme, ": scoring the ", N_RANDOM_SETS, " shared random sets of size ", sz)
    s_sz <- as.character(sz)
    rbindlist(lapply(seq_len(N_RANDOM_SETS), function(b) {
      g  <- DRAW_GENES[[s_sz]][[b]]
      st <- null_stat(draw_score(g, scheme, DRAW_SIGNS[[s_sz]][[b]]))
      data.table(row_type = "random_set", sign_scheme = scheme,
                 signature = NA_character_, weights = NA_character_,
                 set_size = sz, draw = b, gene_set_key = DRAW_KEYS[[s_sz]][b],
                 n = nrow(E_OBS), events = sum(cohort$os_event),
                 rho_noFeature = round(st[["rho_noFeature"]], 4),
                 p_noFeature = signif(st[["p_noFeature"]], 3),
                 rho_lnc_axis = round(st[["rho_lnc_axis"]], 4),
                 n_lnc_axis = as.integer(st[["n_lnc_axis"]]),
                 rho_multimapping = round(st[["rho_multimapping"]], 4),
                 HR = round(st[["HR"]], 4), HR_lo = round(st[["HR_lo"]], 4),
                 HR_hi = round(st[["HR_hi"]], 4), HR_p = signif(st[["HR_p"]], 3),
                 C = round(st[["C"]], 4))
    }))
  }))))
# Check that both schemes share their draws.
shared_keys <- dcast(NULL_SETS[row_type == "random_set"], set_size + draw ~ sign_scheme,
                     value.var = "gene_set_key")
stopifnot(nrow(shared_keys) == length(SET_SIZES) * N_RANDOM_SETS,
          all(shared_keys$outcome_signed == shared_keys$sign_randomised))
msg("Verified: all ", nrow(shared_keys), " (size, draw) pairs carry the same gene set under both ",
    "weightings, so the two nulls are one set of draws re-weighted, not two samples")

null_summary <- NULL_SETS[, .(row_type = "null_summary", signature = NA_character_,
                              weights = NA_character_, draw = NA_integer_,
                              n = n[1], events = events[1], n_lnc_axis = n_lnc_axis[1],
                              null_n_sets = .N,
                              null_rho_noFeature_median = round(median(rho_noFeature), 4),
                              null_rho_noFeature_p2.5   = round(unname(quantile(rho_noFeature, 0.025)), 4),
                              null_rho_noFeature_p97.5  = round(unname(quantile(rho_noFeature, 0.975)), 4),
                              null_abs_rho_noFeature_median = round(median(abs(rho_noFeature)), 4),
                              null_frac_rho_noFeature_p_lt_0.05 = round(mean(p_noFeature < 0.05), 3),
                              null_HR_median = round(median(HR), 4),
                              null_HR_p2.5   = round(unname(quantile(HR, 0.025)), 4),
                              null_HR_p97.5  = round(unname(quantile(HR, 0.975)), 4),
                              null_frac_HR_p_lt_0.05 = round(mean(HR_p < 0.05), 3),
                              null_C_median = round(median(C), 4)),
                          by = .(sign_scheme, set_size)]

# Paired comparison per draw: outcome-signed minus sign-randomised |rho| with
# the non-feature fraction. Pairing removes gene-set variance.
pw <- dcast(NULL_SETS[row_type == "random_set"], set_size + draw ~ sign_scheme,
            value.var = c("rho_noFeature", "C", "HR"))
pw[, d_abs_rho := abs(rho_noFeature_outcome_signed) - abs(rho_noFeature_sign_randomised)]
null_paired <- pw[, {
  w <- suppressWarnings(wilcox.test(abs(rho_noFeature_outcome_signed),
                                    abs(rho_noFeature_sign_randomised), paired = TRUE))
  .(row_type = "null_paired_comparison",
    sign_scheme = "outcome_signed minus sign_randomised (paired on the shared draw)",
    signature = NA_character_, weights = NA_character_, draw = NA_integer_,
    null_n_sets = .N,
    paired_median_diff_abs_rho_noFeature = round(median(d_abs_rho), 4),
    paired_diff_abs_rho_p2.5  = round(unname(quantile(d_abs_rho, 0.025)), 4),
    paired_diff_abs_rho_p97.5 = round(unname(quantile(d_abs_rho, 0.975)), 4),
    paired_frac_outcome_signed_higher = round(mean(d_abs_rho > 0), 4),
    paired_wilcoxon_p = signif(w$p.value, 3))
}, by = set_size]
msg("Paired difference in |rho| with the non-feature fraction (outcome-signed minus ",
    "sign-randomised, same draws), by set size: ",
    paste(sprintf("%d: %+.3f (%.1f%% of draws higher, signed-rank p %.3g)",
                  null_paired$set_size, null_paired$paired_median_diff_abs_rho_noFeature,
                  100 * null_paired$paired_frac_outcome_signed_higher,
                  null_paired$paired_wilcoxon_p), collapse = "; "))

# Percentile of each signature in both nulls of its own size. The HR
# percentile uses the outcome-signed null. For |rho| with quality, the
# sign-randomised null matches arms whose weights are not outcome-signed
# (Jiang, Yu, Hong published arms). The outcome-signed null matches the rest.
pct_in <- function(v, x) mean(v <= x)
SIG_PCT <- rbindlist(lapply(SIG_KEYS, function(k) rbindlist(lapply(c("published", "refit_KIRC"), function(w) {
  sz  <- nrow(USABLE(k))
  nls <- NULL_SETS[sign_scheme == "outcome_signed"  & set_size == sz]
  nlr <- NULL_SETS[sign_scheme == "sign_randomised" & set_size == sz]
  st  <- null_stat(SC[[k]][[w]][["observed"]])
  data.table(row_type = "signature_percentile", sign_scheme = NA_character_,
             signature = k, weights = w,
             weight_provenance = prov_class_of(k, w),
             weight_provenance_detail = prov_detail_of(k, w),
             in_sample_weights = in_sample_of(k, w),
             weights_are_outcome_signed = in_sample_of(k, w),
             reporting_tier = tier_of(k), set_size = sz,
             draw = NA_integer_, n = nrow(E_OBS), events = sum(cohort$os_event),
             rho_noFeature = round(st[["rho_noFeature"]], 4),
             p_noFeature = signif(st[["p_noFeature"]], 3),
             rho_lnc_axis = round(st[["rho_lnc_axis"]], 4),
             n_lnc_axis = as.integer(st[["n_lnc_axis"]]),
             rho_multimapping = round(st[["rho_multimapping"]], 4),
             HR = round(st[["HR"]], 4), HR_lo = round(st[["HR_lo"]], 4),
             HR_hi = round(st[["HR_hi"]], 4), HR_p = signif(st[["HR_p"]], 3),
             C = round(st[["C"]], 4),
             null_n_sets_outcome_signed = nrow(nls), null_n_sets_sign_randomised = nrow(nlr),
             percentile_rho_noFeature_vs_outcome_signed_null =
               round(100 * pct_in(nls$rho_noFeature, st[["rho_noFeature"]]), 1),
             percentile_abs_rho_noFeature_vs_outcome_signed_null =
               round(100 * pct_in(abs(nls$rho_noFeature), abs(st[["rho_noFeature"]])), 1),
             percentile_rho_noFeature_vs_sign_randomised_null =
               round(100 * pct_in(nlr$rho_noFeature, st[["rho_noFeature"]]), 1),
             percentile_abs_rho_noFeature_vs_sign_randomised_null =
               round(100 * pct_in(abs(nlr$rho_noFeature), abs(st[["rho_noFeature"]])), 1),
             percentile_HR_vs_outcome_signed_null = round(100 * pct_in(nls$HR, st[["HR"]]), 1),
             percentile_C_vs_outcome_signed_null  = round(100 * pct_in(nls$C, st[["C"]]), 1))
}))))

NULL_OUT <- rbindlist(list(NULL_SETS, null_summary, null_paired, SIG_PCT),
                      use.names = TRUE, fill = TRUE)
setcolorder(NULL_OUT, c("row_type", "sign_scheme", "signature", "weights", "weight_provenance",
                        "set_size", "draw", "gene_set_key", "n", "events"))
save_tsv(NULL_OUT, "26_random_signature_null.tsv")
print(null_summary[, .(sign_scheme, set_size, null_n_sets, null_rho_noFeature_median,
                       null_abs_rho_noFeature_median, null_HR_median, null_C_median)],
      row.names = FALSE)
print(SIG_PCT[, .(signature, weights, weights_are_outcome_signed, reporting_tier, set_size,
                  rho_noFeature,
                  pct_absrho_signed_null = percentile_abs_rho_noFeature_vs_outcome_signed_null,
                  pct_absrho_random_null = percentile_abs_rho_noFeature_vs_sign_randomised_null,
                  HR, HR_p, percentile_HR_vs_outcome_signed_null,
                  C, percentile_C_vs_outcome_signed_null)], row.names = FALSE)
# Median |rho| in the two nulls, over the same gene sets per size.
cmp <- null_summary[, .(sign_scheme, set_size, m = null_abs_rho_noFeature_median)]
cmp <- dcast(cmp, set_size ~ sign_scheme, value.var = "m")
print(cmp, row.names = FALSE)
msg("Median |rho| with the non-feature fraction, outcome-signed null vs sign-randomised null ",
    "(same draws, both weightings), by set size: ",
    paste(sprintf("%d: %.3f vs %.3f", cmp$set_size, cmp$outcome_signed,
                  cmp$sign_randomised), collapse = "; "))
print(null_paired, row.names = FALSE)

# ---- the percentiles, one weight provenance class at a time -----------------
for (w in c("published", "refit_KIRC")) {
  d <- SIG_PCT[weights == w]
  add_range("Percentile of |rho| with the non-feature fraction, outcome-signed null",
            paste0(w, " weights"), d, "percentile_abs_rho_noFeature_vs_outcome_signed_null", 1)
  add_range("Percentile of |rho| with the non-feature fraction, sign-randomised null",
            paste0(w, " weights"), d, "percentile_abs_rho_noFeature_vs_sign_randomised_null", 1)
  add_range("Percentile of the univariate HR, outcome-signed null",
            paste0(w, " weights"), d, "percentile_HR_vs_outcome_signed_null", 1)
  for (cl in intersect(PROV_CLASSES, unique(d$weight_provenance))) {
    x <- d[weight_provenance == cl]
    msg("percentiles | ", w, " | ", cl, " (", paste(sort(sub("_.*", "", x$signature)), collapse = ", "),
        "): |rho| vs outcome-signed null ",
        rng_txt(x, "percentile_abs_rho_noFeature_vs_outcome_signed_null", 1),
        "; |rho| vs sign-randomised null ",
        rng_txt(x, "percentile_abs_rho_noFeature_vs_sign_randomised_null", 1),
        "; HR vs outcome-signed null ",
        rng_txt(x, "percentile_HR_vs_outcome_signed_null", 1))
  }
}

# ---- the headline ranges, written out -------------------------------------
RANGES <- rbindlist(RANGE_ROWS)
RANGES[, `:=`(minimum = signif(minimum, 6), maximum = signif(maximum, 6))]
setcolorder(RANGES, c("quantity", "stratum", "weight_provenance", "n_arms", "signatures",
                      "range_text", "minimum", "maximum"))
stopifnot(!anyNA(RANGES$weight_provenance), all(RANGES$weight_provenance %in% PROV_CLASSES))
save_tsv(RANGES, "26_headline_ranges_by_provenance.tsv")
msg("Headline ranges written: ", nrow(RANGES), " rows, each confined to one of the ",
    uniqueN(RANGES$weight_provenance), " weight provenance classes present (",
    paste(sort(unique(RANGES$weight_provenance)), collapse = ", "), ")")

# ---- 8. Per-sample scores ----
banner("8 | Per-sample score table")
PS <- data.table(sample_barcode = cohort$sample_barcode, patient = cohort$patient,
                 os_time = cohort$os_time, os_event = cohort$os_event,
                 pct_noFeature = qc$pct_noFeature, pct_multimapping = qc$pct_multimapping,
                 log10_assigned_reads = log10(qc$assigned_reads), lnc_axis = ax)
for (k in SIG_KEYS) for (w in c("published", "refit_KIRC")) for (e in c("observed", "residualised"))
  set(PS, j = paste(k, w, e, sep = "__"),
      value = round(zc(SC[[k]][[w]][[e]][PS$sample_barcode]), 5))
save_tsv(PS, "26_signature_scores_per_sample.tsv")
msg("Per-sample scores: ", nrow(PS), " rows, ", ncol(PS), " columns (each score standardised",
    " over the ", nrow(PS), " discovery tumours)")

# The per-sample table names signatures in its column names, so this key
# carries their provenance, one row per signature x weighting.
KEY <- rbindlist(lapply(SIG_KEYS, function(k) rbindlist(lapply(c("published", "refit_KIRC"), function(w)
  data.table(signature = k, label = SIGNATURES[[k]]$label,
             weight_provenance = prov_class_of(k, w),
             weight_provenance_detail = prov_detail_of(k, w),
             reporting_tier = tier_of(k), weights = w,
             in_sample_weights = in_sample_of(k, w),
             n_genes_published = nrow(SIGNATURES[[k]]$genes),
             n_genes_scored = nrow(USABLE(k)),
             n_genes_in_production_3442 = sum(USABLE(k)$in_production_3442),
             published_coefficients_available = sum(!is.na(SIGNATURES[[k]]$genes$published_coef)),
             per_sample_column_observed     = paste(k, w, "observed", sep = "__"),
             per_sample_column_residualised = paste(k, w, "residualised", sep = "__"),
             doi = SIGNATURES[[k]]$doi,
             weight_provenance_note = PROV_CLASS_NOTE)))))
stopifnot(all(KEY$per_sample_column_observed %in% names(PS)),
          all(KEY$per_sample_column_residualised %in% names(PS)),
          all(KEY$weight_provenance %in% PROV_CLASSES))
save_tsv(KEY, "26_signature_provenance_key.tsv")
print(KEY[, .(signature, weights, weight_provenance, reporting_tier,
              n_genes_scored, n_genes_published)], row.names = FALSE)

banner("26 | Done")
msg("Elapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
write_session_info("26")
