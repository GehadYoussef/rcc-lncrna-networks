# ora_hypergeometric.R: hypergeometric GO over-representation test
# Used by the sign-stratified enrichment (18). Reimplements clusterProfiler::enrichGO on
# the same GO biological-process annotation and gene-set size limits, and reproduces
# the pooled clusterProfiler results.
if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))

suppressPackageStartupMessages({library(data.table); library(AnnotationDbi); library(org.Hs.eg.db)})

# Settings match the enrichGO call in 10_paper_analyses.R: ont = "BP",
# keyType = "ENSEMBL", BH adjustment, gene-set size 10 to 500. GOALL gives the
# same ancestor propagation as clusterProfiler.
build_go_map <- function(universe) {
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = universe,
        keytype = "ENSEMBL", columns = c("GOALL", "ONTOLOGYALL")))
  m <- as.data.table(m)[ONTOLOGYALL == "BP" & !is.na(GOALL)]
  unique(m[, .(ENSEMBL, GO = GOALL)])
}
ora <- function(query, universe, gomap, minGSSize = 10, maxGSSize = 500) {
  gomap <- gomap[ENSEMBL %in% universe]
  # As in DOSE, annotated genes are counted on the full map before terms are
  # filtered by size. Filtering first would shrink the denominators.
  N  <- length(unique(gomap$ENSEMBL))          # annotated universe
  q  <- intersect(query, unique(gomap$ENSEMBL))
  K  <- length(q)                              # annotated query
  sizes <- gomap[, .N, by = GO][N >= minGSSize & N <= maxGSSize]
  hit <- gomap[ENSEMBL %in% q, .(k = .N), by = GO][GO %in% sizes$GO]
  res <- merge(hit, sizes[, .(GO, m = N)], by = "GO")
  res[, pvalue := phyper(k - 1, m, N - m, K, lower.tail = FALSE)]
  res[, p.adjust := p.adjust(pvalue, "BH")]
  res[, `:=`(GeneRatio = paste0(k, "/", K), BgRatio = paste0(m, "/", N))]
  trm <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = res$GO,
           keytype = "GOID", columns = "TERM"))
  res[, Description := trm$TERM[match(GO, trm$GOID)]]
  setorder(res, pvalue)
  res[, .(GO, Description, GeneRatio, BgRatio, Count = k, pvalue, p.adjust,
          annotated_query = K, annotated_universe = N)]
}
