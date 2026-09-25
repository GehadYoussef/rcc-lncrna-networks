# 00_config_sc.R: bootstrap sourced by every single-cell stage script.
# Resolves paths from this file's location, reads all scientific settings and
# thresholds from src/config/singlecell.yml, and defines logging and output
# helpers.

options(stringsAsFactors = FALSE, warn = 1)

if (!exists("SC_DIR") || is.null(SC_DIR)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .f    <- sub("^--file=", "", .args[grep("^--file=", .args)])
  SC_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/"))
            else normalizePath(file.path(getwd()), winslash = "/")
}
SUB_ROOT <- normalizePath(file.path(SC_DIR, "..", ".."), winslash = "/", mustWork = TRUE)

source(file.path(SC_DIR, "lib_sc.R"))

CFG_FILE <- file.path(SUB_ROOT, "src", "config", "singlecell.yml")
CFG <- yaml::read_yaml(CFG_FILE)
CFG_HASH <- digest::digest(file = CFG_FILE, algo = "sha256")
set.seed(CFG$seed)

p <- function(x) normalizePath(file.path(SUB_ROOT, x), winslash = "/", mustWork = FALSE)
MANIFEST     <- p(CFG$paths$manifest)
SAMPLE_SHEET <- p(CFG$paths$sample_sheet)
EXTRA_FILES  <- p(CFG$paths$extra_files)
RAW_DIR      <- p(CFG$paths$raw_dir)
DERIVED_DIR  <- p(CFG$paths$derived_dir)
RESULTS_DIR  <- p(CFG$paths$results_dir)
FIG_DIR      <- p(CFG$paths$figures_dir)
LOG_DIR      <- p(CFG$paths$logs_dir)
REPORT_DIR   <- p(CFG$paths$reports_dir)
GENCODE_GTF  <- p(CFG$paths$gencode_gtf)
LEGACY_GTF   <- p(CFG$paths$legacy_gtf)
BULK_TCGA    <- bulk_path(CFG$paths$bulk_tcga, must_exist = FALSE)
BULK_CPTAC   <- bulk_path(CFG$paths$bulk_cptac, must_exist = FALSE)

for (d in c(RAW_DIR, DERIVED_DIR, RESULTS_DIR, FIG_DIR, LOG_DIR, REPORT_DIR))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- logging and output helpers ----
.LOG_FILE <- NULL
msg <- function(...) {
  line <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., collapse = "")
  cat(line, "\n", sep = "")
  if (!is.null(.LOG_FILE)) cat(line, "\n", sep = "", file = .LOG_FILE, append = TRUE)
}
banner <- function(x) {
  cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
}
start_log <- function(tag) {
  .LOG_FILE <<- file.path(LOG_DIR, paste0(tag, ".log"))
  cat("", file = .LOG_FILE)
  msg("start ", tag, " | config sha256 ", substr(CFG_HASH, 1, 12), NULL)
}
save_tsv <- function(x, file) {
  f <- file.path(RESULTS_DIR, file)
  fwrite(as.data.table(x), f, sep = "\t", na = "", quote = FALSE)
  invisible(f)
}
read_tsv <- function(file) fread(file.path(RESULTS_DIR, file), sep = "\t", na.strings = "")
write_session_info <- function(tag) {
  f <- file.path(LOG_DIR, paste0("sessionInfo_", tag, ".txt"))
  con <- file(f, open = "wt"); sink(con); print(sessionInfo()); sink(); close(con)
  invisible(f)
}
FIG_DPI <- 600
save_fig <- function(plot, name, width = 7, height = 5) {
  for (fmt in c("svg", "png")) {
    f <- file.path(FIG_DIR, paste0(name, ".", fmt))
    ggplot2::ggsave(f, plot, width = width, height = height, units = "in", dpi = FIG_DPI,
                    device = if (fmt == "svg") svglite::svglite else ragg::agg_png)
  }
  invisible(name)
}

read_manifest <- function() {
  man <- fread(MANIFEST, colClasses = "character", na.strings = NULL)
  smp <- fread(SAMPLE_SHEET, colClasses = "character", na.strings = NULL)
  smp[, include := toupper(include) == "TRUE"]
  list(manifest = man, samples = smp)
}
