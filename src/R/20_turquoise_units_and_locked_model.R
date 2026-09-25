# 20_turquoise_units_and_locked_model.R: unit-matched replication and drop test of the locked model in CPTAC-3.
#
# Reads cache/locked_model.rds and results/09_module_replication.tsv (run after
# 09). Nothing is refitted. Writes:
#   20_replication_unitmatched.tsv: per-module CPTAC-3 hazard ratios under
#     cohort and discovery-fixed standardisation, one row per module.
#   20_score_spread_by_standardisation.tsv: spread of the discovery-standardised
#     scores in CPTAC-3 (SD 1.0 = the same spread as in discovery).
#   20_bootstrap_delta_cindex.tsv: drop test. Each retained module coefficient
#     is set to zero in turn, then all lncRNA modules at once, and the change in
#     concordance over the comparator-only predictor is bootstrapped.

if (!exists("R_DIR")) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  R_DIR <- if (length(.f)) dirname(normalizePath(.f[1], winslash = "/")) else getwd()
}
source(file.path(R_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(data.table); library(survival) })
banner("20 | Locked model: unit-matched replication and drop test")

L <- readRDS(file.path(CACHE_DIR, "locked_model.rds"))
stopifnot(identical(L$version, "v9"))
yv <- L$yv
msg("CPTAC-3 evaluated set: n = ", length(L$vsamp), ", events = ", sum(yv[, 2]))

# ---- unit-matched replication and score spread ----
rep09 <- fread(file.path(RESULTS_DIR, "09_module_replication.tsv"))
rep_tbl <- rep09[, .(module = paste(ifelse(biotype == "mRNA", "protein-coding", "lncRNA"), module),
                     HR_cohort = HR_cptac, lo_c = lo, hi_c = hi, p_c = p_cptac,
                     HR_disc = HR_cptac_disc, lo_d = lo_disc, hi_d = hi_disc, p_d = p_disc,
                     n = n, events = events)]
rep_tbl[, `:=`(fdr_c = signif(p.adjust(p_c, "BH"), 3), fdr_d = signif(p.adjust(p_d, "BH"), 3))]
setcolorder(rep_tbl, c("module", "HR_cohort", "lo_c", "hi_c", "p_c",
                       "HR_disc", "lo_d", "hi_d", "p_d", "fdr_c", "fdr_d", "n", "events"))
print(rep_tbl, row.names = FALSE)
save_tsv(rep_tbl, "20_replication_unitmatched.tsv")

# Spread of the discovery-standardised scores within CPTAC-3, and agreement
# between the two standardisations of the same module.
common <- intersect(colnames(L$Ev), colnames(L$Ev_disc))
spread <- rbindlist(lapply(common, function(m) data.table(
  module = m,
  r_cohort_vs_discovery   = round(cor(L$Ev[, m], L$Ev_disc[, m]), 4),
  rho_cohort_vs_discovery = round(cor(L$Ev[, m], L$Ev_disc[, m], method = "spearman"), 4),
  mean_discovery_standardised = round(mean(L$Ev_disc[, m]), 3),
  sd_discovery_standardised   = round(sd(L$Ev_disc[, m]), 3))))
print(spread, row.names = FALSE)
save_tsv(spread, "20_score_spread_by_standardisation.tsv")

# ---- drop test on the locked combined model ----
lp_of <- function(X, b) as.numeric(X[, names(b), drop = FALSE] %*% b)

drop_test <- function(comparator, standardisation, Xv, Ev) {
  M      <- L$models[[comparator]]
  XE     <- cbind(Xv, Ev)
  lp_ref <- lp_of(Xv, M$b_clin)
  b_full <- M$b_full
  sel    <- names(b_full)[names(b_full) %in% colnames(Ev)]
  arms <- list(`+ modules (all retained)` = b_full)
  for (m in sel) { b <- b_full; b[m] <- 0; arms[[paste0("+ modules minus ", m)]] <- b }
  lnc_sel <- grep("^lnc_ME", sel, value = TRUE)
  if (length(lnc_sel)) {
    b <- b_full; b[lnc_sel] <- 0
    arms[["+ modules minus all lncRNA modules"]] <- b
  }
  C_ref <- cindex(yv, lp_ref)
  rbindlist(lapply(names(arms), function(nm) {
    lp <- lp_of(XE, arms[[nm]])
    r  <- paired_boot_delta_c(yv, lp, lp_ref)
    data.table(comparator = comparator, standardisation = standardisation, model = nm,
               n = length(lp), events = sum(yv[, 2]),
               C_comparator = round(C_ref, 4), C_model = round(cindex(yv, lp), 4),
               delta_C = round(r[["delta"]], 4), lo = round(r[["lo"]], 4), hi = round(r[["hi"]], 4),
               p_bootstrap = signif(r[["p_boot"]], 3),
               prop_positive = round(r[["prop_positive"]], 3), n_boot = r[["n_boot"]])
  }))
}

boot_tbl <- rbindlist(list(
  drop_test("clinical",  "cohort",    L$X$clinical$Xv,       L$Ev),
  drop_test("clinical",  "discovery", L$X$clinical$Xv_disc,  L$Ev_disc),
  drop_test("augmented", "cohort",    L$X$augmented$Xv,      L$Ev)))
for (cmp in unique(boot_tbl$comparator)) {
  msg("Comparator: ", cmp, " -- retained modules: ",
      paste(L$models[[cmp]]$sel_me, collapse = ", "))
}
print(boot_tbl, row.names = FALSE)
save_tsv(boot_tbl, "20_bootstrap_delta_cindex.tsv")

write_session_info("20_turquoise_units_and_locked_model")
banner("20 | done")
