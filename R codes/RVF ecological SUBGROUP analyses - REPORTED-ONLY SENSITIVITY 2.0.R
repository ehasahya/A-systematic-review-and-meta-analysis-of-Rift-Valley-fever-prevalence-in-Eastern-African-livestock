# ============================================================================
# RVF IgG seroprevalence: ecological subgroup analyses (reported-only sensitivity)
# ----------------------------------------------------------------------------
# This script reproduces the ecological subgroup analyses of RVF IgG
# seroprevalence in livestock. For each contextual factor (flooding, border,
# national park, outbreak, urban/rural) and each species, it:
#
#   1. Loads the factor-specific dataset and keeps IgG records only.
#   2. Checks that per-study cases and sample sizes are consistent across the
#      five factor files and against the main IgG pooling dataset.
#   3. Pools prevalence on the logit scale, using a three-level model
#      (estimates nested within studies) where a study contributes more than
#      one estimate, and a random-effects model otherwise.
#   4. Estimates pooled prevalence, 95% CI, 95% prediction interval and
#      heterogeneity (Q, I2, tau2) for each subgroup level and overall.
#   5. Tests for differences between subgroups using the factor as a
#      moderator. The test is repeated with robust variance estimation (CR2),
#      a logit GLMM and the Freeman-Tukey transformation.
#   6. Draws a grouped forest plot per factor and species and writes summary
#      tables.
#
# Requirements
#   R packages: metafor, clubSandwich, dplyr, readxl, writexl
#   install.packages(c("metafor", "clubSandwich", "dplyr", "readxl", "writexl"))
#
# Input (Excel files in DATA_DIR)
#   floods reported.xlsx, border reported.xlsx, parks reported.xlsx,
#   outbreaks reported.xlsx, urban or rural reported.xlsx
#   Required columns: Author, Species, Cases, SampleSize, and the factor column
#   (see `datasets` below). Optional column: Test (filtered to IgG).
#   MAIN_PATH: main IgG pooling dataset, used only for the consistency check
#   (set to NA to skip).
#
# Output (DATA_DIR/out)
#   <factor>_<species>.png   grouped forest plots
#   level_estimates.xlsx     pooled estimates per subgroup level
#   between_tests.xlsx       between-subgroup tests and overall estimates
#   reconciliation.xlsx      data consistency check
#
# Usage
#   Set DATA_DIR and MAIN_PATH to the local folder/file locations, then run
#   the whole script. CEX controls font size; BOXMAX controls maximum box size.
# ============================================================================

## ---- settings ---------------------------------------------------------------
CEX    <- 1.3
BOXMAX <- 1.0

ASSAY        <- "IgG"
AX_MAX       <- 100
AX_STEP      <- 20

library(metafor)
library(clubSandwich)
library(dplyr)
library(readxl)
library(writexl)

## ---- paths ------------------------------------------------------------------
DATA_DIR <- "C:/Users/eha/meta"

MAIN_PATH <- "C:/Users/eha/meta"

out_dir <- file.path(DATA_DIR, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# subgroup labels used in plot headers
DISPLAY <- c(
  floods_reported         = "Flooding reported",
  Border_reported         = "Border reported",
  Parks_reported          = "National park reported",
  outbreak_reported       = "Outbreak reported",
  urban_or_rural_reported = "Urban/rural reported"
)

## ---- helper functions -------------------------------------------------------
capfirst <- function(x) { x <- as.character(x); paste0(toupper(substring(x, 1, 1)), substring(x, 2)) }
tf    <- function(x) 100 * plogis(x)
DASH  <- c(" [", " \u2013 ", "]")
fmt_p <- function(p) ifelse(is.na(p), "= NA", ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)))

safe_read <- function(f) tryCatch(readxl::read_excel(f),
                                  error = function(e) { message("Could not read: ", f, "  (", conditionMessage(e), ")"); NULL })

# clean columns, keep the selected assay, drop rows without usable counts
prep <- function(d) {
  names(d) <- trimws(names(d))
  d$Species    <- trimws(tolower(as.character(d$Species)))
  d$Cases      <- suppressWarnings(as.numeric(d$Cases))
  d$SampleSize <- suppressWarnings(as.numeric(d$SampleSize))
  d$Author     <- trimws(as.character(d$Author))
  if ("Test" %in% names(d)) d <- d[tolower(trimws(as.character(d$Test))) == tolower(ASSAY), , drop = FALSE]
  d[!is.na(d$Cases) & !is.na(d$SampleSize) & d$SampleSize > 0, , drop = FALSE]
}

## ---- step 1: load factor datasets -------------------------------------------
datasets <- list(
  floods        = list(file = "floods reported.xlsx",         group = "floods_reported"),
  border        = list(file = "border reported.xlsx",         group = "Border_reported"),
  national_park = list(file = "parks reported.xlsx",          group = "Parks_reported"),
  outbreak      = list(file = "outbreaks reported.xlsx",      group = "outbreak_reported"),
  urban         = list(file = "urban or rural reported.xlsx", group = "urban_or_rural_reported")
)
for (nm in names(datasets)) {
  raw <- safe_read(file.path(DATA_DIR, datasets[[nm]]$file))
  datasets[[nm]]$data <- if (is.null(raw)) NULL else prep(raw)
}

## ---- step 2: data consistency check -----------------------------------------
study_totals <- function(d, label) {
  d %>% group_by(Species, Author) %>%
    summarise(cases = sum(Cases), n = sum(SampleSize), .groups = "drop") %>%
    mutate(file = label)
}
recon <- bind_rows(lapply(names(datasets), function(nm)
  if (is.null(datasets[[nm]]$data)) NULL else study_totals(datasets[[nm]]$data, nm)))

if (nrow(recon)) {
  n_files <- n_distinct(recon$file)
  # (a) studies whose totals differ between factor files, or absent from a file
  mism_files <- recon %>% group_by(Species, Author) %>%
    filter(n_distinct(cases) > 1 | n_distinct(n) > 1 | n() < n_files) %>%
    ungroup() %>% arrange(Species, Author, file)
  
  # (b) per-study totals compared with the main IgG dataset
  mism_main <- tibble()
  if (!is.na(MAIN_PATH)) {
    main_raw <- safe_read(MAIN_PATH)
    if (!is.null(main_raw)) {
      main <- prep(main_raw)
      main <- main[main$Species %in% unique(recon$Species), , drop = FALSE]
      main_tot <- main %>% group_by(Species, Author) %>%
        summarise(cases_main = sum(Cases), n_main = sum(SampleSize), .groups = "drop")
      mism_main <- recon %>% left_join(main_tot, by = c("Species", "Author")) %>%
        mutate(status = case_when(is.na(n_main) ~ "not in main IgG file (check name/year)",
                                  cases != cases_main | n != n_main ~ "totals differ from main IgG file",
                                  TRUE ~ "match")) %>%
        filter(status != "match") %>% arrange(Species, Author, file)
    }
  }
  write_xlsx(list(across_factor_files = mism_files, versus_main_IgG = mism_main),
             file.path(out_dir, "reconciliation.xlsx"))
  cat("\n==== RECONCILIATION ====\n")
  cat("Studies inconsistent ACROSS factor files:", n_distinct(paste(mism_files$Species, mism_files$Author)), "\n")
  if (nrow(mism_files)) print(as.data.frame(mism_files), row.names = FALSE)
  cat("\nRows differing from the MAIN IgG file:", nrow(mism_main), "\n")
  if (nrow(mism_main)) print(as.data.frame(mism_main), row.names = FALSE)
  cat("(full lists in out/reconciliation.xlsx)\n\n")
}

## ---- step 3-4: model functions ----------------------------------------------
# multilevel I2 (Viechtbauer) for rma.mv
ml_I2 <- function(res) tryCatch({
  W <- diag(1 / res$vi, nrow = length(res$vi)); X <- res$X
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  typ <- (res$k - res$p) / sum(diag(P))
  100 * sum(res$sigma2) / (sum(res$sigma2) + typ)
}, error = function(e) NA_real_)

# pool a set of estimates: three-level model if a study gives >1 estimate,
# random-effects model otherwise; results on the logit scale
fit_set <- function(sub) {
  k <- nrow(sub); if (k == 0) return(NULL)
  sub$esid <- seq_len(k)
  es <- escalc(measure = "PLO", xi = Cases, ni = SampleSize, data = sub, add = 1/2, to = "only0")
  base <- list(k = k, ns = n_distinct(sub$StudyID), cases = sum(sub$Cases), n = sum(sub$SampleSize),
               mu = NA, ci.lb = NA, ci.ub = NA, pi.lb = NA, pi.ub = NA,
               I2 = NA, tau2 = NA, QE = NA, QEdf = NA, QEp = NA, repl = FALSE, res = NULL)
  if (k == 1) {
    se <- sqrt(es$vi[1])
    base[c("mu", "ci.lb", "ci.ub")] <- list(es$yi[1], es$yi[1] - qnorm(0.975) * se, es$yi[1] + qnorm(0.975) * se)
    return(base)
  }
  repl <- any(table(sub$StudyID) > 1)
  res <- tryCatch(if (repl) rma.mv(yi, vi, random = ~1 | StudyID/esid, data = es, method = "REML")
                  else       rma(yi, vi, data = es, method = "REML"), error = function(e) NULL)
  if (is.null(res)) return(base)
  mu <- as.numeric(coef(res)); se <- as.numeric(res$se)
  t2 <- if (repl) sum(res$sigma2) else res$tau2
  # prediction interval on the t distribution (df = k - p)
  tcr <- qt(0.975, df = max(1, res$k - res$p))
  base$mu <- mu; base$ci.lb <- as.numeric(res$ci.lb); base$ci.ub <- as.numeric(res$ci.ub)
  base$pi.lb <- mu - tcr * sqrt(se^2 + t2); base$pi.ub <- mu + tcr * sqrt(se^2 + t2)
  base$I2   <- if (repl) ml_I2(res) else res$I2
  base$tau2 <- t2
  base$QE <- res$QE; base$QEdf <- res$k - res$p; base$QEp <- res$QEp
  base$repl <- repl; base$res <- res
  base
}

## ---- step 6: forest plot functions ------------------------------------------
# row layout per subgroup: header, studies, subtotal, heterogeneity, PI, gap
layout_rows <- function(ks) {
  L <- length(ks); top <- sum(ks) + sum(ifelse(ks >= 2, 5.5, 3.5)); cur <- top
  studyrows <- numeric(0); head <- summ <- het <- pir <- numeric(L)
  for (i in seq_len(L)) {
    head[i] <- cur; cur <- cur - 1.5
    studyrows <- c(studyrows, seq(cur, by = -1, length.out = ks[i])); cur <- cur - ks[i]
    summ[i] <- cur; het[i] <- cur - 1; pir[i] <- cur - 2
    cur <- cur - if (ks[i] >= 2) 4 else 2
  }
  list(top = top, studyrows = studyrows, head = head, summ = summ, het = het, pir = pir)
}

draw_forest <- function(d, g, lev_order, ks, overall, per_level, qm_line,
                        n_multi, nm, sp, png_path) {
  lay <- layout_rows(ks)
  res_all <- overall$res
  k <- nrow(d)
  wt <- tryCatch({ w <- weights(res_all); 100 * w / sum(w) }, error = function(e) rep(NA_real_, k))
  wt_lab <- ifelse(is.na(wt), "", sprintf("%.1f", wt))
  
  W <- AX_MAX; at <- seq(0, AX_MAX, by = AX_STEP)
  xlim  <- c(-2.05 * W, 1.75 * W)
  ilabx <- c(-1.00 * W, -0.66 * W, -0.32 * W)
  
  n_foot <- 3
  foot_line0 <- 4.0
  foot_step  <- 1.35 * CEX
  mar_bottom <- foot_line0 + foot_step * n_foot + 1
  width_in  <- max(12, 4 + 7.5 * CEX)
  height_in <- max(7, 1 + 0.30 * CEX * (lay$top + 8) + 0.2 * mar_bottom)
  
  glab <- if (g %in% names(DISPLAY)) unname(DISPLAY[g]) else g
  wt_def <- if (overall$repl) "weight in the three-level model" else "weight in the random-effects model"
  
  png(png_path, width = width_in, height = height_in, units = "in", res = 300)
  op0 <- par(mar = c(mar_bottom, 3, 3, 2))
  tryCatch({
    pi_bar <- function(lo, hi, row, lab, lwd) {
      cap <- 0.18
      segments(tf(lo), row, tf(hi), row, col = "red", lwd = lwd)
      segments(tf(c(lo, hi)), row - cap, tf(c(lo, hi)), row + cap, col = "red", lwd = 1.2 * CEX)
      text(xlim[1], row, lab, pos = 4, cex = CEX)
      text(xlim[2], row, sprintf("[%.2f \u2013 %.2f]", tf(lo), tf(hi)), pos = 2, cex = CEX)
    }
    het_txt <- function(x) sprintf("Heterogeneity: Q = %.1f, df = %d, p %s;  I\u00b2 = %.1f%%,  \u03c4\u00b2 = %.3f",
                                   x$QE, as.integer(x$QEdf), fmt_p(x$QEp), x$I2, x$tau2)
    
    forest(res_all, transf = tf, olim = c(0, 100), refline = NA, rows = lay$studyrows,
           slab = d$Author,
           ilab = data.frame(Cases = as.character(d$Cases), N = as.character(d$SampleSize), Weight = wt_lab),
           ilab.xpos = ilabx, ilab.pos = 2,
           header = c("Study", "Prevalence (%) [95% CI]"),
           xlim = xlim, alim = c(0, AX_MAX), at = at, ylim = c(-3.5, lay$top + 3),
           digits = c(2L, 0L),
           xlab = paste0("RVF ", ASSAY, " prevalence (%) \u2014 ", sp),
           col = "maroon", cex = CEX,
           plim = c(0.4, BOXMAX), annosym = DASH, mlab = "Pooled")
    
    op <- par(xpd = TRUE)
    # column headers
    text(ilabx, lay$top + 2, c("Cases", "N", "Wt %"), pos = 2, cex = CEX, font = 2)
    # subgroup headers, subtotals, heterogeneity and prediction intervals
    for (i in seq_along(lev_order)) {
      pl <- per_level[[i]]
      text(xlim[1], lay$head[i], paste0(glab, " = ", capfirst(lev_order[i])), pos = 4, cex = CEX, font = 2)
      if (!is.null(pl) && !is.na(pl$mu)) {
        addpoly(pl$mu, ci.lb = pl$ci.lb, ci.ub = pl$ci.ub, rows = lay$summ[i], transf = tf,
                col = "navy", annosym = DASH, cex = CEX, mlab = sprintf("Subtotal (k = %d)", pl$k))
        text(ilabx[1], lay$summ[i], pl$cases, pos = 2, cex = CEX, font = 2)
        text(ilabx[2], lay$summ[i], pl$n,     pos = 2, cex = CEX, font = 2)
        if (pl$k >= 2 && !is.na(pl$QE)) {
          text(xlim[1], lay$het[i], het_txt(pl), pos = 4, cex = CEX * 0.85)
          if (!is.na(pl$pi.lb)) pi_bar(pl$pi.lb, pl$pi.ub, lay$pir[i], "95% prediction interval", 1.5 * CEX)
        }
      }
    }
    # overall totals and overall prediction interval
    text(ilabx[1], -1, sum(d$Cases),      pos = 2, cex = CEX, font = 2)
    text(ilabx[2], -1, sum(d$SampleSize), pos = 2, cex = CEX, font = 2)
    text(ilabx[3], -1, "100.0",           pos = 2, cex = CEX, font = 2)
    pi_bar(overall$pi.lb, overall$pi.ub, -2, "95% prediction interval (overall)", 2.5 * CEX)
    
    # footer below the axis
    foot <- function(txt, j, cx = CEX, fnt = 1)
      mtext(txt, side = 1, line = foot_line0 + foot_step * (j - 1), at = xlim[1], adj = 0, cex = cx, font = fnt)
    oh <- het_txt(overall)
    foot(paste0("Overall h", substring(oh, 2)), 1)
    foot(qm_line, 2)
    multi_txt <- if (n_multi == 0) "" else sprintf("  %d stud%s contribute%s to more than one subgroup.",
                                                   n_multi, ifelse(n_multi == 1, "y", "ies"), ifelse(n_multi == 1, "s", ""))
    foot(paste0("N = number sampled;  Wt % = ", wt_def, ".", multi_txt), 3, CEX * 0.85, 3)
    par(op)
  }, error = function(e) message("forest failed ", nm, "/", sp, ": ", conditionMessage(e)))
  par(op0); dev.off()
}

## ---- steps 3-6: analysis per factor x species -------------------------------
level_rows <- list(); between_rows <- list()

analyse <- function(d, nm, g, sp) {
  d <- d[d$Species == sp, , drop = FALSE]
  d <- d[!is.na(d[[g]]) & trimws(as.character(d[[g]])) != "", , drop = FALSE]
  if (nrow(d) < 2) { message("skip ", nm, "/", sp, " (k<2)"); return(invisible()) }
  d$StudyID <- d$Author
  d[[g]] <- trimws(as.character(d[[g]]))
  
  # level order: alphabetical
  lev_order <- sort(unique(d[[g]]))
  d <- d[order(match(d[[g]], lev_order), d$Author), , drop = FALSE]
  d$grp <- factor(d[[g]], levels = lev_order)
  ks <- as.integer(table(d$grp))
  
  # overall pooled estimate
  overall <- fit_set(d)
  if (is.null(overall$res)) { message("overall model failed ", nm, "/", sp); return(invisible()) }
  
  # number of studies contributing to more than one level
  n_multi <- sum(tapply(d[[g]], d$StudyID, function(z) n_distinct(z) > 1))
  
  # step 5: between-subgroup test
  d_t <- d
  d_t$esid <- seq_len(nrow(d_t))
  tested <- paste(levels(d_t$grp), collapse = " vs ")
  QM <- QMdf <- QMp <- rveP <- rveDF <- glmmP <- ftP <- min_st <- NA
  repl_t <- any(table(d_t$StudyID) > 1)
  if (nlevels(d_t$grp) >= 2) {
    min_st <- min(tapply(d_t$StudyID, d_t$grp, n_distinct))
    es_t <- escalc(measure = "PLO", xi = Cases, ni = SampleSize, data = d_t, add = 1/2, to = "only0")
    res_mod <- tryCatch(if (repl_t) rma.mv(yi, vi, mods = ~ grp, random = ~1 | StudyID/esid, data = es_t, method = "REML")
                        else        rma(yi, vi, mods = ~ grp, data = es_t, method = "REML"),
                        error = function(e) NULL)
    if (!is.null(res_mod)) {
      QM <- res_mod$QM; QMdf <- res_mod$QMdf[1]; QMp <- res_mod$QMp
      rob <- tryCatch(robust(res_mod, cluster = es_t$StudyID, clubSandwich = TRUE), error = function(e) NULL)
      if (!is.null(rob)) { rveP <- rob$QMp; rveDF <- rob$QMdf[2] }
    }
    glmmP <- tryCatch(suppressWarnings(rma.glmm(measure = "PLO", xi = Cases, ni = SampleSize,
                                                mods = ~ grp, data = d_t))$QMp, error = function(e) NA)
    es_ft <- escalc(measure = "PFT", xi = Cases, ni = SampleSize, data = d_t)
    ftP <- tryCatch((if (repl_t) rma.mv(yi, vi, mods = ~ grp, random = ~1 | StudyID/esid, data = es_ft, method = "REML")
                     else        rma(yi, vi, mods = ~ grp, data = es_ft, method = "REML"))$QMp,
                    error = function(e) NA)
  }
  # RVE flagged as unreliable (reported in table only) when a tested level has
  # <2 studies or CR2 Satterthwaite df < 4
  rve_flag <- !is.na(QM) && ((!is.na(rveDF) && rveDF < 4) || (!is.na(min_st) && min_st < 2))
  qm_line <- if (is.na(QM)) {
    "Between-subgroup test: not estimable (fewer than two subgroups)"
  } else sprintf("Between-subgroup (%s): Q = %.1f, df = %d, p %s   |   RVE p %s   |   GLMM p %s   |   FT p %s",
                 ifelse(repl_t, "3-level", "RE"), QM, as.integer(QMdf), fmt_p(QMp),
                 fmt_p(rveP), fmt_p(glmmP), fmt_p(ftP))
  
  # step 4: per-level estimates
  per_level <- lapply(lev_order, function(L) fit_set(d[d[[g]] == L, , drop = FALSE]))
  
  # step 6: forest plot
  draw_forest(d, g, lev_order, ks, overall, per_level, qm_line, n_multi, nm, sp,
              file.path(out_dir, paste0(nm, "_", ifelse(sp == "small ruminants", "sr", gsub(" ", "_", sp)), ".png")))
  
  # summary tables
  for (i in seq_along(lev_order)) {
    pl <- per_level[[i]]; if (is.null(pl)) next
    level_rows[[length(level_rows) + 1]] <<- tibble(
      dataset = nm, factor = g, species = sp, level = lev_order[i],
      k = pl$k, n_studies = pl$ns, cases = pl$cases, n = pl$n,
      prev_pct = tf(pl$mu), ci_lo = tf(pl$ci.lb), ci_hi = tf(pl$ci.ub),
      pi_lo = tf(pl$pi.lb), pi_hi = tf(pl$pi.ub),
      QE = pl$QE, QE_df = pl$QEdf, QE_p = pl$QEp, I2_pct = pl$I2, tau2 = pl$tau2,
      model = if (pl$k == 1) "single study" else if (pl$repl) "3-level logit" else "RE logit")
  }
  between_rows[[length(between_rows) + 1]] <<- tibble(
    dataset = nm, factor = g, species = sp, k_total = nrow(d), n_studies = n_distinct(d$StudyID),
    n_levels = length(lev_order), levels_tested = tested,
    n_studies_multi_level = n_multi, min_studies_per_level = min_st,
    overall_prev_pct = tf(overall$mu), overall_ci_lo = tf(overall$ci.lb), overall_ci_hi = tf(overall$ci.ub),
    overall_pi_lo = tf(overall$pi.lb), overall_pi_hi = tf(overall$pi.ub),
    overall_QE = overall$QE, overall_QE_df = overall$QEdf, overall_QE_p = overall$QEp,
    overall_I2_pct = overall$I2, overall_tau2 = overall$tau2,
    QM = QM, QM_df = QMdf, QM_p = QMp, RVE_p = rveP, RVE_df2 = rveDF, RVE_unreliable = rve_flag,
    GLMM_p = glmmP, FT_p = ftP)
  invisible()
}

## ---- run --------------------------------------------------------------------
for (nm in names(datasets)) {
  d <- datasets[[nm]]$data; g <- datasets[[nm]]$group
  if (is.null(d)) { message("SKIP '", nm, "': not loaded."); next }
  if (!g %in% names(d)) { message("SKIP '", nm, "': column '", g, "' not found."); next }
  for (sp in sort(unique(d$Species))) analyse(d, nm, g, sp)
}

if (length(level_rows))   write_xlsx(bind_rows(level_rows),   file.path(out_dir, "level_estimates.xlsx"))
if (length(between_rows)) write_xlsx(bind_rows(between_rows), file.path(out_dir, "between_tests.xlsx"))
cat("\nDone. Forest plots, tables and reconciliation written to:\n", out_dir, "\n")