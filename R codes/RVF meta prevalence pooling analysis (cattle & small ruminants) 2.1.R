# ============================================================================
# RVF prevalence: pooled estimates for cattle and small ruminants
# ----------------------------------------------------------------------------
# This script reproduces the pooled RVF prevalence estimates in cattle and
# small ruminants, separately for each assay: RT-PCR, anti-RVFV IgM and
# anti-RVFV IgG. For each species x assay combination it:
#
#   1. Loads the dataset, keeps cattle and small ruminant records with usable
#      cases and sample sizes, and checks the assay recorded in the Test column.
#   2. Counts the number of estimates contributed by each study.
#   3. Pools prevalence on the logit scale, using a three-level model
#      (estimates nested within studies) where a study contributes more than
#      one estimate, and a random-effects model otherwise.
#   4. Estimates pooled prevalence, 95% CI, 95% prediction interval and
#      heterogeneity (Q, I2, tau2 between and within studies).
#   5. Runs sensitivity analyses: robust variance estimation (CR2) CI, a logit
#      GLMM and the Freeman-Tukey double arcsine transformation.
#   6. Draws a forest plot of the primary model (fixed 0-100% axis, legend,
#      heterogeneity statistics) and writes a results table.
#
# Requirements
#   R packages: metafor, clubSandwich, dplyr, readxl, writexl
#   install.packages(c("metafor", "clubSandwich", "dplyr", "readxl", "writexl"))
#
# Input
#   in_path: meta RVF overall.xlsx
#   Required columns: Author, Species, Test, Cases, SampleSize
#   Test values: "RT-PCR", "IgM", "IgG"
#
# Output (folder "R3_pooling_results_v11" next to the input file)
#   forest_<species>_<assay>.png    forest plots
#   R3_pooled_results.xlsx          pooled estimates and sensitivity analyses
#   R3_estimates_per_study.xlsx     number of estimates per study
#
# Usage
#   Set in_path to the local file location, then run the whole script.
#   CEX controls font size; BOXMAX controls maximum box size.
# ============================================================================

## ---- settings ---------------------------------------------------------------
CEX    <- 2
BOXMAX <- 1.0

library(metafor)
library(clubSandwich)
library(dplyr)
library(readxl)
library(writexl)

DASH <- c(" [", " \u2013 ", "]")   # CI separator
AX_MAX  <- 100                     # same horizontal axis on every plot
AX_STEP <- 20
ASSAYS  <- c("RT-PCR", "IgM", "IgG")   # expected values of the Test column

## ---- paths ------------------------------------------------------------------
in_path <- "C:/Users/eha/meta/meta RVF overall.xlsx"
out_dir <- file.path(dirname(in_path), "R3_pooling_results_v11")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- step 1: load and prepare data ------------------------------------------
d0 <- read_excel(in_path)
names(d0) <- trimws(names(d0))
d0 <- d0 %>%
  mutate(Species = trimws(tolower(as.character(Species))),
         Test    = trimws(as.character(Test)),
         Cases      = suppressWarnings(as.numeric(Cases)),
         SampleSize = suppressWarnings(as.numeric(SampleSize)),
         Author     = trimws(as.character(Author))) %>%
  filter(Species %in% c("cattle", "small ruminants")) %>%
  filter(!is.na(Cases), !is.na(SampleSize), SampleSize > 0)

# stop if the Test column contains anything other than the three assays
bad <- setdiff(unique(d0$Test), ASSAYS)
if (length(bad)) stop("Unexpected Test value(s): ", paste0("'", bad, "'", collapse = ", "),
                      ". Expected: ", paste(ASSAYS, collapse = ", "))
d0$Assay   <- d0$Test
d0$StudyID <- d0$Author

## ---- step 2: estimates per study --------------------------------------------
write_xlsx(d0 %>% count(StudyID, Species, Assay, name = "n_estimates") %>% arrange(desc(n_estimates)),
           file.path(out_dir, "R3_estimates_per_study.xlsx"))

## ---- helper functions -------------------------------------------------------
tf <- function(x) 100 * plogis(x)

# multilevel I2 (Viechtbauer) for rma.mv
ml_I2 <- function(res, vi) tryCatch({
  W <- diag(1 / vi); X <- model.matrix(res)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  typ <- (res$k - res$p) / sum(diag(P))
  100 * sum(res$sigma2) / (sum(res$sigma2) + typ)
}, error = function(e) NA_real_)

## ---- step 6: forest plot function -------------------------------------------
draw_forest <- function(res, dat, es, repl, I2, tau_b, tau_w, sp, assay, png_path) {
  k <- res$k
  wt <- tryCatch({ w <- weights(res); 100 * w / sum(w) }, error = function(e) rep(NA_real_, k))
  wt_lab <- ifelse(is.na(wt), "", sprintf("%.1f", wt))
  
  # prediction interval on the t distribution (df = k - p)
  mu  <- as.numeric(coef(res)); se <- as.numeric(res$se)
  t2  <- if (repl) sum(res$sigma2) else res$tau2
  tcr <- qt(0.975, df = max(1, res$k - res$p))
  pil <- mu - tcr * sqrt(se^2 + t2)
  pih <- mu + tcr * sqrt(se^2 + t2)
  
  # heterogeneity line; tau2 shown as between / within study for the 3-level model
  pstr <- if (is.na(res$QEp)) "" else if (res$QEp < 0.001) " < 0.001" else sprintf(" = %.3f", res$QEp)
  if (repl) {
    modlab <- "Three-level model"
    het <- sprintf("Heterogeneity:  Q = %.0f, df = %d, p%s;  I\u00b2 = %.1f%%;  \u03c4\u00b2 = %.3f / %.3f",
                   res$QE, res$k - res$p, pstr, ifelse(is.na(I2), NA, I2), tau_b, tau_w)
  } else {
    modlab <- "Random-effects model"
    het <- sprintf("Heterogeneity:  Q = %.0f, df = %d, p%s;  I\u00b2 = %.1f%%;  \u03c4\u00b2 = %.3f",
                   res$QE, res$k - res$p, pstr, ifelse(is.na(I2), NA, I2), tau_b)
  }
  
  # fixed x-axis; weight definition for the footnote
  amax <- AX_MAX
  at   <- seq(0, AX_MAX, by = AX_STEP); W <- amax
  wt_def <- if (repl) "weight in the three-level model" else "weight in the random-effects model"
  xlim  <- c(-2.05 * W, 1.75 * W)
  ilabx <- c(-1.00 * W, -0.66 * W, -0.32 * W)          # Cases, N, Wt %
  
  width_in  <- max(11, 4 + 7.0 * CEX)
  height_in <- max(6, 1 + 0.30 * CEX * (k + 10))
  
  png(png_path, width = width_in, height = height_in, units = "in", res = 300)
  op0 <- par(mar = c(17, 3, 3, 2))
  tryCatch({
    forest(res, transf = tf, olim = c(0, 100), refline = NA, slab = dat$Author,
           ilab = data.frame(Cases = as.character(dat$Cases),
                             N     = as.character(dat$SampleSize),
                             Weight = wt_lab),
           ilab.xpos = ilabx, ilab.lab = c("Cases", "N", "Wt %"), ilab.pos = 2,
           header = c("Study", "Prevalence (%) [95% CI]"),
           xlim = xlim, alim = c(0, amax), at = at, ylim = -3.5,
           digits = c(2L, 0L), xlab = "",
           col = "maroon", cex = CEX, cex.lab = CEX, cex.axis = 0.9 * CEX,
           plim = c(0.4, BOXMAX), annosym = DASH, mlab = modlab)
    
    op <- par(xpd = TRUE)
    # totals on the pooled row
    text(ilabx[1], -1, sum(dat$Cases),      pos = 2, cex = CEX, font = 2)
    text(ilabx[2], -1, sum(dat$SampleSize), pos = 2, cex = CEX, font = 2)
    text(ilabx[3], -1, "100.0",             pos = 2, cex = CEX, font = 2)
    # prediction interval
    cap <- 0.18
    segments(tf(pil), -2, tf(pih), -2, col = "red", lwd = 2.5 * CEX)
    segments(tf(pil), -2 - cap, tf(pil), -2 + cap, col = "red", lwd = 1.5 * CEX)
    segments(tf(pih), -2 - cap, tf(pih), -2 + cap, col = "red", lwd = 1.5 * CEX)
    text(xlim[1], -2, "95% prediction interval", pos = 4, cex = CEX)
    text(xlim[2], -2, sprintf("[%.2f \u2013 %.2f]", tf(pil), tf(pih)), pos = 2, cex = CEX)
    # axis label and heterogeneity below the axis
    text(amax / 2, -5.0, paste0("RVF ", assay, " prevalence (%) \u2014 ", sp), cex = CEX, font = 2)
    text(xlim[1], -6.0, het, pos = 4, cex = CEX)
    # legend: square = study estimate (size = weight), diamond = pooled, red bar = PI
    ly <- -7.0; x <- xlim[1]; gap <- 0.03 * diff(xlim); cx <- 0.9 * CEX
    hh <- 0.25
    sq <- 0.012 * diff(xlim)
    rect(x, ly - hh, x + 2 * sq, ly + hh, col = "black", border = NA)
    segments(x - sq, ly, x + 3 * sq, ly, lwd = 1)
    x <- x + 3 * sq + gap / 3
    lab1 <- "Study estimate with 95% CI"
    text(x, ly, lab1, pos = 4, cex = cx, offset = 0); x <- x + strwidth(lab1, cex = cx) + gap
    dw <- 0.035 * diff(xlim)
    polygon(c(x, x + dw / 2, x + dw, x + dw / 2), c(ly, ly + hh, ly, ly - hh), col = "maroon", border = "black")
    x <- x + dw + gap / 3
    lab2 <- "Pooled estimate with 95% CI"
    text(x, ly, lab2, pos = 4, cex = cx, offset = 0); x <- x + strwidth(lab2, cex = cx) + gap
    segments(x, ly, x + dw, ly, col = "red", lwd = 2.5 * CEX)
    x <- x + dw + gap / 3
    text(x, ly, "95% prediction interval", pos = 4, cex = cx, offset = 0)
    # footnote defining the abbreviations
    text(xlim[1], -8.0, paste0("N = number sampled;   Wt % = ", wt_def),
         pos = 4, cex = CEX, font = 3)
    par(op)
  }, error = function(e) message("forest failed ", sp, "/", assay, ": ", conditionMessage(e)))
  par(op0); dev.off()
}

## ---- steps 3-6: analysis per species x assay --------------------------------
analyse <- function(dat, sp, assay) {
  k <- nrow(dat); ns <- dplyr::n_distinct(dat$StudyID)
  if (k < 2) { message("skip ", sp, "/", assay, " (k<2)"); return(NULL) }
  dat$esid <- seq_len(k)
  repl <- any(table(dat$StudyID) > 1)
  
  # step 3-4: primary model (logit)
  es  <- escalc(measure = "PLO", xi = Cases, ni = SampleSize, data = dat, add = 1/2, to = "only0")
  res <- if (repl) rma.mv(yi, vi, random = ~1 | StudyID/esid, data = es, method = "REML")
  else       rma(yi, vi, data = es, method = "REML")
  pr  <- predict(res, transf = tf)
  
  tau_b <- if (repl) res$sigma2[1] else res$tau2
  tau_w <- if (repl) res$sigma2[2] else NA_real_
  I2    <- if (repl) ml_I2(res, es$vi) else res$I2
  
  # step 5: robust variance estimation CI
  rob <- tryCatch(robust(res, cluster = es$StudyID, clubSandwich = TRUE), error = function(e) NULL)
  rve_lo <- if (!is.null(rob)) 100 * plogis(rob$ci.lb) else NA_real_
  rve_hi <- if (!is.null(rob)) 100 * plogis(rob$ci.ub) else NA_real_
  
  # step 5: logit GLMM
  gpr <- tryCatch({ g <- rma.glmm(measure = "PLO", xi = Cases, ni = SampleSize, data = dat)
  predict(g, transf = tf) }, error = function(e) NULL)
  
  # step 5: Freeman-Tukey double arcsine
  esft  <- escalc(measure = "PFT", xi = Cases, ni = SampleSize, data = dat)
  resft <- if (repl) rma.mv(yi, vi, random = ~1 | StudyID/esid, data = esft, method = "REML")
  else       rma(yi, vi, data = esft, method = "REML")
  prft  <- predict(resft, transf = transf.ipft.hm, targ = list(ni = dat$SampleSize))
  
  # step 6: forest plot
  draw_forest(res, dat, es, repl, I2, tau_b, tau_w, sp, assay,
              file.path(out_dir, paste0("forest_", gsub(" ", "_", sp), "_", assay, ".png")))
  
  tibble(species = sp, assay = assay, model = if (repl) "3-level logit" else "RE logit",
         k_estimates = k, n_studies = ns, total_cases = sum(dat$Cases), total_n = sum(dat$SampleSize),
         prev_pct = pr$pred, ci_lo = pr$ci.lb, ci_hi = pr$ci.ub, pi_lo = pr$pi.lb, pi_hi = pr$pi.ub,
         rve_ci_lo = rve_lo, rve_ci_hi = rve_hi,
         tau2_between = tau_b, tau2_within = tau_w, I2_pct = I2,
         glmm_prev = if (!is.null(gpr)) gpr$pred else NA_real_,
         glmm_lo   = if (!is.null(gpr)) gpr$ci.lb else NA_real_,
         glmm_hi   = if (!is.null(gpr)) gpr$ci.ub else NA_real_,
         ft_prev = prft$pred, ft_lo = prft$ci.lb, ft_hi = prft$ci.ub,
         ft_pi_lo = prft$pi.lb, ft_pi_hi = prft$pi.ub)
}

## ---- run --------------------------------------------------------------------
results <- list()
for (sp in c("cattle", "small ruminants"))
  for (assay in ASSAYS)
    results[[paste(sp, assay)]] <- analyse(d0 %>% filter(Species == sp, Assay == assay), sp, assay)

results_tbl <- dplyr::bind_rows(results)
print(as.data.frame(results_tbl))
write_xlsx(results_tbl, file.path(out_dir, "R3_pooled_results.xlsx"))

cat("\nDone. Forest plots and results table written to:\n", out_dir, "\n")