# ============================================================================
# RVF prevalence: publication-bias analyses (supplementary, exploratory)
# ----------------------------------------------------------------------------
# This script reproduces the supplementary publication-bias analyses for the
# pooled RVF prevalence estimates in cattle and small ruminants, separately
# for each assay: RT-PCR, IgM and IgG.
# It uses the same filters, labels and six species x assay strata as the main 
# pooling script (rvf_pooled_prevalence.R). For each stratum it:
#
#   1. Loads the dataset, keeps cattle and small ruminant records with usable
#      cases and sample sizes, and maps the Test column to assay.
#   2. Fits a  random-effects model on the logit scale.
#   3. Runs Egger's regression test.
#   4. Draws a funnel plot with a 95% pseudo-confidence region and Egger's
#      test result below the plot.
#   5. Writes a results table.
#
# Strata with fewer than three estimates are not tested.
#
# Note: Egger's test uses one estimate per row
# and cannot account for the clustering of estimates within studies used in
# the main three-level model. They are unreliable when k is small and I2 is
# high. Results are exploratory and do not alter the pooled estimates.
#
# Requirements
#   R packages: metafor, dplyr, readxl, writexl
#   install.packages(c("metafor", "dplyr", "readxl", "writexl"))
#
# Input
#   in_path: meta RVF overall.xlsx (same file as the main pooling script)
#   Required columns: Author, Species, Test, Cases, SampleSize
#   Test values: "RT-PCR", "IgM", "IgG"
#
# Output (folder "pubbias_supplement_v4" next to the input file)
#   funnel_<species>_<assay>.png    funnel plots
#   S_publication_bias.xlsx         Egger results
#
# Usage
#   Set in_path to the local file location, then run the whole script.
#   CEX controls font size.
# ============================================================================

## ---- settings ---------------------------------------------------------------
CEX <- 1.1
ASSAYS <- c("RT-PCR", "IgM", "IgG")   # expected values of the Test column

library(metafor)
library(dplyr)
library(readxl)
library(writexl)

## ---- paths ------------------------------------------------------------------
in_path <- "C:/Users/eh/meta/meta RVF overall.xlsx"
out_dir <- file.path(dirname(in_path), "pubbias_supplement_v4")
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

## ---- helper functions -------------------------------------------------------
tf <- function(x) 100 * plogis(x)   # logit to prevalence (%)

# x-axis ticks at fixed prevalences (placed on the logit scale)
TICKS_PCT <- c(0.1, 1, 5, 20, 50, 80)

## ---- steps 2-5: analysis per species x assay --------------------------------
pubbias <- function(dat, sp, assay) {
  k <- nrow(dat)
  if (k < 3) { message("skip ", sp, "/", assay, " (k<3; bias tests not run)")
    return(tibble(species = sp, assay = assay, k = k, note = "k<3: not run")) }
  
  # step 2: univariate random-effects model (logit proportion)
  es  <- escalc(measure = "PLO", xi = Cases, ni = SampleSize, data = dat, add = 1/2, to = "only0")
  res <- rma(yi, vi, data = es, method = "REML")
  
  # step 3: Egger, Begg, trim-and-fill
  egg <- tryCatch(regtest(res, model = "lm"),  error = function(e) NULL)
  bgg <- tryCatch(ranktest(res),               error = function(e) NULL)
  tf0 <- tryCatch(trimfill(res),               error = function(e) NULL)
  
  # step 4: funnel plot
  zero <- dat$Cases == 0 | dat$Cases == dat$SampleSize
  # x-range covers every study, the 95% pseudo-CI at the largest SE, and the ticks
  sei_max <- max(sqrt(es$vi))
  mu      <- as.numeric(coef(res))
  xr      <- range(c(es$yi, mu - 1.96 * sei_max, mu + 1.96 * sei_max))
  at      <- qlogis(TICKS_PCT / 100)
  at      <- at[at >= xr[1] - 0.5 & at <= xr[2] + 0.5]
  xlim    <- range(c(xr, at)) + c(-0.3, 0.3)
  
  png(file.path(out_dir, paste0("funnel_", gsub(" ", "_", sp), "_", assay, ".png")),
      width = 6.5, height = 6.0, units = "in", res = 300)
  par(cex = CEX, mar = c(5.6, 4.5, 3.5, 1.5))
  tryCatch({
    funnel(res, atransf = tf, at = at, xlim = xlim, xlab = "Prevalence (%)",
           level = 95, legend = FALSE, refline = mu,
           pch = 19, cex = 0.5,
           main = paste0(sp, " \u2014 ", assay))
    egg_txt <- if (is.null(egg)) sprintf("k = %d; Egger's test: not estimable", k) else
      sprintf("k = %d; Egger's test: t = %.2f (df = %d), p %s", k, egg$zval, as.integer(egg$dfs),
              ifelse(egg$pval < 0.001, "< 0.001", sprintf("= %.3f", egg$pval)))
    mtext(egg_txt, side = 1, line = 4.3, cex = 0.85 * CEX, adj = 0.5)
  }, error = function(e) { message("funnel failed ", sp, "/", assay, ": ", conditionMessage(e)); plot.new() })
  dev.off()
  
  # step 5: results row
  obs_pct <- tf(mu)
  adj_pct <- if (!is.null(tf0)) tf(as.numeric(coef(tf0))) else NA_real_
  
  tibble(
    species = sp, assay = assay, k = k,
    I2_pct  = res$I2, tau2 = res$tau2,
    egger_t     = if (!is.null(egg)) egg$zval else NA_real_,
    egger_p     = if (!is.null(egg)) egg$pval else NA_real_,
    begg_p      = if (!is.null(bgg)) bgg$pval else NA_real_,
    n_imputed   = if (!is.null(tf0)) tf0$k0    else NA_real_,
    fill_side   = if (!is.null(tf0)) tf0$side  else NA_character_,
    prev_observed_pct = obs_pct,
    prev_trimfill_pct = adj_pct,
    n_zero_event = sum(zero),
    note = "exploratory; univariate, ignores 3-level clustering"
  )
}

## ---- run --------------------------------------------------------------------
out <- list()
for (sp in c("cattle", "small ruminants"))
  for (assay in ASSAYS)
    out[[paste(sp, assay)]] <- pubbias(d0 %>% filter(Species == sp, Assay == assay), sp, assay)

tbl <- dplyr::bind_rows(out)
print(as.data.frame(tbl))
write_xlsx(tbl, file.path(out_dir, "S_publication_bias.xlsx"))

cat("\nDone. Funnel plots and results table written to:\n", out_dir, "\n")