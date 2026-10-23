# ============================================================================
# RVF prevalence: symptom subgroup analysis
# ----------------------------------------------------------------------------
# This script reproduces the symptom subgroup analyses of RVF prevalence in
# cattle and small ruminants, separately for each assay (RT-PCR, IgM, IgG),
# comparing symptomatic and asymptomatic animals. For each species x assay
# combination it:
#
#   1. Loads the dataset and keeps records with symptom status ("yes"/"no")
#      and usable cases and sample sizes.
#   2. Pools prevalence on the logit scale, using a three-level model
#      (estimates nested within studies) where a study contributes more than
#      one estimate, and a random-effects model otherwise.
#   3. Estimates pooled prevalence, 95% CI, 95% prediction interval and
#      heterogeneity (Q, I2, tau2) overall and for each symptom group with two
#      or more estimates.
#   4. Tests for a difference between symptomatic and asymptomatic animals
#      using symptom status as a moderator. The test is repeated with robust
#      variance estimation (CR2), a logit GLMM and the Freeman-Tukey
#      transformation.
#   5. Draws a grouped forest plot and writes summary tables.
#
# The overall pooled estimate combines both symptom groups; studies reporting
# both groups contribute to each.
#
# Requirements
#   R packages: metafor, clubSandwich, dplyr, readxl, writexl
#   install.packages(c("metafor", "clubSandwich", "dplyr", "readxl", "writexl"))
#
# Input
#   IN_PATH: symptoms.xlsx
#   Required columns: Author, Species, Symptoms ("yes"/"no"), Test, Cases,
#   SampleSize
#
# Output (folder "symptom_R3_out" next to the input file)
#   symptom_<species>_<assay>.png    grouped forest plots
#   symptom_level_estimates.xlsx     pooled estimates per symptom group
#   symptom_between_tests.xlsx       between-group tests and overall estimates
#
# Usage
#   Set IN_PATH to the local file location, then run the whole script.
#   FONT_PT sets the font size; WIDTH_SCALE and WIDTH_CAP set the plot width;
#   BOXMIN and BOXMAX set the range of study box sizes.
# ============================================================================

## ---- paths and settings -----------------------------------------------------
IN_PATH <- "C:/Users/eha/meta/symptoms.xlsx"
OUT_DIR <- file.path(dirname(IN_PATH), "symptom_R3_output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
SYM_LAB     <- c("yes" = "Symptomatic", "no" = "Asymptomatic")
disp        <- function(x) unname(SYM_LAB[x])
OVERALL_LAB <- "Pooled"
FOOTNOTE    <- "N = number sampled;  Wt % = inverse-variance weight.  'Pooled' mixes both symptom groups; both-arm studies counted twice."
CEX         <- 1.3     # relative text size for plot elements
FONT_PT     <- 24      # base font size (points)
WIDTH_SCALE <- 1.55    # canvas width multiplier
WIDTH_CAP   <- 18      # maximum canvas width (inches)
BOXMIN <- 0.6
BOXMAX <- 1.5

library(metafor); library(clubSandwich)
library(dplyr);   library(readxl); library(writexl)

## ---- helper functions -------------------------------------------------------
DASH <- c(" [", " \u2013 ", "]")
tf   <- function(x) 100 * plogis(x)
pfmt <- function(p) ifelse(is.na(p), "= NA", ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)))
# returns a scalar, or a fallback when the value is missing
sc   <- function(x, fb = NA_real_) if (is.null(x) || length(x) == 0 || !is.finite(x[1])) fb else x[1]

# multilevel I2 (Viechtbauer) for rma.mv
ml_I2 <- function(res, vi) tryCatch({
  W <- diag(1 / vi); X <- model.matrix(res)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  typ <- (res$k - res$p) / sum(diag(P))
  100 * sum(res$sigma2) / (sum(res$sigma2) + typ)
}, error = function(e) NA_real_)

# row layout per group: header, studies, subtotal, heterogeneity, PI, gap
layout_rows <- function(ks, has_sub) {
  L <- length(ks); extra <- ifelse(has_sub, 3L, 0L)
  per <- 1 + ks + extra + 1
  content <- sum(per); top <- content + 1; cur <- content
  studyrows <- integer(0); head <- numeric(L)
  summ <- rep(NA_real_, L); hetrow <- rep(NA_real_, L); pirow <- rep(NA_real_, L)
  for (i in seq_len(L)) {
    head[i] <- cur; cur <- cur - 1
    studyrows <- c(studyrows, seq(cur, cur - ks[i] + 1)); cur <- cur - ks[i]
    if (has_sub[i]) { summ[i] <- cur; cur <- cur - 1
    hetrow[i] <- cur; cur <- cur - 1
    pirow[i] <- cur; cur <- cur - 1 }
    cur <- cur - 1
  }
  list(top = top, studyrows = studyrows, head = head, summ = summ, hetrow = hetrow, pirow = pirow)
}

# Q, df and p from an rma or rma.mv fit
qstats <- function(res) list(
  QE  = sc(res$QE), QEp = sc(res$QEp),
  QEdf = { v <- sc(res$QEdf, NA_real_); if (is.na(v)) sc(res$k - res$p) else v })

## ---- steps 2-3: pooling function --------------------------------------------
# three-level model if a study gives >1 estimate, random-effects model
# otherwise; prediction interval on the t distribution (df = k - p)
fit_level <- function(sub) {
  sub <- sub %>% filter(!is.na(Cases), !is.na(SampleSize), SampleSize > 0)
  k <- nrow(sub); if (k == 0) return(NULL)
  sub$esid <- seq_len(k)
  es <- escalc("PLO", xi = Cases, ni = SampleSize, data = sub, add = 1/2, to = "only0")
  if (k == 1) return(list(mu = es$yi[1], ci.lb = es$yi[1] - 1.96*sqrt(es$vi[1]),
                          ci.ub = es$yi[1] + 1.96*sqrt(es$vi[1]), pi.lb = NA, pi.ub = NA,
                          I2 = NA, tau2 = NA, QE = NA, QEdf = NA, QEp = NA,
                          k = 1, ns = 1, cases = sum(sub$Cases), n = sum(sub$SampleSize)))
  repl <- any(table(sub$StudyID) > 1)
  res <- tryCatch(if (repl) rma.mv(yi, vi, random = ~1 | StudyID/esid, data = es, method = "REML")
                  else       rma(yi, vi, data = es, method = "REML"), error = function(e) NULL)
  if (is.null(res)) return(list(mu = NA, ci.lb = NA, ci.ub = NA, pi.lb = NA, pi.ub = NA,
                                I2 = NA, tau2 = NA, QE = NA, QEdf = NA, QEp = NA,
                                k = k, ns = dplyr::n_distinct(sub$StudyID),
                                cases = sum(sub$Cases), n = sum(sub$SampleSize)))
  mu <- as.numeric(coef(res)); se <- as.numeric(res$se)
  t2 <- if (repl) sum(res$sigma2) else res$tau2
  tcr <- qt(0.975, df = max(1, res$k - res$p)); q <- qstats(res)
  list(mu = mu, ci.lb = res$ci.lb, ci.ub = res$ci.ub,
       pi.lb = mu - tcr*sqrt(se^2 + t2), pi.ub = mu + tcr*sqrt(se^2 + t2),
       I2 = if (repl) ml_I2(res, es$vi) else res$I2, tau2 = t2,
       QE = q$QE, QEdf = q$QEdf, QEp = q$QEp,
       k = k, ns = dplyr::n_distinct(sub$StudyID), cases = sum(sub$Cases), n = sum(sub$SampleSize))
}

het_str <- function(prefix, QE, QEdf, QEp, I2, tau2)
  sprintf("%s Q = %s, df = %s, p %s;  I\u00b2 = %s, \u03c4\u00b2 = %s", prefix,
          ifelse(is.na(QE),  "NA", sprintf("%.1f", QE)),
          ifelse(is.na(QEdf),"NA", as.character(round(QEdf))), pfmt(QEp),
          ifelse(is.na(I2),  "NA", paste0(round(I2, 1), "%")),
          ifelse(is.na(tau2),"NA", sprintf("%.3f", tau2)))

## ---- step 5: forest plot function -------------------------------------------
draw_forest <- function(es, d, lev_order, ks, has_sub, res_all, per_level, overall,
                        qm_line, assay, sp, png_path) {
  lay <- layout_rows(ks, has_sub); k <- nrow(d)
  # box size scaled to study precision
  prec  <- 1 / sqrt(es$vi); rng <- max(prec) - min(prec)
  psz   <- if (rng > 0) BOXMIN + (BOXMAX - BOXMIN) * (prec - min(prec)) / rng else rep((BOXMIN+BOXMAX)/2, k)
  wt_lab <- sprintf("%.1f", 100 * (1/es$vi) / sum(1/es$vi))
  
  # x-axis range set from the largest upper CI / PI limit
  ub  <- tf(es$yi + 1.96 * sqrt(es$vi))
  pis <- tf(c(overall$pi.ub, sapply(per_level, function(p) if (is.null(p)) NA else p$pi.ub)))
  amax <- max(c(ub, pis), na.rm = TRUE)
  step <- if (amax <= 20) 5 else if (amax <= 50) 10 else 20
  amax <- min(100, ceiling(amax / step) * step); amax <- max(amax, 10)
  at <- seq(0, amax, by = step); W <- amax
  xlim  <- c(-2.70 * W, 1.85 * W); ilabx <- c(-1.05 * W, -0.72 * W, -0.38 * W)
  
  # canvas size; row height grows with the square root of FONT_PT
  width_in  <- min(WIDTH_CAP, WIDTH_SCALE * max(11, 4 + 7.0 * CEX))
  height_in <- max(8, 1 + 0.30 * CEX * sqrt(FONT_PT / 12) * (lay$top + 13))
  R_HET <- -3.3; R_BTW1 <- -4.1; R_BTW2 <- -4.9; R_FOOT <- -5.7; Y_AXIS <- -7.2
  # shrink long summary lines to fit the plot width
  avail  <- (xlim[2] - xlim[1]) * 0.99
  fitcex <- function(txt, base) tryCatch({ w <- strwidth(txt, cex = base)
  if (is.finite(w) && w > avail) base * avail / w else base }, error = function(e) base)
  
  png(png_path, width = width_in, height = height_in, units = "in", res = 450, pointsize = FONT_PT)
  op0 <- par(mar = c(6, 3, 3, 2))
  tryCatch({
    forest(res_all, transf = tf, olim = c(0, 100), refline = NA, rows = lay$studyrows,
           slab = d$Author, psize = psz,
           ilab = data.frame(Cases = as.character(d$Cases), N = as.character(d$SampleSize), Weight = wt_lab),
           ilab.xpos = ilabx, ilab.lab = c("Cases", "N", "Wt %"), ilab.pos = 2,
           header = c("Study", "Prevalence (%) [95% CI]"),
           xlim = xlim, alim = c(0, amax), at = at, ylim = c(Y_AXIS, lay$top + 2),
           digits = c(2L, 0L), xlab = paste0("RVF ", assay, " prevalence (%) \u2014 ", sp),
           col = "maroon", cex = CEX, annosym = DASH, mlab = OVERALL_LAB)
    
    op <- par(xpd = TRUE)
    # group headers, subtotals, heterogeneity and prediction intervals
    for (i in seq_along(lev_order)) tryCatch({
      pl <- per_level[[i]]
      text(xlim[1], lay$head[i], disp(lev_order[i]), pos = 4, cex = CEX, font = 2)
      if (has_sub[i] && !is.null(pl) && !is.na(pl$mu)) {
        addpoly(pl$mu, ci.lb = pl$ci.lb, ci.ub = pl$ci.ub, rows = lay$summ[i], transf = tf,
                col = "navy", annosym = DASH, cex = CEX, mlab = sprintf("Subtotal (k=%d)", pl$k))
        text(ilabx[1], lay$summ[i], pl$cases, pos = 2, cex = CEX, font = 2)
        text(ilabx[2], lay$summ[i], pl$n,     pos = 2, cex = CEX, font = 2)
        text(xlim[1], lay$hetrow[i], het_str("Heterogeneity:", pl$QE, pl$QEdf, pl$QEp, pl$I2, pl$tau2),
             pos = 4, cex = CEX * 0.85)
        if (!is.na(pl$pi.lb)) {
          cap <- 0.15
          segments(tf(pl$pi.lb), lay$pirow[i], tf(pl$pi.ub), lay$pirow[i], col = "red", lwd = 1.8 * CEX)
          segments(tf(pl$pi.lb), lay$pirow[i]-cap, tf(pl$pi.lb), lay$pirow[i]+cap, col = "red", lwd = 1.2*CEX)
          segments(tf(pl$pi.ub), lay$pirow[i]-cap, tf(pl$pi.ub), lay$pirow[i]+cap, col = "red", lwd = 1.2*CEX)
          text(xlim[1], lay$pirow[i], "95% prediction interval", pos = 4, cex = CEX * 0.9)
          text(xlim[2], lay$pirow[i], sprintf("[%.2f \u2013 %.2f]", tf(pl$pi.lb), tf(pl$pi.ub)),
               pos = 2, cex = CEX * 0.9)
        }
      }
    }, error = function(e) message("subgroup draw failed ", lev_order[i], ": ", conditionMessage(e)))
    
    # overall totals, overall PI, heterogeneity, between-group test, footnote
    tryCatch({
      text(ilabx[1], -1, sum(d$Cases),      pos = 2, cex = CEX, font = 2)
      text(ilabx[2], -1, sum(d$SampleSize), pos = 2, cex = CEX, font = 2)
      text(ilabx[3], -1, "100.0",           pos = 2, cex = CEX, font = 2)
      cap <- 0.18
      segments(tf(overall$pi.lb), -2, tf(overall$pi.ub), -2, col = "red", lwd = 2.5 * CEX)
      segments(tf(overall$pi.lb), -2-cap, tf(overall$pi.lb), -2+cap, col = "red", lwd = 1.5 * CEX)
      segments(tf(overall$pi.ub), -2-cap, tf(overall$pi.ub), -2+cap, col = "red", lwd = 1.5 * CEX)
      text(xlim[1], -2, "95% prediction interval (overall)", pos = 4, cex = CEX)
      text(xlim[2], -2, sprintf("[%.2f \u2013 %.2f]", tf(overall$pi.lb), tf(overall$pi.ub)), pos = 2, cex = CEX)
      oh_line <- het_str("Overall heterogeneity:", overall$QE, overall$QEdf, overall$QEp, overall$I2, overall$tau2)
      text(xlim[1], R_HET,  oh_line,    pos = 4, cex = fitcex(oh_line, CEX))
      text(xlim[1], R_BTW1, qm_line[1], pos = 4, cex = fitcex(qm_line[1], CEX))
      text(xlim[1], R_BTW2, qm_line[2], pos = 4, cex = fitcex(qm_line[2], CEX))
      text(xlim[1], R_FOOT, FOOTNOTE,   pos = 4, cex = fitcex(FOOTNOTE, CEX * 0.85), font = 3)
    }, error = function(e) message("overall block failed ", assay, "/", sp, ": ", conditionMessage(e)))
    par(op)
  }, error = function(e) message("forest failed ", assay, "/", sp, ": ", conditionMessage(e)))
  par(op0); dev.off()
}

## ---- step 1: load data ------------------------------------------------------
d0 <- read_excel(IN_PATH); names(d0) <- trimws(names(d0))
d0 <- d0 %>%
  mutate(Species = trimws(tolower(as.character(Species))),
         Symptoms = trimws(tolower(as.character(Symptoms))),
         Cases = suppressWarnings(as.numeric(Cases)), SampleSize = suppressWarnings(as.numeric(SampleSize)),
         Author = trimws(as.character(Author)), Assay = trimws(as.character(Test))) %>%
  filter(Species %in% c("cattle","small ruminants")) %>%
  filter(Symptoms %in% c("yes","no")) %>%
  filter(!is.na(Cases), !is.na(SampleSize), SampleSize > 0)
d0$StudyID <- d0$Author
level_rows <- list(); between_rows <- list()

## ---- steps 2-5: analysis per species x assay --------------------------------
analyse <- function(d, assay, sp) {
  d <- d %>% filter(Species == sp, Assay == assay)
  if (nrow(d) < 2) { message("skip ", assay, "/", sp, " (k<2)"); return(invisible()) }
  lev_order <- c("yes","no"); lev_order <- lev_order[lev_order %in% d$Symptoms]
  d <- d[order(match(as.character(d$Symptoms), lev_order)), ]; d$esid <- seq_len(nrow(d))
  es <- escalc("PLO", xi = Cases, ni = SampleSize, data = d, add = 1/2, to = "only0")
  es$grp <- factor(as.character(d$Symptoms), levels = lev_order)
  ks <- as.integer(table(es$grp)); repl <- any(table(d$StudyID) > 1)
  
  # overall pooled estimate
  res_all <- if (repl) rma.mv(yi, vi, random = ~1 | StudyID/esid, data = es, method = "REML")
  else       rma(yi, vi, data = es, method = "REML")
  mu <- as.numeric(coef(res_all)); se <- as.numeric(res_all$se)
  t2a <- if (repl) sum(res_all$sigma2) else res_all$tau2
  tcr <- qt(0.975, df = max(1, res_all$k - res_all$p)); oq <- qstats(res_all)
  overall <- list(mu = mu, ci.lb = res_all$ci.lb, ci.ub = res_all$ci.ub,
                  pi.lb = mu - tcr*sqrt(se^2 + t2a), pi.ub = mu + tcr*sqrt(se^2 + t2a),
                  I2 = if (repl) ml_I2(res_all, es$vi) else res_all$I2, tau2 = t2a,
                  QE = oq$QE, QEdf = oq$QEdf, QEp = oq$QEp)
  
  # per-group estimates; subtotals drawn only for groups with k >= 2
  per_level <- lapply(lev_order, function(L) fit_level(d[as.character(d$Symptoms) == L, ]))
  has_sub <- vapply(per_level, function(p) !is.null(p) && p$k >= 2, logical(1))
  # number of studies reporting both symptom groups
  both <- sum(tapply(d$Symptoms, d$StudyID, function(z) all(c("yes","no") %in% z)))
  
  # step 4: between-group test
  QM <- QMdf <- QMp <- rveP <- glmmP <- ftP <- NA
  if (length(lev_order) >= 2) {
    res_mod <- tryCatch(if (repl) rma.mv(yi, vi, mods = ~grp, random = ~1 | StudyID/esid, data = es)
                        else       rma(yi, vi, mods = ~grp, data = es), error = function(e) NULL)
    if (!is.null(res_mod)) { QM <- res_mod$QM; QMdf <- res_mod$QMdf[1]; QMp <- res_mod$QMp
    rveP <- tryCatch(robust(res_mod, cluster = es$StudyID, clubSandwich = TRUE)$QMp, error = function(e) NA) }
    glmmP <- tryCatch(rma.glmm(measure = "PLO", xi = Cases, ni = SampleSize, mods = ~es$grp, data = d)$QMp, error = function(e) NA)
    esft <- escalc("PFT", xi = Cases, ni = SampleSize, data = d)
    ftP  <- tryCatch((if (repl) rma.mv(yi, vi, mods = ~es$grp, random = ~1 | StudyID/esid, data = esft)
                      else       rma(yi, vi, mods = ~es$grp, data = esft))$QMp, error = function(e) NA)
  }
  qm_line <- c(
    sprintf("Symptomatic vs asymptomatic (%s): Q = %s, df = %s, p %s", ifelse(repl, "3-level", "RE"),
            ifelse(is.na(QM), "NA", sprintf("%.1f", QM)), ifelse(is.na(QMdf), "NA", QMdf), pfmt(QMp)),
    sprintf("Robust-variance p %s   |   GLMM p %s   |   Freeman-Tukey p %s",
            pfmt(rveP), pfmt(glmmP), pfmt(ftP)))
  
  # step 5: forest plot
  draw_forest(es, d, lev_order, ks, has_sub, res_all, per_level, overall, qm_line,
              assay, sp, file.path(OUT_DIR, paste0("symptom_", gsub(" ", "_", sp), "_", assay, ".png")))
  
  # summary tables
  for (i in seq_along(lev_order)) { pl <- per_level[[i]]; if (is.null(pl)) next
  level_rows[[length(level_rows)+1]] <<- tibble(species=sp, assay=assay, symptom=disp(lev_order[i]), symptom_raw=lev_order[i],
                                                k=pl$k, n_studies=pl$ns, cases=pl$cases, n=pl$n, prev_pct=tf(pl$mu), ci_lo=tf(pl$ci.lb), ci_hi=tf(pl$ci.ub),
                                                pi_lo=if(is.na(pl$pi.lb))NA else tf(pl$pi.lb), pi_hi=if(is.na(pl$pi.ub))NA else tf(pl$pi.ub),
                                                QE=pl$QE, QE_df=pl$QEdf, QE_p=pl$QEp, I2_pct=pl$I2, tau2=pl$tau2) }
  between_rows[[length(between_rows)+1]] <<- tibble(species=sp, assay=assay, k_total=nrow(d), n_studies_both_arms=both,
                                                    overall_prev_pct=tf(overall$mu), overall_ci_lo=tf(overall$ci.lb), overall_ci_hi=tf(overall$ci.ub),
                                                    overall_pi_lo=tf(overall$pi.lb), overall_pi_hi=tf(overall$pi.ub),
                                                    overall_QE=overall$QE, overall_QE_df=overall$QEdf, overall_QE_p=overall$QEp,
                                                    QM=QM, QM_df=QMdf, QM_p=QMp, RVE_p=rveP, GLMM_p=glmmP, FT_p=ftP)
  invisible()
}

## ---- run --------------------------------------------------------------------
for (sp in c("cattle","small ruminants")) for (assay in c("RT-PCR","IgM","IgG")) analyse(d0, assay, sp)
if (length(level_rows))   write_xlsx(dplyr::bind_rows(level_rows),   file.path(OUT_DIR, "symptom_level_estimates.xlsx"))
if (length(between_rows)) write_xlsx(dplyr::bind_rows(between_rows), file.path(OUT_DIR, "symptom_between_tests.xlsx"))
cat("\nDone. Symptom forest plots and tables written to:\n", OUT_DIR, "\n")