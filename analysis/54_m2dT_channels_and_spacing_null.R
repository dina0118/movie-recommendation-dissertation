# =====================================================================
# 54_m2dT_channels_and_spacing_null.R
#
# Four follow-ups to 52, all needing the same session and the same
# three probability matrices, so they are done together:
#
#   (1) two-way bootstrap intervals for calibration error and
#       sharpness across M2d -> M2d+T -> M2e            (52 section 4)
#   (2) whether M2d+T is a third instance of a likelihood gain that
#       is not a calibration gain                        (52 section 4)
#   (3) a null for the cut-point spacing agreement       (52 section 5)
#   (4) intervals for the threshold gain by training volume, and a
#       direct test of the 100+ fall-back                (52 section 5)
#
# RUN INSIDE THE SAME CLEAN 06b SESSION, AFTER 52. It reads the cached
# m2dT_frozen_factors.rds rather than retraining. Do not source any
# STARTUP file into this session.
#
# ---------------------------------------------------------------------
# WHY (1) AND (2) MATTER MORE THAN THE LIKELIHOOD DECOMPOSITION
#
# 52 found that 78.4% of the M2d -> M2e likelihood increment comes from
# the cut points. On the thesis's own primary criterion the picture
# looks different: of the 0.01380 improvement in calibration error,
# the cut points supply roughly 0.00389 and the refitted factors
# roughly 0.00991. The likelihood decomposition and the calibration
# decomposition point at different components of the same model.
#
# Neither of those calibration figures has an interval. Section 1
# supplies them. If the intervals overlap heavily, the reversal is not
# established and neither sentence should be written.
#
# Section 2 asks the sharper question. M2d+T's mean predictive standard
# deviation is 0.79148, against M2d's 0.83435 and M2e's 0.82500 -- it
# is markedly the sharpest of the three, and not monotone in the
# ladder. A model that gains likelihood while becoming sharper may be
# gaining on confidence rather than on location. That is the pattern
# already documented for M1 and for M3b(30), and the thesis's central
# methodological claim rests on that pairing. A third instance arising
# inside the main result would strengthen the claim and complicate the
# main result at the same time. It is tested here rather than asserted.
#
# ---------------------------------------------------------------------
# WHY (3) IS A NULL AND NOT A CONFIRMATION
#
# 52 reported Spearman correlations of 0.85, 0.90 and 0.96 between
# M2d+T's cut-point spacings and M2e's, and read them as two fits
# recovering the same structure from different factor matrices.
#
# That reading has no baseline. Both models were fitted to the same
# people's ratings, and a user's cut-point spacings are strongly
# determined by their own rating histogram: someone who rates almost
# everything 4 or 5 will have wide low spacings under any fit. So the
# question is not whether the two agree, but whether they agree by more
# than two independent summaries of the same histogram would.
#
# The null here is the spacing implied directly by the user's empirical
# cumulative rating distribution, with no model involved. If that null
# already correlates at 0.9 with both fits, the agreement in 52 is a
# property of the data and section 6d gains nothing from it.
#
# The null is not defined for users who never use some level, because
# an empirical cumulative probability of 0 or 1 maps to an infinite
# threshold. Only 35.4% of users use all five levels. Both a strict
# version (complete users only, no smoothing) and a smoothed version
# (all users, add-one) are reported, because a result that holds only
# under one of them is a result about the smoothing.
# =====================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

.need <- c("m2e", "dat_fit", "dat_calib", "dat_test", "test", "fit", "calib",
           "ordrec_probs_m3", "tidy_probs", "user_map", "movie_map",
           "LEVELS5", "K_LEVELS", "OUT_DIR", "CACHE_M3", "SEED")
.miss <- .need[!vapply(.need, exists, logical(1))]
if (length(.miss)) stop("Run inside the clean 06b session after 52. Missing: ",
                        paste(.miss, collapse = ", "))
if (!identical(as.integer(test$user_idx), as.integer(dat_test$u)))
  stop("test and dat_test row orders differ. Restart R and re-run 06b.")

N_BOOT_LOCAL <- 2000L
g_nU <- nrow(user_map); g_nI <- nrow(movie_map)
g_y  <- as.integer(test$rating5); g_n <- length(g_y)

g_pm2e <- m2e$par; g_pm2e$w <- NULL; g_pm2e$W <- NULL
g_CACHE7 <- file.path(dirname(normalizePath(CACHE_M3)), "07_M1_M2")
g_m2d  <- readRDS(file.path(g_CACHE7, "m2d_final.rds"))
g_mod  <- readRDS(file.path(CACHE_M3, "m2dT_frozen_factors.rds"))

g_Tmat <- function(par, u) {
  E <- exp(par$beta[u, , drop = FALSE]); t1 <- par$t1[u]
  cbind(t1, t1 + E[, 1], t1 + E[, 1] + E[, 2], t1 + E[, 1] + E[, 2] + E[, 3])
}
g_s_test <- g_m2d$par$bi[dat_test$i] + g_m2d$par$bu[dat_test$u] +
  rowSums(g_m2d$par$P[dat_test$u, , drop = FALSE] *
          g_m2d$par$Q[dat_test$i, , drop = FALSE])

g_P <- list()
g_P[["M2d"]] <- local({
  th <- c(g_m2d$par$t1, g_m2d$par$t1 + cumsum(exp(g_m2d$par$beta)))
  CP <- vapply(th, function(t) plogis(t - g_s_test), numeric(g_n))
  tidy_probs(cbind(CP[, 1], CP[, 2:4] - CP[, 1:3], 1 - CP[, 4])) })
g_P[["M2d+T"]] <- local({
  CP <- plogis(g_Tmat(g_mod$par, dat_test$u) - g_s_test)
  tidy_probs(cbind(CP[, 1], CP[, 2:4] - CP[, 1:3], 1 - CP[, 4])) })
g_P[["M2e"]] <- tidy_probs(ordrec_probs_m3(g_pm2e, dat_test$u, dat_test$i))

g_ll <- lapply(g_P, function(P) log(P[cbind(seq_len(g_n), g_y)]))
g_anchor <- c(M2d = -1.185765, `M2d+T` = -1.161220, M2e = -1.154440)
for (k in names(g_anchor))
  if (abs(mean(g_ll[[k]]) - g_anchor[[k]]) > 1e-4)
    stop(sprintf("%s scores %.6f, expected %.6f. Session or cache is wrong.",
                 k, mean(g_ll[[k]]), g_anchor[[k]]))
cat("Guard passed: all three models reproduce their 52 likelihoods.\n")


# =====================================================================
# 1. Calibration and sharpness, with intervals
# =====================================================================
# Calibration error is a function of column means over the whole test
# set, not a row-level average, so it cannot be bootstrapped by
# resampling a difference vector. It is recomputed inside every
# replicate under that replicate's weights.

g_CP  <- lapply(g_P, function(P) t(apply(P, 1, cumsum)))
g_IND <- outer(g_y, LEVELS5, function(a, b) as.numeric(a <= b))
g_sd  <- lapply(g_P, function(P) {
  E <- as.numeric(P %*% LEVELS5)
  sqrt(pmax(as.numeric(P %*% (LEVELS5^2)) - E^2, 0)) })

g_calib_w <- function(CP, w) {
  sw <- sum(w)
  mean(abs(colSums(CP[, 1:4] * w) / sw - colSums(g_IND[, 1:4] * w) / sw))
}
g_uk <- as.integer(factor(test$userId)); g_ik <- as.integer(factor(test$movieId))
g_nUt <- max(g_uk); g_nIt <- max(g_ik)

set.seed(SEED)
g_W <- matrix(0, N_BOOT_LOCAL, g_n)
for (b in seq_len(N_BOOT_LOCAL)) {
  cu <- tabulate(sample.int(g_nUt, g_nUt, TRUE), nbins = g_nUt)
  ci <- tabulate(sample.int(g_nIt, g_nIt, TRUE), nbins = g_nIt)
  g_W[b, ] <- cu[g_uk] * ci[g_ik]
}
g_wmean <- function(v, w) sum(v * w) / sum(w)

g_lvl <- rbindlist(lapply(names(g_P), function(k) {
  bs_c <- vapply(seq_len(N_BOOT_LOCAL), function(b) g_calib_w(g_CP[[k]], g_W[b, ]), numeric(1))
  bs_s <- vapply(seq_len(N_BOOT_LOCAL), function(b) g_wmean(g_sd[[k]], g_W[b, ]), numeric(1))
  bs_l <- vapply(seq_len(N_BOOT_LOCAL), function(b) g_wmean(g_ll[[k]], g_W[b, ]), numeric(1))
  data.table(model = k,
             calib_err = g_calib_w(g_CP[[k]], rep(1, g_n)),
             calib_lo = quantile(bs_c, .025), calib_hi = quantile(bs_c, .975),
             pred_SD = mean(g_sd[[k]]),
             sd_lo = quantile(bs_s, .025), sd_hi = quantile(bs_s, .975),
             meanLL = mean(g_ll[[k]]),
             ll_lo = quantile(bs_l, .025), ll_hi = quantile(bs_l, .975))
}))
cat("\n=== Levels, two-way bootstrap 95% intervals ===\n")
print(g_lvl[, .(model,
                calib_err = sprintf("%.5f [%.5f, %.5f]", calib_err, calib_lo, calib_hi),
                pred_SD   = sprintf("%.4f [%.4f, %.4f]", pred_SD, sd_lo, sd_hi),
                meanLL    = sprintf("%.5f [%.5f, %.5f]", meanLL, ll_lo, ll_hi))])

# Paired contrasts. The difference is formed inside the replicate, so
# the interval reflects the paired comparison and not two independent
# estimates.
g_pairs <- list(c("M2d", "M2d+T"), c("M2d+T", "M2e"), c("M2d", "M2e"))
g_con <- rbindlist(lapply(g_pairs, function(p) {
  a <- p[1]; z <- p[2]
  bc <- vapply(seq_len(N_BOOT_LOCAL), function(b)
    g_calib_w(g_CP[[z]], g_W[b, ]) - g_calib_w(g_CP[[a]], g_W[b, ]), numeric(1))
  bs <- vapply(seq_len(N_BOOT_LOCAL), function(b)
    g_wmean(g_sd[[z]], g_W[b, ]) - g_wmean(g_sd[[a]], g_W[b, ]), numeric(1))
  bl <- vapply(seq_len(N_BOOT_LOCAL), function(b)
    g_wmean(g_ll[[z]] - g_ll[[a]], g_W[b, ]), numeric(1))
  qq <- function(x) quantile(x, c(.025, .975))
  data.table(step = paste(a, "->", z),
             d_calib = g_calib_w(g_CP[[z]], rep(1, g_n)) - g_calib_w(g_CP[[a]], rep(1, g_n)),
             c_lo = qq(bc)[1], c_hi = qq(bc)[2],
             d_SD = mean(g_sd[[z]]) - mean(g_sd[[a]]),
             s_lo = qq(bs)[1], s_hi = qq(bs)[2],
             d_LL = mean(g_ll[[z]] - g_ll[[a]]),
             l_lo = qq(bl)[1], l_hi = qq(bl)[2])
}))
g_star <- function(lo, hi) ifelse(lo > 0 | hi < 0, "*", "")
cat("\n=== Contrasts. Negative calibration error is better; negative SD is sharper ===\n")
print(g_con[, .(step,
  calibration = sprintf("%+.5f [%+.5f, %+.5f]%s", d_calib, c_lo, c_hi, g_star(c_lo, c_hi)),
  sharpness   = sprintf("%+.5f [%+.5f, %+.5f]%s", d_SD, s_lo, s_hi, g_star(s_lo, s_hi)),
  likelihood  = sprintf("%+.5f [%+.5f, %+.5f]%s", d_LL, l_lo, l_hi, g_star(l_lo, l_hi)))])

cat("\n=== Which component carries which criterion ===\n")
g_tot_c <- g_con[step == "M2d -> M2e", d_calib]
g_tot_l <- g_con[step == "M2d -> M2e", d_LL]
cat(sprintf("  calibration : cut points %.1f%% | refitted factors %.1f%%\n",
            100 * g_con[step == "M2d -> M2d+T", d_calib] / g_tot_c,
            100 * g_con[step == "M2d+T -> M2e", d_calib] / g_tot_c))
cat(sprintf("  likelihood  : cut points %.1f%% | refitted factors %.1f%%\n",
            100 * g_con[step == "M2d -> M2d+T", d_LL] / g_tot_l,
            100 * g_con[step == "M2d+T -> M2e", d_LL] / g_tot_l))
cat("  Write the reversal only if BOTH calibration contrasts have intervals\n")
cat("  excluding zero. Two point estimates in opposite proportions with\n")
cat("  overlapping intervals are not a reversal.\n")
fwrite(g_lvl, file.path(OUT_DIR, "t105_channel_levels.csv"))
fwrite(g_con, file.path(OUT_DIR, "t106_channel_contrasts.csv"))


# =====================================================================
# 2. Is M2d+T a third instance of sharpness without calibration?
# =====================================================================
# The pattern is: likelihood improves, predictive SD falls, calibration
# does not improve. M2d+T is checked against that template explicitly,
# so the answer is a verdict rather than an impression.

g_row <- g_con[step == "M2d -> M2d+T"]
g_ll_up    <- g_row$l_lo > 0
g_sharper  <- g_row$s_hi < 0
g_cal_null <- g_row$c_lo < 0 && g_row$c_hi > 0
cat("\n=== Template check for M2d -> M2d+T ===\n")
cat(sprintf("  likelihood improves, interval excludes zero : %s\n", g_ll_up))
cat(sprintf("  predictive SD falls, interval excludes zero : %s\n", g_sharper))
cat(sprintf("  calibration change indistinguishable from 0 : %s\n", g_cal_null))
if (g_ll_up && g_sharper && g_cal_null) {
  cat("  ALL THREE HOLD. M2d+T is a third instance, and it sits inside the\n")
  cat("  thesis's own main result rather than in an ablation. The +0.0246\n")
  cat("  must not be described as a calibration gain.\n")
} else if (g_ll_up && g_row$c_hi < 0) {
  cat("  Likelihood and calibration both improve. M2d+T is NOT an instance\n")
  cat("  of the pattern; per-user cut points improve the distribution's\n")
  cat("  location, and the sharpening is incidental. This is the reading\n")
  cat("  the chapter currently assumes, and it is now supported.\n")
} else {
  cat("  Mixed. State each channel separately and draw no summary verdict.\n")
}

# PIT, with the anchor checked rather than assumed. The handoff records
# M2e's PIT total variation as 0.0197 but not the binning that produced
# it, so several are tried and the one that reproduces the anchor is
# the one to quote. If none does, PIT is left out for M2d+T.
g_pit <- function(P, y, seed) {
  set.seed(seed)
  CP <- t(apply(P, 1, cumsum))
  lo <- ifelse(y == 1L, 0, CP[cbind(seq_along(y), pmax(y - 1L, 1L))])
  hi <- CP[cbind(seq_along(y), y)]
  lo + runif(length(y)) * (hi - lo)
}
cat("\n=== PIT total variation from uniform, several binnings ===\n")
g_pt <- rbindlist(lapply(c(5L, 10L, 20L), function(nb) {
  as.data.table(c(list(bins = nb), lapply(g_P, function(P) {
    u <- g_pit(P, g_y, SEED)
    0.5 * sum(abs(prop.table(tabulate(pmin(ceiling(u * nb), nb), nb)) - 1 / nb))
  })))
}))
print(g_pt[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
g_hit <- g_pt[abs(M2e - 0.0197) < 5e-4]
if (nrow(g_hit)) {
  cat(sprintf("\n  %d bins reproduces the recorded M2e value of 0.0197.\n", g_hit$bins[1]))
  cat(sprintf("  At that binning M2d+T is %.4f against M2e %.4f and M2d %.4f.\n",
              g_hit[["M2d+T"]][1], g_hit$M2e[1], g_hit$M2d[1]))
} else {
  cat("\n  No binning reproduces 0.0197. The recorded PIT figures were\n")
  cat("  computed some other way. Do not quote a PIT value for M2d+T, and\n")
  cat("  record the binning alongside the existing PIT numbers before they\n")
  cat("  are used again.\n")
}
fwrite(g_pt, file.path(OUT_DIR, "t107_pit_binning.csv"))


# =====================================================================
# 3. The null for cut-point spacing agreement
# =====================================================================

g_sp_m2dT <- exp(g_mod$par$beta)
g_sp_m2e  <- exp(g_pm2e$beta)

g_win <- rbindlist(list(fit[, .(user_idx, rating5 = as.integer(rating5))],
                        calib[, .(user_idx, rating5 = as.integer(rating5))]))
g_cnt <- g_win[, as.list(tabulate(rating5, K_LEVELS)), by = user_idx]
setnames(g_cnt, paste0("V", 1:5), paste0("c", 1:5))
setorder(g_cnt, user_idx)
g_C <- as.matrix(g_cnt[, paste0("c", 1:5), with = FALSE])
g_full <- rowSums(g_C > 0) == 5L
cat(sprintf("\nUsers using all five levels in the training window: %d of %d (%.1f%%)\n",
            sum(g_full), nrow(g_C), 100 * mean(g_full)))

# Spacings implied by an empirical distribution: the cut points that a
# model with no covariates at all would place for this user.
g_implied <- function(Cm, add) {
  P <- (Cm + add) / rowSums(Cm + add)
  Q <- t(apply(P, 1, cumsum))[, 1:4, drop = FALSE]
  Th <- qlogis(pmin(pmax(Q, 1e-9), 1 - 1e-9))
  cbind(Th[, 2] - Th[, 1], Th[, 3] - Th[, 2], Th[, 4] - Th[, 3])
}
g_strict <- g_implied(g_C[g_full, , drop = FALSE], 0)
g_smooth <- g_implied(g_C, 1)

g_cmp <- function(A, B, rows, lab) data.table(
  comparison = lab, n_users = length(rows),
  s1 = cor(A[rows, 1], B[rows, 1], method = "spearman"),
  s2 = cor(A[rows, 2], B[rows, 2], method = "spearman"),
  s3 = cor(A[rows, 3], B[rows, 3], method = "spearman"))

g_rows_f <- which(g_full); g_rows_a <- seq_len(g_nU)
g_null <- rbindlist(list(
  g_cmp(g_sp_m2dT, g_sp_m2e, g_rows_f, "M2d+T vs M2e            (complete users)"),
  g_cmp(g_sp_m2dT[g_rows_f, , drop = FALSE], g_strict, seq_along(g_rows_f),
        "M2d+T vs empirical NULL (complete users)"),
  g_cmp(g_sp_m2e[g_rows_f, , drop = FALSE],  g_strict, seq_along(g_rows_f),
        "M2e   vs empirical NULL (complete users)"),
  g_cmp(g_sp_m2dT, g_sp_m2e, g_rows_a, "M2d+T vs M2e            (all, smoothed)"),
  g_cmp(g_sp_m2dT, g_smooth, g_rows_a, "M2d+T vs empirical NULL (all, smoothed)"),
  g_cmp(g_sp_m2e,  g_smooth, g_rows_a, "M2e   vs empirical NULL (all, smoothed)")))
cat("\n=== Spacing agreement against the model-free null ===\n")
print(g_null[, .(comparison, n_users, s1 = round(s1, 4),
                 s2 = round(s2, 4), s3 = round(s3, 4))])
cat("\n  The two model rows are the claim; the four null rows are the bar.\n")
cat("  If the null rows sit at or above the model rows, the agreement in\n")
cat("  52 section 5 is a property of the users' own histograms and gives\n")
cat("  section 6d no independent support. Only agreement clearly ABOVE the\n")
cat("  null shows the two fits recovering something the histogram alone\n")
cat("  does not contain.\n")
fwrite(g_null, file.path(OUT_DIR, "t108_spacing_null.csv"))


# =====================================================================
# 4. Threshold gain by training volume, with intervals
# =====================================================================
# 52 banded on rows in `fit` alone, which put 82 users below the
# selection floor of ten. The selection rule counts the whole training
# window, so both are reported and the window count is primary.

g_nf <- fit[, .N, by = user_idx]; g_nfv <- rep(0L, g_nU); g_nfv[g_nf$user_idx] <- g_nf$N
g_nw <- g_win[, .N, by = user_idx]; g_nwv <- rep(0L, g_nU); g_nwv[g_nw$user_idx] <- g_nw$N
cat(sprintf("\nUsers below 10 rows: by fit only %d | by training window %d\n",
            sum(g_nfv < 10), sum(g_nwv < 10)))

g_d <- g_ll[["M2d+T"]] - g_ll[["M2d"]]
g_band <- cut(g_nwv[dat_test$u], c(0, 20, 45, 100, Inf),
              labels = c("10-20", "21-45", "46-100", "100+"))
g_bt <- rbindlist(lapply(levels(g_band), function(bd) {
  rows <- which(g_band == bd)
  bs <- vapply(seq_len(N_BOOT_LOCAL), function(b) {
    w <- g_W[b, rows]; if (sum(w) == 0) NA_real_ else sum(g_d[rows] * w) / sum(w)
  }, numeric(1))
  q <- quantile(bs, c(.025, .975), na.rm = TRUE)
  data.table(band = bd, users = uniqueN(dat_test$u[rows]), rows = length(rows),
             gain = mean(g_d[rows]), lo = q[[1]], hi = q[[2]])
}))
g_bt[, sig := g_star(lo, hi)]
cat("\n=== Cut-point gain by training-window volume ===\n")
print(g_bt[, .(band, users, rows, gain = round(gain, 5),
               lo = round(lo, 5), hi = round(hi, 5), sig)])

# The 100+ fall-back, tested rather than read off. Four estimates
# arriving in a shape is not a shape; the difference is formed inside
# each replicate so the comparison is paired.
g_r3 <- which(g_band == "46-100"); g_r4 <- which(g_band == "100+")
g_diff <- vapply(seq_len(N_BOOT_LOCAL), function(b) {
  w3 <- g_W[b, g_r3]; w4 <- g_W[b, g_r4]
  if (sum(w3) == 0 || sum(w4) == 0) return(NA_real_)
  sum(g_d[g_r4] * w4) / sum(w4) - sum(g_d[g_r3] * w3) / sum(w3)
}, numeric(1))
g_q <- quantile(g_diff, c(.025, .975), na.rm = TRUE)
cat(sprintf("\n  100+ minus 46-100, differenced within replicate: %+.5f [%+.5f, %+.5f]%s\n",
            mean(g_d[g_r4]) - mean(g_d[g_r3]), g_q[[1]], g_q[[2]],
            g_star(g_q[[1]], g_q[[2]])))
cat("  If this spans zero, the fall-back at 100+ is not established and the\n")
cat("  banding should be reported as a step at roughly 45 with no further\n")
cat("  structure, matching section 5.\n")
fwrite(g_bt, file.path(OUT_DIR, "t109_threshold_gain_bands.csv"))

g_bt[, band := factor(band, levels = c("10-20", "21-45", "46-100", "100+"))]
p76 <- ggplot(g_bt, aes(band, gain)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, colour = "#D55E00") +
  geom_point(size = 2.2, colour = "#D55E00") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.margin = margin(8, 16, 8, 8)) +
  labs(title = "Log-likelihood gain from per-user cut points, by training volume",
       x = "Ratings in the training window",
       y = "Mean log-likelihood gain")
print(p76)
ggsave(file.path(OUT_DIR, "fig5_2_threshold_gain_bands.png"), p76,
       width = 6.0, height = 3.9, units = "in", dpi = 300, bg = "white")

cat("\nWritten: t105_channel_levels.csv, t106_channel_contrasts.csv,\n")
cat("         t107_pit_binning.csv, t108_spacing_null.csv,\n")
cat("         t109_threshold_gain_bands.csv, p76_threshold_gain_bands.png\n")
