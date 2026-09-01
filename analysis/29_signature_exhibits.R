# =====================================================================
# 29_signature_exhibits.R
#
# The two exhibits that state a core claim in a single figure and that
# nothing in the existing figure set covers.
#
#   (5)  A user's threshold fingerprint  -> p51
#   (23) The sharpness-calibration frontier -> p52
#
# Run INSIDE the file-11 session, after section 5, like 12-15b and 19-25.
# Needs: P_list, test, fit, m2d, m2e, user_threshold_rows,
#        expand_thresholds, tidy_probs, LEVELS5, K_LEVELS, MODEL_LEVELS,
#        MODEL_COLS, theme_pub, save_fig, OUT_DIR, SEED, nU
#
# ---------------------------------------------------------------------
# WHY THESE TWO, AND WHY THEY ARE NOT ALREADY IN THE SET
#
# (5) p15 shows one case-study user's cut points against the global ones,
#     and that user was selected as the most atypical in a filtered pool
#     -- deliberately a favourable case. p21b splits the spacing panel by
#     comparison. Neither shows several users on one axis, which is the
#     thing Koren & Sill (2011) assert in prose and never plot: that one
#     user's "3 stars" sits where another's "4 stars" does. Linero,
#     Bradley & Desai (2018, AoAS 12(4), Fig. 1) is the published idiom.
#     Section 1 also runs the benchmark that decides whether the exhibit
#     is worth showing at all: if the cut points barely move across
#     users, the figure has nothing to say and should be dropped.
#
# (23) The ladder is currently argued one metric at a time. Gneiting,
#     Balabdaoui & Raftery (2007) framed the whole problem as maximising
#     sharpness subject to calibration, and that is a plane, not a line.
#     Sharpness alone is worthless -- a point mass is perfectly sharp --
#     so the axes are only meaningful together, which is exactly why one
#     figure can carry the argument. p04 already does the RMSE/meanLL
#     version of this move; this is the distributional one.
#
# ON WHAT THESE FIGURES MAY NOT CLAIM. Section 2 prints a dominance
# check. If M2e is not lowest on both axes, the caption says it sits on
# the frontier, not that it dominates. Do not write the stronger sentence
# without the printed check supporting it.
# =====================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

.need <- c("P_list", "test", "fit", "m2d", "m2e", "user_threshold_rows",
           "expand_thresholds", "tidy_probs", "LEVELS5", "K_LEVELS",
           "MODEL_LEVELS", "MODEL_COLS", "theme_pub", "save_fig",
           "OUT_DIR", "SEED", "nU")
.miss <- .need[!vapply(.need, exists, logical(1))]
if (length(.miss))
  stop("Run this after section 5 of file 11. Missing: ",
       paste(.miss, collapse = ", "))
stopifnot(all(MODEL_LEVELS %in% names(P_list)))

y_test <- as.integer(test$rating5)
n_test <- nrow(test)
PT     <- lapply(P_list[MODEL_LEVELS], tidy_probs)

# Per-user cut points under M2e, and the one shared set under M2d.
Tall_e <- user_threshold_rows(m2e$par, seq_len(nU))
th_d_g <- expand_thresholds(m2d$par$t1, m2d$par$beta)
stopifnot(ncol(Tall_e) == 4L, length(th_d_g) == 4L)


# =====================================================================
# 1. Exhibit 5 -- a user's threshold fingerprint
# =====================================================================

# ---- 1.1 The benchmark that decides whether to show this at all ------
# The cut points live on the score scale. If they hardly move across
# users, M2e's extra parameters are doing nothing visible and the figure
# should be replaced by the single-pair shift-versus-shape contrast.
spread <- data.table(
  cut  = sprintf("%d|%d", 1:4, 2:5),
  med  = apply(Tall_e, 2, median),
  iqr  = apply(Tall_e, 2, IQR),
  p05  = apply(Tall_e, 2, quantile, 0.05),
  p95  = apply(Tall_e, 2, quantile, 0.95),
  global = th_d_g)
spread[, range_90 := p95 - p05]

cat("=== How far do the cut points move across users? ===\n")
print(spread[, .(cut, global = round(global, 3), median = round(med, 3),
                 IQR = round(iqr, 3), `5%` = round(p05, 3),
                 `95%` = round(p95, 3), `90% range` = round(range_90, 3))])
cat(sprintf("\nDecision rule: show the fingerprint only if some IQR exceeds 0.10.\n"))
cat(sprintf("Largest IQR is %.3f -> %s\n\n", max(spread$iqr),
            if (max(spread$iqr) > 0.10) "SHOW IT" else
              "DROP IT, lead with the single-pair contrast instead"))

# For scale: how large is that spread next to the gaps between the
# shared cut points? A spread comparable to a whole rating band means
# one user's boundary sits where another user's neighbouring one does.
cat(sprintf("Widths between the shared M2d cut points: %s\n",
            paste(sprintf("%.3f", diff(th_d_g)), collapse = "  ")))
cat(sprintf("Median width of the per-user 90%% ranges:  %.3f\n\n",
            median(spread$range_90)))

# ---- 1.2 Choosing the users, by rule, before looking at them ---------
# Only users whose cut points are identified can be read: four cut
# points need the user to have used the levels. The pool rule is the one
# used in files 14/15 so the selection is consistent across the chapter.
u_stats <- fit[, .(n_fit = .N, fit_mean = mean(rating5),
                   n_lv = uniqueN(rating5)), by = user_idx]
u_stats <- merge(u_stats, test[, .(n_test = .N), by = user_idx],
                 by = "user_idx", all.x = TRUE)
u_stats[is.na(n_test), n_test := 0L]

pool <- u_stats[n_fit >= 40L & n_lv == 5L & n_test >= 5L]
stopifnot(nrow(pool) >= 20L)
pool[, centre := rowMeans(Tall_e[user_idx, , drop = FALSE])]
pool[, w3 := Tall_e[user_idx, 3] - Tall_e[user_idx, 2]]

# Percentiles rather than the extremes: the most extreme user of 5,611
# is an outlier and would overstate the effect. The 5th and 95th are
# still ordinary members of the pool.
nearest <- function(x, target) pool$user_idx[which.min(abs(x - target))]
u_gen   <- nearest(pool$centre, quantile(pool$centre, 0.05))
u_typ   <- nearest(pool$centre, median(pool$centre))
u_hsh   <- nearest(pool$centre, quantile(pool$centre, 0.95))
u_clm   <- pool$user_idx[which.max(pool$w3)]

picks <- data.table(
  role = c("generous (5th pct)", "typical (median)",
           "clumps on 3 (widest band)", "harsh (95th pct)"),
  user_idx = c(u_gen, u_typ, u_clm, u_hsh))
if (anyDuplicated(picks$user_idx))
  warning("Two roles selected the same user; widen the pool rule.")
picks <- merge(picks, pool, by = "user_idx", sort = FALSE)
setorder(picks, centre)

# What each picked user gains from having their own cut points.
ll_row <- function(m) log(PT[[m]][cbind(seq_len(n_test), y_test)])
ll_u   <- data.table(user_idx = test$user_idx,
                     M2d = ll_row("M2d"), M2e = ll_row("M2e"))[
                     , .(M2d = mean(M2d), M2e = mean(M2e)), by = user_idx]
picks  <- merge(picks, ll_u, by = "user_idx", all.x = TRUE, sort = FALSE)
picks[, thresholds_gain := M2e - M2d]
setorder(picks, centre)

picks[, lab := factor(sprintf("%s\nuser %d  \u2014  %d ratings, mean %.2f",
                              role, user_idx, n_fit, fit_mean),
                      levels = sprintf("%s\nuser %d  \u2014  %d ratings, mean %.2f",
                                       role, user_idx, n_fit, fit_mean))]

cat("=== The four users, chosen by rule from a pool of",
    nrow(pool), "===\n")
print(picks[, .(role, user_idx, n_fit, fit_mean = round(fit_mean, 2),
                n_test, centre = round(centre, 3),
                thresholds = sapply(user_idx, function(u)
                  paste(sprintf("%+.2f", Tall_e[u, ]), collapse = " ")),
                `M2e-M2d` = round(thresholds_gain, 3))])
fwrite(picks[, .(role, user_idx, n_fit, fit_mean, n_test, centre, w3,
                 t1 = Tall_e[user_idx, 1], t2 = Tall_e[user_idx, 2],
                 t3 = Tall_e[user_idx, 3], t4 = Tall_e[user_idx, 4],
                 M2d, M2e, thresholds_gain)],
       file.path(OUT_DIR, "t51_fingerprint_users.csv"))

# The sentence the figure is built to make concrete: find a pair where
# one user's cut point sits where another user's neighbouring one does.
cross <- CJ(a = picks$user_idx, b = picks$user_idx)[a != b]
cross[, gap := mapply(function(i, j) Tall_e[i, 3] - Tall_e[j, 2], a, b)]
cat(sprintf("\nClosest crossing among the four: user %d's 3|4 boundary sits %.3f from user %d's 2|3.\n",
            cross[which.min(abs(gap)), a], cross[which.min(abs(gap)), gap],
            cross[which.min(abs(gap)), b]))
cat("If that number is near zero, the same score is a 3 for one and a 4 for the other.\n\n")

# ---- 1.3 Panel A: the fingerprints on one axis -----------------------
band <- data.table(k = 1:4, lo = spread$p05, hi = spread$p95)

pts <- rbindlist(lapply(seq_len(nrow(picks)), function(j)
  data.table(lab = picks$lab[j], k = 1:4,
             t = Tall_e[picks$user_idx[j], ], src = "this user (M2e)")))
glb <- data.table(lab = factor("shared cut points (M2d)\nthe same for all 5,611 users",
                               levels = "shared cut points (M2d)\nthe same for all 5,611 users"),
                  k = 1:4, t = th_d_g, src = "everyone (M2d)")
allp <- rbind(pts, glb)
allp[, lab := factor(as.character(lab),
                     levels = c(levels(picks$lab), as.character(glb$lab[1])))]

p51a <- ggplot(allp, aes(t, lab)) +
  geom_rect(data = band, inherit.aes = FALSE,
            aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            fill = "grey88", alpha = 0.55) +
  geom_line(aes(group = lab), colour = "grey80", linewidth = 2.4) +
  geom_point(aes(colour = src), size = 3.1) +
  geom_text(aes(label = sprintf("%d|%d", k, k + 1)), vjust = -1.35,
            size = 2.6, colour = "grey25") +
  scale_colour_manual(values = c("this user (M2e)" = MODEL_COLS[["M2e"]],
                                 "everyone (M2d)"  = MODEL_COLS[["M2d"]]),
                      name = NULL) +
  scale_y_discrete(limits = rev) +
  theme_pub() +
  labs(subtitle = paste0(
         "A.  Where each user's four cut points fall on the score scale. Grey bands span the 5th to ",
         "95th percentile\n     of that cut point across all users, so each fingerprint can be read ",
         "against the population it came from."),
       x = "Position on the score scale", y = NULL)

# ---- 1.4 Panel B: the behaviour the fingerprint came from ------------
hist_u <- rbindlist(lapply(seq_len(nrow(picks)), function(j) {
  u <- picks$user_idx[j]
  h <- merge(data.table(rating5 = LEVELS5),
             fit[user_idx == u, .N, by = rating5], by = "rating5", all.x = TRUE)
  h[is.na(N), N := 0L]
  data.table(lab = picks$lab[j], rating = h$rating5, share = h$N / sum(h$N))
}))
hist_u[, lab := factor(as.character(lab), levels = levels(picks$lab))]

p51b <- ggplot(hist_u, aes(factor(rating), share)) +
  geom_col(fill = MODEL_COLS[["M2e"]], width = 0.72) +
  facet_wrap(~ lab, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_pub() +
  labs(subtitle = paste0(
         "B.  How each of them actually rated in the six-month fit window \u2014 the data the cut ",
         "points above were learned from."),
       x = "Rating given", y = "Share of that user's ratings")

if (requireNamespace("patchwork", quietly = TRUE)) {
  p51 <- patchwork::wrap_plots(p51a, p51b, ncol = 1, heights = c(1.15, 0.85)) +
    patchwork::plot_annotation(
      title = "One user's 4 is another user's 3",
      subtitle = paste0(
        "Four users drawn at stated positions in the distribution of cut-point location, not chosen ",
        "for the story they tell.\nM2d gives every user the same four boundaries; M2e learns them per ",
        "user. Koren & Sill (2011) motivate the design in\nprose \u2014 one user reads 3 stars as close ",
        "to 4, another as close to 1 \u2014 without plotting it; this is that sentence,\nmeasured."),
      theme = theme_pub())
  print(p51); save_fig("p51_threshold_fingerprint.png", p51, w = 9.0, h = 7.4)
} else {
  print(p51a); save_fig("p51a_threshold_fingerprint.png", p51a, w = 8.4, h = 4.0)
  print(p51b); save_fig("p51b_fingerprint_histograms.png", p51b, w = 8.4, h = 2.8)
}


# =====================================================================
# 2. Exhibit 23 -- the sharpness-calibration frontier
# =====================================================================
# Two axes, both computed on the test month, both with the two-way
# bootstrap used everywhere else in the chapter (Owen & Eckles 2012).
#
#   calibration : mean absolute cumulative-calibration gap, the scalar
#                 already defined in file 23 -- mean over the four
#                 thresholds of |mean predicted P(r<=k) - observed|.
#   sharpness   : mean standard deviation of the predictive distribution.
#                 Lower is sharper. On its own this rewards confident
#                 nonsense, which is the whole point of plotting it
#                 against calibration rather than reporting it alone.

BOOT_23 <- 500L

dist_sd <- function(P) {
  mu <- as.numeric(P %*% LEVELS5)
  d  <- outer(rep(1, nrow(P)), LEVELS5) - mu
  sqrt(pmax(rowSums(P * d^2), 1e-12))
}

Fcum <- lapply(PT, function(P)
  t(apply(P, 1, cumsum))[, seq_len(K_LEVELS - 1), drop = FALSE])
Obs  <- vapply(seq_len(K_LEVELS - 1), function(k) as.numeric(y_test <= k),
               numeric(n_test))
SDm  <- vapply(MODEL_LEVELS, function(m) dist_sd(PT[[m]]), numeric(n_test))
LLm  <- vapply(MODEL_LEVELS, function(m)
  log(PT[[m]][cbind(seq_len(n_test), y_test)]), numeric(n_test))

macl_point <- vapply(MODEL_LEVELS, function(m)
  mean(abs(colMeans(Fcum[[m]]) - colMeans(Obs))), numeric(1))
shrp_point <- colMeans(SDm)
ll_point   <- colMeans(LLm)

u_key <- as.integer(factor(test$userId))
i_key <- as.integer(factor(test$movieId))
nU_t  <- max(u_key); nI_t <- max(i_key)

set.seed(SEED)
macl_b <- shrp_b <- matrix(NA_real_, BOOT_23, length(MODEL_LEVELS),
                           dimnames = list(NULL, MODEL_LEVELS))
for (b in seq_len(BOOT_23)) {
  cu <- tabulate(sample.int(nU_t, nU_t, replace = TRUE), nbins = nU_t)
  ci <- tabulate(sample.int(nI_t, nI_t, replace = TRUE), nbins = nI_t)
  w  <- cu[u_key] * ci[i_key]; sw <- sum(w)
  ob <- colSums(Obs * w) / sw
  for (m in MODEL_LEVELS) {
    macl_b[b, m] <- mean(abs(colSums(Fcum[[m]] * w) / sw - ob))
    shrp_b[b, m] <- sum(SDm[, m] * w) / sw
  }
}

front <- data.table(
  model = factor(MODEL_LEVELS, levels = MODEL_LEVELS),
  macl = macl_point[MODEL_LEVELS],
  macl_lo = apply(macl_b, 2, quantile, 0.025)[MODEL_LEVELS],
  macl_hi = apply(macl_b, 2, quantile, 0.975)[MODEL_LEVELS],
  sharp = shrp_point[MODEL_LEVELS],
  sharp_lo = apply(shrp_b, 2, quantile, 0.025)[MODEL_LEVELS],
  sharp_hi = apply(shrp_b, 2, quantile, 0.975)[MODEL_LEVELS],
  meanLL = ll_point[MODEL_LEVELS])

cat("\n=== Sharpness against calibration, test month ===\n")
print(front[, .(model,
                calibration_error = round(macl, 4),
                cal_ci = sprintf("[%.4f, %.4f]", macl_lo, macl_hi),
                sharpness_sd = round(sharp, 4),
                sharp_ci = sprintf("[%.4f, %.4f]", sharp_lo, sharp_hi),
                meanLL = round(meanLL, 4))])
fwrite(front, file.path(OUT_DIR, "t52_sharpness_calibration.csv"))

# ---- the dominance check that governs what the caption may say -------
best <- as.character(front$model[which.min(front$macl)])
cat(sprintf("\nLowest calibration error: %s. Sharpest: %s.\n",
            best, as.character(front$model[which.min(front$sharp)])))
dom <- front[model == "M2e"]
beaten_on_both <- front[model != "M2e" &
                        macl <= dom$macl & sharp <= dom$sharp]
if (nrow(beaten_on_both) == 0L) {
  cat("M2e is not beaten on both axes by any other rung: it is on the frontier.\n")
  strict <- front[model != "M2e" & (macl < dom$macl | sharp < dom$sharp)]
  if (nrow(strict) == 0L)
    cat("It is also strictly better than every other rung on both axes -> 'dominates' is fair.\n")
  else
    cat(sprintf("But %s beat it on one axis, so write 'sits on the frontier', NOT 'dominates'.\n",
                paste(strict$model, collapse = ", ")))
} else {
  cat(sprintf("WARNING: %s is at least as good on both axes. Do not claim M2e is best here.\n",
              paste(beaten_on_both$model, collapse = ", ")))
}
cat("Sharpness is not a virtue on its own \u2014 a point mass is perfectly sharp and\n")
cat("usually wrong. Read the vertical axis only in company with the horizontal one.\n\n")

p52 <- ggplot(front, aes(macl, sharp, colour = model)) +
  geom_path(aes(group = 1), colour = "grey72", linewidth = 0.5,
            arrow = arrow(length = unit(5, "pt"), type = "closed")) +
  geom_linerange(aes(xmin = macl_lo, xmax = macl_hi),
                 linewidth = 0.5, alpha = 0.75) +
  geom_linerange(aes(ymin = sharp_lo, ymax = sharp_hi),
                 linewidth = 0.5, alpha = 0.75) +
  geom_point(size = 3.0) +
  geom_text(aes(label = sprintf("%s   LL %.3f", model, meanLL)),
            vjust = -1.15, size = 3.0, show.legend = FALSE) +
  scale_colour_manual(values = MODEL_COLS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.16)) +
  scale_y_continuous(expand = expansion(mult = 0.14)) +
  theme_pub() +
  labs(title = "Sharper distributions, without paying for it in calibration",
       subtitle = paste0(
         "Each rung of the ladder as a point in the plane Gneiting, Balabdaoui & Raftery (2007) argue ",
         "the problem lives in:\nmaximise sharpness subject to calibration. Bars are two-way bootstrap ",
         "intervals over users and items. Bottom left is\nbetter on both axes. The grey path follows ",
         "the ladder in order."),
       x = "Calibration error  (mean |predicted \u2212 observed| over the four cumulative thresholds)",
       y = "Mean predictive SD\n(down = sharper)")
print(p52)
save_fig("p52_sharpness_calibration.png", p52, w = 7.4, h = 5.2)

cat("Written: p51_threshold_fingerprint.png, p52_sharpness_calibration.png,\n")
cat("         t51_fingerprint_users.csv, t52_sharpness_calibration.csv\n")
