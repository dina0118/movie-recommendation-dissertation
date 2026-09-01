# =====================================================================
# 19_shape_all_models.R  --  completing Seppo's "simple illustration:
#                            mean, standard deviation"
#
# Run INSIDE the file-11 session, after section 5, like 12-15b.
# Needs: P_list, test, tidy_probs, LEVELS5, MODEL_LEVELS, MODEL_COLS,
#        theme_pub, save_fig, OUT_DIR
#
# ---------------------------------------------------------------------
# WHY THIS REPLACES p13 / p13b RATHER THAN ADDING TO THEM
#
# Three problems with the existing pair.
#
# (1) They show M1, M2d and M2e only. M2c carries a per-user random
#     intercept b_u, so a reader will reasonably ask whether M2c already
#     gives each user their own shape. It does not: b_u shifts the linear
#     predictor eta = beta*y + b_u, and the five probabilities are
#     sigma(Theta_k - eta), so the shape is still a function of one
#     scalar. Leaving M2a-M2c out of the picture leaves that question
#     open. Including them turns a three-way contrast into the complete
#     statement: five of the six rungs are one-curve models, and only
#     M2e is not.
#
# (2) They condition on expected ratings 3.0, 3.5 and 4.0 -- the middle
#     of the scale, which is exactly where M1's discretisation skew is
#     smallest. Read from those panels M1 looks symmetric, which appears
#     to confirm the supervisor's premise; read over the whole test month
#     it is not (median skew -0.155, 45.3% of rows above 0.25 in
#     magnitude). Plotting against the point prediction over the full
#     range shows why both readings arise and settles it.
#
# (3) The request was for a *simple* illustration. A jitter strip of
#     16,694 rows is not one. Section 3 gives the plain version: the five
#     bars, for one prediction, model by model.
# =====================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

.need <- c("P_list", "test", "tidy_probs", "LEVELS5", "MODEL_LEVELS",
           "MODEL_COLS", "theme_pub", "save_fig", "OUT_DIR")
.miss <- .need[!vapply(.need, exists, logical(1))]
if (length(.miss))
  stop("Run this after section 5 of file 11. Missing: ",
       paste(.miss, collapse = ", "))
stopifnot(all(MODEL_LEVELS %in% names(P_list)))

# The five rungs whose predictive shape is a deterministic function of
# the point prediction, and the one rung where it is not.
ONE_CURVE <- c("M1", "M2a", "M2b", "M2c", "M2d")
FREE      <- "M2e"


# =====================================================================
# 1. Mean, SD and skew of every predicted distribution, all six models
# =====================================================================

dist_moments <- function(P) {
  P  <- tidy_probs(P)
  mu <- as.numeric(P %*% LEVELS5)
  d  <- outer(rep(1, nrow(P)), LEVELS5) - mu
  v  <- rowSums(P * d^2)
  s  <- sqrt(pmax(v, 1e-12))
  data.table(mean = mu, sd = s, skew = rowSums(P * d^3) / s^3)
}

shape <- rbindlist(lapply(MODEL_LEVELS, function(m)
  data.table(model = m, dist_moments(P_list[[m]]))))
shape[, model := factor(model, levels = MODEL_LEVELS)]

# ---- the claim, as a number rather than a picture --------------------
# If the shape follows from the point prediction alone, then among rows
# whose expected rating agrees to two decimals the SD and the skew must
# agree too, up to the width of the bin. Binning at 0.01 makes any
# residual within-bin variation negligible for a one-curve model and
# leaves M2e's freedom untouched.
shape[, bin := round(mean, 2)]
collapse <- shape[, .(n = .N, sd_within = sd(sd), skew_within = sd(skew)),
                  by = .(model, bin)][n >= 5L,
                  .(bins = .N,
                    max_sd_within   = max(sd_within,   na.rm = TRUE),
                    max_skew_within = max(skew_within, na.rm = TRUE)),
                  by = model]
cat("=== Within a 0.01-wide band of expected rating, how much can the shape still vary? ===\n")
print(collapse[order(match(model, MODEL_LEVELS))])
cat("\nA one-curve model must sit at essentially zero on both columns:\n")
cat("its whole distribution is fixed once the point prediction is fixed.\n")
fwrite(collapse, file.path(OUT_DIR, "t5_shape_collapse.csv"))

# ---- the skew question, answered over the whole test month -----------
skew_tbl <- shape[, .(sd_min = min(sd), sd_max = max(sd),
                      skew_median = median(skew),
                      pct_skew_gt_025 = 100 * mean(abs(skew) > 0.25),
                      skew_q05 = quantile(skew, 0.05),
                      skew_q95 = quantile(skew, 0.95)), by = model]
cat("\n=== Skewness across the whole test month ===\n")
print(skew_tbl[order(match(model, MODEL_LEVELS))])

# The same numbers restricted to the middle of the scale, which is where
# the earlier figure looked. Reporting both is what reconciles them.
mid <- shape[mean %between% c(2.9, 4.1)]
cat("\n=== The same, restricted to expected ratings 2.9-4.1 ===\n")
print(mid[, .(skew_median = median(skew),
              pct_skew_gt_025 = 100 * mean(abs(skew) > 0.25)),
          by = model][order(match(model, MODEL_LEVELS))])
cat("\nM1's skew is near zero in the middle of the scale and grows towards\n")
cat("both ends: it is the discretisation pressing a symmetric Gaussian\n")
cat("against the boundaries at 1 and 5, not a fitted asymmetry.\n")
fwrite(skew_tbl, file.path(OUT_DIR, "t5_skew_summary.csv"))


# =====================================================================
# 2. P33 -- the complete picture, six models, full range
# =====================================================================
# The five one-curve models are drawn as lines because their rows lie
# exactly on a curve; drawing them as clouds would invite the reader to
# see scatter that is not there. M2e is drawn as points because its rows
# genuinely do not lie on a curve. The difference in geometry IS the
# result, and the subtitle says so.

curve_dat <- shape[model %in% ONE_CURVE][order(model, mean)]
free_dat  <- shape[model == FREE]

SKEW_LIM <- c(-4, 3)
n_out <- shape[!(skew %between% SKEW_LIM), .N, by = model]
cat("\n=== Rows outside the plotted skew range ===\n")
print(merge(shape[, .(total = .N), by = model], n_out, by = "model", all.x = TRUE)[
  , .(model, total, outside = fifelse(is.na(N), 0L, N),
      pct = round(100 * fifelse(is.na(N), 0L, N) / total, 2))][
        order(match(model, MODEL_LEVELS))])

panel_lab <- c(sd = "Spread\nSD of the predicted distribution",
               skew = "Asymmetry\nSkewness of the predicted distribution")
mk_long <- function(d) rbind(
  d[, .(model, mean, value = sd,   panel = panel_lab[["sd"]])],
  d[, .(model, mean, value = skew, panel = panel_lab[["skew"]])])
cl <- mk_long(curve_dat); fr <- mk_long(free_dat)
cl[, panel := factor(panel, levels = panel_lab)]
fr[, panel := factor(panel, levels = panel_lab)]

p33 <- ggplot(mapping = aes(mean, value, colour = model)) +
  geom_hline(data = data.table(panel = factor(panel_lab[["skew"]], levels = panel_lab)),
             aes(yintercept = 0), linetype = "dashed",
             colour = "grey55", linewidth = 0.4, inherit.aes = FALSE) +
  geom_point(data = fr[panel == panel_lab[["sd"]] |
                         (panel == panel_lab[["skew"]] & value %between% SKEW_LIM)],
             size = 0.30, alpha = 0.09) +
  geom_line(data = cl, linewidth = 0.75) +
  facet_wrap(~ panel, scales = "free_y") +
  scale_colour_manual(values = MODEL_COLS, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1,
                                                   linewidth = 1.2))) +
  theme_pub() +
  labs(title = "What each model can say about a rating, beyond where it expects it to fall",
       subtitle = paste0(
         "One value per test rating (n = ", format(nrow(test), big.mark = ","),
         "). Five of the six rungs are drawn as lines and one as a cloud, and that\n",
         "is the finding: for M1, M2a, M2b, M2c and M2d the whole predicted ",
         "distribution follows from the point\nprediction, so their rows lie exactly on a ",
         "curve. M2c's per-user intercept shifts that curve along but does not\nleave it. ",
         "Only M2e's per-user thresholds let two rows with the same expected rating carry ",
         "different\nconfidence and hedge in different directions."),
       x = "Expected rating under the model", y = NULL)
print(p33); save_fig("p33_shape_all_models.png", p33, w = 9.0, h = 4.4)

# The skew panel needs clipping to stay readable; drawn separately so the
# clipped fraction can be stated rather than silently dropped.
p33b <- ggplot(mapping = aes(mean, skew, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_point(data = free_dat[skew %between% SKEW_LIM], size = 0.30, alpha = 0.09) +
  geom_line(data = curve_dat, linewidth = 0.8) +
  scale_colour_manual(values = MODEL_COLS, name = NULL) +
  coord_cartesian(ylim = SKEW_LIM) +
  guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1,
                                                   linewidth = 1.2))) +
  theme_pub() +
  labs(title = "Where the asymmetry in each model comes from",
       subtitle = paste0(
         "M1 is a symmetric Gaussian, but it is scored after discretisation at fixed cuts, ",
         "so its predicted\ndistribution is symmetric only near the middle of the scale and ",
         "leans harder the closer the point\nprediction moves to 1 or 5. That asymmetry is ",
         "not fitted to anyone: it is a function of the score alone.\nSkewness divides by ",
         "the cubed spread, so it grows without bound as a distribution concentrates; the\n",
         "axis is clipped and the excluded fraction is reported in the text."),
       x = "Expected rating under the model", y = "Skewness")
print(p33b); save_fig("p33b_skew_all_models.png", p33b, w = 7.4, h = 4.4)


# =====================================================================
# 3. P34 -- the simple illustration
# =====================================================================
# The plain version of the same point: hold the expected rating at 3.5
# and draw the five probabilities. One bar set per model. Because the
# five one-curve models give the same answer to every row in the band,
# any row will do for them. For M2e three rows are shown -- the most
# confident, the median, and the most hedged -- because for M2e the row
# matters and that is the whole claim.

TARGET <- 3.5; HALFW <- 0.02
pick <- rbindlist(lapply(MODEL_LEVELS, function(m) {
  idx <- which(abs(shape[model == m, mean] - TARGET) <= HALFW)
  if (!length(idx)) return(NULL)
  sm <- shape[model == m][idx]
  if (m == FREE) {
    sel <- idx[c(which.min(sm$sd), which.min(abs(sm$sd - median(sm$sd))),
                 which.max(sm$sd))]
    lab <- c("M2e, most confident user", "M2e, typical user",
             "M2e, most hedged user")
  } else {
    sel <- idx[1]; lab <- m
  }
  P <- tidy_probs(P_list[[m]])[sel, , drop = FALSE]
  rbindlist(lapply(seq_along(sel), function(j)
    data.table(facet = lab[j], model = m, rating = 1:5, p = P[j, ],
               sd = shape[model == m][sel[j], sd],
               skew = shape[model == m][sel[j], skew])))
}))
stopifnot(nrow(pick) > 0)
pick[, facet := factor(facet, levels = unique(facet))]

lab_dat <- unique(pick[, .(facet, model, sd, skew)])
lab_dat[, lab := sprintf("SD %.2f   skew %+.2f", sd, skew)]

p34 <- ggplot(pick, aes(factor(rating), p, fill = model)) +
  geom_col(width = 0.72, alpha = 0.9) +
  geom_text(data = lab_dat, aes(x = 3, y = Inf, label = lab),
            vjust = 1.6, size = 2.9, colour = "grey25", inherit.aes = FALSE) +
  facet_wrap(~ facet, nrow = 2) +
  scale_fill_manual(values = MODEL_COLS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  theme_pub() +
  labs(title = sprintf("Every one of these predicts the same rating: %.1f stars", TARGET),
       subtitle = paste0(
         "The five probabilities each model assigns to the five rating categories, for test ",
         "rows whose expected\nrating is ", sprintf("%.2f", TARGET), " \u00b1 ",
         sprintf("%.2f", HALFW), ". The point prediction is the same in every panel; what ",
         "differs is how much\nthe model commits to it and which way it leans. For the ",
         "first five models the panel would look the same\nfor any row with this expected ",
         "rating. For M2e it would not, which is why three of its rows are shown."),
       x = "Rating category", y = "Predicted probability")
print(p34); save_fig("p34_simple_illustration.png", p34, w = 8.6, h = 5.0)

cat("\nWritten to", OUT_DIR, ":\n",
    "  t5_shape_collapse.csv  t5_skew_summary.csv\n",
    "  figs: p33 p33b p34\n")
cat("\nThese supersede p13_shape_freedom.png and p13b_shape_freedom_skew.png.\n")
cat("Keep the old pair only if the banded view is wanted as a supplement;\n")
cat("if both go in, say in the caption that they show the same quantity.\n")
