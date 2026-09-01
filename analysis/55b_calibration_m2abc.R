# =====================================================================
# 55b_calibration_m2abc.R
#
# 55 computed the replacement calibration measures for M1, M2d, M2d+T,
# M2e and the three M3 fits, but not for M2a, M2b and M2c. Every
# calibration claim in the thesis that says "all six models" or "the
# only model in the ladder" therefore rests on a table with three rungs
# missing. This file fills them in.
#
# RUN INSIDE THE SAME CLEAN 06b SESSION, IMMEDIATELY AFTER 55.
# It reuses 55's bootstrap weight matrix h_W, so the new intervals are
# paired with the existing ones rather than drawn independently.
#
# ---------------------------------------------------------------------
# WHY THE MODELS ARE REFITTED RATHER THAN LOADED
#
# m1_m2_results.rds stores the metrics table and the threshold table but
# not the clm/clmm2 objects, so P_m2a, P_m2b and P_m2c cannot be
# reconstructed from cache alone. They are refitted here on the same
# calibration set from the same frozen SVD factors, and then checked
# against the cache in two independent ways before anything is computed:
#
#   (1) the refitted thresholds against h_r7$thresholds
#   (2) the refitted test mean log-likelihood against h_r7$metrics
#
# If either guard fires, the refit is not the model the thesis reports
# and the calibration numbers below would be numbers for a different
# model. Do not widen the tolerances. Diagnose in this order: the row
# count of calib, whether the SVD score is the unclipped one, and the
# random-intercept coverage printed for M2c.
#
# TWO DETAILS THAT ARE EASY TO GET WRONG, BOTH TAKEN FROM 07
#
#   - M1 scores the CLIPPED SVD score, M2a-M2c score the UNCLIPPED one.
#   - M2c evaluates test users with no calibration history at a random
#     intercept of zero.
# =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ordinal)
})

REFIT_M2C <- TRUE     # clmm2 on ~98k rows with ~5.5k intercepts: minutes
TOL_LL    <- 1e-4     # refit vs cached test mean log-likelihood
TOL_TH    <- 1e-3     # refit vs cached thresholds

.need <- c("h_W", "h_P", "h_CP", "h_IND", "h_dk", "h_irls", "h_lg",
           "h_folded", "h_svd", "h_r7", "h_sgn", "h_cox", "h_bias",
           "h_n", "B_BOOT", "B_COX", "SEED", "OUT_DIR",
           "calib", "test", "tidy_probs", "LEVELS5")
.miss <- .need[!vapply(.need, exists, logical(1))]
if (length(.miss))
  stop("Run inside the clean 06b session after 55. Missing: ",
       paste(.miss, collapse = ", "))
if (is.null(h_svd))
  stop("h_svd is NULL: the 07 cache was not found, so the frozen factors ",
       "that M2a-M2c are built on are unavailable.")
stopifnot(nrow(test) == h_n,
          identical(dim(h_W), c(as.integer(B_BOOT), as.integer(h_n))))

MODEL_ORDER <- c("M1", "M2a", "M2b", "M2c", "M2d", "M2d+T", "M2e",
                 "M3a", "M3b(1)", "M3b(30)")
CUTS <- c("1|2", "2|3", "3|4", "4|5")


# =====================================================================
# 1. The frozen score, on the calibration set and on the test month
# =====================================================================
# The same inner product 55 uses for M1, but without the clipping: the
# cumulative link models were fitted on the raw score.

j_score <- function(d)
  rowSums(h_svd$P[d$user_idx, , drop = FALSE] *
          h_svd$Q[d$movie_idx, , drop = FALSE])

j_s_cal <- j_score(calib)
j_s_tst <- j_score(test)

cat(sprintf("calib rows %s | test rows %s\n",
            format(nrow(calib), big.mark = ","),
            format(nrow(test), big.mark = ",")))
if (!is.null(test$svd_score)) {
  j_d <- max(abs(j_s_tst - test$svd_score))
  cat(sprintf("Test score reproduced from factors: max |diff| = %.2e\n", j_d))
  if (j_d > 1e-6)
    stop("The frozen score does not reproduce test$svd_score. The factor ",
         "cache and the index maps are not the pair 07 used.")
}


# =====================================================================
# 2. Refit M2a, M2b, M2c and check them against the 07 cache
# =====================================================================

j_cal <- data.table(
  rating_f  = factor(calib$rating5, levels = LEVELS5, ordered = TRUE),
  svd_score = j_s_cal,
  user_f    = factor(calib$userId))

cat("\nFitting M2a (equidistant) and M2b (flexible) ...\n")
j_m2a <- clm(rating_f ~ svd_score, data = j_cal, link = "logit",
             threshold = "equidistant")
j_m2b <- clm(rating_f ~ svd_score, data = j_cal, link = "logit",
             threshold = "flexible")

j_a  <- as.numeric(j_m2a$alpha)
th_a <- j_a[1] + (0:3) * j_a[2]          # expand [first, spacing] to four
b_a  <- as.numeric(j_m2a$beta)
th_b <- as.numeric(j_m2b$alpha)
b_b  <- as.numeric(j_m2b$beta)

th_c <- b_c <- j_bu <- NULL
if (REFIT_M2C) {
  cat("Fitting M2c (clmm2, user random intercept). This is the slow one.\n")
  j_t0  <- Sys.time()
  j_m2c <- clmm2(rating_f ~ svd_score, random = user_f, data = j_cal,
                 link = "logistic", Hess = TRUE, nAGQ = 1)
  cat(sprintf("  fitted in %.1f minutes | random intercept SD %.4f\n",
              as.numeric(difftime(Sys.time(), j_t0, units = "mins")),
              j_m2c$stDev))
  th_c <- as.numeric(j_m2c$Theta)
  b_c  <- as.numeric(j_m2c$beta)[1]
  j_re <- j_m2c$ranef
  if (is.null(names(j_re))) names(j_re) <- levels(j_cal$user_f)
  j_bu <- j_re[as.character(test$userId)]
  cat(sprintf("  test rows whose user has no calibration history: %d (%.1f%%)\n",
              sum(is.na(j_bu)), 100 * mean(is.na(j_bu))))
  j_bu[is.na(j_bu)] <- 0
  j_bu <- as.numeric(j_bu)
}

# Guard 1: thresholds against the cached threshold table.
j_thc <- as.data.table(h_r7$thresholds)
j_chk_th <- function(m, th) {
  ref <- j_thc[model == m][order(k), threshold]
  if (length(ref) != 4L)
    stop(sprintf("No cached thresholds for %s in h_r7$thresholds.", m))
  d <- max(abs(th - ref))
  cat(sprintf("  %-4s thresholds: max |refit - cached| = %.2e\n", m, d))
  if (d > TOL_TH)
    stop(sprintf("%s thresholds do not reproduce the cached fit.", m))
}
cat("\n=== Guard 1: refitted thresholds against the 07 cache ===\n")
j_chk_th("M2a", th_a); j_chk_th("M2b", th_b)
if (REFIT_M2C) j_chk_th("M2c", th_c)


# =====================================================================
# 3. Predicted distributions on the test month
# =====================================================================
# Identical construction to 55: cumulative logits, differenced, then
# passed through tidy_probs so the floor and the row normalisation match
# every other model in the comparison.

j_probs <- function(th, eta) {
  CP <- vapply(th, function(t) plogis(t - eta), numeric(length(eta)))
  tidy_probs(cbind(CP[, 1], CP[, 2:4] - CP[, 1:3], 1 - CP[, 4]))
}

j_P <- list()
j_P[["M2a"]] <- j_probs(th_a, b_a * j_s_tst)
j_P[["M2b"]] <- j_probs(th_b, b_b * j_s_tst)
if (REFIT_M2C) j_P[["M2c"]] <- j_probs(th_c, b_c * j_s_tst + j_bu)

# Guard 2: test mean log-likelihood against the cached metrics table.
j_y   <- as.integer(test$rating5)
j_ll  <- vapply(j_P, function(P) mean(log(P[cbind(seq_len(h_n), j_y)])),
                numeric(1))
j_met <- as.data.table(h_r7$metrics)
cat("\n=== Guard 2: test mean log-likelihood against the 07 cache ===\n")
for (k in names(j_ll)) {
  ref <- j_met[model == k, meanLL]
  if (!length(ref)) stop(sprintf("No cached meanLL for %s.", k))
  cat(sprintf("  %-4s refit %.6f | cached %.6f | diff %.2e\n",
              k, j_ll[[k]], ref, abs(j_ll[[k]] - ref)))
  if (abs(j_ll[[k]] - ref) > TOL_LL)
    stop(sprintf("%s does not reproduce the reported model.", k))
}
# M1 is already in the session; checking it costs nothing and confirms
# that this session's test rows are the ones 07 scored.
if (!is.null(h_P[["M1"]])) {
  j_m1 <- mean(log(h_P[["M1"]][cbind(seq_len(h_n), j_y)]))
  cat(sprintf("  M1   session %.6f | cached %.6f | diff %.2e\n",
              j_m1, j_met[model == "M1", meanLL],
              abs(j_m1 - j_met[model == "M1", meanLL])))
}

j_CP <- lapply(j_P, function(P) t(apply(P, 1, cumsum))[, 1:4, drop = FALSE])


# =====================================================================
# 4. The three calibration measures, on 55's own bootstrap weights
# =====================================================================

# (a) Folded statistic and its upward bias, for Table B.4.
j_bias <- rbindlist(lapply(names(j_CP), function(k) {
  obs <- h_folded(j_CP[[k]], h_IND, rep(1, h_n))
  bs  <- vapply(seq_len(B_BOOT),
                function(b) h_folded(j_CP[[k]], h_IND, h_W[b, ]), numeric(1))
  data.table(model = k, observed = obs, boot_mean = mean(bs),
             upward_bias = mean(bs) - obs, ratio = mean(bs) / obs)
}))

# (b) Signed deviation at each cut, for Table 5.4 and Table B.1.
#     d_k = mean(CP_k) - mean(IND_k). Positive means more probability is
#     placed at or below cut k than was observed.
j_sgn <- rbindlist(lapply(names(j_CP), function(k) {
  pt <- h_dk(j_CP[[k]], rep(1, h_n))
  bs <- vapply(seq_len(B_BOOT), function(b) h_dk(j_CP[[k]], h_W[b, ]),
               numeric(4))
  rbindlist(lapply(1:4, function(j) {
    q <- quantile(bs[j, ], c(.025, .975))
    data.table(model = k, cut = CUTS[j], d = pt[j],
               lo = q[[1]], hi = q[[2]],
               sig = if (q[[1]] > 0 || q[[2]] < 0) "*" else "")
  }))
}))

# (c) Cox intercept and slope at each cut, for Table B.2 and Figure 5.6.
cat("\nCox slopes for the three new models. 4 cuts x", B_COX,
    "replicates each.\n")
j_cox <- rbindlist(lapply(names(j_CP), function(k) {
  rbindlist(lapply(1:4, function(j) {
    X  <- cbind(1, h_lg(j_CP[[k]][, j])); yv <- h_IND[, j]
    pt <- h_irls(X, yv, rep(1, h_n))
    bs <- vapply(seq_len(B_COX),
                 function(b) h_irls(X, yv, as.numeric(h_W[b, ])), numeric(2))
    qi <- quantile(bs[1, ], c(.025, .975), na.rm = TRUE)
    qs <- quantile(bs[2, ], c(.025, .975), na.rm = TRUE)
    data.table(model = k, cut = CUTS[j],
               intercept = pt[1], i_lo = qi[[1]], i_hi = qi[[2]],
               slope = pt[2], s_lo = qs[[1]], s_hi = qs[[2]])
  }))
}))
j_cox[, `:=`(i_sig = ifelse(i_lo > 0 | i_hi < 0, "*", ""),
             s_sig = ifelse(s_lo > 1 | s_hi < 1, "*", ""))]


# =====================================================================
# 5. Merge with 55's output and print in the shape the tables need
# =====================================================================

j_ord <- function(d) {
  d <- copy(d)
  d[, model := factor(model, levels = MODEL_ORDER)]
  setorder(d, model, na.last = TRUE)
  d[, model := as.character(model)][]
}
all_bias <- j_ord(rbind(as.data.table(h_bias), j_bias, fill = TRUE))
all_sgn  <- j_ord(rbind(as.data.table(h_sgn),  j_sgn,  fill = TRUE))
all_cox  <- j_ord(rbind(as.data.table(h_cox),  j_cox,  fill = TRUE))

cat("\n\n=== Signed deviation, new rows for Table 5.4 / B.1 ===\n")
print(j_sgn[, .(model, cut,
                deviation = sprintf("%+.4f%s", d, sig),
                interval  = sprintf("[%+.4f, %+.4f]", lo, hi))])

cat("\n=== Signed deviation, all models, one row per model ===\n")
print(dcast(all_sgn[, .(model, cut, s = sprintf("%+.4f%s", d, sig))],
            model ~ cut, value.var = "s")[
              order(match(model, MODEL_ORDER))])

cat("\n=== Cox slope, all models. A star marks a departure from 1 ===\n")
print(dcast(all_cox[, .(model, cut, s = sprintf("%.3f%s", slope, s_sig))],
            model ~ cut, value.var = "s")[
              order(match(model, MODEL_ORDER))])

cat("\n=== Cox intercept, all models. A star marks a departure from 0 ===\n")
print(dcast(all_cox[, .(model, cut, s = sprintf("%+.3f%s", intercept, i_sig))],
            model ~ cut, value.var = "s")[
              order(match(model, MODEL_ORDER))])

cat("\n=== Folded statistic and its upward bias, all models ===\n")
print(all_bias[, .(model, observed = round(observed, 5),
                   boot_mean = round(boot_mean, 5),
                   ratio = round(ratio, 2))])

# The claim this file exists to settle.
j_clean <- all_sgn[, .(any_sig = any(sig == "*")), by = model][any_sig == FALSE, model]
cat("\n=== Models with no significant signed deviation at any cut ===\n")
cat("  ", paste(j_clean, collapse = ", "), "\n", sep = "")
cat("  Ladder models only: ",
    paste(intersect(j_clean, c("M1", "M2a", "M2b", "M2c", "M2d", "M2e")),
          collapse = ", "), "\n", sep = "")
cat("  Write the summary sentence from this line, not from memory.\n")

fwrite(all_bias, file.path(OUT_DIR, "t111b_folding_bias_all.csv"))
fwrite(all_sgn,  file.path(OUT_DIR, "t112b_signed_calibration_all.csv"))
fwrite(all_cox,  file.path(OUT_DIR, "t114b_cox_calibration_all.csv"))


# =====================================================================
# 6. Figure 5.6, redrawn with the ladder complete
# =====================================================================

j_fig <- copy(all_cox)
j_fig[, model := factor(model, levels = MODEL_ORDER)]
j_fig[, cut := factor(cut, levels = CUTS)]

p77b <- ggplot(j_fig, aes(cut, slope, group = model)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
  geom_line(linewidth = 0.6, colour = "#D55E00") +
  geom_point(size = 1.8, colour = "#D55E00") +
  facet_wrap(~model, nrow = 2) +
  theme_minimal(base_size = 10) +
  labs(title = "Cox calibration slope at each cut point, by model",
       subtitle = paste0("Dashed line is perfect calibration. Below it the ",
                         "predicted probabilities are too extreme."),
       x = "Cut point", y = "Cox slope")
print(p77b)
ggsave(file.path(OUT_DIR, "p77b_cox_slopes_all.png"), p77b,
       width = 9.5, height = 5.2, dpi = 150)

cat("\nWritten: t111b_folding_bias_all.csv, t112b_signed_calibration_all.csv,\n")
cat("         t114b_cox_calibration_all.csv, p77b_cox_slopes_all.png\n")
cat("The 55 outputs are left untouched; these are new files.\n")
