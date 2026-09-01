# =====================================================================
# 55_calibration_inference.R

# RUN INSIDE THE SAME CLEAN 06b SESSION, AFTER 52 AND 54.
# =====================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

.need <- c("m2e", "m3a", "m3b", "Cmat", "dat_test", "test", "fit", "calib",
           "ordrec_probs_m3", "tidy_probs", "user_map", "movie_map",
           "LEVELS5", "K_LEVELS", "OUT_DIR", "CACHE_M3", "SEED")
.miss <- .need[!vapply(.need, exists, logical(1))]
if (length(.miss)) stop("Run inside the clean 06b session after 54. Missing: ",
                        paste(.miss, collapse = ", "))
if (!identical(as.integer(test$user_idx), as.integer(dat_test$u)))
  stop("test and dat_test row orders differ. Restart R and re-run 06b.")

B_BOOT <- 2000L
B_NULL <- 1000L
B_COX  <- 800L      # fewer: each replicate refits eight logistic models
h_y <- as.integer(test$rating5); h_n <- length(h_y)


# =====================================================================
# 1. Assemble every model that exists in this session
# =====================================================================

h_pm2e <- m2e$par; h_pm2e$w <- NULL; h_pm2e$W <- NULL
h_C7   <- file.path(dirname(normalizePath(CACHE_M3)), "07_M1_M2")
h_rd   <- function(f) { p <- file.path(h_C7, f); if (file.exists(p)) readRDS(p) else NULL }
h_m2d  <- h_rd("m2d_final.rds")
h_svd  <- h_rd("svd_factors.rds")
h_r7   <- h_rd("m1_m2_results.rds")
h_sig  <- if (!is.null(h_r7) && !is.null(h_r7$sigma_m1)) h_r7$sigma_m1 else NA_real_
h_mod  <- readRDS(file.path(CACHE_M3, "m2dT_frozen_factors.rds"))
h_b30  <- readRDS(file.path(CACHE_M3, "m3b_own_lambda_final.rds"))

h_Tm <- function(par, u) {
  E <- exp(par$beta[u, , drop = FALSE]); t1 <- par$t1[u]
  cbind(t1, t1 + E[, 1], t1 + E[, 1] + E[, 2], t1 + E[, 1] + E[, 2] + E[, 3])
}
h_s_frozen <- h_m2d$par$bi[dat_test$i] + h_m2d$par$bu[dat_test$u] +
  rowSums(h_m2d$par$P[dat_test$u, , drop = FALSE] *
          h_m2d$par$Q[dat_test$i, , drop = FALSE])

h_P <- list()
if (!is.null(h_svd) && !is.na(h_sig)) {
  sc <- pmin(pmax(rowSums(h_svd$P[dat_test$u, , drop = FALSE] *
                          h_svd$Q[dat_test$i, , drop = FALSE]), 1), 5)
  CP <- vapply(c(1.5, 2.5, 3.5, 4.5), function(c) pnorm((c - sc) / h_sig), numeric(h_n))
  h_P[["M1"]] <- tidy_probs(cbind(CP[, 1], CP[, 2:4] - CP[, 1:3], 1 - CP[, 4]))
}
h_P[["M2d"]] <- local({
  th <- c(h_m2d$par$t1, h_m2d$par$t1 + cumsum(exp(h_m2d$par$beta)))
  CP <- vapply(th, function(t) plogis(t - h_s_frozen), numeric(h_n))
  tidy_probs(cbind(CP[, 1], CP[, 2:4] - CP[, 1:3], 1 - CP[, 4])) })
h_P[["M2d+T"]] <- local({
  CP <- plogis(h_Tm(h_mod$par, dat_test$u) - h_s_frozen)
  tidy_probs(cbind(CP[, 1], CP[, 2:4] - CP[, 1:3], 1 - CP[, 4])) })
h_P[["M2e"]]     <- tidy_probs(ordrec_probs_m3(h_pm2e,     dat_test$u, dat_test$i))
h_P[["M3a"]]     <- tidy_probs(ordrec_probs_m3(m3a$par,    dat_test$u, dat_test$i))
h_P[["M3b(1)"]]  <- tidy_probs(ordrec_probs_m3(m3b$par,    dat_test$u, dat_test$i))
h_P[["M3b(30)"]] <- tidy_probs(ordrec_probs_m3(h_b30$par,  dat_test$u, dat_test$i))

h_ll <- vapply(h_P, function(P) mean(log(P[cbind(seq_len(h_n), h_y)])), numeric(1))
h_anch <- c(M2e = -1.154440, M2d = -1.185765, `M2d+T` = -1.161220)
for (k in names(h_anch))
  if (abs(h_ll[[k]] - h_anch[[k]]) > 1e-4)
    stop(sprintf("%s scores %.6f, expected %.6f. Session is wrong.", k, h_ll[[k]], h_anch[[k]]))
cat("Guard passed. Models:", paste(names(h_P), collapse = ", "), "\n")

h_CP  <- lapply(h_P, function(P) t(apply(P, 1, cumsum))[, 1:4, drop = FALSE])
h_IND <- outer(h_y, 1:4, function(a, b) as.numeric(a <= b))
h_folded <- function(CPm, IND, w) {
  sw <- sum(w)
  mean(abs(colSums(CPm * w) / sw - colSums(IND * w) / sw))
}


# =====================================================================
# 2. The noise floor: is any model distinguishable from perfect?
# =====================================================================
# Outcomes are drawn from each model's own predicted distribution, so
# the simulated world is one in which that model is exactly right. The
# calibration error observed there is what perfection costs on 16,694
# rows with these predictions. Nothing about the model's quality enters
# except through the shape of its predictions, which is the point: a
# diffuse model and a sharp one do not face the same floor.

h_draw <- function(P, u01) {
  CPfull <- t(apply(P, 1, cumsum))
  1L + rowSums(u01 > CPfull[, 1:4, drop = FALSE])
}
set.seed(SEED)
h_U <- matrix(runif(B_NULL * h_n), B_NULL, h_n)

h_floor <- rbindlist(lapply(names(h_P), function(k) {
  P <- h_P[[k]]; CPm <- h_CP[[k]]
  obs <- h_folded(CPm, h_IND, rep(1, h_n))
  sim <- vapply(seq_len(B_NULL), function(b) {
    ys <- h_draw(P, h_U[b, ])
    h_folded(CPm, outer(ys, 1:4, function(a, z) as.numeric(a <= z)), rep(1, h_n))
  }, numeric(1))
  data.table(model = k, observed = obs,
             floor_median = median(sim), floor_q95 = quantile(sim, .95),
             excess = obs - median(sim),
             p_value = mean(sim >= obs))
}))
cat("\n=== Calibration error against its own perfect-calibration floor ===\n")
print(h_floor[, .(model, observed = round(observed, 5),
                  floor_median = round(floor_median, 5),
                  floor_q95 = round(floor_q95, 5),
                  excess = round(excess, 5),
                  p = round(p_value, 3),
                  verdict = ifelse(p_value > 0.05, "not distinguishable", "miscalibrated"))])
cat("\n  p is the share of perfectly-calibrated worlds producing a value at\n")
cat("  least as large as the one observed. A model whose p exceeds 0.05 is\n")
cat("  not distinguishable from perfectly calibrated ON THIS SAMPLE, which\n")
cat("  is a much weaker statement than being well calibrated and must be\n")
cat("  written as such.\n")
cat("  'excess' is the effect size the thesis's own rule asks for:\n")
cat("  observed minus null, not observed alone.\n")
fwrite(h_floor, file.path(OUT_DIR, "t110_calibration_floor.csv"))

# The folding bias, made explicit. The bootstrap mean of a folded
# statistic overstates it, and by more when the truth is nearer zero.
h_uk <- as.integer(factor(test$userId)); h_ik <- as.integer(factor(test$movieId))
h_nU <- max(h_uk); h_nI <- max(h_ik)
set.seed(SEED)
h_W <- matrix(0L, B_BOOT, h_n)
for (b in seq_len(B_BOOT)) {
  cu <- tabulate(sample.int(h_nU, h_nU, TRUE), nbins = h_nU)
  ci <- tabulate(sample.int(h_nI, h_nI, TRUE), nbins = h_nI)
  h_W[b, ] <- as.integer(cu[h_uk] * ci[h_ik])
}
h_bias <- rbindlist(lapply(names(h_P), function(k) {
  obs <- h_folded(h_CP[[k]], h_IND, rep(1, h_n))
  bs  <- vapply(seq_len(B_BOOT), function(b) h_folded(h_CP[[k]], h_IND, h_W[b, ]), numeric(1))
  data.table(model = k, observed = obs, boot_mean = mean(bs),
             upward_bias = mean(bs) - obs, ratio = mean(bs) / obs)
}))
setorder(h_bias, observed)
cat("\n=== Upward bias of the bootstrap, ordered by how well calibrated ===\n")
print(h_bias[, .(model, observed = round(observed, 5), boot_mean = round(boot_mean, 5),
                 upward_bias = round(upward_bias, 5), ratio = round(ratio, 2))])
cat("\n  If the ratio rises as the observed value falls, the statistic\n")
cat("  penalises the best-calibrated models hardest and no interval built\n")
cat("  on it is usable. Retire the folded statistic from all inference;\n")
cat("  it may remain in the summary table as a description with no CI.\n")
fwrite(h_bias, file.path(OUT_DIR, "t111_folding_bias.csv"))


# =====================================================================
# 3. Signed calibration-in-the-large, per threshold
# =====================================================================
# d_k = mean(CP_k) - mean(IND_k). Positive means the model puts more
# probability on "at most k" than actually occurred: too pessimistic at
# that cut. These are means, so they bootstrap without difficulty and
# their contrasts are meaningful.

h_dk <- function(CPm, w) { sw <- sum(w); colSums(CPm * w) / sw - colSums(h_IND * w) / sw }
h_sgn <- rbindlist(lapply(names(h_P), function(k) {
  pt <- h_dk(h_CP[[k]], rep(1, h_n))
  bs <- vapply(seq_len(B_BOOT), function(b) h_dk(h_CP[[k]], h_W[b, ]), numeric(4))
  rbindlist(lapply(1:4, function(j) {
    q <- quantile(bs[j, ], c(.025, .975))
    data.table(model = k, cut = c("1|2", "2|3", "3|4", "4|5")[j],
               d = pt[j], lo = q[[1]], hi = q[[2]],
               sig = if (q[[1]] > 0 || q[[2]] < 0) "*" else "")
  }))
}))
cat("\n=== Signed deviation at each cut point, two-way bootstrap ===\n")
print(dcast(h_sgn, model ~ cut,
            value.var = "d")[, lapply(.SD, function(x)
              if (is.numeric(x)) round(x, 5) else x)])
cat("\n  With significance:\n")
print(dcast(h_sgn[, .(model, cut, s = sprintf("%+.4f%s", d, sig))],
            model ~ cut, value.var = "s"))
cat("\n  A model with deviations of opposite sign at different cuts is not\n")
cat("  well calibrated with a small average; it is miscalibrated in two\n")
cat("  directions that cancel. Read this table before quoting any single\n")
cat("  calibration number.\n")
fwrite(h_sgn, file.path(OUT_DIR, "t112_signed_calibration.csv"))

# Paired contrasts on the signed components, for the steps that matter.
h_steps <- list(c("M2d", "M2d+T"), c("M2d+T", "M2e"), c("M2d", "M2e"),
                c("M2e", "M3a"), c("M2e", "M3b(30)"))
h_steps <- Filter(function(p) all(p %in% names(h_P)), h_steps)
h_ct <- rbindlist(lapply(h_steps, function(p) {
  bs <- vapply(seq_len(B_BOOT), function(b)
    abs(h_dk(h_CP[[p[2]]], h_W[b, ])) - abs(h_dk(h_CP[[p[1]]], h_W[b, ])), numeric(4))
  pt <- abs(h_dk(h_CP[[p[2]]], rep(1, h_n))) - abs(h_dk(h_CP[[p[1]]], rep(1, h_n)))
  rbindlist(lapply(1:4, function(j) {
    q <- quantile(bs[j, ], c(.025, .975))
    data.table(step = paste(p[1], "->", p[2]), cut = c("1|2", "2|3", "3|4", "4|5")[j],
               d_abs = pt[j], lo = q[[1]], hi = q[[2]],
               sig = if (q[[1]] > 0 || q[[2]] < 0) "*" else "")
  }))
}))
cat("\n=== Change in absolute deviation, per cut. Negative is better ===\n")
print(dcast(h_ct[, .(step, cut, s = sprintf("%+.4f%s", d_abs, sig))],
            step ~ cut, value.var = "s"))
cat("\n  The folding is applied AFTER the difference is formed inside each\n")
cat("  replicate, so this is a paired comparison of two folded quantities\n")
cat("  rather than a folded comparison, and the bias in section 2 largely\n")
cat("  cancels. This is the table that should replace the calibration\n")
cat("  column of the 7d contrasts.\n")
fwrite(h_ct, file.path(OUT_DIR, "t113_calibration_contrasts.csv"))


# =====================================================================
# 4. Weak calibration: Cox intercept and slope
# =====================================================================
# For each cut k, regress the binary outcome I(y <= k) on the model's
# own logit for that event. Perfect calibration gives intercept 0 and
# slope 1. A slope below 1 means the predicted probabilities are too
# extreme -- too confident -- which is precisely the channel the thesis
# argues M1 and M3b(30) exploit. The slope separates that from a level
# error in a way the deviation d_k cannot.

h_irls <- function(X, y, w, iter = 25L) {
  b <- c(0, 1)
  for (it in seq_len(iter)) {
    eta <- as.numeric(X %*% b)
    mu  <- plogis(eta)
    v   <- pmax(mu * (1 - mu), 1e-9)
    z   <- eta + (y - mu) / v
    Wv  <- w * v
    A   <- crossprod(X, Wv * X); rhs <- crossprod(X, Wv * z)
    nb  <- tryCatch(solve(A, rhs), error = function(e) NULL)
    if (is.null(nb) || any(!is.finite(nb))) return(c(NA_real_, NA_real_))
    if (max(abs(nb - b)) < 1e-8) { b <- as.numeric(nb); break }
    b <- as.numeric(nb)
  }
  b
}
h_lg <- function(p) qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))

set.seed(SEED)
h_cox <- rbindlist(lapply(names(h_P), function(k) {
  rbindlist(lapply(1:4, function(j) {
    X <- cbind(1, h_lg(h_CP[[k]][, j])); yv <- h_IND[, j]
    pt <- h_irls(X, yv, rep(1, h_n))
    bs <- vapply(seq_len(B_COX), function(b) h_irls(X, yv, as.numeric(h_W[b, ])), numeric(2))
    qi <- quantile(bs[1, ], c(.025, .975), na.rm = TRUE)
    qs <- quantile(bs[2, ], c(.025, .975), na.rm = TRUE)
    data.table(model = k, cut = c("1|2", "2|3", "3|4", "4|5")[j],
               intercept = pt[1], i_lo = qi[[1]], i_hi = qi[[2]],
               slope = pt[2], s_lo = qs[[1]], s_hi = qs[[2]])
  }))
}))
h_cox[, `:=`(i_sig = ifelse(i_lo > 0 | i_hi < 0, "*", ""),
             s_sig = ifelse(s_lo > 1 | s_hi < 1, "*", ""))]
cat("\n=== Cox slope per cut. 1 is perfect; below 1 means too extreme ===\n")
print(dcast(h_cox[, .(model, cut, s = sprintf("%.3f%s", slope, s_sig))],
            model ~ cut, value.var = "s"))
cat("\n=== Cox intercept per cut. 0 is perfect ===\n")
print(dcast(h_cox[, .(model, cut, s = sprintf("%+.3f%s", intercept, i_sig))],
            model ~ cut, value.var = "s"))
cat("\n  A star on the slope marks a departure from 1, not from 0.\n")
cat("  Slopes systematically below 1 are the signature of overconfident\n")
cat("  predictions. Compare M1 and M3b(30) against M2e here: if the thesis\n")
cat("  is right that their likelihood gains are sharpness, their slopes\n")
cat("  should sit further below 1.\n")
h_summ <- h_cox[, .(mean_slope = mean(slope), mean_abs_intercept = mean(abs(intercept)),
                    cuts_slope_below_1 = sum(s_hi < 1)), by = model]
print(h_summ[, .(model, mean_slope = round(mean_slope, 4),
                 mean_abs_intercept = round(mean_abs_intercept, 4),
                 cuts_slope_below_1)])
fwrite(h_cox, file.path(OUT_DIR, "t114_cox_calibration.csv"))

p77 <- ggplot(h_cox, aes(cut, slope, group = model)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45", linewidth = 0.4) +
  geom_line(colour = "#D55E00", linewidth = 0.6) +
  geom_point(colour = "#D55E00", size = 1.8) +
  facet_wrap(~ model, nrow = 2) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.margin = margin(8, 16, 8, 8)) +
  labs(title = "Cox calibration slope at each cut point, by model",
       x = "Cut point", y = "Cox slope")
print(p77)
ggsave(file.path(OUT_DIR, "fig5_4_cox_slopes.png"), p77,
       width = 6.5, height = 4.5, units = "in", dpi = 300, bg = "white")


# =====================================================================
# 5. The paired test 54 section 3 was missing
# =====================================================================
# 54 compared a model-model spacing correlation with a model-null one
# by reading their magnitudes. Both are computed on the same users and
# are strongly dependent, so the difference needs an interval. Users
# are resampled; the correlations are recomputed inside the replicate.

h_win <- rbindlist(list(fit[, .(user_idx, rating5 = as.integer(rating5))],
                        calib[, .(user_idx, rating5 = as.integer(rating5))]))
h_cnt <- h_win[, as.list(tabulate(rating5, K_LEVELS)), by = user_idx]
setorder(h_cnt, user_idx)
h_Cm <- as.matrix(h_cnt[, paste0("V", 1:5), with = FALSE])
h_ok <- which(rowSums(h_Cm > 0) == 5L)
h_impl <- local({
  P <- h_Cm[h_ok, , drop = FALSE] / rowSums(h_Cm[h_ok, , drop = FALSE])
  Q <- t(apply(P, 1, cumsum))[, 1:4, drop = FALSE]
  Th <- qlogis(pmin(pmax(Q, 1e-9), 1 - 1e-9))
  cbind(Th[, 2] - Th[, 1], Th[, 3] - Th[, 2], Th[, 4] - Th[, 3]) })
h_A <- exp(h_mod$par$beta)[h_ok, , drop = FALSE]
h_B <- exp(h_pm2e$beta)[h_ok, , drop = FALSE]

set.seed(SEED)
h_nk <- length(h_ok)
h_gap <- rbindlist(lapply(1:3, function(j) {
  f <- function(idx) c(cor(h_A[idx, j], h_B[idx, j], method = "spearman"),
                       cor(h_B[idx, j], h_impl[idx, j], method = "spearman"))
  pt <- f(seq_len(h_nk))
  bs <- vapply(seq_len(B_BOOT), function(b) { idx <- sample.int(h_nk, h_nk, TRUE)
                                              d <- f(idx); d[1] - d[2] }, numeric(1))
  q <- quantile(bs, c(.025, .975), na.rm = TRUE)
  data.table(spacing = c("1|2 to 2|3", "2|3 to 3|4", "3|4 to 4|5")[j],
             model_model = pt[1], model_null = pt[2], difference = pt[1] - pt[2],
             lo = q[[1]], hi = q[[2]],
             sig = if (q[[1]] > 0 || q[[2]] < 0) "*" else "")
}))
cat(sprintf("\n=== Spacing agreement above the null, %d complete users, paired ===\n", h_nk))
print(h_gap[, .(spacing, model_model = round(model_model, 4),
                model_null = round(model_null, 4),
                difference = round(difference, 4),
                lo = round(lo, 4), hi = round(hi, 4), sig)])
cat("\n  A difference whose interval excludes zero means the two fits agree\n")
cat("  by more than two summaries of the same histogram would, which is\n")
cat("  what section 6d needs. A difference spanning zero means the\n")
cat("  agreement reported in 52 is the histogram and nothing more.\n")
fwrite(h_gap, file.path(OUT_DIR, "t115_spacing_paired.csv"))

cat("\nWritten: t110_calibration_floor.csv, t111_folding_bias.csv,\n")
cat("         t112_signed_calibration.csv, t113_calibration_contrasts.csv,\n")
cat("         t114_cox_calibration.csv, t115_spacing_paired.csv,\n")
cat("         p77_cox_slopes.png\n")
