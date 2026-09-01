# =====================================================================
# 16_ten_month_evaluation.R  
# Three things are measured, and they answer different questions:
#   (A) all eligible rows per month      -- the deployment-realistic view
#   (B) balanced panel (users in all 10) -- isolates model staleness from
#                                           a moving test population
#   (C) coverage                         -- the share of raw activity the
#                                           frozen model can score at all,
#                                           which no likelihood metric sees
# =====================================================================

theme_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          axis.title  = element_text(size = base_size),
          axis.text   = element_text(size = base_size - 1),
          strip.text  = element_text(face = "bold", size = base_size),
          legend.text = element_text(size = base_size - 1),
          legend.position = "bottom",
          plot.margin = margin(t = 8, r = 16, b = 8, l = 8))
}

save_fig <- function(name, plot, w = 6.5, h = 4.0) {
  plot <- plot + labs(title = NULL, subtitle = NULL, caption = NULL)
  ggsave(file.path(FIG_DIR, name), plot,
         width = w, height = h, units = "in", dpi = 300, bg = "white")
  invisible(plot)
}

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

# ---- paths ------------------------------------------------------------
DIR_CACHE00 <- "D:/桌面/movie-recommendation-dissertation/06_models/cache/00_data_selection"
RAW_FILE    <- file.path(DIR_CACHE00, "raw_study_period.rds")
ANOM_FILE   <- file.path(DIR_CACHE00, "anomalous_users.rds")
MOVIES_CSV  <- "D:/桌面/movie-recommendation-dissertation/data/ml-32m/movies.csv"

FIRST_MONTH <- "2018-10"
N_MONTHS    <- 10L

stopifnot(file.exists(RAW_FILE), file.exists(ANOM_FILE))

# Fail immediately, not thirty lines in, if this was run standalone.
.need <- c("user_map", "movie_map", "test", "Pmat", "Qmat", "sigma_m1",
           "th_a", "b_a", "th_b", "b_b", "m2c_par", "m2d", "m2e",
           "gauss_probs", "cum_probs", "ordrec_probs_d", "ordrec_probs_e",
           "metric_set", "tidy_probs", "ndcg_per_user", "fcp_stat",
           "LEVELS5", "SEED", "N_BOOT", "MODEL_LEVELS", "MODEL_COLS",
           "theme_pub", "save_fig", "OUT_DIR")
.miss <- .need[!vapply(.need, exists, logical(1))]
if (length(.miss))
  stop("Run this after section 5 of file 11. Missing: ",
       paste(.miss, collapse = ", "))

# month helpers, copied verbatim from 00_data_selection_eda.R so the
# month indexing is identical to the one that built the selection
month_index <- function(ym) {
  y <- as.integer(substr(ym, 1, 4)); m <- as.integer(substr(ym, 6, 7))
  y * 12L + (m - 1L)
}
index_month <- function(idx) sprintf("%04d-%02d", idx %/% 12L, idx %% 12L + 1L)

MI0    <- month_index(FIRST_MONTH)
MONTHS <- index_month(MI0 + seq_len(N_MONTHS) - 1L)
cat("Evaluation months:", paste(MONTHS, collapse = " "), "\n")


# =====================================================================
# 1. Re-slice the ten months
# =====================================================================
# Rule, verified against 00_data_selection_eda.R lines 903-914:
#   (i)   drop the flagged batch-import accounts, as `dat` does
#   (ii)  rating5 := ceiling(rating), so half-stars round UP
#   (iii) keep rows whose user AND movie are in the training maps
#   (iv)  do NOT re-apply Nu >= 10 / Nm >= 10 per month. That threshold
#         defines the training population. Re-applying it monthly would
#         make the denominator move with the window and would confound
#         model staleness with a changing test population.

raw_all <- readRDS(RAW_FILE)
anom    <- readRDS(ANOM_FILE)
setDT(raw_all)

ev <- raw_all[!userId %in% anom][mi >= MI0 & mi < MI0 + N_MONTHS]
rm(raw_all); gc()

ev[, rating5 := as.integer(ceiling(rating))]
stopifnot(all(ev$rating5 %in% 1:5))

ev <- merge(ev, user_map,  by = "userId",  all.x = TRUE, sort = FALSE)
ev <- merge(ev, movie_map, by = "movieId", all.x = TRUE, sort = FALSE)
setDT(ev)
ev[, scoreable := !is.na(user_idx) & !is.na(movie_idx)]
setorder(ev, mi)

cat(sprintf("Raw rows across the ten months: %s | scoreable: %s (%.1f%%)\n",
            format(nrow(ev), big.mark = ","),
            format(sum(ev$scoreable), big.mark = ","),
            100 * mean(ev$scoreable)))


# =====================================================================
# 2. Replication guard: month 1 must reproduce the master table
# =====================================================================
# If the re-slicing rule is right, 2018-10 rebuilt from raw is the same
# set of rows as the session's `test`, and must give the master-table
# numbers to the last decimal. This is the single check that decides
# whether the other nine months mean anything.

oct <- ev[ym == FIRST_MONTH & scoreable == TRUE]

if (nrow(oct) != nrow(test)) {
  cat("\n!! Month-1 row count differs from the session's test object.\n")
  cat(sprintf("   rebuilt: %d | session: %d\n", nrow(oct), nrow(test)))
  a <- oct[,  .(userId, movieId)]; b <- test[, .(userId, movieId)]
  cat("   in rebuilt but not in test:\n"); print(head(fsetdiff(a, b), 10))
  cat("   in test but not in rebuilt:\n"); print(head(fsetdiff(b, a), 10))
  cat("   Likely cause: eligibility built from `fit` (user_map) vs from\n")
  cat("   train_sel (fit + calib). Find the cause, do not widen this.\n")
}
stopifnot(nrow(oct) == nrow(test))


# =====================================================================
# 3. Scoring: six frozen models, one month at a time
# =====================================================================
# The reconstruction is copied from file 11 section 2 expression for
# expression, so a month scored here and the same month scored there
# cannot diverge.

build_probs <- function(d) {
  # gauss_probs()/cum_probs() are only correct for multi-row input
  # (the drop = TRUE bug documented in the handoff). Never call on 1 row.
  stopifnot(nrow(d) >= 2L)
  y_svd  <- rowSums(Pmat[d$user_idx, , drop = FALSE] *
                      Qmat[d$movie_idx, , drop = FALSE])
  y_clip <- pmin(pmax(y_svd, 1), 5)
  bu <- m2c_par$ranef[as.character(d$userId)]
  bu[is.na(bu)] <- 0
  list(M1  = gauss_probs(y_clip, sigma_m1),
       M2a = cum_probs(th_a, b_a * y_svd),
       M2b = cum_probs(th_b, b_b * y_svd),
       M2c = cum_probs(m2c_par$Theta, m2c_par$beta * y_svd + as.numeric(bu)),
       M2d = ordrec_probs_d(m2d$par, d$user_idx, d$movie_idx),
       M2e = ordrec_probs_e(m2e$par, d$user_idx, d$movie_idx))
}

# NDCG@10 and FCP both use min_items = 5L, so they are computed on a
# smaller user base than meanLL, and that base thins as monthly activity
# falls. n_users_rank is carried alongside so the two curves are never
# read as if they described the same sample.
eval_month <- function(d, panel_label) {
  P  <- build_probs(d)
  ms <- lapply(P, metric_set, y = d$rating5)
  out <- rbindlist(lapply(names(P), function(m) {
    E  <- as.numeric(tidy_probs(P[[m]]) %*% LEVELS5)
    nd <- ndcg_per_user(E, d$rating5, d$userId)
    f  <- fcp_stat(E, d$rating5, d$userId)
    cbind(data.table(model = m), ms[[m]]$tbl,
          data.table(NDCG10       = if (nrow(nd)) mean(nd$ndcg, na.rm = TRUE) else NA_real_,
                     FCP          = if (is.finite(f$FCP)) f$FCP else NA_real_,
                     n_users_rank = nrow(nd)))
  }))
  out[, `:=`(ym = d$ym[1], panel = panel_label,
             n_rows = nrow(d), n_users = uniqueN(d$userId))]
  # row-level log-likelihoods, carried for the monthly cluster bootstrap
  L <- do.call(cbind, lapply(ms, `[[`, "row_ll"))
  colnames(L) <- names(P)
  attr(out, "ll") <- L
  out
}

res_oct <- eval_month(oct, "all")
cat("\n--- Month 1 (2018-10) rebuilt from raw ---\n")
print(res_oct[, .(model, meanLL, RMSE, RPS, CatAcc, NDCG10, FCP)])

stopifnot(abs(res_oct[model == "M2e", meanLL] - (-1.154440)) < 1e-4,
          abs(res_oct[model == "M2d", meanLL] - (-1.185765)) < 1e-4,
          abs(res_oct[model == "M1",  meanLL] - (-1.248970)) < 1e-4)
cat("Guard passed: the re-sliced month reproduces the master table.\n")


# =====================================================================
# 4. Panel A -- all eligible rows per month
# =====================================================================

scored <- ev[scoreable == TRUE]
ll_by_month <- vector("list", N_MONTHS); names(ll_by_month) <- MONTHS

panelA <- rbindlist(lapply(MONTHS, function(mth) {
  d <- scored[ym == mth]
  if (nrow(d) < 50L) { cat(sprintf("  %s: only %d rows, skipped\n", mth, nrow(d))); return(NULL) }
  r <- eval_month(d, "all")
  ll_by_month[[mth]] <<- list(ll = attr(r, "ll"), uid = d$userId)
  cat(sprintf("  %s: %6d rows, %4d users, M2e meanLL %.5f\n",
              mth, nrow(d), uniqueN(d$userId), r[model == "M2e", meanLL]))
  r
}))


# =====================================================================
# 5. Panel B -- balanced panel
# =====================================================================
# Users with at least one scoreable row in every one of the ten months.
# Composition is then constant by construction, so any movement is the
# model going stale rather than the population changing. The price is
# survivorship bias towards heavy users: this panel is not a sample of
# the user base and must never be described as one.

u_months <- scored[, .(k = uniqueN(ym)), by = userId]
panel_users <- u_months[k == N_MONTHS, userId]
cat(sprintf("\nBalanced panel: %d users of %d ever active (%.1f%%)\n",
            length(panel_users), nrow(u_months),
            100 * length(panel_users) / nrow(u_months)))

panelB <- if (length(panel_users) >= 30L) {
  bal <- scored[userId %in% panel_users]
  cat(sprintf("  panel rows: %s (%.1f%% of scoreable rows)\n",
              format(nrow(bal), big.mark = ","),
              100 * nrow(bal) / nrow(scored)))
  rbindlist(lapply(MONTHS, function(mth) {
    d <- bal[ym == mth]
    if (nrow(d) < 50L) return(NULL)
    eval_month(d, "balanced")
  }))
} else {
  cat("  Fewer than 30 users survive all ten months. The balanced panel\n")
  cat("  is not estimable here; report panel A and say why.\n")
  NULL
}

monthly <- rbindlist(list(panelA, panelB), use.names = TRUE)
monthly[, month_k := match(ym, MONTHS)]
monthly[, model := factor(model, levels = MODEL_LEVELS)]
setcolorder(monthly, c("panel", "ym", "month_k", "model", "n_rows",
                       "n_users", "n_users_rank"))
setorder(monthly, panel, month_k, model)
fwrite(monthly, file.path(OUT_DIR, "t4_monthly_metrics.csv"))


# =====================================================================
# 6. Coverage -- what the frozen model can no longer score
# =====================================================================
# Films released after the training window are absent from movie_map, so
# an increasing share of raw activity falls outside the model's reach as
# the window advances. A model that stays accurate on what it can score
# while scoring less and less is degrading in a way meanLL cannot show.

covg <- ev[, .(n_raw       = .N,
              n_scoreable = sum(scoreable),
              lost_user   = sum(is.na(user_idx)),
              lost_item   = sum(!is.na(user_idx) & is.na(movie_idx)),
              n_items_raw = uniqueN(movieId),
              n_items_new = uniqueN(movieId[is.na(movie_idx)])),
          by = ym][order(ym)]
covg[, retained := n_scoreable / n_raw]
covg[, month_k := match(ym, MONTHS)]

# Decompose the unscoreable-item rows: films released after the training
# window versus films that simply never appeared in it. Only the first
# is genuine catalogue drift; the second is a sparsity artefact.
if (file.exists(MOVIES_CSV)) {
  mv <- fread(MOVIES_CSV)
  mv[, release_year := as.integer(sub(".*\\((\\d{4})\\)\\s*$", "\\1", title))]
  ev <- merge(ev, mv[, .(movieId, release_year)], by = "movieId",
              all.x = TRUE, sort = FALSE)
  newrel <- ev[is.na(movie_idx),
               .(rows_post2018 = sum(release_year >= 2018L, na.rm = TRUE),
                 rows_unseen   = sum(release_year <  2018L, na.rm = TRUE)),
               by = ym]
  covg <- merge(covg, newrel, by = "ym", all.x = TRUE)
  setorder(covg, month_k)
}
fwrite(covg, file.path(OUT_DIR, "t4_monthly_coverage.csv"))
cat("\n--- Coverage by month ---\n")
print(covg[, .(ym, n_raw, n_scoreable, retained = round(retained, 3),
              n_items_new)])


# =====================================================================
# 7. Does M2e win every month?
# =====================================================================
# The decay curve is the visible result, but the more useful one is
# whether the ordering of the ladder survives ten independent months.
# Ten out of ten is a stronger statement about robustness than one
# bootstrap interval on one month.

MET <- c("meanLL", "RMSE", "RPS", "CatAcc", "NDCG10", "FCP")
LOWER_BETTER <- c(RMSE = TRUE, RPS = TRUE, meanLL = FALSE,
                  CatAcc = FALSE, NDCG10 = FALSE, FCP = FALSE)

winners <- rbindlist(lapply(MET, function(m) {
  monthly[!is.na(get(m)), .(metric = m,
                            winner = if (LOWER_BETTER[[m]])
                              as.character(model[which.min(get(m))])
                            else as.character(model[which.max(get(m))])),
          by = .(panel, ym)]
}))
win_tab <- dcast(winners[, .N, by = .(panel, metric, winner)],
                 panel + metric ~ winner, value.var = "N", fill = 0L)
fwrite(win_tab, file.path(OUT_DIR, "t4_monthly_winners.csv"))
cat("\n--- Months won, out of the months evaluated ---\n")
print(win_tab)


# =====================================================================
# 8. Decay slope and the two contrasts that matter
# =====================================================================
# Slope in nats per month, from a straight line through ten points. It
# is a description of the ten months observed, not a forecast, and with
# n = 10 the interval is wide -- report it with the standard error.

slopes <- monthly[panel == "all", {
  f <- lm(meanLL ~ month_k)
  .(slope = coef(f)[2], se = summary(f)$coefficients[2, 2],
    r2 = summary(f)$r.squared)
}, by = model]
cat("\n--- meanLL trend, nats per month (panel A) ---\n")
print(slopes[, .(model, slope = round(slope, 5), se = round(se, 5),
                 r2 = round(r2, 3))])
fwrite(slopes, file.path(OUT_DIR, "t4_decay_slopes.csv"))

# The contrast that isolates per-user thresholds is M2e - M2d, not
# M2e - M1: same architecture, same protocol, differing only in global
# versus per-user cut points. M2e - M1 spans the whole ladder and cannot
# answer a question about thresholds.
gaps <- dcast(monthly[panel == "all"], ym + month_k ~ model,
              value.var = "meanLL")
gaps[, `:=`(g_e_1 = M2e - M1, g_e_d = M2e - M2d, g_d_c = M2d - M2c)]

set.seed(SEED)
boot_gap <- rbindlist(lapply(MONTHS, function(mth) {
  z <- ll_by_month[[mth]]
  if (is.null(z)) return(NULL)
  rb <- split(seq_len(nrow(z$ll)), z$uid); un <- names(rb)
  B <- matrix(NA_real_, N_BOOT, 2)
  for (b in seq_len(N_BOOT)) {
    idx <- unlist(rb[sample(un, length(un), replace = TRUE)], use.names = FALSE)
    mm  <- colMeans(z$ll[idx, , drop = FALSE])
    B[b, ] <- c(mm["M2e"] - mm["M1"], mm["M2e"] - mm["M2d"])
  }
  data.table(ym = mth,
             contrast = c("M2e - M1", "M2e - M2d"),
             lo = c(quantile(B[, 1], 0.025), quantile(B[, 2], 0.025)),
             hi = c(quantile(B[, 1], 0.975), quantile(B[, 2], 0.975)))
}))
gap_long <- melt(gaps[, .(ym, month_k, `M2e - M1` = g_e_1, `M2e - M2d` = g_e_d)],
                 id.vars = c("ym", "month_k"), variable.name = "contrast",
                 value.name = "diff")
gap_long <- merge(gap_long, boot_gap, by = c("ym", "contrast"), all.x = TRUE)
setorder(gap_long, contrast, month_k)
fwrite(gap_long, file.path(OUT_DIR, "t4_monthly_gaps.csv"))
cat("\n--- Monthly contrasts with 95% cluster bootstrap intervals ---\n")
print(gap_long[, .(ym, contrast, diff = round(diff, 4),
                   lo = round(lo, 4), hi = round(hi, 4))])


# =====================================================================
# 9. Figures
# =====================================================================

lab_panel <- c(all = "All eligible rows (composition changes)",
               balanced = "Balanced panel (same users throughout)")

# ---- P22: the six ladders over ten months ---------------------------
p22 <- ggplot(monthly, aes(month_k, meanLL, colour = model, group = model)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y",
             labeller = labeller(panel = lab_panel)) +
  scale_colour_manual(values = MODEL_COLS, name = NULL) +
  scale_x_continuous(breaks = seq(1, N_MONTHS, by = 2),      # 每两个月一个刻度
                     labels = MONTHS[seq(1, N_MONTHS, by = 2)]) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # 去掉 size = 8
  labs(title = "Test mean log-likelihood by month, October 2018 to July 2019",
       x = NULL, y = "Test mean log-likelihood")
print(p22); save_fig("p22_ten_month_meanLL.png", p22, w = 6.5, h = 6.3)

# ---- P23: the contrasts, with intervals -----------------------------
p23 <- ggplot(gap_long, aes(month_k, diff, colour = contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = contrast),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_colour_manual(values = c("M2e - M1" = "#D55E00",
                                 "M2e - M2d" = "#0072B2"), name = NULL) +
  scale_fill_manual(values = c("M2e - M1" = "#D55E00",
                               "M2e - M2d" = "#0072B2"), guide = "none") +
  scale_x_continuous(breaks = seq_len(N_MONTHS), labels = MONTHS) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Does the advantage of the ordinal model decay with distance from training?",
       subtitle = "M2e - M2d isolates per-user thresholds; M2e - M1 spans the whole ladder.\n95% cluster bootstrap intervals, users resampled within each month.",
       x = NULL, y = expression(Delta * " mean log-likelihood"))
print(p23); save_fig("p23_ten_month_contrasts.png", p23, w = 6.8, h = 4.2)

# ---- P24: coverage decay --------------------------------------------
covg_long <- melt(covg[, .(month_k, ym,
                         `scoreable share of raw rows` = retained,
                         `rows scored (thousands)` = n_scoreable / 1000)],
                 id.vars = c("month_k", "ym"))
p24 <- ggplot(covg_long, aes(month_k, value)) +
  geom_col(fill = "#0072B2", alpha = 0.85, width = 0.65) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq_len(N_MONTHS), labels = MONTHS) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "What the frozen model can still score",
       subtitle = "Films released after the training window never enter the item map.\nThis loss is invisible to every likelihood metric, which conditions on the rows that survive.",
       x = NULL, y = NULL)
print(p24); save_fig("p24_coverage_decay.png", p24, w = 6.8, h = 4.6)

# ---- P25: ranking metric, on its own (smaller) user base ------------
p25 <- ggplot(monthly[panel == "all"],
              aes(month_k, NDCG10, colour = model, group = model)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
  scale_colour_manual(values = MODEL_COLS, name = NULL) +
  scale_x_continuous(breaks = seq_len(N_MONTHS), labels = MONTHS) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "NDCG@10 over the same ten months",
       subtitle = sprintf("Users with at least five rated items in the month: %s. This is a different sample from the likelihood curve.",
                          paste(monthly[panel == "all" & model == "M2e", n_users_rank],
                                collapse = ", ")),
       x = NULL, y = "NDCG@10")
print(p25); save_fig("p25_ten_month_ndcg.png", p25, w = 6.8, h = 4.0)

cat("\nWritten to", OUT_DIR, ":\n",
    " t4_monthly_metrics.csv  t4_monthly_coverage.csv\n",
    " t4_monthly_winners.csv  t4_decay_slopes.csv  t4_monthly_gaps.csv\n",
    " figs: p22 p23 p24 p25\n")
