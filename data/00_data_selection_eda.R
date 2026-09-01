# =====================================================================
# 00_data_selection_eda.R
# Data selection for the v7 pipeline: temporal window choice, anomalous
# account removal, user/movie activity thresholds, stratified splits.
#
# Everything downstream (M1, M2e, M2f, M2g, M2h, M3+) consumes the
# artefacts written at the end of this script. Nothing here touches the
# test month except to describe it; all eligibility rules are decided on
# the training window only.
#
# Run order:
#   Stage A (sections 1-5)  : global profiling, anomaly detection
#   Stage B (sections 6-9)  : window + threshold selection
#   Stage C (sections 10-14): final split construction + audits
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

setDTthreads(0)
options(datatable.print.nrows = 30)


# =====================================================================
# 0. Configuration
# =====================================================================

CFG <- list(

  # ---- paths --------------------------------------------------------
  dir_ml       = "D:/桌面/movie-recommendation-dissertation/data/ml-32m/",
  dir_tmdb     = "D:/桌面/movie-recommendation-dissertation/data/TMBD/",
  dir_work     = "D:/桌面/movie-recommendation-dissertation/06_models/",
  dir_out      = "output/00_data_selection/",
  dir_cache    = "cache/00_data_selection/",

  # ---- reproducibility ----------------------------------------------
  seed         = 2026L,

  # ---- global study period ------------------------------------------
  # Restrict all profiling to a period where MovieLens is in its modern
  # steady state. Pre-2010 volumes and rating conventions differ enough
  # that pooling them distorts every month-level statistic below.
  study_start  = "2013-10-12",
  study_end    = "2023-10-12",

  # ---- month-level anomaly detection --------------------------------
  # Robust z on log10 monthly volume. MAD-based so that the anomalous
  # months do not inflate their own reference scale.
  month_z_warn = 2.0,
  month_z_flag = 2.5,

  # ---- user-level anomaly (batch-import) detection -------------------
  # A "batch import" account dumps a large rating history in a very
  # short wall-clock window: high volume, near-zero span, one dominant
  # day, and often degenerate rating variance. Each rule is deliberately
  # conservative on its own; a user is flagged if ANY rule fires.
  ba_min_n              = 300L,   # rules only apply above this volume
  ba_span_days          = 1.0,   # rule 1: n >= ba_min_n within <= 1 day
  ba_frac_max_day       = 0.90,  # rule 2: >=90% of ratings on one day
  ba_rate_per_active_day= 600,   # rule 3: sustained impossible throughput
  ba_const_min_n        = 100L,   # rule 4: >= 100 ratings and sd(rating) == 0
  ba_median_gap_sec     = 2,     # rule 5: median inter-arrival <= 2s

  # ---- candidate windows --------------------------------------------
  win_lengths      = c(3L, 6L, 9L),          # training window, months
  test_length      = 1L,                     # test window, months
  cand_first_month = "2015-01",              # earliest training start
  cand_last_month  = "2023-08",              # latest test month

  # ---- window quality gates -----------------------------------------
  win_max_cv       = 0.25,   # CV of monthly volume inside the window
  win_max_absz     = 2.0,    # no month in window/test may exceed this
  win_max_step     = 0.50,   # max |month-on-month| relative change

  # ---- activity thresholds to scan ----------------------------------
  Nu_grid          = c(5L, 10L, 20L, 30L, 50L),
  Nm_grid          = c(5L, 10L, 20L, 50L),
  prune_mode       = "iterative",  # "single" or "iterative"
  prune_max_iter   = 20L,

  # ---- final choices (set after inspecting Stage B output) ----------
  final_win_start  = "2018-04",  # e.g. "2018-01"
  final_win_len    = 6L,    # e.g. 6L
  final_Nu         = 10L,    # e.g. 20L
  final_Nm         = 10L,    # e.g. 10L
  

  # ---- split construction -------------------------------------------
  calib_frac       = 0.20,
  stratify_by      = "rating5",   # "rating10" or "rating5"
  drop_cold_items  = TRUE,        # test rows on unseen movies removed
  drop_cold_users  = TRUE,        # test rows from ineligible users removed

  # ---- M2h feasibility ----------------------------------------------
  m2h_min_calib_per_user = 10L,   # rows needed to estimate 4 thresholds

  # ---- plotting -----------------------------------------------------
  col_a = "steelblue",
  col_b = "coral",
  col_c = "grey40",
  fig_w = 9, fig_h = 5, fig_dpi = 150
)

setwd(CFG$dir_work)
dir.create(CFG$dir_out,   showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$dir_cache, showWarnings = FALSE, recursive = TRUE)
set.seed(CFG$seed)


theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  legend.position = "bottom",
                  plot.margin = margin(8, 16, 8, 8)))

save_fig <- function(p, name, w = 6.5, h = 4.0) {
  p <- p + labs(subtitle = NULL, caption = NULL) +
    theme(plot.title = element_text(size = 12, face = "plain", hjust = 0))
  ggsave(file.path(CFG$dir_out, paste0(name, ".png")), p,
         width = w, height = h, units = "in", dpi = 300, bg = "white")
  invisible(p)
}

write_tab <- function(dt, name) {
  fwrite(dt, file.path(CFG$dir_out, paste0(name, ".csv")))
  invisible(dt)
}

# Month arithmetic helpers. Months are handled as "YYYY-MM" strings and
# as integer month indices so that window enumeration is pure integer
# arithmetic rather than date arithmetic.
month_index <- function(ym) {
  y <- as.integer(substr(ym, 1, 4)); m <- as.integer(substr(ym, 6, 7))
  y * 12L + (m - 1L)
}
index_month <- function(idx) {
  sprintf("%04d-%02d", idx %/% 12L, idx %% 12L + 1L)
}
month_start_ts <- function(ym) {
  as.numeric(as.POSIXct(paste0(ym, "-01"), tz = "UTC"))
}


# =====================================================================
# 1. Load and integrity-check the raw data
# =====================================================================

cache_raw <- file.path(CFG$dir_cache, "raw_study_period.rds")

if (file.exists(cache_raw)) {
  raw <- readRDS(cache_raw)
} else {
  ratings <- fread(file.path(CFG$dir_ml, "ratings.csv"))
  setnames(ratings, tolower(names(ratings)),
           skip_absent = TRUE)   # tolerate userId/userid variants
  setnames(ratings,
           old = names(ratings),
           new = c("userId", "movieId", "rating", "timestamp")[seq_along(names(ratings))])

  # Integrity audit on the FULL file, before any subsetting.
  integrity_full <- data.table(
    check = c("rows", "distinct_users", "distinct_movies",
              "na_rows", "duplicate_user_movie",
              "rating_min", "rating_max", "off_grid_ratings",
              "ts_min", "ts_max"),
    value = c(
      format(nrow(ratings), big.mark = ","),
      format(uniqueN(ratings$userId), big.mark = ","),
      format(uniqueN(ratings$movieId), big.mark = ","),
      format(ratings[, sum(is.na(userId) | is.na(movieId) |
                             is.na(rating) | is.na(timestamp))], big.mark = ","),
      format(nrow(ratings) - uniqueN(ratings, by = c("userId", "movieId")),
             big.mark = ","),
      as.character(min(ratings$rating)),
      as.character(max(ratings$rating)),
      format(ratings[, sum(abs(rating * 2 - round(rating * 2)) > 1e-9)],
             big.mark = ","),
      as.character(as.Date(as.POSIXct(min(ratings$timestamp), origin = "1970-01-01", tz = "UTC"))),
      as.character(as.Date(as.POSIXct(max(ratings$timestamp), origin = "1970-01-01", tz = "UTC")))
    ))
  write_tab(integrity_full, "01_integrity_full_file")
  print(integrity_full)

  # Full-history monthly series is needed for context plots, so keep a
  # lightweight aggregate before subsetting to the study period.
  ratings[, ym := format(as.POSIXct(timestamp, origin = "1970-01-01",
                                    tz = "UTC"), "%Y-%m")]
  monthly_full <- ratings[, .(n = .N,
                              n_users  = uniqueN(userId),
                              n_movies = uniqueN(movieId),
                              mean_rating = mean(rating)), by = ym][order(ym)]
  saveRDS(monthly_full, file.path(CFG$dir_cache, "monthly_full_history.rds"))

  raw <- ratings[timestamp >= as.numeric(as.POSIXct(CFG$study_start, tz = "UTC")) &
                   timestamp <  as.numeric(as.POSIXct(CFG$study_end,   tz = "UTC"))]
  raw[, mi := month_index(ym)]
  setkey(raw, userId, timestamp)
  saveRDS(raw, cache_raw)
  rm(ratings); gc()
}

monthly_full <- readRDS(file.path(CFG$dir_cache, "monthly_full_history.rds"))

movies <- fread(file.path(CFG$dir_ml, "movies.csv"))
movies[, release_year := as.integer(sub(".*\\((\\d{4})\\)\\s*$", "\\1", title))]

cat(sprintf("\nStudy period rows: %s | users: %s | movies: %s\n",
            format(nrow(raw), big.mark = ","),
            format(uniqueN(raw$userId), big.mark = ","),
            format(uniqueN(raw$movieId), big.mark = ",")))


# =====================================================================
# 2. Global temporal profile
# =====================================================================

# ---- 2a. Full-history context --------------------------------------
monthly_full[, date := as.Date(paste0(ym, "-01"))]

p_hist <- ggplot(monthly_full, aes(date, n)) +
  geom_line(colour = CFG$col_c, linewidth = 0.4) +
  annotate("rect",
           xmin = as.Date(CFG$study_start), xmax = as.Date(CFG$study_end),
           ymin = -Inf, ymax = Inf, fill = CFG$col_a, alpha = 0.12) +
  scale_y_continuous(labels = comma) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "Monthly rating volume in MovieLens-32M, 1995-2023",
       x = NULL, y = "Ratings per month")
save_fig(p_hist, "02a_monthly_volume_full_history", w = 6.5, h = 3.6)

# ---- 2b. Study-period monthly aggregates ---------------------------
monthly <- raw[, .(n            = .N,
                   n_users      = uniqueN(userId),
                   n_movies     = uniqueN(movieId),
                   mean_rating  = mean(rating),
                   sd_rating    = sd(rating),
                   p_half_star  = mean(rating * 2 %% 2 == 1),
                   p_extreme_lo = mean(rating <= 1.0),
                   p_extreme_hi = mean(rating >= 4.5)),
               by = .(ym, mi)][order(mi)]

monthly[, rows_per_user := n / n_users]
monthly[, date := as.Date(paste0(ym, "-01"))]

write_tab(monthly, "02b_monthly_profile")

# ---- 2c. Monthly rating-category composition -----------------------
month_cat <- raw[, .N, by = .(ym, rating)][order(ym, rating)]
month_cat[, share := N / sum(N), by = ym]
month_cat[, date := as.Date(paste0(ym, "-01"))]

p_cat_drift <- ggplot(month_cat, aes(date, share, group = factor(rating),
                                     colour = factor(rating))) +
  geom_line(linewidth = 0.5) +
  scale_y_continuous(labels = percent) +
  scale_colour_viridis_d(option = "D", end = 0.9) +
  labs(title = "Rating-category composition over time",
       subtitle = "Category drift is the mechanism behind distribution-shift effects in M1 vs M2",
       x = NULL, y = "share of month's ratings", colour = "rating")
save_fig(p_cat_drift, "02c_category_composition_drift")


# =====================================================================
# 3. Month-level anomaly detection
# =====================================================================

# Robust z on log10 volume: median and MAD are computed across the study
# period so a handful of extreme months cannot inflate their own scale.
monthly[, log_n := log10(n)]
monthly[, `:=`(med_log = median(log_n),
               mad_log = mad(log_n, constant = 1.4826))]
monthly[, z_robust := (log_n - med_log) / mad_log]

# Month-on-month relative step, which catches level shifts that a
# level-based z can miss when the shift is sustained.
monthly[, step := (n - shift(n)) / shift(n)]

monthly[, flag_month := fifelse(abs(z_robust) >= CFG$month_z_flag, "flag",
                        fifelse(abs(z_robust) >= CFG$month_z_warn, "warn", "ok"))]

write_tab(monthly[, .(ym, n, n_users, n_movies, mean_rating,
                      z_robust, step, flag_month)],
          "03a_monthly_anomaly_scores")

p_month_z <- ggplot(monthly, aes(date, n)) +
  geom_col(aes(fill = flag_month), width = 25) +
  geom_hline(yintercept = 10^monthly$med_log[1], linetype = 2,
             colour = CFG$col_c) +
  scale_fill_manual(values = c(ok = CFG$col_a, warn = "goldenrod",
                               flag = CFG$col_b)) +
  scale_y_continuous(labels = comma) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Monthly volume with robust anomaly flags",
       subtitle = sprintf("Robust z on log10 volume; warn |z|>=%.1f, flag |z|>=%.1f; dashed line is the median month",
                          CFG$month_z_warn, CFG$month_z_flag),
       x = NULL, y = "ratings per month", fill = NULL)
save_fig(p_month_z, "03b_monthly_volume_flagged")

p_month_mean <- ggplot(monthly, aes(date, mean_rating)) +
  geom_line(colour = CFG$col_a) +
  geom_point(aes(colour = flag_month), size = 1.4) +
  scale_colour_manual(values = c(ok = CFG$col_a, warn = "goldenrod",
                                 flag = CFG$col_b)) +
  labs(title = "Mean rating level by month",
       subtitle = "A volume anomaly that also moves the rating level is a different problem from one that does not",
       x = NULL, y = "mean rating", colour = NULL)
save_fig(p_month_mean, "03c_monthly_mean_rating", h = 4)

cat("\n=== Flagged months ===\n")
print(monthly[flag_month != "ok",
              .(ym, n, z_robust = round(z_robust, 2),
                step = round(step, 3), mean_rating = round(mean_rating, 3),
                flag_month)])


# =====================================================================
# 4. User-level anomaly detection (batch-import accounts)
# =====================================================================

# ---- 4a. Behavioural features, computed once over the study period --
# All aggregation is by-reference within data.table; no row loops.
cache_uf <- file.path(CFG$dir_cache, "user_features.rds")

if (file.exists(cache_uf)) {
  user_feat <- readRDS(cache_uf)
} else {
  raw[, day := as.integer(timestamp %/% 86400L)]

  # per-user per-day counts, used for the concentration features
  ud <- raw[, .N, by = .(userId, day)]
  ud_stats <- ud[, .(n_days = .N, max_day_n = max(N)), by = userId]

  # inter-arrival statistics; shift() within user is vectorised
  setkey(raw, userId, timestamp)
  raw[, gap := timestamp - shift(timestamp), by = userId]

  user_feat <- raw[, .(
    n            = .N,
    ts_min       = min(timestamp),
    ts_max       = max(timestamp),
    span_days    = (max(timestamp) - min(timestamp)) / 86400,
    sd_rating    = sd(rating),
    mean_rating  = mean(rating),
    n_distinct_r = uniqueN(rating),
    p_extreme    = mean(rating <= 1.0 | rating >= 4.5),
    p_half_star  = mean(rating * 2 %% 2 == 1),
    median_gap   = as.numeric(median(gap, na.rm = TRUE)),
    n_months     = uniqueN(mi)
  ), by = userId]

  user_feat <- merge(user_feat, ud_stats, by = "userId", all.x = TRUE)

  user_feat[, frac_max_day := max_day_n / n]
  user_feat[, rate_per_active_day := n / pmax(n_days, 1L)]
  user_feat[is.na(sd_rating), sd_rating := 0]

  # Rating entropy: a degenerate profile (one or two categories only) is
  # a secondary signal of non-human input.
  ent <- raw[, .N, by = .(userId, rating)]
  ent[, p := N / sum(N), by = userId]
  ent <- ent[, .(entropy = -sum(p * log(p))), by = userId]
  user_feat <- merge(user_feat, ent, by = "userId", all.x = TRUE)

  raw[, gap := NULL]
  saveRDS(user_feat, cache_uf)
  rm(ud, ud_stats, ent); gc()
}

# ---- 4b. Rule battery ----------------------------------------------
user_feat[, `:=`(
  r1_dump_1day   = n >= CFG$ba_min_n & span_days <= CFG$ba_span_days,
  r2_one_day_hog = n >= CFG$ba_min_n & frac_max_day >= CFG$ba_frac_max_day,
  r3_throughput  = n >= CFG$ba_min_n & rate_per_active_day >= CFG$ba_rate_per_active_day,
  r4_constant    = n >= CFG$ba_const_min_n & sd_rating == 0,
  r5_machine_gap = n >= CFG$ba_min_n & !is.na(median_gap) &
                     median_gap <= CFG$ba_median_gap_sec
)]
rule_cols <- c("r1_dump_1day", "r2_one_day_hog", "r3_throughput",
               "r4_constant", "r5_machine_gap")
user_feat[, n_rules := rowSums(.SD), .SDcols = rule_cols]
user_feat[, anomalous := n_rules >= 1L]

rule_summary <- rbindlist(lapply(rule_cols, function(rc) {
  data.table(rule       = rc,
             n_users    = user_feat[get(rc) == TRUE, .N],
             n_ratings  = user_feat[get(rc) == TRUE, sum(n)],
             pct_ratings= 100 * user_feat[get(rc) == TRUE, sum(n)] /
                            user_feat[, sum(n)])
}))
rule_summary <- rbind(rule_summary,
  data.table(rule = "ANY (union)",
             n_users     = user_feat[anomalous == TRUE, .N],
             n_ratings   = user_feat[anomalous == TRUE, sum(n)],
             pct_ratings = 100 * user_feat[anomalous == TRUE, sum(n)] /
                             user_feat[, sum(n)]))
write_tab(rule_summary, "04a_anomaly_rule_summary")
cat("\n=== Batch-import rule battery ===\n"); print(rule_summary)

# ---- 4c. Diagnostic plots ------------------------------------------
p_anom_scatter <- ggplot(user_feat[n >= 5],
                         aes(pmax(span_days, 0.01), n, colour = anomalous)) +
  geom_point(alpha = 0.25, size = 0.6) +
  scale_x_log10(labels = comma) + scale_y_log10(labels = comma) +
  scale_colour_manual(values = c(`FALSE` = CFG$col_a, `TRUE` = CFG$col_b)) +
  labs(title = "Volume against activity span, by user",
       subtitle = "Batch-import accounts occupy the high-volume / short-span corner",
       x = "activity span (days, log)", y = "ratings (log)",
       colour = "flagged")
save_fig(p_anom_scatter, "04b_user_volume_vs_span")

p_anom_conc <- ggplot(user_feat[n >= CFG$ba_min_n],
                      aes(frac_max_day, fill = anomalous)) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.7) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = c(`FALSE` = CFG$col_a, `TRUE` = CFG$col_b)) +
  labs(title = "Concentration of a user's ratings on their busiest day",
       subtitle = sprintf("Users with >= %d ratings only", CFG$ba_min_n),
       x = "share of ratings on the single busiest day",
       y = "users (log)", fill = "flagged")
save_fig(p_anom_conc, "04c_user_day_concentration")

# ---- 4d. Monthly volume before and after removal --------------------
anom_users <- user_feat[anomalous == TRUE, userId]
saveRDS(anom_users, file.path(CFG$dir_cache, "anomalous_users.rds"))

monthly_clean <- raw[!userId %in% anom_users, .(n_clean = .N), by = .(ym, mi)]
month_cmp <- merge(monthly[, .(ym, mi, date, n, z_robust)],
                   monthly_clean, by = c("ym", "mi"), all.x = TRUE)
month_cmp[is.na(n_clean), n_clean := 0L]
month_cmp[, `:=`(log_c = log10(pmax(n_clean, 1)))]
month_cmp[, z_clean := (log_c - median(log_c)) / mad(log_c, constant = 1.4826)]
month_cmp[, removed_pct := 100 * (n - n_clean) / n]

write_tab(month_cmp[, .(ym, n, n_clean, removed_pct,
                        z_raw = z_robust, z_clean)],
          "04d_monthly_before_after_cleaning")

p_clean_cmp <- ggplot(melt(month_cmp[, .(date, raw = z_robust, cleaned = z_clean)],
                           id.vars = "date", variable.name = "series",
                           value.name = "z"),
                      aes(date, z, colour = series)) +
  geom_hline(yintercept = c(-CFG$month_z_flag, CFG$month_z_flag),
             linetype = 2, colour = CFG$col_c) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c(raw = CFG$col_b, cleaned = CFG$col_a)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Monthly volume anomaly before and after removing flagged accounts",
       subtitle = "Months whose anomaly survives cleaning are genuine demand shifts, not import artefacts",
       x = NULL, y = "robust z on log10 volume", colour = NULL)
save_fig(p_clean_cmp, "04e_anomaly_before_after_cleaning")

# ---- 4e. Working dataset -------------------------------------------
dat <- raw[!userId %in% anom_users]
cat(sprintf("\nAfter anomaly removal: %s rows (-%.2f%%), %s users (-%.2f%%)\n",
            format(nrow(dat), big.mark = ","),
            100 * (1 - nrow(dat) / nrow(raw)),
            format(uniqueN(dat$userId), big.mark = ","),
            100 * (1 - uniqueN(dat$userId) / uniqueN(raw$userId))))


# =====================================================================
# 5. Rating category schemes
# =====================================================================

# Seppo's suggestion: collapse the 10-point half-star scale to 5 whole
# stars. This is carried as a parallel column so every downstream table
# can be produced under either scheme without re-running the pipeline.
dat[, rating10 := factor(rating, levels = seq(0.5, 5.0, by = 0.5),
                         ordered = TRUE)]
dat[, rating5  := factor(ceiling(rating), levels = 1:5, ordered = TRUE)]

cat_cmp <- rbind(
  dat[, .(scheme = "10-category", level = as.character(rating10))][, .N, by = .(scheme, level)],
  dat[, .(scheme = "5-category",  level = as.character(rating5 ))][, .N, by = .(scheme, level)]
)
cat_cmp[, share := N / sum(N), by = scheme]

write_tab(cat_cmp, "05a_category_scheme_comparison")

p_cat <- ggplot(cat_cmp, aes(level, share, fill = scheme)) +
  geom_col(position = "dodge", width = 0.75) +
  facet_wrap(~scheme, scales = "free_x") +  # free_x 已经让两个 scheme 用各自的 levels
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c(`10-category` = CFG$col_a,
                               `5-category`  = CFG$col_b)) +
  labs(title = "Rating scale before and after collapsing half-stars",
       subtitle = "Collapsing removes the sparsest categories, which is what makes per-user thresholds estimable",
       x = NULL, y = "share of ratings") +
  theme(legend.position = "none")
save_fig(p_cat, "05b_category_schemes")

# Sparsity of the rarest category is the binding constraint for M2h.
rare_cat <- rbind(
  dat[, .N, by = rating10][, .(scheme = "10-category",
                               rarest_share = min(N) / sum(N),
                               rarest_level = as.character(rating10[which.min(N)]))],
  dat[, .N, by = rating5 ][, .(scheme = "5-category",
                               rarest_share = min(N) / sum(N),
                               rarest_level = as.character(rating5[which.min(N)]))]
)
write_tab(rare_cat, "05c_rarest_category")
print(rare_cat)


# =====================================================================
# 6. Candidate window enumeration
# =====================================================================

# For every (start month, length) pair, describe the training window and
# the month immediately following it. The loop is over a few hundred
# candidates, not over rows.

mi_first <- month_index(CFG$cand_first_month)
mi_last  <- month_index(CFG$cand_last_month)

monthly_clean_full <- dat[, .(n = .N,
                              n_users = uniqueN(userId),
                              n_movies = uniqueN(movieId),
                              mean_rating = mean(rating)), by = mi][order(mi)]
setkey(monthly_clean_full, mi)

mstats <- copy(monthly_clean_full)
mstats[, log_n := log10(n)]
mstats[, z := (log_n - median(log_n)) / mad(log_n, constant = 1.4826)]

cand_grid <- CJ(start = mi_first:(mi_last), len = CFG$win_lengths)
cand_grid <- cand_grid[start + len + CFG$test_length - 1L <= mi_last]

window_scan <- rbindlist(lapply(seq_len(nrow(cand_grid)), function(i) {
  s   <- cand_grid$start[i]; L <- cand_grid$len[i]
  tr  <- s:(s + L - 1L)
  te  <- (s + L):(s + L + CFG$test_length - 1L)
  mtr <- mstats[.(tr)]; mte <- mstats[.(te)]
  if (anyNA(mtr$n) || anyNA(mte$n)) return(NULL)

  data.table(
    win_start   = index_month(s),
    win_len     = L,
    test_month  = index_month(te[1]),
    n_train     = sum(mtr$n),
    n_test      = sum(mte$n),
    cv_monthly  = sd(mtr$n) / mean(mtr$n),
    max_absz    = max(abs(c(mtr$z, mte$z))),
    max_step    = max(abs(diff(mtr$n) / head(mtr$n, -1L))),
    train_mean_rating = weighted.mean(mtr$mean_rating, mtr$n),
    test_mean_rating  = weighted.mean(mte$mean_rating, mte$n),
    test_train_ratio  = sum(mte$n) / sum(mtr$n)
  )
}))
# 每个窗长内部用自己的分位数作门槛，让三个面板可比
window_scan[, cv_gate := quantile(cv_monthly, 0.60, na.rm = TRUE), by = win_len]
window_scan[, passes := cv_monthly <= cv_gate &
              max_absz   <= CFG$win_max_absz &
              max_step   <= CFG$win_max_step]

write_tab(window_scan, "06a1_window_scan_all")

cat(sprintf("\n%d of %d candidate windows pass the balance gates.\n",
            sum(window_scan$passes), nrow(window_scan)))

p_win <- ggplot(window_scan, aes(n_train, cv_monthly, colour = passes)) +
  geom_point(alpha = 0.8, size = 1.6) +
  geom_hline(data = unique(window_scan[, .(win_len, cv_gate)]),
             aes(yintercept = cv_gate), linetype = 2, colour = CFG$col_c) +
  facet_wrap(~paste0(win_len, "-month window")) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6, "M")) +
  scale_colour_manual(values = c(`FALSE` = CFG$col_c, `TRUE` = CFG$col_a),
                      labels = c("No", "Yes")) +
  labs(title = "Training rows against monthly volume variation, by window length",
       x = "Training rows (millions)", y = "CV of monthly volume",
       colour = "Passes gates")
save_fig(p_win, "06b1_window_scan", h = 4.3)

cat("\n=== Best-balanced passing windows, by length ===\n")
print(window_scan[passes == TRUE][order(win_len, cv_monthly),
                                  head(.SD, 5), by = win_len,
                                  .SDcols = c("win_start", "test_month",
                                              "n_train", "n_test",
                                              "cv_monthly", "max_absz",
                                              "cv_gate")])


# =====================================================================
# 7. Activity distributions inside a candidate window
# =====================================================================

# Eligibility is decided on the TRAINING window only. This function is
# the single definition of "who is in the study", reused by the
# threshold grid, the final split, and the downstream pipeline.

window_rows <- function(d, win_start, win_len, test_len = CFG$test_length) {
  s  <- month_index(win_start)
  list(train = d[mi >= s & mi < s + win_len],
       test  = d[mi >= s + win_len & mi < s + win_len + test_len])
}

# Reference window used for the distribution plots. Defaults to the
# best-balanced 6-month candidate until final_win_start is set.
ref_win <- if (!is.na(CFG$final_win_start)) {
  list(start = CFG$final_win_start, len = CFG$final_win_len)
} else {
  b <- window_scan[passes == TRUE & win_len == 6L][order(cv_monthly)][1]
  list(start = b$win_start, len = b$win_len)
}
cat(sprintf("\nReference window for activity profiling: %s + %d months\n",
            ref_win$start, ref_win$len))

W  <- window_rows(dat, ref_win$start, ref_win$len)
tr <- W$train; te <- W$test

user_act  <- tr[, .(n_train = .N), by = userId][order(-n_train)]
movie_act <- tr[, .(n_train = .N), by = movieId][order(-n_train)]

# ---- 7a. Histograms with candidate thresholds marked ---------------
p_user_hist <- ggplot(user_act, aes(n_train)) +
  geom_histogram(bins = 60, fill = CFG$col_a, alpha = 0.85) +
  geom_vline(xintercept = CFG$Nu_grid, linetype = 2, colour = CFG$col_b) +
  scale_x_log10(labels = comma) + scale_y_continuous(labels = comma) +
  labs(title = "User activity within the training window",
       subtitle = paste0("Dashed lines are the candidate thresholds: ",
                         paste(CFG$Nu_grid, collapse = ", ")),
       x = "ratings in training window (log)", y = "users")
save_fig(p_user_hist, "07a_user_activity_hist")

p_movie_hist <- ggplot(movie_act, aes(n_train)) +
  geom_histogram(bins = 60, fill = CFG$col_a, alpha = 0.85) +
  geom_vline(xintercept = CFG$Nm_grid, linetype = 2, colour = CFG$col_b) +
  scale_x_log10(labels = comma) + scale_y_continuous(labels = comma) +
  labs(title = "Movie support within the training window",
       subtitle = paste0("Dashed lines are the candidate thresholds: ",
                         paste(CFG$Nm_grid, collapse = ", ")),
       x = "ratings in training window (log)", y = "movies")
save_fig(p_movie_hist, "07b_movie_support_hist")

# ---- 7b. ECDFs: users kept vs ratings kept -------------------------
ecdf_dt <- function(act, grid, label) {
  tot_e <- nrow(act); tot_r <- sum(act$n_train)
  rbindlist(lapply(grid, function(k) {
    keep <- act[n_train >= k]
    data.table(entity = label, threshold = k,
               pct_entities = 100 * nrow(keep) / tot_e,
               pct_ratings  = 100 * sum(keep$n_train) / tot_r)
  }))
}
retention <- rbind(
  ecdf_dt(user_act,  CFG$Nu_grid, "users"),
  ecdf_dt(movie_act, CFG$Nm_grid, "movies")
)
write_tab(retention, "07c_threshold_retention")
cat("\n=== Retention at each candidate threshold ===\n"); print(retention)

p_ret <- ggplot(melt(retention, id.vars = c("entity", "threshold"),
                     variable.name = "measure", value.name = "pct"),
                aes(threshold, pct, colour = measure)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_wrap(~entity, scales = "free_x") +
  scale_x_log10() +
  scale_colour_manual(values = c(pct_entities = CFG$col_b,
                                 pct_ratings  = CFG$col_a),
                      labels = c("entities kept", "ratings kept")) +
  labs(title = "What a threshold actually removes",
       subtitle = "Pruning removes many entities but few ratings; the gap between the two curves is the whole argument for pruning",
       x = "minimum ratings in training window (log)", y = "percent retained",
       colour = NULL)
save_fig(p_ret, "07d_retention_curves")

# ---- 7c. Lorenz / concentration ------------------------------------
lorenz <- function(act, label) {
  a <- act[order(n_train)]
  data.table(entity = label,
             p_entities = seq_len(nrow(a)) / nrow(a),
             p_ratings  = cumsum(a$n_train) / sum(a$n_train))
}
lz <- rbind(lorenz(user_act, "users"), lorenz(movie_act, "movies"))
p_lorenz <- ggplot(lz, aes(p_entities, p_ratings, colour = entity)) +
  geom_abline(slope = 1, linetype = 2, colour = CFG$col_c) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(labels = percent) + scale_y_continuous(labels = percent) +
  scale_colour_manual(values = c(users = CFG$col_a, movies = CFG$col_b)) +
  labs(title = "Concentration of ratings across users and movies",
       subtitle = "Distance from the diagonal is the inequality that makes the long tail uninformative",
       x = "cumulative share of entities", y = "cumulative share of ratings",
       colour = NULL)
save_fig(p_lorenz, "07e_lorenz_concentration")

# ---- 7d. Who the low-activity users are ----------------------------
# If pruning removed a distinct kind of user rather than a thin slice of
# ordinary ones, the target population has changed and must be described.
user_act[, band := cut(n_train, c(0, 4, 9, 19, 29, 49, Inf),
                       labels = c("1-4", "5-9", "10-19", "20-29", "30-49", "50+"))]
band_profile <- merge(tr, user_act[, .(userId, band)], by = "userId")[
  , .(n_users        = uniqueN(userId),
      n_ratings      = .N,
      mean_rating    = mean(rating),
      sd_rating      = sd(rating),
      p_extreme      = mean(rating <= 1 | rating >= 4.5),
      p_half_star    = mean(rating * 2 %% 2 == 1),
      median_movie_pop = as.numeric(median(.N)) ), by = band][order(band)]
write_tab(band_profile, "07f_activity_band_profile")
print(band_profile)

p_band <- ggplot(band_profile, aes(band, mean_rating)) +
  geom_col(fill = CFG$col_a, width = 0.7) +
  geom_hline(yintercept = tr[, mean(rating)], linetype = 2, colour = CFG$col_b) +
  labs(title = "Do low-activity users rate differently?",
       subtitle = "Dashed line is the window mean; a flat profile means pruning does not shift the rating distribution",
       x = "ratings in training window", y = "mean rating")
save_fig(p_band, "07g_activity_band_mean_rating", h = 4)


# =====================================================================
# 8. Threshold grid: single-pass and iterative pruning
# =====================================================================

# Single pass applies both thresholds once. Iterative k-core repeats
# until the joint condition holds, which removes strictly more data. The
# two are reported side by side because the literature uses both and the
# choice materially changes what "20-core" means.

apply_prune <- function(train, Nu, Nm, mode = "single",
                        max_iter = CFG$prune_max_iter) {
  d <- train
  iter <- 0L
  repeat {
    iter <- iter + 1L
    uc <- d[, .N, by = userId][N >= Nu, userId]
    mc <- d[, .N, by = movieId][N >= Nm, movieId]
    d2 <- d[userId %in% uc & movieId %in% mc]
    if (mode == "single" || nrow(d2) == nrow(d) || iter >= max_iter) {
      d <- d2; break
    }
    d <- d2
  }
  list(data = d, iters = iter)
}

grid <- CJ(Nu = CFG$Nu_grid, Nm = CFG$Nm_grid, mode = c("single", "iterative"),
           sorted = FALSE)

threshold_grid <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i]
  pr <- apply_prune(tr, g$Nu, g$Nm, g$mode)
  d  <- pr$data
  if (nrow(d) == 0L) return(NULL)

  keep_u <- unique(d$userId); keep_m <- unique(d$movieId)
  te_k   <- te[userId %in% keep_u & movieId %in% keep_m]

  data.table(
    Nu = g$Nu, Nm = g$Nm, mode = g$mode, iters = pr$iters,
    train_rows   = nrow(d),
    train_users  = length(keep_u),
    train_movies = length(keep_m),
    pct_rows     = 100 * nrow(d) / nrow(tr),
    pct_users    = 100 * length(keep_u) / uniqueN(tr$userId),
    pct_movies   = 100 * length(keep_m) / uniqueN(tr$movieId),
    density      = nrow(d) / (length(keep_u) * length(keep_m)),
    rows_per_user  = nrow(d) / length(keep_u),
    rows_per_movie = nrow(d) / length(keep_m),
    test_rows_kept = nrow(te_k),
    pct_test_kept  = 100 * nrow(te_k) / nrow(te),
    test_users_kept = uniqueN(te_k$userId),
    mean_rating_train = d[, mean(rating)],
    mean_rating_test  = if (nrow(te_k)) te_k[, mean(rating)] else NA_real_
  )
}))
threshold_grid[, rating_shift := mean_rating_test - mean_rating_train]
write_tab(threshold_grid, "08a_threshold_grid")

cat("\n=== Threshold grid (iterative) ===\n")
print(threshold_grid[mode == "iterative",
                     .(Nu, Nm, train_rows, train_users, train_movies,
                       pct_rows = round(pct_rows, 1),
                       pct_users = round(pct_users, 1),
                       density = signif(density, 3),
                       pct_test_kept = round(pct_test_kept, 1))])

p_grid_rows <- ggplot(threshold_grid,
                      aes(factor(Nm), factor(Nu), fill = pct_rows)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%", pct_rows)), size = 3) +
  facet_wrap(~mode) +
  scale_fill_gradient(low = "white", high = CFG$col_a) +
  labs(title = "Share of training ratings retained",
       x = "minimum ratings per movie", y = "minimum ratings per user",
       fill = "% rows")
save_fig(p_grid_rows, "08b_grid_rows_retained")

p_grid_test <- ggplot(threshold_grid,
                      aes(factor(Nm), factor(Nu), fill = pct_test_kept)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%", pct_test_kept)), size = 3) +
  facet_wrap(~mode) +
  scale_fill_gradient(low = "white", high = CFG$col_b) +
  labs(title = "Share of the test month that survives the eligibility rules",
       subtitle = "This is the evaluation coverage: what fraction of next month's traffic the study actually speaks to",
       x = "minimum ratings per movie", y = "minimum ratings per user",
       fill = "% test rows")
save_fig(p_grid_test, "08c_grid_test_coverage")

p_grid_dens <- ggplot(threshold_grid[mode == "iterative"],
                      aes(train_rows, density, colour = factor(Nu),
                          shape = factor(Nm))) +
  geom_point(size = 2.5) +
  scale_x_log10(labels = comma) + scale_y_log10() +
  scale_colour_viridis_d(end = 0.85) +
  labs(title = "The pruning trade-off: size against density",
       subtitle = "Density is the quantity that drives collaborative-filtering performance; size drives runtime and variance",
       x = "training rows (log)", y = "matrix density (log)",
       colour = "min/user", shape = "min/movie")
save_fig(p_grid_dens, "08d_size_vs_density")

# ---- 8b. Does pruning change the rating distribution? --------------
# A pruning rule that also shifts the response distribution is not a
# neutral quality filter; it silently redefines the prediction target.
prune_dist <- rbindlist(lapply(CFG$Nu_grid, function(k) {
  pr <- apply_prune(tr, k, 1L, "single")$data
  pr[, .(Nu = k, share = .N / nrow(pr)), by = .(level = as.character(rating5))]
}))
prune_dist <- rbind(prune_dist,
                    tr[, .(Nu = 0L, share = .N / nrow(tr)),
                       by = .(level = as.character(rating5))])
write_tab(prune_dist, "08e_prune_rating_distribution")

p_prune_dist <- ggplot(prune_dist, aes(level, share, fill = factor(Nu))) +
  geom_col(position = "dodge", width = 0.8) +
  scale_y_continuous(labels = percent) +
  scale_fill_viridis_d(end = 0.9) +
  labs(title = "Rating distribution under increasing user thresholds",
       subtitle = "Bars that move with the threshold mean pruning is changing the target, not just the sample size",
       x = "rating (5-category)", y = "share", fill = "min ratings/user")
save_fig(p_prune_dist, "08f_prune_rating_distribution")


# =====================================================================
# 9. Sensitivity of the window choice to the threshold choice
# =====================================================================

# The window and the thresholds are not independent: a longer window
# lets more users clear a given bar. This crosses the top candidate
# windows with the threshold grid so the final pair is chosen jointly.

top_windows <- window_scan[passes == TRUE][order(win_len, cv_monthly),
                                           head(.SD, 3), by = win_len]

joint_scan <- rbindlist(lapply(seq_len(nrow(top_windows)), function(i) {
  w  <- top_windows[i]
  Wi <- window_rows(dat, w$win_start, w$win_len)
  rbindlist(lapply(seq_len(nrow(grid[mode == "iterative"])), function(j) {
    g  <- grid[mode == "iterative"][j]
    pr <- apply_prune(Wi$train, g$Nu, g$Nm, "iterative")$data
    if (nrow(pr) == 0L) return(NULL)
    ku <- unique(pr$userId); km <- unique(pr$movieId)
    tk <- Wi$test[userId %in% ku & movieId %in% km]
    data.table(win_start = w$win_start, win_len = w$win_len,
               test_month = w$test_month, Nu = g$Nu, Nm = g$Nm,
               train_rows = nrow(pr), train_users = length(ku),
               train_movies = length(km),
               density = nrow(pr) / (length(ku) * length(km)),
               test_rows = nrow(tk), pct_test = 100 * nrow(tk) / nrow(Wi$test),
               test_users = uniqueN(tk$userId))
  }))
}))
write_tab(joint_scan, "09a_window_threshold_joint_scan")

p_joint <- ggplot(joint_scan, aes(train_rows, pct_test,
                                  colour = paste0(win_len, "m"),
                                  shape = factor(Nu))) +
  geom_point(size = 2.4, alpha = 0.9) +
  scale_x_log10(labels = comma) +
  scale_colour_manual(values = c(`3m` = CFG$col_c, `6m` = CFG$col_a,
                                 `9m` = CFG$col_b)) +
  labs(title = "Joint choice of window length and activity thresholds",
       subtitle = "Upper-left is the goal: small training set, high evaluation coverage",
       x = "training rows (log)", y = "% of test month retained",
       colour = "window", shape = "min/user")
save_fig(p_joint, "09b_joint_window_threshold")


# =====================================================================
# 10. Build the final selected dataset
# =====================================================================

# Set the four finals in CFG after inspecting sections 6-9, then rerun
# from here. Defaults fall back to the best-balanced 6-month window and
# the mid grid point so the script is runnable end to end on first pass.

FIN <- list(
  win_start = "2018-04",
  win_len   = 6L,
  Nu        = 10L,
  Nm        = 10L
)
cat(sprintf("\n=== FINAL SELECTION: window %s +%dm | Nu=%d | Nm=%d | %s pruning ===\n",
            FIN$win_start, FIN$win_len, FIN$Nu, FIN$Nm, CFG$prune_mode))

WF <- window_rows(dat, FIN$win_start, FIN$win_len)
train_raw <- WF$train
test_raw  <- WF$test

pr <- apply_prune(train_raw, FIN$Nu, FIN$Nm, CFG$prune_mode)
train_sel <- pr$data
elig_users  <- unique(train_sel$userId)
elig_movies <- unique(train_sel$movieId)

# Test-side eligibility. Cold rows are separated rather than silently
# dropped so the audit in section 12 can quantify what was excluded.
test_raw[, cold_user := !userId  %in% elig_users]
test_raw[, cold_item := !movieId %in% elig_movies]
test_sel <- test_raw[
  (!CFG$drop_cold_users | cold_user == FALSE) &
  (!CFG$drop_cold_items | cold_item == FALSE)]

# ---- 10a. Stratified fit / calib split -----------------------------
# Stratification is applied WITHIN the training window only, so the
# global timeline is untouched: every training row precedes every test
# row. Rank-within-stratum assignment is fully vectorised.
strat_col <- CFG$stratify_by
set.seed(CFG$seed)
train_sel[, .u := runif(.N)]
train_sel[, .r := frank(.u, ties.method = "first") / .N, by = c(strat_col)]
train_sel[, part := fifelse(.r <= CFG$calib_frac, "calib", "fit")]
train_sel[, c(".u", ".r") := NULL]

fit_set   <- train_sel[part == "fit"]
calib_set <- train_sel[part == "calib"]

cat(sprintf("fit: %s | calib: %s | test: %s\n",
            format(nrow(fit_set),   big.mark = ","),
            format(nrow(calib_set), big.mark = ","),
            format(nrow(test_sel),  big.mark = ",")))

# ---- 10b. Alternative user-disjoint calib (diagnostic only) --------
# Row-level calibration lets a user appear in both fit and calib. That
# is correct for calibrating a score the SVD produced from OTHER rows,
# but the user-disjoint version is built here so the sensitivity can be
# reported rather than assumed away.
set.seed(CFG$seed + 1L)
u_all  <- unique(train_sel$userId)
u_cal  <- sample(u_all, floor(CFG$calib_frac * length(u_all)))
train_sel[, part_userdisjoint := fifelse(userId %in% u_cal, "calib", "fit")]


# =====================================================================
# 11. Split diagnostics
# =====================================================================

# ---- 11a. Category distributions across the three parts ------------
split_dist <- rbindlist(list(
  fit_set  [, .(part = "fit",   N = .N), by = .(level = get(strat_col))],
  calib_set[, .(part = "calib", N = .N), by = .(level = get(strat_col))],
  test_sel [, .(part = "test",  N = .N), by = .(level = get(strat_col))]
))
split_dist[, share := N / sum(N), by = part]
split_dist <- dcast(split_dist, level ~ part, value.var = "share")
setnafill(split_dist, fill = 0, cols = c("fit", "calib", "test"))
write_tab(split_dist, "11a_split_category_distribution")
print(split_dist)

kl <- function(p, q) { i <- p > 0 & q > 0; sum(p[i] * log(p[i] / q[i])) }
kl_tab <- data.table(
  comparison = c("calib || fit", "test || fit", "test || calib"),
  KL = c(kl(split_dist$calib, split_dist$fit),
         kl(split_dist$test,  split_dist$fit),
         kl(split_dist$test,  split_dist$calib)),
  TV = c(0.5 * sum(abs(split_dist$calib - split_dist$fit)),
         0.5 * sum(abs(split_dist$test  - split_dist$fit)),
         0.5 * sum(abs(split_dist$test  - split_dist$calib)))
)
write_tab(kl_tab, "11b_split_divergence")
cat("\n=== Divergence between split parts ===\n"); print(kl_tab)

p_split_dist <- ggplot(melt(split_dist, id.vars = "level",
                            variable.name = "part", value.name = "share"),
                       aes(level, share, fill = part)) +
  geom_col(position = "dodge", width = 0.8) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c(fit = CFG$col_a, calib = CFG$col_b,
                               test = CFG$col_c)) +
  labs(title = "Rating distribution across fit, calibration and test",
       subtitle = "Stratification equalises fit and calib by construction; the fit-to-test gap is genuine temporal drift",
       x = "rating category", y = "share", fill = NULL)
save_fig(p_split_dist, "11c_split_category_distribution")

# ---- 11b. Per-month composition of the selected window -------------
sel_monthly <- rbind(
  train_sel[, .(part = "train", ym, rating)],
  test_sel [, .(part = "test",  ym, rating)]
)[, .(n = .N, mean_rating = mean(rating)), by = .(part, ym)][order(ym)]
write_tab(sel_monthly, "11d_selected_monthly_profile")

p_sel_month <- ggplot(sel_monthly, aes(ym, n, fill = part)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c(train = CFG$col_a, test = CFG$col_b)) +
  labs(title = "Monthly volume in the selected window after all filtering",
       subtitle = "Even months are what stop a single month from dominating the fitted thresholds",
       x = NULL, y = "ratings", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_sel_month, "11e_selected_window_monthly", h = 4)


# =====================================================================
# 12. Cold-start and coverage audit
# =====================================================================

cold_audit <- data.table(
  stage = c("test rows, raw window",
            "…from ineligible users",
            "…on unseen/ineligible movies",
            "…both",
            "test rows retained"),
  n = c(nrow(test_raw),
        test_raw[cold_user == TRUE, .N],
        test_raw[cold_item == TRUE, .N],
        test_raw[cold_user == TRUE & cold_item == TRUE, .N],
        nrow(test_sel))
)
cold_audit[, pct := 100 * n / nrow(test_raw)]
write_tab(cold_audit, "12a_cold_audit")
cat("\n=== Test-month coverage audit ===\n"); print(cold_audit)

# Comparison against the previous annual design, which is the number
# that motivated this whole redesign.
coverage_cmp <- data.table(
  design = c("v6: calendar-year split, no eligibility rules",
             sprintf("v7: %s +%dm, Nu=%d, Nm=%d", FIN$win_start,
                     FIN$win_len, FIN$Nu, FIN$Nm)),
  pct_test_warm = c(17.0, 100 * nrow(test_sel) / nrow(test_raw))
)
write_tab(coverage_cmp, "12b_coverage_vs_v6")

p_cold <- ggplot(coverage_cmp, aes(reorder(design, pct_test_warm),
                                   pct_test_warm)) +
  geom_col(fill = CFG$col_a, width = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", pct_test_warm)), hjust = -0.15) +
  coord_flip(clip = "off") + ylim(0, 115) +
  labs(title = "Share of evaluated rows the model can actually score",
       subtitle = "Under the annual design most of the test set fell back to a global mean",
       x = NULL, y = "% of test rows in the modelled population")
save_fig(p_cold, "12c_coverage_comparison", h = 3.2)

# Are the retained test users representative of the window's users?
rep_audit <- rbind(
  train_raw[, .(group = "all users in window",
                users = uniqueN(userId), mean_rating = mean(rating),
                sd_rating = sd(rating))],
  train_sel[, .(group = "eligible users",
                users = uniqueN(userId), mean_rating = mean(rating),
                sd_rating = sd(rating))],
  test_sel [, .(group = "evaluated users (test)",
                users = uniqueN(userId), mean_rating = mean(rating),
                sd_rating = sd(rating))]
)
write_tab(rep_audit, "12d_representativeness")
print(rep_audit)


# =====================================================================
# 13. Feasibility audit for the model ladder
# =====================================================================

# M2g (user random intercept) and M2h (per-user free thresholds) impose
# per-user sample-size requirements that the data selection has to meet.
# Checking this here, before the pipeline runs, is cheaper than
# discovering non-convergence downstream.

calib_per_user <- calib_set[, .(n_calib = .N), by = userId]
test_per_user  <- test_sel [, .(n_test  = .N), by = userId]

m2h_feas <- data.table(
  metric = c("users in calib",
             sprintf("…with >= %d calib rows", CFG$m2h_min_calib_per_user),
             "…with all categories present in calib",
             "median calib rows per user",
             "users in test", "median test rows per user"),
  value = c(
    nrow(calib_per_user),
    calib_per_user[n_calib >= CFG$m2h_min_calib_per_user, .N],
    calib_set[, uniqueN(get(strat_col)), by = userId][
      V1 == nlevels(calib_set[[strat_col]]), .N],
    as.numeric(median(calib_per_user$n_calib)),
    nrow(test_per_user),
    as.numeric(median(test_per_user$n_test))
  ))
write_tab(m2h_feas, "13a_m2h_feasibility")
cat("\n=== Per-user sample sizes for M2g / M2h ===\n"); print(m2h_feas)

p_m2h <- ggplot(calib_per_user, aes(n_calib)) +
  geom_histogram(bins = 50, fill = CFG$col_a, alpha = 0.85) +
  geom_vline(xintercept = CFG$m2h_min_calib_per_user, linetype = 2,
             colour = CFG$col_b) +
  scale_x_log10(labels = comma) + scale_y_continuous(labels = comma) +
  labs(title = "Calibration rows available per user",
       subtitle = "Per-user thresholds need enough rows per user; the dashed line is the minimum used for M2h",
       x = "calibration rows (log)", y = "users")
save_fig(p_m2h, "13b_calib_rows_per_user")

# Empty-category incidence, the known failure mode for per-user
# thresholds, reported under both category schemes.
empty_cat <- rbindlist(lapply(c("rating10", "rating5"), function(sc) {
  k <- nlevels(calib_set[[sc]])
  cc <- calib_set[, .(n_present = uniqueN(get(sc))), by = userId]
  data.table(scheme = sc,
             users = nrow(cc),
             pct_all_categories = 100 * cc[n_present == k, .N] / nrow(cc),
             median_categories  = as.numeric(median(cc$n_present)))
}))
write_tab(empty_cat, "13c_empty_category_incidence")
print(empty_cat)


# =====================================================================
# 14. Provenance table and export
# =====================================================================

provenance <- data.table(
  step = c("raw study period",
           "after anomalous-account removal",
           "training window rows (pre-eligibility)",
           "after user threshold",
           "after movie threshold / k-core",
           "  of which fit",
           "  of which calib",
           "test month rows (pre-eligibility)",
           "test month rows retained"),
  rows = c(nrow(raw), nrow(dat), nrow(train_raw),
           nrow(train_raw[userId %in% train_raw[, .N, by = userId][N >= FIN$Nu, userId]]),
           nrow(train_sel), nrow(fit_set), nrow(calib_set),
           nrow(test_raw), nrow(test_sel)),
  users = c(uniqueN(raw$userId), uniqueN(dat$userId), uniqueN(train_raw$userId),
            NA_integer_, uniqueN(train_sel$userId), uniqueN(fit_set$userId),
            uniqueN(calib_set$userId), uniqueN(test_raw$userId),
            uniqueN(test_sel$userId)),
  movies = c(uniqueN(raw$movieId), uniqueN(dat$movieId), uniqueN(train_raw$movieId),
             NA_integer_, uniqueN(train_sel$movieId), uniqueN(fit_set$movieId),
             uniqueN(calib_set$movieId), uniqueN(test_raw$movieId),
             uniqueN(test_sel$movieId))
)
provenance[, pct_of_raw := 100 * rows / nrow(raw)]
write_tab(provenance, "14a_provenance")
cat("\n=== Provenance ===\n"); print(provenance)

p_prov <- ggplot(provenance[!grepl("^  ", step)],
                 aes(reorder(step, -rows), rows)) +
  geom_col(fill = CFG$col_a, width = 0.65) +
  geom_text(aes(label = comma(rows)), hjust = -0.1, size = 3) +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.25))) +
  labs(title = "From raw ratings to the modelled dataset",
       subtitle = "Every reduction below is a decision that has to be defended in the write-up",
       x = NULL, y = "rows")
save_fig(p_prov, "14b_provenance_waterfall", h = 4.5)

selection <- list(
  cfg          = CFG,
  final        = FIN,
  window_scan  = window_scan,
  threshold_grid = threshold_grid,
  joint_scan   = joint_scan,
  anom_users   = anom_users,
  elig_users   = elig_users,
  elig_movies  = elig_movies,
  provenance   = provenance,
  fit          = fit_set,
  calib        = calib_set,
  test         = test_sel,
  test_raw     = test_raw,
  built_at     = Sys.time()
)
saveRDS(selection, file.path(CFG$dir_cache,
                             sprintf("selection_v7_%s_%dm_Nu%d_Nm%d.rds",
                                     FIN$win_start, FIN$win_len,
                                     FIN$Nu, FIN$Nm)))

cat("\nWritten to", CFG$dir_out, "and", CFG$dir_cache, "\n")
sessionInfo()
