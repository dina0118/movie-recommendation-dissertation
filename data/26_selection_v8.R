# =====================================================================
# 26_selection_v8.R
# Phase B, step 1: rebuild the selection at a lower movie threshold.
#
# Only Nm changes. The window (2018-04 + 6m, test 2018-10), the user
# threshold (Nu = 10), the anomaly exclusions, the pruning mode and the
# fit/calib split are all held at v7's values, so the user population is
# unchanged and the only manipulated quantity is the length of the item
# tail.
#
# SELF-VALIDATION. The same function is first run at Nm = 10 and its
# output compared row-for-row against the stored v7 selection. If that
# comparison fails, the rebuild is not faithful and the Nm = 1 output
# must not be used. This is the whole reason the script is structured
# this way rather than just pruning harder.
#
# Output: cache/00_data_selection/selection_v8_2018-04_6m_Nu10_Nm1.rds
#         (same list structure as v7, so 05c/06c can read either)
#
# Runtime: ~2-4 minutes.
# =====================================================================

NM_NEW <- 1L          # the new movie threshold
NU     <- 10L         # unchanged from v7
SEED   <- 2024L

# --- paths -----------------------------------------------------------
BASE      <- "D:/桌面/movie-recommendation-dissertation"
SEL_DIR   <- file.path(BASE, "/06_models/cache/00_data_selection")
SEL_V7    <- file.path(SEL_DIR, "/selection_v7_2018-04_6m_Nu10_Nm10.rds")
RAW_CACHE <- file.path(SEL_DIR, "/raw_study_period.rds")
OUT_DIR   <- file.path(BASE, "/06_models/output/08_M3")
# ---------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

CUTS <- c(1.5, 2.5, 3.5, 4.5)


# =====================================================================
# 1. Inputs
# =====================================================================
stopifnot("raw_study_period.rds not found - run 00 first" = file.exists(RAW_CACHE),
          "v7 selection not found" = file.exists(SEL_V7))

raw <- readRDS(RAW_CACHE)
v7  <- readRDS(SEL_V7)
setDT(raw)

CFG <- v7$cfg; FIN <- v7$final
cat(sprintf("v7 config: window %s +%dm | Nu=%d | Nm=%d | pruning=%s | calib_frac=%.2f\n",
            FIN$win_start, FIN$win_len, FIN$Nu, FIN$Nm,
            CFG$prune_mode, CFG$calib_frac))
stopifnot(FIN$Nu == NU)

month_index <- function(ym) {
  y <- as.integer(substr(ym, 1, 4)); m <- as.integer(substr(ym, 6, 7))
  y * 12L + (m - 1L)
}
if (!"mi" %in% names(raw)) raw[, mi := month_index(ym)]

# The 5-category collapse. Derived here rather than assumed present;
# the v7 comparison in Section 3 checks it is the same mapping.
if (!"rating5" %in% names(raw))
  raw[, rating5 := findInterval(rating, CUTS) + 1L]

dat <- raw[!userId %in% v7$anom_users]
cat(sprintf("raw %s -> after anomaly exclusion %s rows\n",
            format(nrow(raw), big.mark = ","),
            format(nrow(dat), big.mark = ",")))


# =====================================================================
# 2. The pipeline, as a function of Nm
# =====================================================================
# Verbatim logic from 00 sections 8 and 10: iterative k-core pruning,
# then the same rank-within-stratum fit/calib assignment under the same
# seed. Nothing here is new; it is factored so the Nm = 10 call can be
# checked against v7.

apply_prune <- function(train, Nu, Nm, mode = CFG$prune_mode,
                        max_iter = CFG$prune_max_iter) {
  d <- train; iter <- 0L
  repeat {
    iter <- iter + 1L
    uc <- d[, .N, by = userId][N >= Nu, userId]
    mc <- d[, .N, by = movieId][N >= Nm, movieId]
    d2 <- d[userId %in% uc & movieId %in% mc]
    if (mode == "single" || nrow(d2) == nrow(d) || iter >= max_iter) { d <- d2; break }
    d <- d2
  }
  list(data = d, iters = iter)
}

build_selection <- function(Nm) {
  s <- month_index(FIN$win_start)
  train_raw <- dat[mi >= s & mi < s + FIN$win_len]
  test_raw  <- dat[mi >= s + FIN$win_len & mi < s + FIN$win_len + CFG$test_length]

  pr <- apply_prune(train_raw, NU, Nm)
  train_sel <- copy(pr$data)
  elig_users  <- unique(train_sel$userId)
  elig_movies <- unique(train_sel$movieId)

  test_raw <- copy(test_raw)
  test_raw[, cold_user := !userId  %in% elig_users]
  test_raw[, cold_item := !movieId %in% elig_movies]
  test_sel <- test_raw[
    (!CFG$drop_cold_users | cold_user == FALSE) &
    (!CFG$drop_cold_items | cold_item == FALSE)]

  # Stratified fit/calib split, same seed and same rank-within-stratum
  # rule as 00 section 10a.
  set.seed(CFG$seed)
  train_sel[, .u := runif(.N)]
  train_sel[, .r := frank(.u, ties.method = "first") / .N,
            by = c(CFG$stratify_by)]
  train_sel[, part := fifelse(.r <= CFG$calib_frac, "calib", "fit")]
  train_sel[, c(".u", ".r") := NULL]

  list(cfg = CFG,
       final = list(win_start = FIN$win_start, win_len = FIN$win_len,
                    Nu = NU, Nm = Nm),
       anom_users  = v7$anom_users,
       elig_users  = elig_users,
       elig_movies = elig_movies,
       fit   = train_sel[part == "fit"],
       calib = train_sel[part == "calib"],
       test  = test_sel,
       test_raw = test_raw,
       prune_iters = pr$iters,
       built_at = Sys.time())
}


# =====================================================================
# 3. SELF-VALIDATION: rebuild v7 and compare
# =====================================================================
cat("\n=== Reproducing v7 at Nm = 10 ===\n")
chk <- build_selection(10L)

key_of <- function(d) sort(paste(d$userId, d$movieId, d$timestamp, sep = "_"))
ok_fit   <- identical(key_of(chk$fit),   key_of(as.data.table(v7$fit)))
ok_calib <- identical(key_of(chk$calib), key_of(as.data.table(v7$calib)))
ok_test  <- identical(key_of(chk$test),  key_of(as.data.table(v7$test)))

cat(sprintf("fit   rows %d vs %d | identical: %s\n",
            nrow(chk$fit), nrow(v7$fit), ok_fit))
cat(sprintf("calib rows %d vs %d | identical: %s\n",
            nrow(chk$calib), nrow(v7$calib), ok_calib))
cat(sprintf("test  rows %d vs %d | identical: %s\n",
            nrow(chk$test), nrow(v7$test), ok_test))

# rating5 mapping must also agree, not just the row sets
d5_new <- chk$fit[, .N, by = rating5][order(rating5)]$N
d5_old <- as.data.table(v7$fit)[, .N, by = .(rating5 = as.integer(as.character(rating5)))][order(rating5)]$N
cat("rating5 counts new:", d5_new, "\n")
cat("rating5 counts v7 :", d5_old, "\n")

stopifnot("Rebuild does not reproduce v7 - do NOT use the v8 output" =
            all(ok_fit, ok_calib, ok_test),
          identical(d5_new, d5_old))
cat("REPRODUCTION CHECK: PASS\n")


# =====================================================================
# 4. Build v8
# =====================================================================
cat(sprintf("\n=== Building v8 at Nm = %d ===\n", NM_NEW))
v8 <- build_selection(NM_NEW)

out_file <- file.path(SEL_DIR,
                      sprintf("selection_v8_%s_%dm_Nu%d_Nm%d.rds",
                              FIN$win_start, FIN$win_len, NU, NM_NEW))
saveRDS(v8, out_file)
cat(sprintf("Saved %s\n", out_file))


# =====================================================================
# 5. What actually changed
# =====================================================================
summarise <- function(s, label) {
  tr <- rbind(s$fit[, .(movieId, userId)], s$calib[, .(movieId, userId)])
  data.table(
    selection    = label,
    Nm           = s$final$Nm,
    prune_iters  = s$prune_iters,
    fit_rows     = nrow(s$fit),
    calib_rows   = nrow(s$calib),
    train_rows   = nrow(tr),
    users        = uniqueN(tr$userId),
    movies       = uniqueN(tr$movieId),
    density      = nrow(tr) / (uniqueN(tr$userId) * uniqueN(tr$movieId)),
    test_rows    = nrow(s$test),
    test_users   = uniqueN(s$test$userId),
    test_movies  = uniqueN(s$test$movieId),
    pct_test_kept = 100 * nrow(s$test) / nrow(s$test_raw))
}
cmp <- rbind(summarise(chk, "v7 (Nm=10)"), summarise(v8, sprintf("v8 (Nm=%d)", NM_NEW)))
cat("\n=== v7 vs v8 ===\n"); print(cmp)
fwrite(cmp, file.path(OUT_DIR, "m3g_v7_vs_v8_selection.csv"))

# The user population must be untouched: that is the point of holding
# Nu fixed. If it moved, iterative pruning cascaded and the comparison
# is no longer a clean item-side manipulation.
u7 <- sort(unique(c(chk$fit$userId, chk$calib$userId)))
u8 <- sort(unique(c(v8$fit$userId,  v8$calib$userId)))
cat(sprintf("\nUsers v7 %d | v8 %d | v7 subset of v8: %s | identical: %s\n",
            length(u7), length(u8), all(u7 %in% u8), identical(u7, u8)))
if (!identical(u7, u8))
  cat("NOTE: the user set changed. With Nm lowered this can only ADD users",
      "(more of their ratings survive), which is expected under iterative",
      "pruning; report it rather than treating the comparison as user-fixed.\n")

# ---- item support distribution, the variable Phase B is about --------
sup8 <- rbind(v8$fit[, .(movieId)], v8$calib[, .(movieId)])[, .(n_train = .N), by = movieId]
sup8[, band := cut(n_train, c(0, 2, 5, 9, 20, 50, 200, Inf),
                   labels = c("1-2", "3-5", "6-9", "10-20", "21-50", "51-200", "200+"))]

test8 <- merge(v8$test[, .(movieId)], sup8, by = "movieId", all.x = TRUE)
band_tbl <- merge(
  sup8[, .(films = .N, train_ratings = sum(n_train)), by = band],
  test8[, .(test_rows = .N), by = band], by = "band", all = TRUE)[order(band)]
cat("\n=== v8 support strata ===\n"); print(band_tbl)
fwrite(band_tbl, file.path(OUT_DIR, "m3h_v8_support_strata.csv"))

# The three new bands are the entire point of Phase B. If they carry too
# few test rows the stratified contrast will have no power there, and
# that has to be known before the models are trained, not after.
new_rows <- band_tbl[band %in% c("1-2", "3-5", "6-9"), sum(test_rows, na.rm = TRUE)]
cat(sprintf("\nTest rows in the three NEW bands (1-9 ratings): %s\n",
            format(new_rows, big.mark = ",")))
if (is.na(new_rows) || new_rows < 500)
  warning("Fewer than 500 test rows below 10 ratings. Phase B will be ",
          "underpowered exactly where it matters; consider a longer test ",
          "window or accept that the sub-10 bands are descriptive only.")

p_sup <- ggplot(sup8, aes(n_train)) +
  geom_histogram(bins = 60, fill = "steelblue", alpha = 0.85) +
  geom_vline(xintercept = 10, linetype = 2, colour = "coral") +
  scale_x_log10() + theme_minimal() +
  labs(title = sprintf("Item support under v8 (Nm = %d)", NM_NEW),
       subtitle = "Dashed line is v7's threshold: everything left of it is newly visible",
       x = "training ratings per film (log)", y = "films")
ggsave(file.path(OUT_DIR, "m3i_v8_support_hist.png"), p_sup,
       width = 8, height = 4, dpi = 150)

cat("\n26 COMPLETE.\n")
