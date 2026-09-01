# =====================================================================
# 27_selection_v9.R
# Phase C, step 1: rebuild the selection at a higher USER threshold.
#
# Nu 10 -> 30. Nm stays at 10, the window, anomaly exclusions, pruning
# mode and split rule stay at v7's values. This is the mirror image of
# Phase B: there the item tail was extended, here the user tail is cut.
# Output: cache/00_data_selection/selection_v9_2018-04_6m_Nu30_Nm10.rds
# =====================================================================

NU_NEW <- 30L         # the new user threshold
NM     <- 10L         # unchanged from v7

BASE    <- "D:/桌面/movie-recommendation-dissertation"
SEL_DIR <- file.path(BASE, "06_models/cache/00_data_selection")
SEL_V7  <- file.path(SEL_DIR, "selection_v7_2018-04_6m_Nu10_Nm10.rds")
OUT_DIR <- file.path(BASE, "06_models/output/09_phaseC")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
CUTS <- c(1.5, 2.5, 3.5, 4.5)


# =====================================================================
# 1. Reassemble the window from v7
# =====================================================================
stopifnot("v7 selection not found" = file.exists(SEL_V7))
v7  <- readRDS(SEL_V7)
CFG <- v7$cfg; FIN <- v7$final
stopifnot(FIN$Nu == 10L, FIN$Nm == NM, NU_NEW > FIN$Nu)

month_index <- function(ym) {
  y <- as.integer(substr(ym, 1, 4)); m <- as.integer(substr(ym, 6, 7))
  y * 12L + (m - 1L)
}

train_v7 <- rbind(as.data.table(v7$fit), as.data.table(v7$calib), fill = TRUE)
test_raw <- as.data.table(v7$test_raw)
for (d in list(train_v7, test_raw)) {
  if (!"mi" %in% names(d))      d[, mi := month_index(ym)]
  if (!"rating5" %in% names(d)) d[, rating5 := findInterval(rating, CUTS) + 1L]
}
train_v7[, rating5 := as.integer(as.character(rating5))]
test_raw[, rating5 := as.integer(as.character(rating5))]

cat(sprintf("v7 training rows %s | users %d | movies %d\n",
            format(nrow(train_v7), big.mark = ","),
            uniqueN(train_v7$userId), uniqueN(train_v7$movieId)))


# =====================================================================
# 2. The pipeline, as a function of Nu
# =====================================================================
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

build_selection <- function(Nu) {
  pr <- apply_prune(train_v7, Nu, NM)
  train_sel <- copy(pr$data)
  elig_users  <- unique(train_sel$userId)
  elig_movies <- unique(train_sel$movieId)

  ts <- copy(test_raw)
  ts[, cold_user := !userId  %in% elig_users]
  ts[, cold_item := !movieId %in% elig_movies]
  test_sel <- ts[(!CFG$drop_cold_users | cold_user == FALSE) &
                 (!CFG$drop_cold_items | cold_item == FALSE)]

  set.seed(CFG$seed)
  train_sel[, .u := runif(.N)]
  train_sel[, .r := frank(.u, ties.method = "first") / .N, by = c(CFG$stratify_by)]
  train_sel[, part := fifelse(.r <= CFG$calib_frac, "calib", "fit")]
  train_sel[, c(".u", ".r") := NULL]

  list(cfg = CFG,
       final = list(win_start = FIN$win_start, win_len = FIN$win_len,
                    Nu = Nu, Nm = NM),
       anom_users = v7$anom_users,
       elig_users = elig_users, elig_movies = elig_movies,
       fit = train_sel[part == "fit"], calib = train_sel[part == "calib"],
       test = test_sel, test_raw = ts,
       prune_iters = pr$iters, built_at = Sys.time())
}


# =====================================================================
# 3. Guard: re-pruning at (10, 10) must be a no-op
# =====================================================================
cat("\n=== Idempotence check at Nu = 10 ===\n")
chk <- build_selection(10L)
key_of <- function(d) sort(paste(d$userId, d$movieId, d$timestamp, sep = "_"))
same_train <- identical(key_of(rbind(chk$fit, chk$calib, fill = TRUE)),
                        key_of(train_v7))
same_test  <- identical(key_of(chk$test), key_of(as.data.table(v7$test)))
cat(sprintf("training set unchanged: %s | test set matches v7: %s\n",
            same_train, same_test))
stopifnot("Re-pruning v7 at its own thresholds changed it - investigate" =
            same_train, same_test)
cat("IDEMPOTENCE CHECK: PASS (see header for what this does not prove)\n")


# =====================================================================
# 4. Build v9
# =====================================================================
cat(sprintf("\n=== Building v9 at Nu = %d ===\n", NU_NEW))
v9 <- build_selection(NU_NEW)
out_file <- file.path(SEL_DIR, sprintf("selection_v9_%s_%dm_Nu%d_Nm%d.rds",
                                       FIN$win_start, FIN$win_len, NU_NEW, NM))
saveRDS(v9, out_file)
cat(sprintf("Saved %s\n", out_file))

# The item set is expected to move even though Nm is fixed: dropping
# low-activity users removes ratings, which can push films below Nm and
# cascade. If prune_iters > 1 that cascade actually happened.
cat(sprintf("Pruning iterations: %d\n", v9$prune_iters))


# =====================================================================
# 5. What changed, and is there enough left to test anything
# =====================================================================
summarise <- function(s, label) {
  tr <- rbind(s$fit, s$calib, fill = TRUE)
  data.table(selection = label, Nu = s$final$Nu,
             fit_rows = nrow(s$fit), calib_rows = nrow(s$calib),
             users = uniqueN(tr$userId), movies = uniqueN(tr$movieId),
             median_user_n = as.numeric(median(tr[, .N, by = userId]$N)),
             density = nrow(tr) / (uniqueN(tr$userId) * uniqueN(tr$movieId)),
             test_rows = nrow(s$test), test_users = uniqueN(s$test$userId),
             median_test_n = as.numeric(median(s$test[, .N, by = userId]$N)))
}
cmp <- rbind(summarise(chk, "v7 (Nu=10)"), summarise(v9, sprintf("v9 (Nu=%d)", NU_NEW)))
cat("\n=== v7 vs v9 ===\n"); print(cmp)
fwrite(cmp, file.path(OUT_DIR, "pc1_v7_vs_v9_selection.csv"))

# POWER. This is the number that decides whether Phase C can answer its
# question. The ladder increments are small (M2d->M2e was ~0.01 nats on
# v7); resolving them needs users, and Nu=30 spends users to buy
# per-user data. Check before training, not after.
cat(sprintf("\nTest users: v7 %d -> v9 %d (%.0f%% retained)\n",
            cmp$test_users[1], cmp$test_users[2],
            100 * cmp$test_users[2] / cmp$test_users[1]))
cat(sprintf("Test rows:  v7 %s -> v9 %s (%.0f%% retained)\n",
            format(cmp$test_rows[1], big.mark = ","),
            format(cmp$test_rows[2], big.mark = ","),
            100 * cmp$test_rows[2] / cmp$test_rows[1]))
if (cmp$test_users[2] < 300)
  warning("Fewer than 300 test users. The two-way bootstrap resamples ",
          "users, so its intervals will be wide enough that small ladder ",
          "increments become untestable. Phase C may only be able to ",
          "confirm the ORDERING, not the individual increments.")

# The mechanism Phase C is really about: per-user thresholds need
# per-user data, and this is how much more of it v9 users have.
ua <- rbind(chk$fit, chk$calib, fill = TRUE)[, .N, by = userId][, sel := "v7 (Nu=10)"]
ub <- rbind(v9$fit,  v9$calib,  fill = TRUE)[, .N, by = userId][, sel := "v9 (Nu=30)"]
p_u <- ggplot(rbind(ua, ub), aes(N, fill = sel)) +
  geom_density(alpha = 0.45, colour = NA) + scale_x_log10() + theme_minimal() +
  labs(title = "Training ratings per user",
       subtitle = "Phase C buys per-user data at the cost of users. Both axes matter.",
       x = "ratings per user (log)", y = "density", fill = NULL)
ggsave(file.path(OUT_DIR, "pc2_user_activity_shift.png"), p_u,
       width = 8, height = 4, dpi = 150)

# Calibration rows per user drive the per-user threshold estimates
# directly, so report them separately from training volume.
cu <- v9$calib[, .N, by = userId]$N
cat(sprintf("\nCalibration rows per user (v9): min %d | median %.0f | mean %.1f\n",
            min(cu), median(cu), mean(cu)))
cat(sprintf("Users with < 4 calibration rows: %d (%.1f%%) - four thresholds",
            sum(cu < 4), 100 * mean(cu < 4)),
    "cannot be identified from fewer than four observations.\n")

cat("\n27 COMPLETE.\n")
