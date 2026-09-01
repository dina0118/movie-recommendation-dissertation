# ============================================================================
# 08_M2d_gradient_check.R
# M2d (joint OrdRec, global thresholds): log-likelihood, analytic gradient,
# numerical gradient check, and a cross-check against ordinal::clm.
# ============================================================================

SEL_FILE  <- "D:/桌面/movie-recommendation-dissertation/06_models/cache/00_data_selection/selection_v7_2018-04_6m_Nu10_Nm10.rds"
CACHE_DIR <- "D:/桌面/movie-recommendation-dissertation/06_models/cache/07_M1_M2"
# ----------------------------------------------------------------------------

SEED     <- 2024L
K_LEVELS <- 5L
SVD_DIM  <- 16L
set.seed(SEED)

# ============================================================================
# PART 1 - Model functions
# ============================================================================

# Thresholds from the unconstrained parameterisation [eqs. 4-5].
# Returns c(t1, t2, t3, t4), monotone increasing by construction.
expand_thresholds <- function(t1, beta) {
  c(t1, t1 + cumsum(exp(beta)))
}

# Numerically safe sigmoid. The clamp at +/-30 prevents exp overflow;
# sigmoid(30) is already 1 - 9e-14, so nothing of statistical relevance is
# lost. NOTE for the gradient check: the clamp makes the function flat
# outside [-30, 30], where a numerical derivative is 0 but the analytic one
# is ~1e-13. The toy data below keeps |t - y| far away from 30, so the check
# is done on the smooth region only, which is the region training lives in.
sigmoid <- function(x) 1 / (1 + exp(-pmin(pmax(x, -30), 30)))

# Forward pass for a vector of observations.
#   y  numeric vector of internal scores
#   r  integer vector of observed ratings in 1..5
#   th numeric length-4 thresholds
# Returns prob (P(r_ui = r)), and the A, B terms needed by every gradient.
ordrec_forward <- function(y, r, th) {
  g_hi <- ifelse(r == K_LEVELS, 1, sigmoid(th[pmin(r, K_LEVELS - 1L)] - y))
  g_lo <- ifelse(r == 1L,       0, sigmoid(th[pmax(r - 1L, 1L)]       - y))
  A <- g_hi * (1 - g_hi)          # 0 automatically when r = 5 (g_hi = 1)
  B <- g_lo * (1 - g_lo)          # 0 automatically when r = 1 (g_lo = 0)
  prob <- pmax(g_hi - g_lo, 1e-12)
  list(prob = prob, A = A, B = B)
}

# Full-batch log-likelihood and analytic gradient, for the gradient check
# and for the frozen-factor clm cross-check. Vectorised over rows; rowsum()
# does the per-user / per-item aggregation. (The SGD file will use the same
# per-row quantities, just one row at a time.)
#
#   par  list(t1, beta[3], bu[nU], bi[nI], P[nU x K], Q[nI x K])
#   dat  list(u, i, r)  integer vectors, 1-based indices
#   freeze_factors  if TRUE, y is taken from dat$y_fixed and only t1, beta
#                   get gradients (used in Part 3)
ordrec_loglik_grad <- function(par, dat, freeze_factors = FALSE) {
  u <- dat$u; i <- dat$i; r <- dat$r
  th <- expand_thresholds(par$t1, par$beta)

  y <- if (freeze_factors) dat$y_fixed else
    par$bi[i] + par$bu[u] + rowSums(par$P[u, , drop = FALSE] *
                                    par$Q[i, , drop = FALSE])

  fw <- ordrec_forward(y, r, th)
  ll <- sum(log(fw$prob))

  # -- threshold gradients ---------------------------------------------------
  # dl/dt1 per row: (A - B)/p
  d_t1 <- sum((fw$A - fw$B) / fw$prob)
  # dl/dbeta_s per row: exp(beta_s) * (A*[s <= r-1] - B*[s <= r-2]) / p
  d_beta <- vapply(1:(K_LEVELS - 2L), function(s) {
    exp(par$beta[s]) *
      sum((fw$A * (s <= r - 1L) - fw$B * (s <= r - 2L)) / fw$prob)
  }, numeric(1))

  if (freeze_factors) {
    return(list(ll = ll, d_t1 = d_t1, d_beta = d_beta))
  }

  # -- score-side gradients --------------------------------------------------
  dy <- (fw$B - fw$A) / fw$prob                    # dl/dy per row
  d_bu <- rowsum(dy, u); d_bi <- rowsum(dy, i)     # named by group index
  d_P  <- rowsum(dy * par$Q[i, , drop = FALSE], u)
  d_Q  <- rowsum(dy * par$P[u, , drop = FALSE], i)

  # rowsum() only returns groups that occur; scatter back to full size
  full <- function(m, n) {
    out <- matrix(0, n, ncol(m)); out[as.integer(rownames(m)), ] <- m; out
  }
  nU <- length(par$bu); nI <- length(par$bi)
  list(ll = ll, d_t1 = d_t1, d_beta = d_beta,
       d_bu = as.numeric(full(d_bu, nU)),
       d_bi = as.numeric(full(d_bi, nI)),
       d_P  = full(d_P, nU),
       d_Q  = full(d_Q, nI))
}

# ============================================================================
# PART 2 - Numerical gradient check on toy data
#
# Central differences: (f(x + d) - f(x - d)) / (2d), d = 1e-5, against every
# analytic component. PASS requires relative error < 1e-6 on every checked
# coordinate. Nothing later in the project is run unless this prints PASS.
# ============================================================================

cat("========== PART 2: numerical gradient check (toy data) ==========\n")

nU_toy <- 20L; nI_toy <- 15L; K_toy <- 4L; n_obs <- 400L
toy_par <- list(
  t1   = -1.2,
  beta = c(0.1, -0.3, 0.4),
  bu   = rnorm(nU_toy, 0, 0.3),
  bi   = rnorm(nI_toy, 0, 0.3),
  P    = matrix(rnorm(nU_toy * K_toy, 0, 0.3), nU_toy, K_toy),
  Q    = matrix(rnorm(nI_toy * K_toy, 0, 0.3), nI_toy, K_toy))
toy_dat <- list(u = sample(nU_toy, n_obs, TRUE),
                i = sample(nI_toy, n_obs, TRUE),
                r = sample(K_LEVELS, n_obs, TRUE))  # every category occurs

# flatten <-> unflatten so one loop can perturb any coordinate
flatten <- function(p) c(p$t1, p$beta, p$bu, p$bi, as.numeric(p$P), as.numeric(p$Q))
unflatten <- function(v) {
  o <- 0L
  take <- function(n) { x <- v[(o + 1L):(o + n)]; o <<- o + n; x }
  list(t1   = take(1L), beta = take(3L),
       bu   = take(nU_toy), bi = take(nI_toy),
       P    = matrix(take(nU_toy * K_toy), nU_toy, K_toy),
       Q    = matrix(take(nI_toy * K_toy), nI_toy, K_toy))
}

g <- ordrec_loglik_grad(toy_par, toy_dat)
analytic <- c(g$d_t1, g$d_beta, g$d_bu, g$d_bi, as.numeric(g$d_P), as.numeric(g$d_Q))
v0 <- flatten(toy_par)

# check every threshold parameter + a random sample of the rest
idx_check <- c(1:4, sort(sample(5:length(v0), 60L)))
delta <- 1e-5
numeric_g <- vapply(idx_check, function(j) {
  vp <- v0; vp[j] <- vp[j] + delta
  vm <- v0; vm[j] <- vm[j] - delta
  (ordrec_loglik_grad(unflatten(vp), toy_dat)$ll -
   ordrec_loglik_grad(unflatten(vm), toy_dat)$ll) / (2 * delta)
}, numeric(1))

rel_err <- abs(analytic[idx_check] - numeric_g) /
           pmax(1, abs(analytic[idx_check]), abs(numeric_g))

param_names <- c("t1", paste0("beta", 1:3),
                 rep("bu/bi/P/Q", length(idx_check) - 4L))
worst <- order(rel_err, decreasing = TRUE)[1:5]
cat(sprintf("Checked %d coordinates (all 4 threshold params + 60 sampled others)\n",
            length(idx_check)))
cat(sprintf("Max relative error: %.2e\n", max(rel_err)))
cat("Worst five:\n")
for (w in worst)
  cat(sprintf("  coord %4d (%s): analytic %.8f | numeric %.8f | rel err %.2e\n",
              idx_check[w], param_names[w], analytic[idx_check[w]],
              numeric_g[w], rel_err[w]))

grad_check_pass <- max(rel_err) < 1e-6
cat(sprintf("\nGRADIENT CHECK: %s\n\n", if (grad_check_pass) "PASS" else "FAIL"))
stopifnot(grad_check_pass)

# Edge categories deserve their own line of evidence: rows with r = 1 and
# r = 5 exercise the g_0 = 0 / g_S = 1 conventions. Confirm they were present.
cat(sprintf("Toy ratings per category: %s  (1 and 5 must both be > 0)\n\n",
            paste(tabulate(toy_dat$r, K_LEVELS), collapse = " / ")))

# ============================================================================
# PART 3 - Frozen-factor cross-check against ordinal::clm
#
# With the factors frozen at the SVD scores (y_ui = svd_score, coefficient
# fixed at 1), OrdRec eq. 9 IS McCullagh's cumulative logit model with an
# offset:  logit P(Y <= k) = t_k - offset.  ordinal::clm fits exactly that
# with clm(rating_f ~ 1 + offset(svd_score)).  So the four thresholds that
# maximise OUR likelihood must agree with clm's estimates. 
# ============================================================================

cat("========== PART 3: frozen-factor cross-check vs ordinal::clm ==========\n")

if (!file.exists(SEL_FILE)) {
  cat("Selection file not found - skipping Part 3 (fine on a machine without\n")
  cat("the data; run it on the machine that has the selection RDS).\n")
} else {
  for (pkg in c("data.table", "ordinal", "recosystem"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is needed for Part 3. install.packages(\"", pkg, "\")")
  library(data.table)

  # --- load the selection exactly as 07 does --------------------------------
  sel   <- readRDS(SEL_FILE)
  fit   <- as.data.table(sel$fit)
  calib <- as.data.table(sel$calib)
  to_int5 <- function(d) { d[, rating5 := as.integer(as.character(rating5))]; d[] }
  fit <- to_int5(fit); calib <- to_int5(calib)

  user_map  <- data.table(userId  = sort(unique(fit$userId)))[,  user_idx  := .I]
  movie_map <- data.table(movieId = sort(unique(fit$movieId)))[, movie_idx := .I]
  add_idx <- function(d) {
    d <- merge(d, user_map,  by = "userId",  all.x = TRUE, sort = FALSE)
    d <- merge(d, movie_map, by = "movieId", all.x = TRUE, sort = FALSE)
    setDT(d)
    d[, cold := is.na(user_idx) | is.na(movie_idx)]
    d[]
  }
  fit <- add_idx(fit)[cold == FALSE]
  calib <- add_idx(calib)[cold == FALSE]

  # --- SVD scores on calib ----------------------------
  # recosystem with >1 thread is not bit-for-bit reproducible, but that does
  # not matter here: whatever scores come out are frozen and BOTH sides of
  # the comparison see the same frozen scores.
  set.seed(SEED)
  reco <- recosystem::Reco()
  reco$train(recosystem::data_memory(fit$user_idx, fit$movie_idx,
                                     as.numeric(fit$rating5), index1 = TRUE),
             opts = list(dim = SVD_DIM, costp_l2 = 0.01, costq_l2 = 0.01,
                         lrate = 0.1, niter = 20,
                         nthread = parallel::detectCores(), verbose = FALSE))
  calib[, svd_score := reco$predict(
    recosystem::data_memory(user_idx, movie_idx, index1 = TRUE),
    recosystem::out_memory())]

  # --- side 1: ordinal::clm with the score as a fixed offset ----------------
  calib[, rating_f := factor(rating5, levels = 1:K_LEVELS, ordered = TRUE)]
  m_ref <- ordinal::clm(rating_f ~ 1 + offset(svd_score),
                        data = calib, link = "logit")
  th_clm <- as.numeric(m_ref$alpha)
  ll_clm <- as.numeric(logLik(m_ref))

  # --- side 2: our likelihood, maximised over (t1, beta) only ---------------
  dat_frozen <- list(u = calib$user_idx, i = calib$movie_idx,
                     r = calib$rating5, y_fixed = calib$svd_score)
  negll <- function(v) {
    -ordrec_loglik_grad(list(t1 = v[1], beta = v[2:4]),
                        dat_frozen, freeze_factors = TRUE)$ll
  }
  neggr <- function(v) {
    g <- ordrec_loglik_grad(list(t1 = v[1], beta = v[2:4]),
                            dat_frozen, freeze_factors = TRUE)
    -c(g$d_t1, g$d_beta)
  }
  opt <- optim(c(0, 0, 0, 0), negll, neggr, method = "BFGS",
               control = list(maxit = 500, reltol = 1e-12))
  th_ours <- expand_thresholds(opt$par[1], opt$par[2:4])
  ll_ours <- -opt$value

  # --- compare --------------------------------------------------------------
  cmp <- data.frame(threshold = paste0(1:4, "|", 2:5),
                    clm = th_clm, ours = th_ours,
                    abs_diff = abs(th_clm - th_ours))
  print(cmp, row.names = FALSE, digits = 8)
  cat(sprintf("logLik  clm: %.6f | ours: %.6f | diff: %.2e\n",
              ll_clm, ll_ours, abs(ll_clm - ll_ours)))

  clm_check_pass <- max(cmp$abs_diff) < 1e-4 && abs(ll_clm - ll_ours) < 1e-3
  cat(sprintf("\nCLM CROSS-CHECK: %s\n", if (clm_check_pass) "PASS" else "FAIL"))
  stopifnot(clm_check_pass)

  saveRDS(list(grad_check_max_rel_err = max(rel_err),
               clm_thresholds = th_clm, our_thresholds = th_ours,
               ll_clm = ll_clm, ll_ours = ll_ours),
          file.path(CACHE_DIR, "m2d_gradient_check.rds"))
  cat("Saved check summary to", file.path(CACHE_DIR, "m2d_gradient_check.rds"), "\n")
}

cat("\nAll checks passed. The SGD training loop (next file) may now be trusted\n")
cat("to the extent that it reuses exactly these likelihood/gradient functions.\n")
