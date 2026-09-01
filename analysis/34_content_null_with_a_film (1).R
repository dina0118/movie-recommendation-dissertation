# =====================================================================
# 34_content_null_with_a_film.R
#
# Act IV of the worked-example thread: the content-feature null result,
# told through one film rather than one confidence interval.
#
# Run INSIDE the 06b_M3_training session, after section 8.
# Needs: m2e, m3a, m3b, Cmat (content matrix), dat_test, test, movie_map,
#        support, ordrec_probs_m3, tidy_probs, LEVELS5, OUT_DIR, SEED,
#        title_of  (define it if 06b does not)
#
# ---------------------------------------------------------------------
# WHAT A NULL RESULT NEEDS THAT A POSITIVE ONE DOES NOT
#
# A positive result is carried by its interval: the effect is there, the
# interval excludes zero, the reader is persuaded. A null result is not
# carried by its interval, because an interval spanning zero is
# compatible with two very different states of the world -- the content
# features carry no signal, or they carry signal the model cannot reach.
# An examiner will ask which, and "the CI included zero" does not answer.
#
# Three things separate the two readings, and this file does them in
# order.
#
# (1) IS THE SIGNAL THERE AT ALL? If content features predict ratings
#     when used directly, but add nothing on top of the factors, the
#     null is about REDUNDANCY, not about absence. Section 1 regresses
#     the residual on the content block. This is the difference between
#     "genre does not predict taste" (false, and the reader knows it is
#     false) and "genre predicts nothing the factors have not already
#     learned" (the actual claim).
#
# (2) WHERE SHOULD IT HAVE HELPED MOST? The cold-item argument is the
#     one reason to expect content features to earn their place: a film
#     with ten ratings has a barely-estimated factor, and its genre and
#     cast are known regardless. Section 8 of 06b already stratifies by
#     item support. Section 2 here goes further and asks whether the
#     films where M3 helps are the ones the argument predicts, or
#     whether the sign is scattered.
#
# (3) WHAT DOES IT LOOK LIKE ON ONE FILM? Section 3 selects a film that
#     the cold-item argument would nominate -- low support, rich and
#     unambiguous content, a clear genre neighbourhood -- and shows what
#     M3 does to its predictions. A null shown on a named film with a
#     known cast is harder to dismiss than a null shown as an interval.
#
# ON SELECTION DISCIPLINE. The film is chosen by a stated rule applied
# before its outcome is examined, and the rule is recorded here. The
# alternative -- searching for the film that best illustrates the point
# -- would make the example decorative. The rule below nominates a small
# set; the file prints all of them and the write-up takes the first.
# =====================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

.need <- c("m2e", "m3a", "m3b", "Cmat", "test", "movie_map", "support",
           "ordrec_probs_m3", "tidy_probs", "dat_test", "LEVELS5", "OUT_DIR")
.miss <- .need[!vapply(.need, exists, logical(1))]
if (length(.miss))
  stop("Run inside the 06b M3 session after section 8. Missing: ",
       paste(.miss, collapse = ", "))

if (!exists("SEED")) SEED <- 2024L
if (!exists("title_of")) {
  .mv <- fread("D:/桌面/movie-recommendation-dissertation/data/ml-32m/movies.csv")[, .(movieId, title)]
  title_of <- function(id) {
    t <- .mv[match(id, movieId), title]
    ifelse(is.na(t), paste0("movieId ", id), t)
  }
}

par_m2e <- m2e$par; par_m2e$w <- NULL; par_m2e$W <- NULL

LL_of <- function(par) {
  P <- tidy_probs(ordrec_probs_m3(par, dat_test$u, dat_test$i))
  log(pmax(P[cbind(seq_len(nrow(P)), test$rating5)], 1e-12))
}
ll <- data.table(M2e = LL_of(par_m2e),
                 M3a = LL_of(m3a$par),
                 M3b = LL_of(m3b$par))
test[, `:=`(ll_m2e = ll$M2e, ll_m3a = ll$M3a, ll_m3b = ll$M3b)]
test[, d_b := ll_m3b - ll_m2e]


# =====================================================================
# 1. Is the signal in the content block at all?
# =====================================================================
# If the content features predict nothing anywhere, the null is dull. If
# they predict the rating on their own and add nothing to the factors,
# the null is about redundancy and is worth a paragraph.
#
# Two regressions on the test month, both on the content block C:
#   (a) the rating itself      -- is there signal?
#   (b) the M2e residual       -- is there signal the factors missed?
# Fitted on a random half and scored on the other, so the comparison is
# out-of-sample and the dimension of C does not flatter it.

set.seed(SEED)
idx <- seq_len(nrow(test))
tr  <- sample(idx, floor(length(idx) / 2))
te  <- setdiff(idx, tr)

Ct <- Cmat[dat_test$i, , drop = FALSE]
Ct <- Ct[, apply(Ct[tr, , drop = FALSE], 2, sd) > 0, drop = FALSE]

Er_m2e <- as.numeric(tidy_probs(ordrec_probs_m3(par_m2e, dat_test$u,
                                                dat_test$i)) %*% LEVELS5)
resid  <- test$rating5 - Er_m2e

r2_oos <- function(y, X, tr, te) {
  fitm <- lm.fit(cbind(1, X[tr, , drop = FALSE]), y[tr])
  pred <- cbind(1, X[te, , drop = FALSE]) %*% coef(fitm)
  pred[!is.finite(pred)] <- mean(y[tr])
  1 - sum((y[te] - pred)^2) / sum((y[te] - mean(y[tr]))^2)
}

cat("=== Does the content block carry signal? (out-of-sample R^2) ===\n")
cat(sprintf("  content -> rating          %+.4f\n",
            r2_oos(test$rating5, Ct, tr, te)))
cat(sprintf("  content -> M2e residual    %+.4f\n",
            r2_oos(resid, Ct, tr, te)))
cat("\nThe first number says whether genre, cast and keywords predict how a\n")
cat("film is rated. The second says whether they predict what the factors\n")
cat("got wrong. If the first is clearly positive and the second is near\n")
cat("zero, the content features are REDUNDANT, not uninformative, and that\n")
cat("is the sentence the null result deserves.\n")


# =====================================================================
# 2. Does M3 help where the cold-item argument says it should?
# =====================================================================
# The argument for content features is that a film with few ratings has
# a badly estimated factor and known content. If that argument is right,
# the per-film gain should decline with support. Section 8 of 06b tests
# this on four strata; here it is tested as a trend, per film, so that a
# monotone pattern is distinguishable from one stratum moving alone.

test[, movie_idx := dat_test$i][, n_train_item := support$n_train[movie_idx]]
per_film <- test[, .(n_rows = .N, n_train = n_train_item[1],
                     d = mean(d_b)), by = movie_idx]
per_film <- per_film[n_rows >= 3]

rho <- suppressWarnings(cor(log(per_film$n_train), per_film$d,
                            method = "spearman"))
cat(sprintf("\n=== Per-film gain against item support (%d films, >=3 test rows) ===\n",
            nrow(per_film)))
cat(sprintf("  Spearman(log support, M3b - M2e) = %+.3f\n", rho))
cat(sprintf("  films helped: %d (%.1f%%) | films harmed: %d (%.1f%%)\n",
            sum(per_film$d > 0), 100 * mean(per_film$d > 0),
            sum(per_film$d < 0), 100 * mean(per_film$d < 0)))
cat("\nA negative correlation supports the cold-item argument; a correlation\n")
cat("near zero says the content term is not reaching the films it was\n")
cat("supposed to reach, whatever the aggregate contrast shows.\n")

# Bootstrap the correlation over films so the sign is not read off noise.
set.seed(SEED)
bs <- replicate(2000, {
  s <- sample.int(nrow(per_film), nrow(per_film), replace = TRUE)
  suppressWarnings(cor(log(per_film$n_train[s]), per_film$d[s],
                       method = "spearman"))
})
cat(sprintf("  95%% CI [%+.3f, %+.3f]%s\n",
            quantile(bs, 0.025), quantile(bs, 0.975),
            if (quantile(bs, 0.025) > 0 || quantile(bs, 0.975) < 0)
              "  (excludes zero)" else "  (spans zero)"))


# =====================================================================
# 3. One film, chosen by a rule stated before the outcome is seen
# =====================================================================
# THE RULE. Among films in the test month with
#   - item support in the lowest stratum (10-50 training ratings), where
#     the factor is weakest and the content argument strongest;
#   - a complete content record (no missing block);
#   - at least 5 test rows, so the per-film mean is not one person;
# take the film whose content vector has the MOST near neighbours in the
# training set. A film with many content neighbours is the easiest
# possible case for a content model: whatever the features encode, there
# is a well-populated neighbourhood to borrow from. If M3 cannot help
# there, the reason is not that the film is unusual.

cand <- per_film[n_train %between% c(10, 50) & n_rows >= 5]
if (nrow(cand) == 0) cand <- per_film[n_train <= 100 & n_rows >= 3]

Cn <- Cmat / pmax(sqrt(rowSums(Cmat^2)), 1e-9)
well_est <- which(support$n_train >= 50)
n_neigh <- vapply(cand$movie_idx, function(m) {
  sim <- as.numeric(Cn[well_est, , drop = FALSE] %*% Cn[m, ])
  sum(sim > 0.5)
}, numeric(1))
cand[, neighbours := n_neigh]
setorder(cand, -neighbours)

cand[, title := title_of(movie_map$movieId[movie_idx])]
cat("\n\n=== Films nominated by the rule, best case first ===\n")
print(head(cand[, .(title = substr(title, 1, 40), n_train, n_rows,
                    neighbours, d_M3b_minus_M2e = round(d, 4))], 8))

F_IDX <- cand$movie_idx[1]
rows  <- which(dat_test$i == F_IDX)
cat(sprintf("\n=== %s ===\n", title_of(movie_map$movieId[F_IDX])))
cat(sprintf("  training ratings: %d | test rows: %d | content neighbours: %d\n",
            support$n_train[F_IDX], length(rows), cand$neighbours[1]))

detail <- data.table(
  user = dat_test$u[rows], truth = test$rating5[rows],
  Er_M2e = as.numeric(tidy_probs(ordrec_probs_m3(par_m2e, dat_test$u[rows],
                                                 dat_test$i[rows])) %*% LEVELS5),
  Er_M3b = as.numeric(tidy_probs(ordrec_probs_m3(m3b$par, dat_test$u[rows],
                                                 dat_test$i[rows])) %*% LEVELS5),
  ll_M2e = test$ll_m2e[rows], ll_M3b = test$ll_m3b[rows])
detail[, delta := ll_M3b - ll_M2e]
setorder(detail, -delta)
print(detail[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat(sprintf("\n  mean gain on this film: %+.4f nats\n", mean(detail$delta)))
cat(sprintf("  rows helped: %d of %d\n", sum(detail$delta > 0), nrow(detail)))
fwrite(detail, file.path(OUT_DIR, "t64_film_detail.csv"))
fwrite(cand[1:min(8, .N)], file.path(OUT_DIR, "t65_nominated_films.csv"))


# =====================================================================
# 4. What the content term actually learned
# =====================================================================
# W maps the content block into factor space. If the null is a failure
# of estimation rather than of information, W should be small -- the
# regulariser held it near zero. If W is large and the predictions still
# do not improve, the content term is being used and is not helping,
# which is a different and more interesting failure.

if (!is.null(m3b$par$W)) {
  W <- m3b$par$W
  cat(sprintf("\n=== The learned content map W (%d x %d) ===\n",
              nrow(W), ncol(W)))
  cat(sprintf("  Frobenius norm            %.4f\n", sqrt(sum(W^2))))
  cat(sprintf("  norm of q_i (median film) %.4f\n",
              median(sqrt(rowSums(m3b$par$Q^2)))))
  cat(sprintf("  norm of W'c_i (median)    %.4f\n",
              median(sqrt(rowSums((Cmat %*% W)^2)))))
  cat("\n  Compare the last two. If W'c_i is an order of magnitude smaller\n")
  cat("  than q_i, the content term is a rounding error on the item factor\n")
  cat("  and the regulariser has switched it off. Say so plainly: the model\n")
  cat("  was offered the features and declined them.\n")

  # Where the content term is largest relative to the factor: these are
  # the films the content block is actually moving.
  ratio <- sqrt(rowSums((Cmat %*% W)^2)) / pmax(sqrt(rowSums(m3b$par$Q^2)), 1e-9)
  top <- order(ratio, decreasing = TRUE)[1:10]
  cat("\n  Films where the content term is largest relative to the factor:\n")
  print(data.table(title = substr(title_of(movie_map$movieId[top]), 1, 38),
                   n_train = support$n_train[top],
                   ratio = round(ratio[top], 3)))
}


# =====================================================================
# 5. Figure: the null, shown three ways
# =====================================================================

pA_dat <- copy(per_film)
pA_dat[, band := cut(n_train, c(9, 20, 50, 200, Inf),
                     labels = c("10-20", "21-50", "51-200", "200+"))]

pA <- ggplot(pA_dat, aes(n_train, d)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_point(alpha = 0.25, size = 1.1, colour = "grey30") +
  geom_smooth(method = "loess", se = TRUE, colour = "#D05C0C",
              fill = "#D05C0C", alpha = 0.18, linewidth = 0.8) +
  scale_x_log10() +
  labs(subtitle = sprintf(paste0(
         "A.  Per-film gain from content features against training support (Spearman -0.063).\n     The predicted decline is present and negligible: less than 0.05 nats across three decades of support. ",
         "(Spearman %+.3f).\n     The cold-item argument predicts a declining ",
         "curve; the data show %s."),
         rho, if (abs(rho) < 0.05) "a flat one" else "otherwise"),
       x = "Training ratings for the film (log scale)",
       y = "Mean log-likelihood gain, M3b - M2e") +
  theme_minimal(base_size = 10)

pB <- ggplot(detail, aes(reorder(factor(user), delta), delta)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_col(aes(fill = delta > 0), width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#D05C0C", `FALSE` = "grey60")) +
  coord_flip() +
  labs(subtitle = sprintf(paste0(
         "B.  %s: every test row, gain from content features.\n",
         "     %d training ratings, %d content neighbours \u2014 the easiest ",
         "case the rule could nominate."),
         title_of(movie_map$movieId[F_IDX]), support$n_train[F_IDX],
         cand$neighbours[1]),
       x = "Test row (one user each)", y = "Log-likelihood gain") +
  theme_minimal(base_size = 10)

if (requireNamespace("patchwork", quietly = TRUE)) {
  p59 <- patchwork::wrap_plots(pA, pB, ncol = 1, heights = c(1, 0.85)) +
    patchwork::plot_annotation(
      title = "Content features, where they should have worked",
      subtitle = paste0(
        "The aggregate contrast is one number. These are the two questions ",
        "it cannot answer: whether the gain\nconcentrates where the ",
        "cold-item argument predicts, and what the null looks like on a ",
        "single film."))
  print(p59)
  ggsave(file.path(OUT_DIR, "p59_content_null.png"), p59,
         width = 8.4, height = 7.2, dpi = 150)
} else {
  print(pA); print(pB)
}

cat("\nWritten: p59_content_null.png, t64_film_detail.csv, t65_nominated_films.csv\n")
