# =====================================================================
# 05c_M3_features_v8.R
# Phase B, step 2: rebuild the content features for the v8 movie set.
# =====================================================================

FEATURE_VERSION <- "v8"
OVERVIEW_DIM    <- 32L
MIN_DOC_COUNT   <- 20L   # absolute, as in 05b, so the vocabulary rule is
                         # identical rather than rescaled to the new size
SEED            <- 2024L

BASE      <- "D:/桌面/movie-recommendation-dissertation"
SEL_FILE  <- file.path(BASE, "06_models/cache/00_data_selection/selection_v8_2018-04_6m_Nu10_Nm1.rds")
ML_DIR    <- file.path(BASE, "data/ml-32m")
TMDB_CSV  <- file.path(BASE, "data/TMDB/tmdb_meta_full.csv")
CACHE_DIR <- file.path(BASE, "06_models/cache/08_M3")
OUT_DIR   <- file.path(BASE, "06_models/output/08_M3")

dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(ggplot2)
})
library(text2vec); library(irlba)

set.seed(SEED)
n_cores <- max(1L, parallel::detectCores() - 1L)
is_missing_text <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))
BANNED <- c("popularity", "tmdb_rating", "tmdb_votes", "revenue", "budget")


# =====================================================================
# 1. Movie universe and support
# =====================================================================
sel   <- readRDS(SEL_FILE)
fit   <- as.data.table(sel$fit)
calib <- as.data.table(sel$calib)
stopifnot(sel$final$Nm < 10L)   # this file is for the lowered-Nm selection

movie_map <- data.table(movieId = sort(unique(fit$movieId)))[, movie_idx := .I]
nI <- nrow(movie_map)
cat(sprintf("v8 movies in movie_map: %d\n", nI))

support <- rbind(fit[, .(movieId)], calib[, .(movieId)])[, .(n_train = .N), by = movieId]
support <- merge(movie_map, support, by = "movieId", all.x = TRUE)
support[is.na(n_train), n_train := 0L]
setorder(support, movie_idx)
cat(sprintf("support: min %d | median %d | max %d | below 10: %d films\n",
            min(support$n_train), median(support$n_train), max(support$n_train),
            sum(support$n_train < 10)))

# A film can sit in movie_map with zero fit rows if all its training
# ratings landed in calib. Its factor gets no gradient at all, which is
# a legitimate state here but must be visible.
n_nofit <- nI - uniqueN(fit$movieId)
if (n_nofit > 0)
  cat(sprintf("NOTE: %d films have calib rows but no fit rows.\n", n_nofit))


# =====================================================================
# 2. Metadata
# =====================================================================
movies <- fread(file.path(ML_DIR, "movies.csv"), nThread = n_cores)
links  <- fread(file.path(ML_DIR, "links.csv"),  nThread = n_cores)
tmdb   <- unique(fread(TMDB_CSV, nThread = n_cores), by = "tmdbId")

meta <- merge(movie_map, movies[, .(movieId, title, genres)], by = "movieId", all.x = TRUE)
meta <- merge(meta, links[!is.na(tmdbId), .(movieId, tmdbId)], by = "movieId", all.x = TRUE)
meta <- merge(meta, tmdb[, .(tmdbId, overview, runtime, release_year)],
              by = "tmdbId", all.x = TRUE)
setorder(meta, movie_idx)
stopifnot(identical(meta$movie_idx, seq_len(nI)))

cat(sprintf("TMDb match rate: %.1f%% (v7 was 99.x%% - the tail matches worse)\n",
            100 * mean(!is.na(meta$tmdbId))))


# =====================================================================
# 3. Feature blocks
# =====================================================================
cat("\n[1/3] genre multi-hot...\n")
gen_long <- meta[!grepl("no genres listed", genres, ignore.case = TRUE),
                 .(movie_idx, genre = trimws(unlist(strsplit(genres, "|", fixed = TRUE)))),
                 by = .I][, .(movie_idx, genre)]
gen_levels <- sort(unique(gen_long$genre))
G <- matrix(0, nI, length(gen_levels),
            dimnames = list(NULL, paste0("genre_", gsub("[^A-Za-z0-9]", "_", gen_levels))))
G[cbind(gen_long$movie_idx, match(gen_long$genre, gen_levels))] <- 1
cat(sprintf("  %d genres | films with none: %d\n", ncol(G), sum(rowSums(G) == 0)))

cat("\n[2/3] overview TF-IDF + SVD...\n")
ov_cache <- file.path(CACHE_DIR, sprintf("overview_svd_%s.rds", FEATURE_VERSION))
if (file.exists(ov_cache)) {
  ov <- readRDS(ov_cache); cat("  loaded from cache\n")
} else {
  has_ov <- !is_missing_text(meta$overview)
  it  <- itoken(tolower(meta$overview[has_ov]), tokenizer = word_tokenizer,
                progressbar = FALSE)
  vc  <- prune_vocabulary(
    create_vocabulary(it, stopwords = stopwords::stopwords("en")),
    doc_count_min = MIN_DOC_COUNT, doc_proportion_max = 0.5)
  X <- fit_transform(create_dtm(it, vocab_vectorizer(vc)), TfIdf$new())
  set.seed(SEED)
  sv <- irlba::irlba(X, nv = OVERVIEW_DIM)
  Z <- matrix(0, nI, OVERVIEW_DIM,
              dimnames = list(NULL, sprintf("ov%02d", seq_len(OVERVIEW_DIM))))
  Z[has_ov, ] <- sv$u %*% diag(sv$d)
  ov <- list(Z = Z, has_ov = has_ov, vocab_size = nrow(vc),
             var_explained = sv$d^2 / sum(sv$d^2))
  saveRDS(ov, ov_cache)
}
cat(sprintf("  vocab %d | with overview %d (%.1f%%)\n",
            ov$vocab_size, sum(ov$has_ov), 100 * mean(ov$has_ov)))

cat("\n[3/3] runtime and release year...\n")
rt <- as.numeric(meta$runtime); rt[!is.finite(rt) | rt < 40 | rt > 300] <- NA_real_
has_rt <- !is.na(rt); rt[!has_rt] <- median(rt, na.rm = TRUE)
ry <- as.numeric(meta$release_year); ry[!is.finite(ry) | ry < 1900 | ry > 2018] <- NA_real_
has_ry <- !is.na(ry); ry[!has_ry] <- median(ry, na.rm = TRUE)
film_age <- 2018 - ry
cat(sprintf("  runtime missing %d | year missing %d\n", sum(!has_rt), sum(!has_ry)))


# =====================================================================
# 4. Assemble and standardise
# =====================================================================
C_raw <- cbind(G, ov$Z, runtime = rt, film_age = film_age,
               has_overview = as.numeric(ov$has_ov),
               has_runtime  = as.numeric(has_rt),
               has_year     = as.numeric(has_ry))
block <- c(rep("genre", ncol(G)), rep("overview", ncol(ov$Z)),
           "numeric", "numeric", "missing", "missing", "missing")

keep <- apply(C_raw, 2, function(x) sd(x) > 1e-8)
if (any(!keep))
  cat(sprintf("  dropping %d zero-variance columns\n", sum(!keep)))
C_raw <- C_raw[, keep, drop = FALSE]; block <- block[keep]

ctr <- colMeans(C_raw); scl <- apply(C_raw, 2, sd)
C <- as.matrix(scale(C_raw, center = ctr, scale = scl))
attr(C, "scaled:center") <- NULL; attr(C, "scaled:scale") <- NULL
stopifnot(all(is.finite(C)), nrow(C) == nI, !any(BANNED %in% colnames(C)))

support[, nm_band := cut(n_train, c(0, 2, 5, 9, 20, 50, 200, Inf),
                         labels = c("1-2", "3-5", "6-9", "10-20",
                                    "21-50", "51-200", "200+"))]

saveRDS(list(C = C, block = block, movie_map = movie_map,
             center = ctr, scale = scl, support = support,
             meta = list(version = FEATURE_VERSION, built = Sys.time(),
                         sel_file = SEL_FILE, overview_dim = OVERVIEW_DIM,
                         min_doc_count = MIN_DOC_COUNT,
                         excluded = BANNED, seed = SEED)),
        file.path(CACHE_DIR, sprintf("m3_features_%s.rds", FEATURE_VERSION)))
cat(sprintf("\nSaved m3_features_%s.rds  |  C: %d x %d\n",
            FEATURE_VERSION, nrow(C), ncol(C)))
print(table(block))


# =====================================================================
# 5. The confound check, now across the full support range
# =====================================================================
# In v7 the strongest correlate of support was runtime (r = 0.25) and
# film_age (r = -0.10). Whether that steepens once the sub-10 tail is
# included decides how much of any stratified gain can be attributed to
# sparsity rather than to content composition.
cor_sup <- apply(C, 2, function(x) cor(x, log(support$n_train)))
cat("\nTop |correlation| with log support:\n")
print(round(sort(cor_sup[order(-abs(cor_sup))][1:10]), 3))

d_cor <- data.table(feature = names(cor_sup), r = cor_sup, block = block)[order(-abs(r))][1:25]
p_cor <- ggplot(d_cor, aes(reorder(feature, r), r, fill = block)) +
  geom_col() + coord_flip() + theme_minimal() +
  labs(title = "v8: content features vs item support",
       subtitle = "Compare against the v7 version: a steeper profile means a stronger composition confound",
       x = NULL, y = "correlation with log(training ratings)")
ggsave(file.path(OUT_DIR, "m3j_v8_feature_support_correlation.png"),
       p_cor, width = 8, height = 6, dpi = 150)

cat("\n05c COMPLETE.\n")
