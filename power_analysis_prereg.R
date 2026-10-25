# ===========================================================================
# power_analysis_prereg.R
#
# Power analysis for every R-side analysis planned in the replication
# preregistration (preregistration_draft.html): RQ1 (profile replication,
# both Wasserstein+permutation and Pearson correlation), RQ2 (sign-reversal
# replication), and RQ3 (PCA out-of-sample reconstruction + noise ceiling,
# and Tucker's congruence on a freshly-fit solution). RQ4 (frozen neural
# network model transfer) is covered separately in power_analysis_prereg.py,
# since it depends on the Python-side Brier skill results.
#
# All effect sizes are drawn from the original study's own data
# (coefficients.Rdata, mcoef.Rdata) -- this is a retrospective power
# analysis, not a prospective one with independently-justified priors; see
# power_analysis_report.html for that caveat spelled out in more detail.
# ===========================================================================
# setwd("PATH/TO/graph_friendship")  # uncomment and edit to point at this repo, or run from this directory
if(!require(psych)) install.packages("psych"); require(psych)
if(!require(pwr))   install.packages("pwr");   require(pwr)

load("coefficients.Rdata")   # zcoef, ucoef
load("mcoef.Rdata")          # mcoef
ns <- nrow(zcoef)

predictor_registry <- data.frame(
    col = c("degree","closeness","betweenness","eigenvector","samecom","edgbet",
           "chgmod","chgast","chapl","udists","chub","chuc","chue",
           "mutfr","clust","uclust","chuclust","constr","uconstr","chuconstr",
           "vdist","ecc","uecc","chuecc","vecc","uvecc","dhom","structsim",
           "esize","uesize","chuesize","aa","katz","srcos","comsize","ucomsize",
           "netcent","chgcent","netvar","chgvar"),
    name = c("Target degree","Target closeness","Target betweenness","Target eigenvector",
            "Same community","Edge betweenness",
            "Change modularity","Change assortativity","Change average path length",
            "Target distance","Change own betweenness","Change own closeness","Change own eigenvector",
            "Mutual friends","Target clustering coefficient","Self clustering coefficient",
            "Change own clustering coefficient","Target Burt's constraint","Self Burt's constraint",
            "Change own Burt's constraint",
            "Visual distance","Target eccentricity","Self eccentricity","Change own eccentricity",
            "Target visual eccentricity","Self visual eccentricity","Degree homophily",
            "Structural similarity","Target effective size","Self effective size",
            "Change own effective size","Adamic-Adar","Katz index",
            "Successor representation similarity","Target community size","Self community size",
            "Network centralization","Change centralization","Network degree variance","Change degree variance"),
    stringsAsFactors = FALSE
)
visual_cols  <- c("vdist", "vecc", "uvecc")
context_cols <- c("uclust", "uconstr", "uecc", "uvecc", "uesize", "ucomsize", "netcent", "netvar")
predictor_registry$context <- predictor_registry$col %in% context_cols
predictor_registry$visual  <- predictor_registry$col %in% visual_cols
predictor_registry <- predictor_registry[!predictor_registry$context, ]
visualidx <- predictor_registry$visual
colnames(zcoef) <- predictor_registry$name
colnames(ucoef) <- predictor_registry$name
colnames(mcoef) <- predictor_registry$name

Ns <- c(100, 150, 200, 300, 400, 500, 750, 1000, 1500, 2000)
N_REPS <- 50
TARGET_N <- 300   # the replication's actual planned sample size

sink("power_prereg_r_results.txt")
cat("Power analysis for the R-side preregistration analyses (RQ1-RQ3)\n")
cat("Original sample N =", ns, "\n\n")

# ===========================================================================
# RQ2 -- sign-reversal replication (two-sided, mutual friends & target
# distance unique/partial coefficients)
# ===========================================================================
cat("=== RQ2: sign-reversal replication (two-sided t-tests) ===\n")
for (pred in c("Mutual friends", "Target distance")) {
    x <- ucoef[, pred]
    d <- t.test(x)$statistic / sqrt(ns)
    cat(sprintf("%s: unique/partial d = %.3f\n", pred, d))
    power_at_target <- pwr.t.test(d = abs(d), n = TARGET_N, sig.level = .05, type = "one.sample")$power
    cat(sprintf("  power at N=%d (two-sided): %.4f\n", TARGET_N, power_at_target))
    for (pw in c(.80, .90, .95)) {
        n_needed <- ceiling(pwr.t.test(d = abs(d), sig.level = .05, power = pw, type = "one.sample")$n)
        cat(sprintf("  N for %.0f%% power (two-sided): %d\n", pw * 100, n_needed))
    }
}

# ===========================================================================
# RQ1 -- profile replication: Pearson correlation (fast, closed-form-
# adjacent via simulation) and Wasserstein distance + permutation (fully
# simulation-based, no closed form)
# ===========================================================================
cat("\n=== RQ1: profile replication ===\n")

manual_r <- function(a, b) {
    am <- a - mean(a); bm <- b - mean(b)
    sum(am * bm) / sqrt(sum(am^2) * sum(bm^2))
}

# ---- Pearson profile correlation, by simulation ----
profile_corr_power <- function(coef_mat, label, Ns, n_reps) {
    ref <- colMeans(coef_mat)
    k <- ncol(coef_mat)
    out <- data.frame(N = Ns, mean_r = NA, power_sig = NA)
    for (ni in seq_along(Ns)) {
        N <- Ns[ni]
        rs <- numeric(n_reps); sig <- 0
        for (r in seq_len(n_reps)) {
            idx <- sample(ns, min(N, ns), replace = (N > ns))
            prof <- colMeans(coef_mat[idx, , drop = FALSE])
            rr <- manual_r(prof, ref)
            rs[r] <- rr
            tstat <- rr * sqrt((k - 2) / max(1e-9, 1 - rr^2))
            p <- 2 * (1 - pt(abs(tstat), df = k - 2))
            if (!is.na(p) && p < .05 && rr > 0) sig <- sig + 1
        }
        out$mean_r[ni] <- mean(rs); out$power_sig[ni] <- sig / n_reps
    }
    cat("\n--", label, "(Pearson profile correlation) --\n")
    print(out, row.names = FALSE)
    write.csv(out, sprintf("power_prereg_rq1_pearson_%s.csv", label), row.names = FALSE)
    out
}
set.seed(1)
rq1_pearson_zero     <- profile_corr_power(zcoef, "zero_order", Ns, N_REPS)
rq1_pearson_unique   <- profile_corr_power(ucoef, "unique", Ns, N_REPS)
rq1_pearson_marginal <- profile_corr_power(mcoef[, !is.na(colSums(mcoef))], "marginal", Ns, N_REPS)

# ---- Wasserstein distance + permutation ----
# Fast 1-D Wasserstein-1 distance via a shared quantile grid (avoids the
# O(n log n) exact transport solve -- accurate to the grid resolution,
# which is more than enough precision for a power curve).
wasserstein1d <- function(x, qy_grid, probs) {
    qx <- quantile(x, probs, type = 7, names = FALSE)
    mean(abs(qx - qy_grid))
}

wasserstein_power <- function(coef_mat, label, Ns, n_reps, n_perm = 200, n_grid = 200) {
    k <- ncol(coef_mat)
    probs <- (seq_len(n_grid) - 0.5) / n_grid
    full_grids <- apply(coef_mat, 2, function(x) quantile(x, probs, type = 7, names = FALSE))  # n_grid x k

    out <- data.frame(N = Ns, power = NA)
    for (ni in seq_along(Ns)) {
        N <- Ns[ni]
        sig <- 0
        for (r in seq_len(n_reps)) {
            idx <- sample(ns, min(N, ns), replace = (N > ns))
            sub <- coef_mat[idx, , drop = FALSE]
            # distance matrix: sub predictor i vs. full-sample predictor j
            D <- matrix(NA, k, k)
            for (i in seq_len(k)) {
                for (j in seq_len(k)) {
                    D[i, j] <- wasserstein1d(sub[, i], full_grids[, j], probs)
                }
            }
            true_stat <- mean(diag(D))
            perm_stats <- numeric(n_perm)
            for (p in seq_len(n_perm)) {
                perm <- sample(k)
                perm_stats[p] <- mean(D[cbind(seq_len(k), perm)])
            }
            p_val <- (sum(perm_stats <= true_stat) + 1) / (n_perm + 1)
            if (p_val < .05) sig <- sig + 1
        }
        out$power[ni] <- sig / n_reps
        cat(sprintf("  %s N=%d power=%.2f\n", label, N, out$power[ni]))
    }
    write.csv(out, sprintf("power_prereg_rq1_wasserstein_%s.csv", label), row.names = FALSE)
    out
}

cat("\n-- Wasserstein distance + permutation (this is the slow part) --\n")
Ns_wass <- c(100, 200, 300, 500, 750, 1000)   # coarser grid -- this test is expensive
set.seed(2)
rq1_wass_zero     <- wasserstein_power(zcoef, "zero_order", Ns_wass, 30)
rq1_wass_unique   <- wasserstein_power(ucoef, "unique", Ns_wass, 30)
rq1_wass_marginal <- wasserstein_power(mcoef[, !is.na(colSums(mcoef))], "marginal", Ns_wass, 30)

# ===========================================================================
# RQ3 -- individual-differences generalization
# (a) out-of-sample reconstruction + noise ceiling
# (b) Tucker's congruence on a freshly-fit solution
# ===========================================================================
cat("\n=== RQ3: PCA out-of-sample reconstruction + congruence ===\n")

partial_out_visual <- function(zc, uc) {
    uc_visual <- uc[, predictor_registry$name[visualidx], drop = FALSE]
    zc_sub    <- zc[, predictor_registry$name[!visualidx], drop = FALSE]
    out <- apply(zc_sub, 2, function(y) resid(lm(y ~ uc_visual)))
    colnames(out) <- predictor_registry$name[!visualidx]
    out
}
zcoef_clean <- partial_out_visual(zcoef, ucoef)
full_fit  <- principal(zcoef_clean, 2)
full_load <- unclass(full_fit$loadings)
full_lambda_lambda_t <- full_load %*% t(full_load)
full_lower <- full_lambda_lambda_t[lower.tri(full_lambda_lambda_t)]
full_cor_full <- cor(zcoef_clean)
full_cor_lower <- full_cor_full[lower.tri(full_cor_full)]

tucker <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))
match_congruence <- function(sub_load, ref_load) {
    best <- -Inf; best_cong <- c(NA, NA)
    for (perm in list(c(1, 2), c(2, 1))) {
        c1 <- tucker(sub_load[, perm[1]], ref_load[, 1])
        c2 <- tucker(sub_load[, perm[2]], ref_load[, 2])
        tot <- abs(c1) + abs(c2)
        if (tot > best) { best <- tot; best_cong <- c(abs(c1), abs(c2)) }
    }
    best_cong
}

Ns_pca <- c(100, 150, 200, 300, 500, 750, 1000, 1500, 2000)
N_REPS_PCA <- 40
rows <- list()
set.seed(3)
for (N in Ns_pca) {
    for (rep in seq_len(N_REPS_PCA)) {
        idx <- sample(ns, min(N, ns), replace = (N > ns))
        sub <- zcoef_clean[idx, , drop = FALSE]
        sub <- sub[, apply(sub, 2, function(x) sd(x) > 0), drop = FALSE]
        if (ncol(sub) < ncol(zcoef_clean)) next
        sub_cor <- tryCatch(cor(sub), error = function(e) NULL)
        if (is.null(sub_cor)) next
        sub_cor_lower <- sub_cor[lower.tri(sub_cor)]

        recon_r      <- manual_r(full_lower, sub_cor_lower)              # (a) model-implied vs. replication empirical
        noise_ceil_r <- manual_r(full_cor_lower, sub_cor_lower)          # (a) noise ceiling
        ratio <- recon_r / noise_ceil_r

        fit2 <- tryCatch(principal(sub, 2), error = function(e) NULL)    # (b) fresh fit + congruence
        cong <- if (!is.null(fit2)) match_congruence(unclass(fit2$loadings), full_load) else c(NA, NA)

        rows[[length(rows) + 1]] <- data.frame(N = N, rep = rep, recon_r = recon_r,
                                               noise_ceiling_r = noise_ceil_r, ratio = ratio,
                                               cong1 = cong[1], cong2 = cong[2])
    }
    cat("  PCA N=", N, "done\n")
}
pca_res <- do.call(rbind, rows)
pca_summary <- aggregate(cbind(recon_r, noise_ceiling_r, ratio, cong1, cong2) ~ N, pca_res, mean, na.rm = TRUE)
pca_summary$power_ratio_ge90 <- sapply(split(pca_res$ratio, pca_res$N), function(x) mean(x >= .90, na.rm = TRUE))[as.character(pca_summary$N)]
pca_summary$power_cong_ge95  <- sapply(split(pmin(pca_res$cong1, pca_res$cong2), pca_res$N), function(x) mean(x >= .95, na.rm = TRUE))[as.character(pca_summary$N)]
write.csv(pca_summary, "power_prereg_rq3_pca_summary.csv", row.names = FALSE)
cat("\n")
print(pca_summary, row.names = FALSE)

# ===========================================================================
# Headline summary at the replication's actual planned sample size
# (TARGET_N = 300), pulled from the grids above rather than re-simulated.
# ===========================================================================
cat("\n\n=== HEADLINE: power at the planned N =", TARGET_N, "===\n")
at_n <- function(df, n = TARGET_N) df[df$N == n, , drop = FALSE]

cat("\nRQ1 -- Pearson profile correlation (descriptive, not a separate test):\n")
cat(sprintf("  zero-order: mean r = %.4f\n", at_n(rq1_pearson_zero)$mean_r))
cat(sprintf("  unique:     mean r = %.4f\n", at_n(rq1_pearson_unique)$mean_r))
cat(sprintf("  marginal:   mean r = %.4f\n", at_n(rq1_pearson_marginal)$mean_r))

cat("\nRQ1 -- Wasserstein distance + permutation (omnibus test):\n")
cat(sprintf("  zero-order: power = %.4f\n", at_n(rq1_wass_zero)$power))
cat(sprintf("  unique:     power = %.4f\n", at_n(rq1_wass_unique)$power))
cat(sprintf("  marginal:   power = %.4f\n", at_n(rq1_wass_marginal)$power))

cat("\nRQ3 -- PCA generalization:\n")
pca_at_n <- at_n(pca_summary)
cat(sprintf("  reconstruction ratio: mean = %.4f, power(ratio>=.90) = %.4f\n",
            pca_at_n$ratio, pca_at_n$power_ratio_ge90))
cat(sprintf("  congruence:           mean(min component) = %.4f, power(cong>=.95) = %.4f\n",
            min(pca_at_n$cong1, pca_at_n$cong2), pca_at_n$power_cong_ge95))

cat("\n(RQ2's power at N =", TARGET_N, "is printed above, in the RQ2 section.)\n")

sink()
cat("done -- see power_prereg_r_results.txt and power_prereg_rq*.csv\n")
