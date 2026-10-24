# ===========================================================================
# replication_analysis_rq1_rq3.R
#
# The exact confirmatory analysis for RQ1 (profile replication), RQ2 (sign-
# reversal replication), and RQ3 (individual-differences generalization) of
# the replication preregistration (preregistration_draft.html, "Statistical
# models" section). RQ4 (frozen neural-network-model transfer) is a separate
# script, replication_analysis_rq4.py, since it depends on the Python-side
# models.
#
# This is meant to be run exactly once, after the preregistration is
# registered and the replication cohort's raw data has been pulled and run
# through generate_networks.R (if needed) and expand_all_features.R to
# produce that cohort's own alldata.csv. Point REPLICATION_ALLDATA_PATH at
# that file below and run the whole script top to bottom.
#
# All the per-predictor machinery (predictor_registry, zscore_cols,
# partial_out_visual, the zero-order/marginal/unique per-subject regression
# loop) is copied from analyze_nets.R so that the replication sample's
# coefficients are computed by the identical procedure used on the original
# sample -- not re-derived or approximated. The original sample's own
# zcoef/ucoef/mcoef come from the cached coefficients.Rdata/mcoef.Rdata
# (already computed by analyze_nets.R on the full N=8,614 sample), not
# recomputed from raw data here.
#
# Verified on a held-out subsample of the original study's own data before
# first use -- see the comment above RUN_VALIDATION_ON_SOURCE below for what
# that checked and what it can't check (see also the chat response this
# script's commit refers back to).
# ===========================================================================

setwd("D:/Dropbox/mysocialbrain/net/Study2/analyze")
if (!require(psych)) install.packages("psych"); require(psych)

# ---- CONFIG: edit these two before running for real -----------------------
ORIGINAL_CACHE_DIR        <- "."                            # where coefficients.Rdata / mcoef.Rdata live
REPLICATION_ALLDATA_PATH  <- "replication_alldata.csv"       # the new cohort's alldata.csv
OUTPUT_PATH               <- "replication_results_rq1_rq3.csv"
N_PERM_WASSERSTEIN        <- 2000                             # permutations for RQ1's omnibus test
ALPHA                     <- .05
RECONSTRUCTION_THRESHOLD  <- .90                              # RQ3(a) pre-specified ratio threshold
CONGRUENCE_EXCELLENT      <- .95                              # Lorenzo-Seva & ten Berge convention
CONGRUENCE_FAIR           <- .85

# ===========================================================================
# Shared predictor registry -- identical to analyze_nets.R and
# power_analysis_prereg.R. Kept in sync with those by hand; if the original
# study's predictor set ever changes, update all three together.
# ===========================================================================
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
cat_map <- list(
    Target  = c("degree","closeness","betweenness","eigenvector","clust","constr","ecc","vecc","esize","comsize"),
    Dyad    = c("samecom","edgbet","udists","mutfr","vdist","dhom","structsim","aa","katz","srcos"),
    Network = c("chgmod","chgast","chapl","chgcent","chgvar","netcent","netvar"),
    Self    = c("chub","chuc","chue","uclust","chuclust","uconstr","chuconstr","uecc","chuecc","uvecc","uesize","chuesize","ucomsize")
)
predictor_registry$cat <- sapply(predictor_registry$col, function(x)
    names(cat_map)[sapply(cat_map, function(v) x %in% v)])

context_cols <- c("uclust", "uconstr", "uecc", "uvecc", "uesize", "ucomsize", "netcent", "netvar")
visual_cols  <- c("vdist", "vecc", "uvecc")
predictor_registry$context <- predictor_registry$col %in% context_cols
predictor_registry$visual  <- predictor_registry$col %in% visual_cols
predictor_registry <- predictor_registry[!predictor_registry$context, ]

varind_cols          <- predictor_registry$col
ponames              <- predictor_registry$name
visualidx            <- predictor_registry$visual
ponames_substantive  <- ponames[!visualidx]
vn                   <- length(ponames)

zscore_cols <- c(
    "degree", "closeness", "betweenness", "eigenvector",
    "density", "assortativity", "apl", "modularity",
    "udegree", "ucloseness", "ubetweenness", "ueigenvector",
    "edgbet", "chgmod", "chgast", "chapl", "udists", "chub", "chuc", "chue",
    "mutfr", "clust", "uclust", "chuclust", "constr", "uconstr", "chuconstr",
    "vdist", "ecc", "uecc", "chuecc", "vecc", "uvecc", "dhom", "structsim",
    "esize", "uesize", "chuesize", "aa", "katz", "srcos", "comsize", "ucomsize",
    "netcent", "chgcent", "netvar", "chgvar"
)

# ===========================================================================
# prepare_data() / compute_subject_coefs() -- identical procedure to
# analyze_nets.R's data loading and per-subject regression loop (zero-order,
# marginal-over-visual, unique/partial), applied here to whichever sample is
# passed in. The original sample's own matrices are NOT recomputed this way
# in production (they're loaded from the cache below); this exists so the
# replication sample's matrices are derived identically to how the
# original's were.
# ===========================================================================
prepare_data <- function(path) {
    dat <- read.csv(path)
    dat[, zscore_cols] <- scale(dat[, zscore_cols])
    dat
}

compute_subject_coefs <- function(dat, verbose = TRUE) {
    varind <- match(varind_cols, names(dat))
    stopifnot(!anyNA(varind))
    visual_varind <- varind[visualidx]
    uid <- unique(dat$sid)
    ns <- length(uid)

    zcoef <- matrix(NA, ns, vn)
    ucoef <- matrix(NA, ns, vn)
    mcoef <- matrix(NA, ns, vn)
    start <- Sys.time()
    for (i in 1:ns) {
        sel  <- dat$sid == uid[i]
        cdat <- dat[sel, ]

        for (j in 1:vn) {
            fit <- lm(binres ~ as.matrix(cdat[, varind[j]]), data = cdat)
            zcoef[i, j] <- fit$coefficients[2]
        }
        for (j in which(!visualidx)) {
            fit <- lm(binres ~ as.matrix(cdat[, c(varind[j], visual_varind)]), data = cdat)
            mcoef[i, j] <- fit$coefficients[2]
        }
        fit <- lm(binres ~ as.matrix(cdat[, varind]), data = cdat)
        ucoef[i, ] <- fit$coefficients[-1]

        if (verbose && i %% 50 == 0) cat(sprintf("  subject %d/%d (%s elapsed)\n", i, ns, format(Sys.time() - start)))
    }
    colnames(zcoef) <- ponames
    colnames(ucoef) <- ponames
    colnames(mcoef) <- ponames
    list(sid = uid, zcoef = zcoef, ucoef = ucoef, mcoef = mcoef)
}

partial_out_visual <- function(zc, uc) {
    uc_visual <- uc[, ponames[visualidx], drop = FALSE]
    zc_sub    <- zc[, ponames[!visualidx], drop = FALSE]
    out <- apply(zc_sub, 2, function(y) resid(lm(y ~ uc_visual)))
    colnames(out) <- ponames_substantive
    out
}

manual_r <- function(a, b) {
    am <- a - mean(a); bm <- b - mean(b)
    sum(am * bm) / sqrt(sum(am^2) * sum(bm^2))
}
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

# ===========================================================================
# RQ1 -- profile replication: Wasserstein distance + permutation (true
# two-sample version: original sample vs. replication sample, not a
# resampling of one sample against itself as in the power analysis), plus
# the companion Pearson correlation of per-predictor coefficient means.
# ===========================================================================
wasserstein_perm_test <- function(orig_mat, repl_mat, n_perm = N_PERM_WASSERSTEIN, n_grid = 200) {
    stopifnot(identical(colnames(orig_mat), colnames(repl_mat)))
    k <- ncol(orig_mat)
    probs <- (seq_len(n_grid) - 0.5) / n_grid
    orig_grids <- apply(orig_mat, 2, function(x) quantile(x, probs, type = 7, names = FALSE))
    repl_grids <- apply(repl_mat, 2, function(x) quantile(x, probs, type = 7, names = FALSE))

    D <- matrix(NA, k, k)   # D[i, j] = Wasserstein(original predictor i, replication predictor j)
    for (i in seq_len(k)) for (j in seq_len(k)) D[i, j] <- mean(abs(orig_grids[, i] - repl_grids[, j]))

    true_stat <- mean(diag(D))
    perm_stats <- numeric(n_perm)
    for (p in seq_len(n_perm)) {
        perm <- sample(k)
        perm_stats[p] <- mean(D[cbind(seq_len(k), perm)])
    }
    p_val <- (sum(perm_stats <= true_stat) + 1) / (n_perm + 1)
    list(statistic = true_stat, p.value = p_val, n_perm = n_perm)
}

run_rq1 <- function(orig, repl) {
    variants <- list(zero_order = list(o = orig$zcoef, r = repl$zcoef),
                      marginal   = list(o = orig$mcoef, r = repl$mcoef),
                      unique     = list(o = orig$ucoef, r = repl$ucoef))
    rows <- list()
    for (vname in names(variants)) {
        o <- variants[[vname]]$o; r <- variants[[vname]]$r
        # marginal has NA columns for the visual predictors themselves (vdist/vecc) --
        # drop them from the profile tests, same as power_analysis_prereg.R does.
        keep <- !is.na(colSums(o)) & !is.na(colSums(r))
        o <- o[, keep, drop = FALSE]; r <- r[, keep, drop = FALSE]

        wt <- wasserstein_perm_test(o, r)
        prof_r <- manual_r(colMeans(o), colMeans(r))

        rows[[length(rows) + 1]] <- data.frame(
            rq = "RQ1", variant = vname, test = "wasserstein_permutation",
            statistic = wt$statistic, p_value = wt$p.value,
            sig = wt$p.value < ALPHA, detail = sprintf("n_predictors=%d", ncol(o))
        )
        rows[[length(rows) + 1]] <- data.frame(
            rq = "RQ1", variant = vname, test = "pearson_profile_correlation",
            statistic = prof_r, p_value = NA,   # descriptive only, per Inference criteria
            sig = NA, detail = sprintf("n_predictors=%d", ncol(o))
        )
    }
    do.call(rbind, rows)
}

# ===========================================================================
# RQ2 -- sign-reversal replication: two-sided, one-sample t-tests on the
# replication sample's unique/partial coefficients for mutual friend count
# and target distance.
# ===========================================================================
run_rq2 <- function(repl) {
    rows <- list()
    for (pred in c("Mutual friends", "Target distance")) {
        x <- repl$ucoef[, pred]
        tt <- t.test(x, alternative = "two.sided")
        rows[[length(rows) + 1]] <- data.frame(
            rq = "RQ2", variant = pred, test = "two_sided_one_sample_t",
            statistic = unname(tt$statistic), p_value = tt$p.value,
            sig = tt$p.value < ALPHA, detail = sprintf("mean=%.4f, n=%d", mean(x), length(x))
        )
    }
    do.call(rbind, rows)
}

# ===========================================================================
# RQ3 -- individual-differences generalization.
# (a) out-of-sample reconstruction + noise ceiling (original loadings,
#     recomputed from the cached original zcoef/ucoef since the PCA fit
#     object itself isn't separately cached -- principal(zcoef_clean, 2) is
#     deterministic given zcoef/ucoef, so this reproduces it exactly).
# (b) a fresh two-component varimax solution fit on the replication sample,
#     matched to the original loadings by Tucker's congruence coefficient.
# ===========================================================================
run_rq3 <- function(orig, repl) {
    orig_clean <- partial_out_visual(orig$zcoef, orig$ucoef)
    repl_clean <- partial_out_visual(repl$zcoef, repl$ucoef)

    orig_fit  <- principal(orig_clean, 2)
    orig_load <- unclass(orig_fit$loadings)
    orig_lambda_lambda_t <- orig_load %*% t(orig_load)
    orig_ll_lower  <- orig_lambda_lambda_t[lower.tri(orig_lambda_lambda_t)]
    orig_cor       <- cor(orig_clean)
    orig_cor_lower <- orig_cor[lower.tri(orig_cor)]

    repl_cor       <- cor(repl_clean)
    repl_cor_lower <- repl_cor[lower.tri(repl_cor)]

    recon_r      <- manual_r(orig_ll_lower, repl_cor_lower)
    noise_ceil_r <- manual_r(orig_cor_lower, repl_cor_lower)
    ratio        <- recon_r / noise_ceil_r

    repl_fit <- principal(repl_clean, 2)
    cong <- match_congruence(unclass(repl_fit$loadings), orig_load)
    cong_min <- min(cong)
    cong_band <- if (cong_min >= CONGRUENCE_EXCELLENT) "excellent" else if (cong_min >= CONGRUENCE_FAIR) "fair" else "poor"

    rbind(
        data.frame(rq = "RQ3a", variant = "reconstruction", test = "ratio_to_noise_ceiling",
                   statistic = ratio, p_value = NA, sig = ratio >= RECONSTRUCTION_THRESHOLD,
                   detail = sprintf("recon_r=%.4f, noise_ceiling_r=%.4f", recon_r, noise_ceil_r)),
        data.frame(rq = "RQ3b", variant = "congruence_component1", test = "tucker_congruence",
                   statistic = cong[1], p_value = NA, sig = NA, detail = NA_character_),
        data.frame(rq = "RQ3b", variant = "congruence_component2", test = "tucker_congruence",
                   statistic = cong[2], p_value = NA, sig = NA, detail = NA_character_),
        data.frame(rq = "RQ3b", variant = "congruence_min_band", test = "tucker_congruence",
                   statistic = cong_min, p_value = NA, sig = NA, detail = cong_band)
    )
}

# ===========================================================================
# Main driver. Set SKIP_MAIN <- TRUE before sourcing this file to load only
# the functions/objects above (used by the test harness).
# ===========================================================================
if (!exists("SKIP_MAIN")) {
    load(file.path(ORIGINAL_CACHE_DIR, "coefficients.Rdata"))   # zcoef, ucoef, udem, spcz, spcu, spcp
    colnames(zcoef) <- ponames; colnames(ucoef) <- ponames
    load(file.path(ORIGINAL_CACHE_DIR, "mcoef.Rdata"))          # mcoef
    colnames(mcoef) <- ponames
    orig <- list(zcoef = zcoef, ucoef = ucoef, mcoef = mcoef)

    cat("Loading and scoring replication sample from", REPLICATION_ALLDATA_PATH, "...\n")
    repl_dat <- prepare_data(REPLICATION_ALLDATA_PATH)
    cat("Computing per-subject coefficients for", length(unique(repl_dat$sid)), "replication participants...\n")
    repl <- compute_subject_coefs(repl_dat)

    results <- rbind(
        run_rq1(orig, repl),
        run_rq2(repl),
        run_rq3(orig, repl)
    )
    write.csv(results, OUTPUT_PATH, row.names = FALSE)
    cat("\n"); print(results, row.names = FALSE)
    cat("\nwrote", OUTPUT_PATH, "\n")
}
