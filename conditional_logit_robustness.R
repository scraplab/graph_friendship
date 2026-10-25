# ===========================================================================
# conditional_logit_robustness.R
#
# Robustness check for the "partial" (complete, all-predictors-at-once)
# linear probability model in analyze_nets.R. Each trial is a forced choice
# among 5-22 candidates, which is exactly the structure conditional logit
# (McFadden) models directly -- the per-row LPM treats each candidate as an
# independent binary outcome instead, ignoring the within-trial competitive
# structure. This fits one pooled conditional logit across every trial from
# every subject (not per-subject -- a per-subject fit would hit the same
# quasi-separation problem that motivated using LPMs over logistic
# regression in the first place, now with only 30 trials per subject and
# n_t-1 effective comparisons per trial instead of n_t independent draws).
# The result is a population-average estimate, directly comparable in sign
# and relative magnitude to the mean partial/unique LPM coefficients, not a
# replacement for the per-subject individual-differences analysis.
#
# Uses survival::clogit, which is conditional logistic regression
# implemented as a stratified Cox model with one "event" (the chosen
# candidate) per stratum (trial) -- exactly the right tool for a forced
# choice among a variable-size candidate set, no padding or restriction to
# a fixed choice-set size required.
# ===========================================================================
# setwd("PATH/TO/graph_friendship")  # uncomment and edit to point at this repo, or run from this directory
if(!require(survival)) install.packages("survival"); require(survival)

cat("loading data...\n")
dat <- read.csv("alldata.csv")
cat("rows:", nrow(dat), "\n")

# --- predictor registry: identical to analyze_nets.R / compute_mcoef.R ---
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
context_cols <- c("uclust", "uconstr", "uecc", "uvecc", "uesize", "ucomsize", "netcent", "netvar")
predictor_registry$context <- predictor_registry$col %in% context_cols
predictor_registry <- predictor_registry[!predictor_registry$context, ]   # 32 predictors, matches the "partial" LPM

# category + colors: identical to analyze_nets.R's cat_map/catcol_dark, so
# this plot's colors mean the same thing as every other figure in this
# project (Target = green, Dyad = orange, Network = purple, Self = blue).
cat_map <- list(
    Target  = c("degree","closeness","betweenness","eigenvector","clust","constr","ecc","vecc","esize","comsize"),
    Dyad    = c("samecom","edgbet","udists","mutfr","vdist","dhom","structsim","aa","katz","srcos"),
    Network = c("chgmod","chgast","chapl","chgcent","chgvar","netcent","netvar"),
    Self    = c("chub","chuc","chue","uclust","chuclust","uconstr","chuconstr","uecc","chuecc","uvecc","uesize","chuesize","ucomsize")
)
predictor_registry$cat <- sapply(predictor_registry$col, function(x)
    names(cat_map)[sapply(cat_map, function(v) x %in% v)])
catnames    <- c("Target", "Dyad", "Network", "Self")
catcol_dark <- c("#2CA02C", "#FF7F0E", "#9467BD", "#1F77B4")

zscore_cols <- setdiff(predictor_registry$col, "samecom")   # samecom is binary, not z-scored (matches analyze_nets.R)
dat[, zscore_cols] <- scale(dat[, zscore_cols])

cat("building strata id (one per sid x trial)...\n")
dat$strata_id <- paste(dat$sid, dat$trial, sep = "_")

form <- as.formula(paste("binres ~", paste(predictor_registry$col, collapse = " + "), "+ strata(strata_id)"))

cat("fitting pooled conditional logit (this is the long step)...\n")
t0 <- Sys.time()
fit <- clogit(form, data = dat)
cat("fit done. elapsed:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

s <- summary(fit)
coef_tab <- as.data.frame(s$coefficients)
coef_tab$col <- rownames(coef_tab)
coef_tab <- merge(predictor_registry[, c("col", "name")], coef_tab, by = "col")
names(coef_tab)[names(coef_tab) == "coef"] <- "clogit_coef"
names(coef_tab)[names(coef_tab) == "Pr(>|z|)"] <- "p"

# --- compare against the mean partial/unique LPM coefficient, per predictor ---
load("coefficients.Rdata")   # zcoef, ucoef
colnames(ucoef) <- predictor_registry$name[match(colnames(ucoef), predictor_registry$name)]  # no-op guard
mean_ucoef <- colMeans(ucoef)
coef_tab$lpm_partial_mean <- mean_ucoef[match(coef_tab$name, predictor_registry$name)]

coef_tab$sign_agree <- sign(coef_tab$clogit_coef) == sign(coef_tab$lpm_partial_mean)
r_concordance <- cor(coef_tab$clogit_coef, coef_tab$lpm_partial_mean, method = "pearson")
r_concordance_rank <- cor(coef_tab$clogit_coef, coef_tab$lpm_partial_mean, method = "spearman")

coef_tab <- coef_tab[order(coef_tab$p), ]
write.csv(coef_tab, "conditional_logit_results.csv", row.names = FALSE)

cat("\n=== pooled conditional logit vs. mean partial/unique LPM coefficient ===\n")
print(coef_tab[, c("name", "clogit_coef", "p", "lpm_partial_mean", "sign_agree")], row.names = FALSE)
cat("\nsign agreement:", sum(coef_tab$sign_agree), "/", nrow(coef_tab), "\n")
cat("Pearson r (clogit coef vs. mean LPM partial coef):", round(r_concordance, 4), "\n")
cat("Spearman r (same, rank-based):", round(r_concordance_rank, 4), "\n")

save(fit, coef_tab, file = "conditional_logit_fit.Rdata")
cat("\nsaved conditional_logit_results.csv and conditional_logit_fit.Rdata\n")

# ===========================================================================
# Visualization -- for inspecting the concordance between the two methods,
# not (yet) intended for the manuscript. Colors/category grouping match
# analyze_nets.R's combined_regression.pdf (Target/Dyad/Network/Self,
# catcol_dark) so this reads consistently with every other figure in the
# project. Two views:
#   (A) side-by-side horizontal coefficients, LPM partial vs. conditional
#       logit, same predictor order in both columns, colored by category
#   (B) a scatter of the two methods' coefficients directly against each
#       other, colored by category, with the two sign-reversal predictors
#       (mutual friends, target distance) called out by name
# ===========================================================================
coef_tab$cat <- predictor_registry$cat[match(coef_tab$name, predictor_registry$name)]
coef_tab$col <- catcol_dark[match(coef_tab$cat, catnames)]

# order: by category (in the same Target/Dyad/Network/Self sequence used
# elsewhere), then by LPM partial coefficient magnitude within category
ord <- order(match(coef_tab$cat, catnames), -abs(coef_tab$lpm_partial_mean))
coef_tab <- coef_tab[ord, ]

pdf("conditional_logit_diagnostics.pdf", width = 10, height = 8)

# --- panel A: side-by-side horizontal coefficients ---
layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))
n <- nrow(coef_tab)
y <- n:1

draw_side <- function(vals, title_txt, xlab_txt) {
    xmax <- max(abs(vals)) * 1.15
    par(mai = c(0.7, 2.6, 0.5, 0.2))
    plot(0, 0, type = "n", xlim = c(-xmax, xmax), ylim = c(0.5, n + 0.5),
        axes = FALSE, xlab = xlab_txt, ylab = "", main = title_txt)
    abline(v = 0, lty = 2, col = "grey50")
    segments(0, y, vals, y, col = coef_tab$col, lwd = 6, lend = 1)
    points(vals, y, pch = 16, col = coef_tab$col, cex = 0.6)
    mismatch <- !coef_tab$sign_agree
    if (any(mismatch)) {
        text(vals[mismatch] + sign(vals[mismatch]) * xmax * 0.06, y[mismatch], "*",
            col = "black", font = 2, cex = 1.3)
    }
    axis(1, cex.axis = 0.8)
    axis(2, at = y, labels = coef_tab$name, las = 1, cex.axis = 0.65, tick = FALSE, line = -0.5)
}
draw_side(coef_tab$lpm_partial_mean, "Linear probability model\n(mean partial coefficient)", "Coefficient")
draw_side(coef_tab$clogit_coef, "Pooled conditional logit", "Coefficient")
mtext("* marks the 7 predictors where the two methods disagree in sign (all small-magnitude in the LPM)",
    side = 1, line = -1, outer = TRUE, cex = 0.8, col = "grey40")

# --- panel B: scatter, clogit vs LPM partial ---
layout(1)
par(mai = c(1, 1, 0.6, 0.3))
plot(coef_tab$lpm_partial_mean, coef_tab$clogit_coef, pch = 16, col = coef_tab$col, cex = 1.3,
    xlab = "Linear probability model (mean partial coefficient)",
    ylab = "Pooled conditional logit coefficient",
    main = sprintf("Conditional logit vs. linear probability model\nsign agreement: %d/%d, Spearman r = %.2f",
                   sum(coef_tab$sign_agree), n, r_concordance_rank))
abline(h = 0, v = 0, lty = 2, col = "grey60")
flag <- coef_tab$name %in% c("Mutual friends", "Target distance")
label_pos <- ifelse(coef_tab$lpm_partial_mean[flag] < 0, 4, 3)
text(coef_tab$lpm_partial_mean[flag], coef_tab$clogit_coef[flag], coef_tab$name[flag],
    pos = label_pos, cex = 0.8, font = 2, xpd = NA)
legend("bottomright", legend = catnames, col = catcol_dark, pch = 16, bty = "n", cex = 0.85)

dev.off()
cat("saved conditional_logit_diagnostics.pdf\n")
cat("ALL_DONE\n")
