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
setwd("D:/Dropbox/mysocialbrain/net/Study2/analyze")
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
cat("ALL_DONE\n")
