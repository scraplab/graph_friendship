# ===========================================================================
# merge_features.R
#
# Merges fulldata.csv (expand_data.R's output) and extra_features.csv
# (expand_extra_features.R's output) into one file, alldata.csv -- the
# single input analyze_nets.R and neural_network_models.ipynb now read,
# instead of each independently loading and merging the two files
# themselves (which they used to do, redundantly, in two different
# languages).
#
# Also fixes the one known data-quality issue in this data: clust/uclust/
# chuclust (local clustering coefficient, and its self/change counterparts)
# are NaN for any node with degree < 2, since igraph's
# transitivity(type="local") divides by the number of possible triangles
# among a node's neighbors -- zero when there are fewer than 2 of them, so
# 0/0 -- which is undefined, not a real "structural" NaN. Filled to 0 below
# (no neighbors, or only one, means no possible closure), the same
# convention expand_extra_features.R already uses for its analogous aa
# deg<=1 guard. Confirmed via df[numeric_features].isna().sum() in the
# notebook that these are the ONLY three columns with any missing values at
# all -- everything else (including constr/uconstr/chuconstr, which IS
# defined down to degree 1, just not degree 0) is complete. The fill below
# is written generically (any numeric, non-identifier column with NA gets
# this treatment, reported as it happens) rather than hardcoding just these
# three names, so a future feature round that reintroduces this same
# failure mode doesn't silently slip through unfilled.
#
# Column order is reorganized into one coherent scheme (identifiers ->
# outcome -> demographics -> network baseline -> self baseline -> target ->
# dyad -> network change -> self change -> visual), matching the category
# structure both analyze_nets.R's predictor_registry and the notebook's
# FEATURE_GROUPS already use, instead of "however merge()/cbind() happened
# to leave things" (the previous arrangement: fulldata.csv's original 38
# columns, un-reordered, followed by two later-appended batches of new
# columns in whatever order each expansion script happened to compute
# them).
#
# Does NOT z-score -- that stays in analyze_nets.R (R) and the
# z-score-normalization cell (Python), same as before; both now just start
# from this one file instead of each redoing the same merge separately.
# ===========================================================================

dat   <- read.csv("fulldata.csv")
extra <- read.csv("extra_features.csv")

shared_cols <- intersect(names(dat), names(extra))
stopifnot(setequal(shared_cols, c("rid", "sid", "trial", "netid", "id")))
alldata <- merge(dat, extra, by = shared_cols, sort = FALSE)
stopifnot(nrow(alldata) == nrow(dat))   # inner join should not drop or duplicate rows -- see
                                        # extra_features.csv's header note: it deliberately
                                        # keeps rows (friends) that fulldata.csv already
                                        # filtered out, precisely so this join drops them for
                                        # free, matching 1:1 on the surviving rid set

# ---- zero-fill any numeric, non-identifier column that has NA -------------
id_cols <- c("rid", "sid", "trial", "netid", "unode", "id")
for (col in names(alldata)) {
    if (col %in% id_cols) next
    x <- alldata[[col]]
    n_na <- sum(is.na(x))
    if (n_na == 0) next
    if (!is.numeric(x)) {
        warning(col, ": ", n_na, " NA in a non-numeric column -- not auto-filled, ",
                "needs its own handling")
        next
    }
    cat(sprintf("%s: filling %d NA -> 0 (of %d rows)\n", col, n_na, nrow(alldata)))
    alldata[[col]][is.na(x)] <- 0
}
# every numeric, non-identifier column should now be complete (non-numeric
# columns, e.g. demographics, are exempt -- a real "declined to answer" NA
# there is legitimate data, not a computation artifact)
numeric_check_cols <- setdiff(names(alldata), id_cols)
stopifnot(all(sapply(alldata[numeric_check_cols], function(x)
    !is.numeric(x) || sum(is.na(x)) == 0)))

# ---- coherent column order -------------------------------------------------
col_order <- c(
    # identifiers / trial structure
    "rid", "sid", "trial", "netid", "unode", "id", "nnode",
    # outcome
    "response", "rt", "binres",
    # demographics
    "age", "race", "ethnicity", "sex", "education", "nationality", "gender",
    # network (baseline / context)
    "density", "assortativity", "apl", "modularity", "netcent", "netvar",
    # self (baseline / context)
    "udegree", "ucloseness", "ubetweenness", "ueigenvector",
    "uclust", "uconstr", "uecc", "uesize", "ucomsize",
    # target
    "degree", "closeness", "betweenness", "eigenvector",
    "clust", "constr", "ecc", "esize", "comsize",
    # dyad
    "samecom", "edgbet", "udists", "mutfr", "dhom", "structsim", "aa", "katz", "srcos",
    # network change
    "chgmod", "chgast", "chapl", "chgcent", "chgvar",
    # self change
    "chub", "chuc", "chue", "chuclust", "chuconstr", "chuecc", "chuesize",
    # visual (layout-derived)
    "vdist", "vecc", "uvecc"
)
stopifnot(setequal(col_order, names(alldata)))   # every column placed, nothing dropped/misspelled
alldata <- alldata[, col_order]

write.csv(alldata, "alldata.csv", row.names = FALSE)
cat(sprintf("\nwrote alldata.csv: %d rows, %d columns\n", nrow(alldata), ncol(alldata)))
