
# load packages -----------------------------------------------------------

if(!require(RColorBrewer)) install.packages("RColorBrewer"); require(RColorBrewer)
if(!require(psych)) install.packages("psych"); require(psych)
if(!require(MASS)) install.packages("MASS"); require(MASS)

# load data ---------------------------------------------------------------

# CHANGED: reads alldata.csv -- merge_features.R's already-merged, already-
# NaN-filled, coherently-column-ordered union of fulldata.csv and
# extra_features.csv -- instead of merging the two files itself. Both this
# script and the Python notebook used to redo that same merge independently;
# merge_features.R now does it once, upstream, for both. It has also already
# 0-filled clust/uclust/chuclust (NaN for any degree<2 node -- see that
# script's header comment), which this script never handled before.
dat <- read.csv("alldata.csv")
head(dat)

# z-score metric variables ------------------------------------------------
# CHANGED: every column referenced BY NAME rather than by position
# (dat[,17:24] etc.). alldata.csv's column order is no longer the same as
# fulldata.csv's -- merge_features.R reorganized it into one coherent
# scheme -- so the old position ranges would silently z-score the wrong
# columns now. This list is every numeric predictor column except the
# binary samecom (not a metric variable) and the identifier/outcome/
# demographic columns handled elsewhere; it covers the exact same 47
# columns the old four-part position/name patchwork did, just name-based
# and reorder-proof.
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
dat[, zscore_cols] <- scale(dat[, zscore_cols])

# ---------------------------------------------------------------------------
# Predictor registry: one row per analyzed predictor -- column name in `dat`,
# display name, category, and two flags (context/visual). Built as a single
# table rather than hand-maintained parallel vectors (ponames/catidx/varind
# used to each be redefined separately in two different figure sections;
# that duplication is exactly how 7 columns got missed a moment ago). Defined
# once here and reused everywhere below.
#
# context: TRUE for anything constant across every candidate within a trial
#   (self-node or whole-network baselines, e.g. udegree). These CANNOT in
#   principle have a linear effect on which candidate is chosen within a
#   trial -- any apparent association with binres in a regression pooled
#   across a subject's trials would be a between-trial artifact, not a real
#   "this candidate over that one" effect. CHANGED: now excluded from every
#   analysis in this script entirely (analyzed separately, at the
#   participant level, in the Python notebook instead) -- previously only
#   the *original* context columns (density, assortativity, apl, modularity,
#   udegree, ucloseness, ubetweenness, ueigenvector -- columns 21-29, never
#   referenced by varind to begin with) were excluded; the new ones added
#   since (self/network baselines from both the expand_data.R extra round
#   and expand_extra_features.R) were not.
# visual: TRUE for the layout/visualization-derived predictors. CHANGED:
#   excluded from the individual-differences PCA itself below, and the
#   substantive predictors going into that PCA have subjects' unique/partial
#   visual coefficients partialed out of their zero-order coefficients
#   first -- so visual predictors still get zero-order/unique coefficients
#   and appear in the main coefficient figure, just not the PCA.
# ---------------------------------------------------------------------------
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

# category assignment by name lookup (not by position) -- deliberately
# avoids a fourth hand-typed parallel vector that has to stay aligned with
# `col` above purely by row order
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

predictor_registry <- predictor_registry[!predictor_registry$context, ]   # CHANGED: drop context entirely

varind   <- match(predictor_registry$col, names(dat))
vn       <- length(varind)
ponames  <- predictor_registry$name
catnames <- c("Target", "Dyad", "Network", "Self")
catidx   <- match(predictor_registry$cat, catnames)

# ADDED: substantive-only (non-visual) versions, for the PCA below -- L/ev
# from a PCA fit on visual-excluded, residualized data have 30 rows/values,
# not vn's 32, so anything indexing against that PCA's output needs these,
# not the full ponames/catidx/vn.
visualidx           <- predictor_registry$visual
ponames_substantive <- ponames[!visualidx]
catidx_substantive  <- catidx[!visualidx]
vn_substantive       <- length(ponames_substantive)

# ADDED: for each substantive predictor's per-subject zero-order coefficient,
# regress out that subject's own unique/partial coefficient(s) for the
# visual control predictors, keeping the residuals -- removes any between-
# subjects variance in a substantive predictor's apparent zero-order effect
# that's really just that subject's sensitivity to the visual/layout
# confounds, before it feeds a PCA. Takes the full zcoef/ucoef (context
# already excluded, visual included) and returns a subjects x
# vn_substantive matrix of residuals.
partial_out_visual <- function(zc, uc) {
    uc_visual <- uc[, visualidx, drop = FALSE]
    zc_sub    <- zc[, !visualidx, drop = FALSE]
    out <- apply(zc_sub, 2, function(y) resid(lm(y ~ uc_visual)))
    colnames(out) <- ponames_substantive
    out
}

# compute summary statistics ----------------------------------------------

uid <- unique(dat$sid)
ns <- length(uid)
# CHANGED: varind/vn now come from predictor_registry above (context
# variables excluded, both new and old); no longer redefined here
zcoef <- matrix(NA,ns,vn)
ucoef <- matrix(NA,ns,vn)
# ADDED: marginal-over-visual coefficient -- for each non-visual predictor,
# a regression of binres on that predictor PLUS only the visual predictors
# already in this (linear-regression) analysis (vdist, vecc -- see
# visualidx/predictor_registry above; NOT uvecc, which is excluded from
# this whole script as a context variable, unlike in the neural-network
# notebook, where it's part of the 'visual' group). Controls for the
# visual-layout confound specifically, not for the other graph-theoretic
# predictors -- parallels the notebook's "marginal benefit over visual"
# analysis (Analysis 1), just at the per-predictor grain this script uses
# instead of per feature-group. Left NA for the visual predictors
# themselves (vdist, vecc): there's no "controlling for visual" version of
# a visual predictor, same as the notebook leaves 'visual' blank in its own
# marginal-over-visual panel.
mcoef <- matrix(NA,ns,vn)
visual_varind <- varind[visualidx]
start <- Sys.time()
for (i in 1:ns){

    # select subject data
    sel <- dat$sid==uid[i]
    cdat <- dat[sel,]

    # zero-order
    for (j in 1:vn){
        fit <- lm(binres~as.matrix(cdat[,varind[j]]),data=cdat)
        zcoef[i,j] <- fit$coefficients[2]
    }

    # marginal over visual
    for (j in which(!visualidx)){
        fit <- lm(binres~as.matrix(cdat[,c(varind[j],visual_varind)]),data=cdat)
        mcoef[i,j] <- fit$coefficients[2]   # predictor j's own coefficient (first column after intercept)
    }

    # unique
    fit <- lm(binres~as.matrix(cdat[,varind]),data=cdat)
    ucoef[i,] <- fit$coefficients[-1]


    if (i %% 10 == 0){
        cat("\014")
        print(i/ns)
        print(Sys.time() -start)
    }
}

# inferential testing -----------------------------------------------------

# mean zero order coefficients
mzc <- apply(zcoef,2,mean)
mzc

# zero-order p-values
zp <- apply(zcoef,2,function(x) t.test(x)$p.value)
zpc <- p.adjust(zp,method = "holm")
zpc

# zero-order cohen's d values
zd <- apply(zcoef,2,function(x) t.test(x)$statistic/sqrt(ns))
zd

# pairwise comparisons
zpmat <- matrix(NA,vn,vn)
for (i in 1:(vn-1)){
    for (j in (i+1):vn){
        zpmat[i,j] <- t.test(zcoef[,i]-zcoef[,j])$p.value
    }
}
zpmatc <- zpmat
zpmatc[!is.na(zpmat)]<-p.adjust(zpmatc[!is.na(zpmat)])


for (i in 1:(vn-1)){
    for (j in (i+1):vn){
        if (zpmatc[i,j]>.05){
            print(paste(ponames[i],ponames[j]))
        }
    }
}
# all pairwise comparisons significant except:
# "Change modularity Change average path length"

# mean unique coefficients
muc <- apply(ucoef,2,mean)
muc

# unique p-values
up <- apply(ucoef,2,function(x) t.test(x)$p.value)
upc <- p.adjust(up,method = "holm")
upc

# unique cohen's d values
ud <- apply(ucoef,2,function(x) t.test(x)$statistic/sqrt(ns))
ud

# pairwise comparisons
upmat <- matrix(NA,vn,vn)
for (i in 1:(vn-1)){
    for (j in (i+1):vn){
        upmat[i,j] <- t.test(ucoef[,i]-ucoef[,j])$p.value
    }
}
upmatc <- upmat
upmatc[!is.na(upmat)]<-p.adjust(upmatc[!is.na(upmat)])
mean(upmatc<.05,na.rm=T)


for (i in 1:(vn-1)){
    for (j in (i+1):vn){
        if (upmatc[i,j]>.05){
            print(paste(ponames[i],ponames[j]))
        }
    }
}

# mean marginal-over-visual coefficients (NA for vdist/vecc themselves)
mmc <- apply(mcoef,2,mean,na.rm=TRUE)
mmc

# marginal-over-visual p-values -- t.test() can't run on an all-NA column
# (vdist/vecc, which mcoef never fills in), so those positions are passed
# through as NA rather than erroring; p.adjust() in turn leaves NA inputs
# as NA in its output and Holm-adjusts only the non-NA entries (confirmed
# via a quick standalone test), so mpc stays length vn with NA exactly at
# the 2 visual-predictor positions.
mp  <- apply(mcoef,2,function(x) if (all(is.na(x))) NA else t.test(x)$p.value)
mpc <- p.adjust(mp,method = "holm")
mpc

# marginal-over-visual cohen's d values
md <- apply(mcoef,2,function(x) if (all(is.na(x))) NA else t.test(x)$statistic/sqrt(ns))
md

# pairwise comparisons -- only between non-visual predictors (mcoef[,j] is
# all-NA for vdist/vecc, which breaks t.test() on the difference; those
# rows/columns of mpmat are left NA, same as the unused upper-triangle
# positions already are)
mpmat <- matrix(NA,vn,vn)
for (i in 1:(vn-1)){
    for (j in (i+1):vn){
        if (visualidx[i] || visualidx[j]) next
        mpmat[i,j] <- t.test(mcoef[,i]-mcoef[,j])$p.value
    }
}
mpmatc <- mpmat
mpmatc[!is.na(mpmat)]<-p.adjust(mpmatc[!is.na(mpmat)])
mean(mpmatc<.05,na.rm=T)


for (i in 1:(vn-1)){
    for (j in (i+1):vn){
        if (!is.na(mpmatc[i,j]) && mpmatc[i,j]>.05){
            print(paste(ponames[i],ponames[j]))
        }
    }
}

# split half reliabilty ---------------------------------------------------

uid <- unique(dat$sid)
ns <- length(uid)
ntr <- length(unique(dat$trial))
# varind/vn: see predictor_registry above
zcoef1 <- matrix(NA,ns,vn)
ucoef1 <- matrix(NA,ns,vn)
zcoef2 <- matrix(NA,ns,vn)
ucoef2 <- matrix(NA,ns,vn)
set.seed(1)
start <- Sys.time()
for (i in 1:ns){
    
    # select subject data
    sel <- dat$sid==uid[i]
    cdat <- dat[sel,]
    sel <- sample(ntr,ntr/2,F)
    sel <- 1:ntr %in% sel
    cdat1 <- cdat[sel,]
    cdat2 <- cdat[!sel,]
    
    # zero-order
    for (j in 1:vn){
        fit <- lm(binres~as.matrix(cdat1[,varind[j]]),data=cdat1)
        zcoef1[i,j] <- fit$coefficients[2]
        fit <- lm(binres~as.matrix(cdat2[,varind[j]]),data=cdat2)
        zcoef2[i,j] <- fit$coefficients[2]
    }
    
    # unique
    fit <- lm(binres~as.matrix(cdat1[,varind]),data=cdat1)
    ucoef1[i,] <- fit$coefficients[-1]
    
    fit <- lm(binres~as.matrix(cdat2[,varind]),data=cdat2)
    ucoef2[i,] <- fit$coefficients[-1]
    
    if (i %% 10 == 0){
        cat("\014")
        print(i/ns)
        print(Sys.time() -start)
    }
}

spcz <- diag(cor(zcoef1,zcoef2))
spcz*2/(1+spcz)

spcu <- diag(cor(ucoef1,ucoef2))
spcu*2/(1+spcu)

# CHANGED: partial_out_visual() before the split-half PCA too, so this
# reliability check reflects the same visual-confound-removed, visual-
# excluded PCA that's actually used everywhere else below, not a plain PCA
# on the raw (unpartialed, visual-included) zero-order coefficients.
spcp <- diag(cor(principal(partial_out_visual(zcoef1, ucoef1), 2)$scores,
               principal(partial_out_visual(zcoef2, ucoef2), 2)$scores))
spcp*2/(1+spcp)


# individual differences --------------------------------------------------

# CHANGED: partial out subjects' unique visual-control coefficients from the
# substantive predictors' zero-order coefficients, and exclude the visual
# predictors themselves, before any of the PCA-related diagnostics/fit below
# -- see predictor_registry's comments and partial_out_visual() above.
zcoef_clean <- partial_out_visual(zcoef, ucoef)

# VSS 1 and 2 favor 2 factor
vss(zcoef_clean)

# fit PCA
p <- principal(zcoef_clean,2)

# compute demographics
# CHANGED: was dat[,8:14] -- a position range that only worked because
# fulldata.csv happened to put these 7 columns there; alldata.csv's
# reorganized column order breaks that, so this is now by name instead.
demo_cols <- c("age", "race", "ethnicity", "sex", "education", "nationality", "gender")
udem <- aggregate(dat[,demo_cols],list(dat$sid),function(x) x[1])

summary(udem$age)
table(udem$gender)
table(udem$sex)
table(udem$race)
table(udem$ethnicity)
table(udem$education)
sort(table(udem$nationality))

save(zcoef,ucoef,udem,spcz,spcu,spcp,file="coefficients.Rdata")






















# =================================================================
# Demographics cleaning + PC/demographic inference + combined figure.
#
# Runs immediately after:
#     udem <- aggregate(dat[,demo_cols], list(dat$sid), function(x) x[1])
#
# Expects in the environment: dat, uid, zcoef, p (psych::principal fit), udem
# Requires RColorBrewer, psych, MASS (loaded at the top of your script).
#
# Produces: udemc (cleaned demographics), dem_audit, lev_stats, res, age_res,
#           pca_scree_loadings.pdf (panels A-D) and pca_demographics.pdf
#           (panels E-G), both 6.5in wide.
# =================================================================

# -----------------------------------------------------------------
# 0. sid column and ROW ALIGNMENT
#
# aggregate(x, list(dat$sid), ...) names the grouping column "Group.1",
# and returns rows sorted by sid. zcoef / p$scores rows follow
# unique(dat$sid), i.e. first-appearance order. If those orders differ,
# every demographic result below is scrambled with no error raised, so
# the reorder is done explicitly and then asserted.
# -----------------------------------------------------------------
names(udem)[1] <- "sid"
udem <- udem[match(uid, udem$sid), ]
rownames(udem) <- NULL
stopifnot(identical(as.character(udem$sid), as.character(uid)))
stopifnot(nrow(udem) == nrow(p$scores))

# =================================================================
# 1. CLEANING
# =================================================================

min_n         <- 10            # smallest group retained
decline_codes <- c(0, "0")     # 0 = declined to answer
age_valid     <- c(18, 100)

# Canonical key: strip non-ASCII (kills both "â€™" and "’"), DROP
# punctuation (so "Master's", "Master’s" and "Masterâ€™s" all land on
# "masters degree"), lowercase, squeeze whitespace. Every recode and
# ordering rule below is keyed on this form.
canon <- function(x) {
    y <- as.character(x)
    y <- gsub("[^\x20-\x7E]", "", y, useBytes = TRUE)
    y <- gsub("[[:punct:]]", "", y)
    tolower(trimws(gsub("[[:space:]]+", " ", y)))
}

# -----------------------------------------------------------------
# SUBSTANTIVE recodes, keyed on canonical labels, applied BEFORE any
# size-based policy.
#
# This is where category merges belong. merge_adjacent only merges
# levels that fall BELOW min_n, which is why "Elementary" (10) and
# "None" (36) previously survived as their own levels: both cleared
# the threshold, so the size-based rule had nothing to do. Deciding
# that Elementary and None belong with "high school, no degree" is a
# substantive judgement, not a statistical one, so it is stated here
# explicitly and does not depend on any cell count.
# -----------------------------------------------------------------
recode_map <- list(
    gender = list(
        # "Other" (34) folded into Non-binary
        "Non-binary" = c("nonbinary", "non binary", "genderqueer",
                         "gender queer", "other")
    ),
    race = list(
        "Other Asian" = c("japanese", "korean", "vietnamese"),
        "Indigenous American, Alaskan, or Pacific Islander" = c(
            "guamanian or chamorro", "native hawaiian",
            "american indian or alaska native", "samoan",
            "other pacific islander"),
        # "Latino" (15) folded into the residual race category
        "Other" = c("latino", "hispanic", "some other race", "other race",
                    "other")
    ),
    education = list(
        # Elementary + None + high-school-no-degree, named as one level
        "Some high school or less" = c("elementary", "none", "no schooling",
                                       "less than high school",
                                       "highschoolnodegree",
                                       "high school no degree",
                                       "some high school"),
        # the two Master's variants both canonicalise to "masters degree"
        "Graduate degree" = c("doctorate degree", "professional school degree",
                              "masters degree")
    )
)

# Levels to DROP outright (set NA), keyed on canonical labels. Distinct
# from a recode (which relabels) and from the size policy (which is a
# threshold rule): "Other" for sex is uninterpretable rather than small.
drop_map <- list(
    sex = c("other")
)

# Education is ordinal and merge_adjacent folds NEIGHBOURING levels, so the
# order must be stated. Canonical keys, lowest to highest.
edu_order <- c(
    "some high school or less",
    "high school graduate", "high school diploma or equivalent",
    "some college no degree", "some college",
    "associate degree", "associates degree",
    "bachelors degree",
    "graduate degree"
)

dem_policy <- list(
    age         = list(type = "numeric"),
    gender      = list(type = "nominal", small = "na"),
    sex         = list(type = "nominal", small = "na"),
    race        = list(type = "nominal", small = "na"),
    ethnicity   = list(type = "nominal", small = "na"),
    education   = list(type = "ordinal", small = "merge_adjacent"),
    nationality = list(type = "nominal", small = "binary_mode")
)

to_na <- function(x, codes) { x[x %in% codes] <- NA; x }

# Collapse labels differing only by encoding/punctuation. Display label for
# each group is the most frequent pure-ASCII variant; if all are mangled,
# the canonical form is used.
dedupe_levels <- function(x) {
    ch  <- as.character(x)
    key <- canon(ch)
    dsp <- tapply(ch, key, function(z) {
        clean <- z[!grepl("[^\x20-\x7E]", z, useBytes = TRUE)]
        pool  <- if (length(clean)) clean else z
        names(sort(table(pool), decreasing = TRUE))[1]
    })
    list(key = key, label = unname(dsp[key]))
}

apply_recode <- function(key, label, map) {
    if (is.null(map)) return(list(key = key, label = label, hits = NULL))
    hits <- list()
    for (target in names(map)) {
        sel <- key %in% map[[target]]
        hits[[target]] <- sort(unique(key[sel]))
        label[sel] <- target
        key[sel]   <- canon(target)
    }
    list(key = key, label = label, hits = hits)
}

small_to_na <- function(f, min_n) {
    tb <- table(f)
    f[f %in% names(tb)[tb < min_n]] <- NA
    droplevels(f)
}

# Repeatedly merge the smallest level into its smaller-count neighbour.
# Requires levels(f) to already be in substantive order. NOTE: only acts on
# levels below min_n -- anything above it is left alone, by design.
merge_adjacent <- function(f, min_n) {
    f   <- droplevels(f)
    lv  <- levels(f)
    grp <- setNames(seq_along(lv), lv)
    repeat {
        tb <- table(factor(grp[as.character(f)], levels = sort(unique(grp))))
        if (length(tb) <= 2 || min(tb) >= min_n) break
        i   <- which.min(tb)
        ids <- as.integer(names(tb))
        nb  <- if (i == 1) 2 else if (i == length(tb)) i - 1 else
            if (tb[i-1] <= tb[i+1]) i - 1 else i + 1
        grp[grp == ids[i]] <- ids[nb]
    }
    newlab <- tapply(names(grp), grp, function(z) {
        z <- z[order(match(z, lv))]
        if (length(z) == 1) z[1] else paste0(z[1], " / ", z[length(z)])
    })
    factor(newlab[as.character(grp[as.character(f)])],
           levels = newlab[order(as.numeric(names(newlab)))])
}

binary_mode <- function(f, min_n) {
    f  <- droplevels(f)
    tb <- sort(table(f), decreasing = TRUE)
    if (length(tb) < 2) return(f)
    modal <- names(tb)[1]
    out <- factor(ifelse(as.character(f) == modal, modal, "Other"),
                  levels = c(modal, "Other"))
    if (min(table(out)) < min_n) out[] <- NA   # contrast not estimable
    droplevels(out)
}

udemc     <- udem
dem_audit <- list()

for (v in names(dem_policy)) {
    if (!v %in% names(udem)) { warning("column not found: ", v); next }
    pol <- dem_policy[[v]]
    raw <- udem[[v]]
    
    if (pol$type == "numeric") {
        x <- to_na(as.numeric(as.character(raw)), decline_codes)
        x[x < age_valid[1] | x > age_valid[2]] <- NA
        udemc[[v]] <- x
        dem_audit[[v]] <- list(type = "numeric", n_na = sum(is.na(x)),
                               summary = summary(x))
        next
    }
    
    ch <- to_na(as.character(raw), decline_codes)
    ch[trimws(ch) == ""] <- NA
    n_declined <- sum(is.na(ch))
    ok <- !is.na(ch)
    
    dd <- dedupe_levels(ch[ok])
    n_dupe_groups <- sum(tapply(ch[ok], dd$key,
                                function(z) length(unique(z))) > 1)
    rc <- apply_recode(dd$key, dd$label, recode_map[[v]])
    
    lab <- key <- rep(NA_character_, length(ch))
    lab[ok] <- rc$label
    key[ok] <- rc$key
    
    # explicit drops, after the recodes so a recode target can still be kept
    n_dropped <- 0L
    if (!is.null(drop_map[[v]])) {
        sel <- !is.na(key) & key %in% drop_map[[v]]
        n_dropped <- sum(sel)
        lab[sel] <- NA_character_; key[sel] <- NA_character_
    }
    
    if (pol$type == "ordinal" && v == "education") {
        obs_keys <- unique(na.omit(key))
        unknown  <- setdiff(obs_keys, edu_order)
        if (length(unknown))
            warning("education: level(s) missing from edu_order, merge order ",
                    "is unreliable for: ", paste(unknown, collapse = " | "))
        ord_keys <- c(intersect(edu_order, obs_keys), unknown)
        lev <- unname(sapply(ord_keys, function(kk) lab[which(key == kk)[1]]))
        f0  <- factor(lab, levels = unique(lev))
    } else {
        f0 <- factor(lab, levels = names(sort(table(lab), decreasing = TRUE)))
    }
    before <- table(f0, useNA = "no")
    
    f1 <- switch(pol$small,
                 na             = small_to_na(f0, min_n),
                 merge_adjacent = merge_adjacent(f0, min_n),
                 binary_mode    = binary_mode(f0, min_n),
                 stop("unknown policy: ", pol$small))
    
    udemc[[v]] <- f1
    dem_audit[[v]] <- list(type = pol$type, policy = pol$small,
                           n_declined = n_declined, n_dupe_groups = n_dupe_groups,
                           n_dropped = n_dropped,
                           recodes = rc$hits, before = before,
                           after = table(f1, useNA = "no"),
                           n_na_total = sum(is.na(f1)),
                           still_below_min_n = names(table(f1))[table(f1) < min_n])
}

# ---- cleaning report ----
for (v in names(dem_audit)) {
    a <- dem_audit[[v]]
    cat("\n==== ", v, " (", a$type,
        if (!is.null(a$policy)) paste0(", ", a$policy) else "", ")\n", sep = "")
    if (a$type == "numeric") { print(a$summary); cat("n NA:", a$n_na, "\n"); next }
    if (a$n_dupe_groups > 0)
        cat("encoding duplicates collapsed:", a$n_dupe_groups, "label group(s)\n")
    for (tg in names(a$recodes))
        if (length(a$recodes[[tg]]))
            cat("recoded -> '", tg, "': ", paste(a$recodes[[tg]], collapse = ", "),
                "\n", sep = "")
    if (a$n_dropped > 0) cat("dropped outright:", a$n_dropped, "participant(s)\n")
    cat("before small-group policy:\n"); print(a$before)
    cat("after:\n");                     print(a$after)
    cat("declined:", a$n_declined, " total NA:", a$n_na_total, "\n")
    if (length(a$still_below_min_n))
        cat("WARNING: still below min_n:",
            paste(a$still_below_min_n, collapse = ", "), "\n")
}
cat("\n---- usable N per variable ----\n")
print(sapply(udemc[names(dem_policy)], function(x) sum(!is.na(x))))

# =================================================================
# 2. FIGURE SETUP: ordering, display labels, sizes
# =================================================================

# CHANGED: ponames/vn/catidx/catnames (and the substantive-only/visual
# variants used by the PCA panels further down) now come from
# predictor_registry, defined once near the top of the script instead of
# being redefined here and again in "data visualization" below -- see that
# block's comments for why (two independently hand-maintained copies is how
# 7 columns got missed in an earlier round).
stopifnot(vn == length(varind))

demvars <- c("gender", "sex", "ethnicity", "race", "nationality", "education")
demlabs <- c("Gender", "Sex", "Ethnicity", "Race", "Nationality", "Education")

tick_cex   <- 0.75   # all tick labels
axlab_cex  <- 0.9    # all axis titles
title_cex  <- 1.0    # panel main titles
group_cex  <- 0.9    # demographic variable names
star_cex   <- 1.3
panel_cex  <- 1.4    # panel letters -- largest text on the page
# Label distances, in INCHES. Only the primitives are set here; the x-title
# distance and the vertical margins are DERIVED from them further down, to
# satisfy three constraints exactly:
#   (1) y tick labels sit ax_in_tick from the y axis (unchanged)
#   (2) x tick labels sit the same physical distance from the x axis -- mgp[2]
#       governs both axes at once, so this holds by construction
#   (3) gap(x tick labels -> x-axis title) = 0.5 * gap(x-axis title -> the
#       panel title of the row below), identically in every panel
# mtext() CENTRES text on `line`, so an item at distance d with height h
# occupies d +/- h/2; the gaps below are edge-to-edge.
ax_in_tick     <- 0.07   # BOTH axes: tick labels to axis (mgp[2], unchanged)
gap_tick_title <- 0.06   # tick labels -> x-axis title, edge to edge
gap_title_next <- 0.12   # x-axis title -> panel title of the row below (= 2x)
ax_in_ytitle_a <- 0.22   # A: eigenvalue ticks are 1 char wide, so this can be close
ax_in_ytitle_b <- 0.28   # B: score ticks are wider
ax_in_ptitle   <- 0.06   # gap from a plot's top edge up to its panel title
# Asterisk ink sits high in its em box, so adj = 0.5 still reads high.
star_dy_frac <- 0.34

baseline <- "zero"   # "zero" = level vs 0; "ref" = level vs first displayed level

# CHANGED: reverted to 4 tab20 colors (Target=green, Dyad=orange,
# Network=purple, Self=blue) -- context predictors (which briefly justified
# a 6-way context/change split, each family in both a pink/teal "context"
# and purple/blue "change" shade) are now excluded from this script
# entirely (see predictor_registry above), so there's no more "context"
# sibling for Network/Self to be told apart from here. The pink/teal pair
# is simply unused in this script now; the Python notebook still analyzes
# context predictors and keeps the full 6-color scheme.
catcol    <- c("#2CA02C", "#FF7F0E", "#9467BD", "#1F77B4")
# col_pc (E/F/G panels) stays on tab20's red pair -- still not one of the
# 4 colors in catcol above, so PC1/PC2 still don't collide with anything.
col_pc    <- c("#D62728", "#FF9896")
dens_ramp <- colorRampPalette(c("white", brewer.pal(9, "YlOrBr")))

# ---- PCA quantities ----
colnames(zcoef) <- ponames   # zcoef itself stays full-width (visual included, for
                              # the main coefficient figure); only the PCA (p, L, S,
                              # ev, fp below) was fit on zcoef_clean (visual-excluded,
                              # residualized) and so only spans vn_substantive columns
L  <- unclass(p$loadings)[, 1:2, drop = FALSE]
S  <- p$scores[, 1:2]
ve <- p$Vaccounted["Proportion Var", 1:2]

# Component signs are arbitrary, but the previous rule ("make the largest
# loading positive") flipped RC1, whose largest loading is
# "Change average path length" at -0.91 -- leaving the displayed component
# oriented as "prefers close, same-community targets". Orient each component
# on a named anchor instead, so the sign is interpretable. Set orient_by to
# NULL to restore the largest-loading rule.
orient_by <- c("Adamic-Adar", "Target degree")
for (k in 1:2) {
    j <- if (is.null(orient_by)) which.max(abs(L[, k])) else
        which(ponames_substantive == orient_by[k])   # CHANGED: was ponames -- L has
                                                       # vn_substantive rows, not vn
    if (L[j, k] < 0) { L[, k] <- -L[, k]; S[, k] <- -S[, k] }
}

ev <- eigen(cor(zcoef_clean), only.values = TRUE)$values   # CHANGED: was zcoef
fp <- fa.parallel(zcoef_clean, fa = "pc", plot = FALSE)     # CHANGED: was zcoef
sim_mean <- if (!is.null(fp$pc.sim)) fp$pc.sim else fp$pc.simr

stopifnot(nrow(S) == nrow(udemc))

# ---- bimodality of PC1 -----------------------------------------------------
# ADDED back: this had apparently been run before, but there's no trace of
# it left in this script, v1/, or .Rhistory, so this is a fresh addition
# (the standard choices for this question) rather than a reconstruction of
# whatever was run previously -- flag if a different test was intended.
# Hartigan's dip test (H0: unimodal) is the standard formal test; the
# bimodality coefficient (Pfister et al., 2013; SAS/JMP convention) is a
# simpler, commonly-reported descriptive companion, with BC > 5/9 (~0.555)
# conventionally read as evidence of bimodality. Both computed on PC1 only,
# as asked -- not PC2.
if (!require(diptest)) install.packages("diptest"); require(diptest)
dip.test(S[,1])

bimodality_coefficient <- function(x) {
    x    <- x[!is.na(x)]
    n    <- length(x)
    m    <- mean(x)
    v    <- sum((x-m)^2)/n
    skew <- (sum((x-m)^3)/n) / v^1.5
    kurt <- (sum((x-m)^4)/n) / v^2   # regular (not excess) kurtosis -- normal = 3
    (skew^2 + 1) / kurt
}
bimodality_coefficient(S[,1])   # > 5/9 ~ 0.556 conventionally read as evidence of bimodality

# ---- level ordering ----
race_groups <- list(
    "Asian" = c("asian", "asian indian", "chinese", "filipino", "japanese",
                "korean", "vietnamese", "south asian", "east asian",
                "southeast asian", "other asian"),
    "Black or African American" = c("black", "black or african american",
                                    "african american"),
    "Indigenous American, Alaskan, or Pacific Islander" = c(
        "american indian or alaska native", "native hawaiian", "samoan",
        "guamanian or chamorro", "other pacific islander",
        "indigenous american alaskan or pacific islander"),
    "Middle Eastern or North African" = c("middle eastern or north african",
                                          "middle eastern", "north african"),
    "White" = c("white", "caucasian")
)
race_group_of <- function(levs) {
    k <- canon(levs); out <- rep(NA_character_, length(levs))
    for (g in names(race_groups)) out[k %in% race_groups[[g]]] <- g
    out[is.na(out)] <- "Other or unspecified"
    out
}

# Keyword ranks, matched with spaces removed so run-together labels
# ("Highschoolnodegree") rank correctly. Test order matters:
# "highschoolgraduate" contains "graduate", so high-school branches come first.
edu_rank <- function(levs) {
    rank1 <- function(s) {
        s2 <- gsub(" ", "", s)
        if (grepl("noschooling|elementary|lessthanhighschool", s2))                       1
        else if (grepl("somehighschool|highschoolnodegree|highschoolnodiploma|none", s2)) 2
        else if (grepl("highschoolgraduate|highschooldiploma|ged|highschool", s2))        3
        else if (grepl("somecollege", s2))                                                4
        else if (grepl("associate|technical|vocational", s2))                             5
        else if (grepl("bachelor", s2))                                                   6
        else if (grepl("graduate|master|doctor|professional", s2))                        7
        else                                                                             99
    }
    sapply(canon(levs), rank1)
}

order_levels <- function(v, levs) {
    if (v == "education") {
        levs[order(edu_rank(levs), canon(levs))]
    } else if (v == "race") {
        g    <- race_group_of(levs)
        gkey <- ifelse(g == "Other or unspecified", "zzz", g)
        # within an umbrella group: alphabetical, but a residual "Other X"
        # sorts last (so South Asian precedes Other Asian)
        levs[order(gkey, grepl("^other", canon(levs)), canon(levs))]
    } else {
        isother <- canon(levs) %in% c("other", "other or unspecified", "some other")
        levs[order(isother, canon(levs))]
    }
}

# ---- display labels: the E/F margin affords ~24 characters ----
short_labels <- c(
    "Indigenous American, Alaskan, or Pacific Islander" = "Indigenous/Pac. Isl.",
    "American Indian or Alaska Native"                  = "Am. Ind./AK Native",
    "Native Hawaiian or Other Pacific Islander"         = "Native Haw./Pac. Isl.",
    "Black or African American"                         = "Black/African Am.",
    "High school diploma or equivalent"                 = "HS diploma/equiv.",
    "Hispanic, Latino, or Spanish origin"               = "Hispanic/Latino",
    "Not Hispanic, Latino, or Spanish origin"           = "Not Hispanic/Latino",
    "Middle Eastern or North African"                   = "MENA",
    "Some college no degree"                            = "Some college"
)
disp <- function(x) {
    x   <- as.character(x)
    hit <- x %in% names(short_labels)
    x[hit] <- short_labels[x[hit]]
    ifelse(nchar(x) > 24, abbreviate(x, minlength = 24), x)
}

star <- function(pv) ifelse(is.na(pv), "",
                            ifelse(pv < .001, "***", ifelse(pv < .01, "**", ifelse(pv < .05, "*", ""))))

# =================================================================
# 3. INFERENCE (analytic throughout)
# =================================================================
lev_order <- list(); lev_stats <- list()
for (i in seq_along(demvars)) {
    v  <- demvars[i]
    g  <- droplevels(factor(udemc[[v]]))
    lv <- order_levels(v, levels(g))
    lev_order[[v]] <- lv
    for (k in 1:2) {
        pv <- est <- setNames(rep(NA_real_, length(lv)), lv)
        for (j in seq_along(lv)) {
            x <- S[!is.na(g) & g == lv[j], k]; x <- x[!is.na(x)]
            if (baseline == "zero") {
                if (length(x) > 2) { tt <- t.test(x)
                pv[j] <- tt$p.value; est[j] <- unname(tt$estimate) }
            } else {
                y <- S[!is.na(g) & g == lv[1], k]; y <- y[!is.na(y)]
                if (j == 1) est[j] <- mean(x)
                else if (length(x) > 2 && length(y) > 2) {
                    tt <- t.test(x, y); pv[j] <- tt$p.value
                    est[j] <- unname(tt$estimate[1] - tt$estimate[2])
                }
            }
        }
        # raw p only; the correction is applied once, pooled, below
        lev_stats[[paste(v, k)]] <- data.frame(level = lv, est = est, p = pv,
                                               row.names = NULL)
    }
}

# Omnibus tests: a SEPARATE family, not the one the asterisks use. aov is
# what you ran; Welch and Kruskal-Wallis are alongside because cells are
# unequal in n and likely in variance, which is where the equal-variance F
# misstates p.
res <- data.frame()
for (k in 1:2) for (i in seq_along(demvars)) {
    v <- demvars[i]; g <- droplevels(factor(udemc[[v]])); y <- S[, k]
    ok <- !is.na(g) & !is.na(y)
    if (length(unique(g[ok])) < 2) next
    a <- summary(aov(y[ok] ~ g[ok]))[[1]]; ss <- a[["Sum Sq"]]
    res <- rbind(res, data.frame(
        pc = k, var = v, n = sum(ok), k_lev = length(unique(g[ok])),
        F_aov = a[["F value"]][1], df1 = a[["Df"]][1], df2 = a[["Df"]][2],
        p_aov = a[["Pr(>F)"]][1], eta2 = ss[1]/sum(ss),
        p_welch = tryCatch(oneway.test(y[ok] ~ g[ok], var.equal = FALSE)$p.value,
                           error = function(e) NA_real_),
        p_kw = kruskal.test(y[ok] ~ g[ok])$p.value))
}
res$p_holm_omnibus <- p.adjust(res$p_aov, method = "holm")

age_res <- data.frame()
for (k in 1:2) {
    ok <- !is.na(udemc$age) & !is.na(S[, k])
    ct <- cor.test(udemc$age[ok], S[ok, k], method = "kendall")
    age_res <- rbind(age_res, data.frame(pc = k, tau = unname(ct$estimate),
                                         p = ct$p.value, n = sum(ok)))
}

# ---- ONE Holm-Bonferroni family: every level test, both components,
#      plus the two age correlations. These adjusted p-values drive the
#      asterisks. Untested levels carry NA; p.adjust drops NAs and takes
#      the family size from the non-NA count. Holm controls FWER under
#      arbitrary dependence, so the correlations among these variables
#      do not invalidate it.
fam_src <- c(unlist(lapply(lev_stats, function(d) d$p), use.names = FALSE),
             age_res$p)
fam_q   <- p.adjust(fam_src, method = "holm")
fam_m   <- sum(!is.na(fam_src))

pos <- 1L
for (nm in names(lev_stats)) {
    ni <- nrow(lev_stats[[nm]])
    lev_stats[[nm]]$p_holm <- fam_q[pos:(pos + ni - 1L)]
    pos <- pos + ni
}
age_res$p_holm <- fam_q[pos:(pos + nrow(age_res) - 1L)]
stopifnot(pos + nrow(age_res) - 1L == length(fam_q))

# Age CI at a level consistent with the family correction. Bonferroni, since
# Holm has no clean simultaneous-interval counterpart: conservative but
# coherent, so an interval excluding zero will not contradict a missing star.
ci_adjust <- TRUE
alpha  <- if (ci_adjust) 0.05 / fam_m else 0.05
zcrit  <- qnorm(1 - alpha/2)
zt     <- atanh(age_res$tau); se_t <- sqrt(0.437 / (age_res$n - 4))
age_res$lo <- tanh(zt - zcrit*se_t)
age_res$hi <- tanh(zt + zcrit*se_t)

cat(sprintf("\nHolm family: %d tests (%d level tests + 2 age taus)\n",
            fam_m, sum(sapply(lev_stats, nrow))))
cat(sprintf("age tau CIs at %.2f%% (%s)\n", 100*(1-alpha),
            if (ci_adjust) "Bonferroni-adjusted to the family" else "nominal"))
print(res, digits = 3); print(age_res, digits = 3)
for (nm in names(lev_stats)) { cat("\n--", nm, "\n"); print(lev_stats[[nm]], digits = 3) }
cat("\nrealised level order:\n"); print(lev_order)

# =================================================================
# 4. FIGURE  (6.5 x 9 in: 8.5 x 11 sheet, 1 in margins)
#
#   row 1  A scree        | B score density
#   row 2  C PC1 loadings | D PC2 loadings
#   row 3  E PC1 scores   | F PC2 scores   (all variables, blocked)
#   row 4  G age (spans all columns)
#
# THREE layout columns, not two. A single layout() grid forces every row
# onto the same column boundaries, which is why C previously inherited
# E's wide label margin: C ended up with ~0.6in of unused space and the
# C/D bar areas were pinned to the E/F split. With a third column of
# 0.54in, row 2 can split at 3.92in (C | D) while row 3 splits at 4.46in
# (E | F), giving each pair equal plot regions of its own.
#
# Margins are specified in INCHES and converted to lines using the
# device's actual par("csi"), read after pdf() opens. Hard-coding lines
# assumes 1/6in per line, and on a device where that is wrong the two
# columns' plot regions come out unequal -- which is what made E
# physically wider than F.
# =================================================================
slot_lab <- character(0); slot_var <- character(0); slot_gap <- logical(0)
for (i in seq_along(demvars)) {
    lv <- lev_order[[demvars[i]]]
    slot_lab <- c(slot_lab, lv)
    slot_var <- c(slot_var, rep(demvars[i], length(lv)))
    slot_gap <- c(slot_gap, rep(FALSE, length(lv)))
    if (i < length(demvars)) {
        slot_lab <- c(slot_lab, ""); slot_var <- c(slot_var, "")
        slot_gap <- c(slot_gap, TRUE)
    }
}
nslot <- length(slot_lab)
ypos  <- nslot - seq_len(nslot) + 1          # slot 1 at the top

panel_letter <- function(lab) {
    text(grconvertX(0.015, from = "nfc", to = "user"),
         grconvertY(0.985, from = "nfc", to = "user"),
         lab, font = 2, cex = panel_cex, adj = c(0, 1), xpd = NA)
}
star_y <- function(y) y - star_dy_frac * strheight("*", cex = star_cex)

fig_width  <- 6.5
fig_height <- 9.0

# margins in INCHES
# Sized to the text they must hold, not padded: level labels need ~1.20in
# (24 chars at tick_cex), the bold variable names ~0.76in ("Nationality" at
# group_cex), predictor names ~1.50in. Everything beyond that was white space
# and is now plot area -- C/D gain 0.075in each, E/F gain 0.29in each.
lab_dem_in  <- 2.25   # E: 1.45in of level labels + 0.80in for the name tier
d_left_in   <- 0.42   # right-hand panels (D, F)
r_in        <- 0.10   # right margin everywhere
grp_in      <- 1.45   # right edge of the variable-name tier, from plot edge
a_left_in   <- 0.62   # A
b_left_in   <- 0.60   # B
g_left_in   <- 0.62   # G: only needs "PC1 (nnnn)", so it spans the full width

# CHANGED: lab_load_in is now measured from the actual (longest) predictor
# name at tick_cex, rather than a hardcoded guess -- the guess (1.60in) was
# sized for the original, shorter predictor names and started clipping the
# longest of the new ones ("Successor representation similarity") once
# ponames_substantive grew to 30 predictors (confirmed by test-rendering
# panel C at the old value). Padded by 0.20in: a plain measured-width pad
# of 0.05in still clipped by about one character in a rendered check --
# likely a small font-metric mismatch between R's built-in metrics and the
# actual PDF renderer's glyph widths -- so the padding is deliberately
# generous rather than exact.
.tmp_lab <- tempfile(fileext = ".pdf")
pdf(.tmp_lab, width = 6.5, height = 9); plot.new()
lab_load_in <- max(strwidth(ponames_substantive, units = "inches", cex = tick_cex)) + 0.20
dev.off(); unlink(.tmp_lab)

# Column widths follow from the inch margins alone -- no font units involved.
W_F <- (fig_width - (lab_dem_in  - d_left_in)) / 2; W_E <- fig_width - W_F
W_D <- (fig_width - (lab_load_in - d_left_in)) / 2; W_C <- fig_width - W_D
col_w <- c(W_C, W_E - W_C, W_F)
stopifnot(all(col_w > 0.05), abs(sum(col_w) - fig_width) < 1e-9)

# ---- derived vertical geometry -------------------------------------------
# Text heights are MEASURED on a throwaway device of the same type and size.
# The axis and panel titles are then drawn with text() at explicit
# coordinates rather than mtext(): mtext(line = L) indexes margin LINES
# counting outward and the glyph occupies the band L..L+1, so treating L as
# the glyph centre puts every label about a line-height too far out -- which
# is what pushed the A/B titles off the page and the x-titles into the row
# below.
.tmp <- tempfile(fileext = ".pdf")
pdf(.tmp, width = fig_width, height = fig_height); plot.new()
h_tick_in  <- strheight("0123456789",          units = "inches", cex = tick_cex)
h_axlab_in <- strheight("Component loading",   units = "inches", cex = axlab_cex)
h_title_in <- strheight("Component retention", units = "inches", cex = title_cex)
dev.off(); unlink(.tmp)

# distances from the plot edge to each label's CENTRE
tick_extent_in <- ax_in_tick + h_tick_in          # outer edge of the tick labels
xtitle_ctr_in  <- tick_extent_in + gap_tick_title + h_axlab_in/2
ptitle_ctr_in  <- ax_in_ptitle   + h_title_in/2

# Each margin holds its label plus half the required onward gap: the space
# between an x-title and the panel title below it is the leftover inside this
# bottom margin plus the leftover above that title inside the next top
# margin, so splitting gap_title_next between them makes the ratio exact.
mai_bot <- xtitle_ctr_in + h_axlab_in/2 + gap_title_next/2
mai_top <- ptitle_ctr_in + h_title_in/2 + gap_title_next/2

# CHANGED: h_load now scales with vn_substantive (predictor count) instead
# of a hardcoded 2.15in -- that guess was sized for far fewer predictors and
# left C/D's loading bars too short to fit a label once vn_substantive grew
# to 30 (confirmed by test-rendering panel C at the old height: labels
# overlapped/were illegible). 1.8x the measured tick-label height per row
# leaves comfortable room for barplot's own inter-bar gaps (space=0.35) plus
# descenders -- confirmed by test-rendering panel C at this height with all
# 30 real (long) predictor names: every label legible, none clipped or
# overlapping.
h_load <- vn_substantive * (h_tick_in * 1.8) + mai_bot + mai_top
h_age <- 1.10
# h_scree chosen so panel B's plot region is exactly square with no padding
h_scree <- (fig_width - (lab_dem_in - d_left_in))/2 - b_left_in - r_in +
    mai_bot + mai_top
# CHANGED: h_dem is pinned to what it would have been under the ORIGINAL
# combined-figure sizing (fig_height=9, h_load=2.15), not derived from the
# new (bigger) h_load above -- A-D and E-G are now two SEPARATE figures (see
# below), and only A-D was asked to grow, so E-G's proportions are
# deliberately held at their previous values rather than being pushed around
# by A-D's own change.
h_dem <- 9.0 - h_scree - 2.15 - h_age
stopifnot(h_dem > 1, mai_bot > 0, mai_top > 0)

cat(sprintf("measured heights: tick %.4f axlab %.4f title %.4f in\n",
            h_tick_in, h_axlab_in, h_title_in))
cat(sprintf("gap tick-labels -> x-title %.4f in (%.2f pt)\n",
            gap_tick_title, gap_tick_title*72))
cat(sprintf("gap x-title -> title below %.4f in (%.2f pt)  ratio %.3f\n",
            gap_title_next, gap_title_next*72, gap_title_next/gap_tick_title))
cat(sprintf("mai bottom %.4f top %.4f; A-D rows %.3f %.3f; E-G rows %.3f %.3f\n",
            mai_bot, mai_top, h_scree, h_load, h_dem, h_age))
if (nslot > 30)
    warning("nslot = ", nslot, ": demographic rows will be very thin; ",
            "consider raising min_n")

ylim_s <- range(S, na.rm = TRUE) + c(-0.20, 0.20) * diff(range(S, na.rm = TRUE))
pad    <- 0.03 * diff(ylim_s)

# CHANGED: split into two figures/PDFs -- panels A-D and E-G used to share
# one 9in page; with 30 substantive predictors now, A-D needed to grow tall
# enough (see h_load above) that keeping E-G on the same page no longer
# made sense. Each figure below validates and lays out only its own panels;
# figure 2 (E-G) is opened later, right before its own drawing code.

# ---- Figure 1: A-D ----
fig_height_AD <- h_scree + h_load
pdf("pca_scree_loadings.pdf", height = fig_height_AD, width = fig_width)
# Check every panel's inch margins fit its own figure region: a panel whose
# margins exceed it fails at plot.new() with "figure margins too large".
reg <- rbind(A = c(col_w[1]+col_w[2], h_scree, a_left_in),
             B = c(col_w[3],          h_scree, b_left_in),
             C = c(col_w[1],          h_load,  lab_load_in),
             D = c(col_w[2]+col_w[3], h_load,  d_left_in))
pw <- reg[,1] - reg[,3] - r_in
ph <- reg[,2] - mai_bot - mai_top
# B is squared by padding, so its plot region is the smaller of the two
pw["B"] <- ph["B"] <- min(pw["B"], ph["B"])
if (any(pw < 0.05) || any(ph < 0.05))
    stop("margins too large for panel(s): ",
         paste(rownames(reg)[pw < 0.05 | ph < 0.05], collapse = ", "))

layout(matrix(c(1,1,2,
                3,4,4), 2, 3, byrow = TRUE),
       widths = col_w, heights = c(h_scree, h_load))

# Margins are set with par("mai"), which is in INCHES. par("mar") is in
# units of csi * cex, and whether layout() applies the multi-panel cex
# reduction the way mfrow does is exactly the assumption that produced the
# unequal column widths and then the margins-too-large error. mai removes
# the font units from the geometry entirely, so the widths below are exact
# by construction.
#
# inl() converts inches to LINES and is still needed for mgp and mtext(line=),
# which have no inch equivalent. It reads par() at call time -- inside a
# panel, after par(mai=) -- rather than from a value captured earlier.
inl <- function(inches) inches / (par("csi") * par("cex"))

# Axis and panel titles placed by explicit coordinate, so that "centred"
# really is the glyph centre. d is the distance in inches from the relevant
# plot edge to the centre of the text.
xtitle <- function(lab, d = xtitle_ctr_in) {
    y <- grconvertY(grconvertY(0, "npc", "inches") - d, "inches", "user")
    text(mean(par("usr")[1:2]), y, lab, cex = axlab_cex, xpd = NA)
}
ytitle <- function(lab, d) {
    x <- grconvertX(grconvertX(0, "npc", "inches") - d, "inches", "user")
    text(x, mean(par("usr")[3:4]), lab, cex = axlab_cex, srt = 90, xpd = NA)
}
ptitle <- function(lab, d = ptitle_ctr_in) {
    y <- grconvertY(grconvertY(1, "npc", "inches") + d, "inches", "user")
    text(mean(par("usr")[1:2]), y, lab, cex = title_cex, font = 2, xpd = NA)
}

mai_of <- function(left) c(mai_bot, left, mai_top, r_in)

# Panel B holds a density over two axes with identical ranges, so it should be
# square on the page as well as nominally. Its figure region is wider than it
# is tall, so the surplus width goes into the right margin and any surplus
# height is split between top and bottom. Computed rather than hand-tuned, so
# it stays square if the row heights or column widths change.
mai_square <- function(region_w, region_h, left) {
    avail_w <- region_w - left     - r_in
    avail_h <- region_h - mai_bot  - mai_top
    side    <- min(avail_w, avail_h)
    c(mai_bot + (avail_h - side)/2, left,
      mai_top + (avail_h - side)/2, r_in + (avail_w - side))
}

cat(sprintf("columns: %.3f %.3f %.3f in\n", col_w[1], col_w[2], col_w[3]))
cat(sprintf("target plot widths: C/D %.3f, E/F %.3f in\n",
            W_C - lab_load_in - r_in, W_E - lab_dem_in - r_in))

# ---- A: scree ----
# CHANGED: vn -> vn_substantive throughout this panel -- ev/sim_mean come
# from zcoef_clean (visual-excluded), which has vn_substantive columns, not
# vn's full count.
par(mai = mai_of(a_left_in), mgp = c(3, inl(ax_in_tick), 0))
# CHANGED: parallel-analysis (sim_mean) and Kaiser (eigenvalue = 1) reference
# lines removed from the plot -- both now favor 4 components while VSS (and
# the retained solution) favors 2, so showing them was confusing rather than
# informative. fp/sim_mean are still computed above (kept as a diagnostic
# available in the console/log), just no longer drawn. ylim no longer needs
# to accommodate sim_mean.
plot(seq_len(vn_substantive), ev, type = "n", xlab = "", ylab = "",
     ylim = range(0, ev)*c(1,1.05), xlim = c(0.6, vn_substantive+0.4), axes = FALSE)
# gap.axis = -1 disables axis()'s overlap check so all 13 numbers draw at the
# same size as every other tick label
axis(1, at = seq_len(vn_substantive), cex.axis = tick_cex, gap.axis = -1)
xtitle("Component")
axis(2, cex.axis = tick_cex, las = 1)
ytitle("Eigenvalue", ax_in_ytitle_a)
box(bty = "l")
lines(seq_len(vn_substantive), ev, col = "black", lwd = 1.8)
points(seq_len(vn_substantive), ev, pch = 21, bg = ifelse(seq_len(vn_substantive) <= 2, "black", "white"),
       col = "black", cex = 0.9)
ptitle("PCA scree plot")
panel_letter("A")

# ---- B: score density ----
par(mai = mai_square(col_w[3], h_scree, b_left_in),
    mgp = c(3, inl(ax_in_tick), 0))
lim <- max(abs(S), na.rm = TRUE) * 1.15
kd  <- kde2d(S[,1], S[,2], n = 120, lims = c(-lim, lim, -lim, lim))
image(kd, col = dens_ramp(128), xlab = "", ylab = "", axes = FALSE, useRaster = TRUE)
axis(1, cex.axis = tick_cex); axis(2, cex.axis = tick_cex, las = 1); box(bty = "l")
xtitle(sprintf("PC1 (%.0f%%)", 100*ve[1]))
ytitle(sprintf("PC2 (%.0f%%)", 100*ve[2]), ax_in_ytitle_b)
contour(kd, add = TRUE, drawlabels = FALSE, col = rgb(0,0,0,0.28), nlevels = 6)
abline(h = 0, v = 0, lty = 3, col = "grey55")
points(S[,1], S[,2], pch = 20, cex = 0.055, col = rgb(0,0,0,0.13))
ptitle("Score density")
cat(sprintf("panel B plot region: %.3f x %.3f in (square: %s)\n",
            par("pin")[1], par("pin")[2],
            abs(diff(par("pin"))) < 1e-6))
panel_letter("B")

# ---- C, D: loadings (own column split, narrower left margin) ----
ordL <- order(L[, 1])
for (k in 1:2) {
    par(mai = mai_of(if (k == 1) lab_load_in else d_left_in),
        mgp = c(3, inl(ax_in_tick), 0))
    bp <- barplot(L[ordL, k], horiz = TRUE, col = catcol[catidx_substantive[ordL]], border = NA,   # CHANGED: catidx -> catidx_substantive
                  xlim = c(-1, 1), names.arg = rep("", vn_substantive), axes = FALSE, space = 0.35)   # CHANGED: vn -> vn_substantive
    axis(1, at = seq(-1, 1, 0.5), cex.axis = tick_cex)
    xtitle("Component loading")
    if (k == 1) axis(2, at = bp, labels = ponames_substantive[ordL], las = 1, tick = FALSE,   # CHANGED: ponames -> ponames_substantive
                     cex.axis = tick_cex, line = -0.4)
    abline(v = 0, col = "black", lwd = 0.8)
    if (k == 1) legend("bottomright", ncol = 1, cex = tick_cex, inset = c(0.02, 0.02),
                       legend = catnames, fill = catcol, border = NA,
                       bty = "o", box.col = NA, bg = rgb(1,1,1,0.85))
    ptitle(sprintf("PC%d loadings", k))
    cat(sprintf("panel %s plot width: %.3f in\n", c("C","D")[k], par("pin")[1]))
    panel_letter(c("C","D")[k])
}

dev.off()
cat(sprintf("\nwrote pca_scree_loadings.pdf: %.2f x %.2f in\n",
            fig_width, fig_height_AD))

# ---- Figure 2: E-G ----
fig_height_EG <- h_dem + h_age
pdf("pca_demographics.pdf", height = fig_height_EG, width = fig_width)
reg_EG <- rbind(E = c(col_w[1]+col_w[2], h_dem, lab_dem_in),
                F = c(col_w[3],          h_dem, d_left_in),
                G = c(sum(col_w),        h_age, g_left_in))
pw <- reg_EG[,1] - reg_EG[,3] - r_in
ph <- reg_EG[,2] - mai_bot - mai_top
if (any(pw < 0.05) || any(ph < 0.05))
    stop("margins too large for panel(s): ",
         paste(rownames(reg_EG)[pw < 0.05 | ph < 0.05], collapse = ", "))

layout(matrix(c(1,1,2,
                3,3,3), 2, 3, byrow = TRUE),
       widths = col_w, heights = c(h_dem, h_age))

# ---- E, F: all demographics, blocked, one component per panel ----
for (k in 1:2) {
    par(mai = mai_of(if (k == 1) lab_dem_in else d_left_in),
        mgp = c(3, inl(ax_in_tick), 0))
    box_at <- ypos[!slot_gap]
    dk <- lapply(which(!slot_gap), function(si) {
        g <- droplevels(factor(udemc[[slot_var[si]]]))
        x <- S[!is.na(g) & g == slot_lab[si], k]
        x[!is.na(x)]
    })
    # notch = 95% CI of the median; R warns when a notch passes the hinges,
    # which flags a cell too small for visual comparison
    boxplot(dk, at = box_at, horizontal = TRUE, boxwex = 0.66, notch = TRUE,
            col = col_pc[k], border = "grey35", outline = FALSE,
            ylim = ylim_s, xlim = c(0.5, nslot + 0.5),
            axes = FALSE, xlab = "", ylab = "")
    abline(v = 0, lty = 3, col = "grey55")
    for (si in which(slot_gap)) abline(h = ypos[si], col = "grey80", lwd = 0.6)
    axis(1, cex.axis = tick_cex)
    xtitle("Component score")
    if (k == 1) {
        # named lvl_n, not ns: `ns` is the participant count from the upstream
        # analysis loop and must not be shadowed here
        lvl_n <- sapply(which(!slot_gap), function(si) {
            g <- droplevels(factor(udemc[[slot_var[si]]]))
            sum(!is.na(g) & g == slot_lab[si]) })
        axis(2, at = box_at, labels = sprintf("%s (%d)", disp(slot_lab[!slot_gap]), lvl_n),
             las = 1, tick = FALSE, cex.axis = tick_cex, line = -0.4)
        # Variable names horizontal, on an outer tier of the same margin,
        # positioned in inch coordinates from the device's own par("csi").
        left_in <- grconvertX(0, from = "npc", to = "inches")
        # clamp: if the tier would fall outside the device it is silently
        # clipped and the labels disappear, which is what happened before
        x_in    <- max(left_in - grp_in, 0.02)
        x_grp   <- grconvertX(x_in, from = "inches", to = "user")
        if (left_in - grp_in < 0.02)
            warning("variable-name tier clamped: grp_in (", grp_in,
                    "in) exceeds the left margin (", round(left_in, 2), "in)")
        for (i in seq_along(demvars)) {
            at_i <- mean(ypos[slot_var == demvars[i]])
            text(x_grp, at_i, demlabs[i], adj = c(1, 0.5), font = 2,
                 cex = group_cex, xpd = NA)
        }
    }
    # asterisks: side of the box indicates direction
    for (si in which(!slot_gap)) {
        st  <- lev_stats[[paste(slot_var[si], k)]]
        sr  <- st[match(slot_lab[si], st$level), ]
        lab <- star(sr$p_holm)
        if (!nzchar(lab)) next
        x <- dk[[which(which(!slot_gap) == si)]]
        if (!length(x)) next
        bsx <- boxplot.stats(x)$stats
        yy  <- star_y(ypos[si])
        if (!is.na(sr$est) && sr$est > 0)
            text(min(bsx[5] + pad, ylim_s[2]), yy, lab, cex = star_cex, adj = c(0, 0.5))
        else
            text(max(bsx[1] - pad, ylim_s[1]), yy, lab, cex = star_cex, adj = c(1, 0.5))
    }
    box(bty = "l")
    ptitle(sprintf("PC%d scores", k))
    cat(sprintf("panel %s plot width: %.3f in\n", c("E","F")[k], par("pin")[1]))
    panel_letter(c("E","F")[k])
}

# ---- G: age ----
# G needs only "PC1 (nnnn)" on the left, so it is NOT aligned to E: a small
# margin lets it use the full figure width.
par(mai = mai_of(g_left_in), mgp = c(3, inl(ax_in_tick), 0))
tl  <- max(abs(c(age_res$tau, age_res$lo, age_res$hi)), na.rm = TRUE) * 1.5
bpa <- barplot(rev(age_res$tau), horiz = TRUE, col = rev(col_pc), border = NA,
               xlim = c(-tl, tl), axes = FALSE, names.arg = rep("", 2),
               space = 0.45, xlab = "", ylab = "")
arrows(rev(age_res$lo), bpa, rev(age_res$hi), bpa, angle = 90, code = 3,
       length = 0.03, col = "grey30", lwd = 1)
abline(v = 0, col = "black", lwd = 0.8)
axis(1, cex.axis = tick_cex)
axis(2, at = bpa, labels = rev(sprintf("PC%d (%d)", age_res$pc, age_res$n)),
     las = 1, tick = FALSE, cex.axis = tick_cex, line = -0.4)
xtitle(expression("Kendall's " * tau))
for (j in 1:2) {
    r   <- age_res[rev(seq_len(2))[j], ]
    lab <- star(r$p_holm); if (!nzchar(lab)) next
    yy  <- star_y(bpa[j])
    if (r$tau > 0) text(r$hi + 0.06*tl, yy, lab, cex = star_cex, adj = c(0, 0.5))
    else           text(r$lo - 0.06*tl, yy, lab, cex = star_cex, adj = c(1, 0.5))
}
box(bty = "l")
ptitle("Age")
panel_letter("G")

dev.off()
cat(sprintf("\nwrote pca_demographics.pdf: %.2f x %.2f in\n",
            fig_width, fig_height_EG))
cat(sprintf("demographic rows: %d slots in %.2f in -> %.1f pt per row\n",
            nslot, h_dem, (h_dem - mai_bot - mai_top)/nslot*72))







































# total variance explained by individual differences ----------------------

p <- principal(zcoef_clean,2)   # CHANGED: was zcoef -- see predictor_registry above

fit <- lm(p$scores[,1]~age+gender+sex+race+ethnicity+education+nationality,data=udemc)
summary(fit)
# Adjusted R-squared:  0.007839  

fit <- lm(p$scores[,2]~age+gender+sex+race+ethnicity+education+nationality,data=udemc)
summary(fit)
# Adjusted R-squared:  0.04067




# data visualization ------------------------------------------------------


# CHANGED: ponames/vn/catidx/catnames now come from predictor_registry near
# the top of the script (context predictors excluded there -- see that
# block's comments) instead of being redefined here, a second time.
stopifnot(vn == length(varind))

# CHANGED: reverted to 4 tab20 colors, same reasoning as catcol in FIGURE
# SETUP above -- context predictors (which briefly justified a 6-way split)
# are excluded from this script entirely now, so Network/Self are back to
# a single (change) shade each, no separate "context" sibling to tell them
# apart from here.
catcol_dark  <- c("#2CA02C", "#FF7F0E", "#9467BD", "#1F77B4")
catcol_light <- c("#98DF8A", "#FFBB78", "#C5B0D5", "#AEC7E8")
# CHANGED: dropped the intermediate "mid" shade that used to distinguish
# the marginal-over-visual column -- per request, keep the color scheme as
# close to the Python notebook's as possible, where every conditioned
# panel (unique, marginal-over-visual, Shapley) uses the SAME light shade
# and only the zero-order panel is dark. The column headers (added below)
# now carry the "which column is which" job the mid shade used to.
col_light <- catcol_light[catidx]   # marginal over visual AND partial
col_dark  <- catcol_dark[catidx]    # zero-order

colnames(zcoef) <- ponames
colnames(ucoef) <- ponames
colnames(mcoef) <- ponames

mzc <- apply(zcoef,2,mean)
muc <- apply(ucoef,2,mean)
mmc <- apply(mcoef,2,mean,na.rm=TRUE)   # NaN at the 2 visual-predictor positions; never plotted

zpsym <- rep(NA,vn); zpsym[zpc<.05]<-"*"; zpsym[zpc<.01]<-"**"; zpsym[zpc<.001]<-"***"
upsym <- rep(NA,vn); upsym[upc<.05]<-"*"; upsym[upc<.01]<-"**"; upsym[upc<.001]<-"***"
mpsym <- rep(NA,vn); mpsym[mpc<.05]<-"*"; mpsym[mpc<.01]<-"**"; mpsym[mpc<.001]<-"***"

ord <- order(mzc, decreasing = FALSE)

# CHANGED: combined_regression.pdf restructured from one stacked column
# (zero-order hist+box directly above partial hist+box, per predictor) into
# three side-by-side panels -- A: zero-order, B: marginal over visual
# (ADDED), C: partial -- sharing one predictor-name column, with the same
# predictor row order (by `ord`, zero-order mean) in all three panels. With
# vn predictors the old single-column version stacked hist+box+hist+box
# rows on one page and had become too tall/dense to read; putting the
# panels side by side instead of stacked (each predictor now costs 2
# stacked rows shared by all panels, not 4+a label row) roughly halves the
# number of stacked rows needed per added panel. This mirrors the panel
# design already used for the same comparison in the Python notebook's
# plot_group_results (zero-order / marginal-over-visual / unique).
fig_height_in <- 8.5     # <- the knob; adjust if rows feel cramped/sparse
fig_width_in  <- 9.5     # CHANGED: was 6.5 -- wider, per request (more room
                          # per plot column now that there are 3 of them)

g <- 0.3

# CHANGED: every text element on this page unified to one size, matching
# the legend's category labels ("Target"/"Dyad"/"Network"/"Self") -- was
# previously a mix of 1.2 (predictor labels, axes), 1.4 (panel letters, via
# the shared panel_cex), and 0.85 (legend headers).
txt_cex <- 0.95

csi         <- 12/72                # inches per margin line at pdf() pointsize 12
axis_mar    <- 3.5                  # bottom margin of the axis row, in lines
axis_mai_in <- axis_mar*csi         # pure bottom margin, in inches
axis_in     <- axis_mai_in + 1.5/72 # axis row TOTAL height: margin + ~1.5pt of plot region
top_pad_in  <- 0.02   # tiny gap between a predictor's stacked hist and box panels
# CHANGED: now reserved on EVERY column boundary, not just the outer page
# edges -- split in half so the total gap between any two adjacent
# columns' tick labels equals side_pad_in, the same amount already
# confirmed sufficient at the page edges. The middle column used to get
# ml=mr=0, leaving its tick labels ("-0.4"/"0.4") with no breathing room
# against columns 2/4's tick labels and visibly overlapping them
# (confirmed by test-rendering, then fixed here).
side_pad_in  <- 0.09
inner_pad_in <- side_pad_in / 2
col2_ml <- side_pad_in;  col2_mr <- inner_pad_in   # zero-order
col3_ml <- inner_pad_in; col3_mr <- inner_pad_in   # marginal over visual
col4_ml <- inner_pad_in; col4_mr <- side_pad_in    # partial

# ---- legend geometry (2 rows: zero-order/partial-shade x 4 columns: category) ----
# CHANGED: legend transposed from 4 rows x 2 (stat-type) columns to 2 rows
# (stat-type) x 4 (category) columns, and moved out of its own dedicated
# bottom row into the SAME row as the two x-axes, in the middle column --
# directly between them, reusing space that would otherwise sit blank
# rather than adding a whole separate row (which left the old legend row's
# left/right thirds wasted, and put daylight between it and everything
# above it). CHANGED again: moved (along with the predictor-name column it
# shares a column with) from the middle to the far left -- see the
# column-width/lmat comments below. CHANGED again: back to 2 rows (was
# briefly 3) -- marginal-over-visual no longer has its own shade (same
# light color as partial now, per request), so a 3rd legend row would just
# show two identically-colored rows; the new column headers below now
# carry the "which column is which" job that 3rd row used to.
leg_sq        <- 0.08   # CHANGED: was 0.11 -- shrunk (see draw_legend for the other fixes)
leg_rgap      <- 0.03   # vertical gap between the 2 stat-type rows -- unchanged, not asked to shrink
leg_cgap      <- 0.05   # ADDED: explicit small lateral gap between category swatches,
                        # replacing the old scheme of spreading them across the full
                        # available width
leg_hdr_h     <- 0.16   # vertical room for the category-name header above the swatches
leg_row_lab_w <- max(strwidth(c("Zero-order","Partial"), units = "inches", cex = txt_cex)) + 0.06
leg_need_h    <- leg_hdr_h + 2*leg_sq + leg_rgap
axis_legend_in <- max(axis_in, leg_need_h)   # the shared bottom row's height

# ADDED: header row, holding both the descriptive column header
# ("Zero-order"/"Marginal over visual"/"Partial") and the bare panel
# letter (A/B/C) at the same height -- mirrors the Python notebook's
# header+letter convention. Replaces the old in-panel draw_panel_letter(),
# which just placed a bare letter partway down each panel with no
# column-identifying text at all.
.tmp_hdr <- tempfile(fileext = ".pdf")
pdf(.tmp_hdr, width = fig_width_in, height = fig_height_in); plot.new()
h_header_in <- strheight("Marginal over visual", units = "inches", cex = txt_cex) + 0.10
dev.off(); unlink(.tmp_hdr)

# CHANGED: each predictor now costs only 2 height units (one hist row, one
# box row) shared across both panels, instead of 4+g (two stacked hist+box
# pairs plus a label row) -- the predictor-name label moved into its own
# middle COLUMN (spanning a predictor's hist+box rows) instead of its own
# row, since it no longer needs to separate a zero-order block from a
# partial block that are now side by side rather than stacked.
bg <- 2*g   # between-predictor gap weight (proportional analog of the old
            # 4*g, scaled down to match the smaller per-predictor block size)

W   <- vn*2 + (vn-1)*bg               # predictor-block weight (page-height independent)
hpu <- (fig_height_in - h_header_in - axis_legend_in)/W    # inches per weight unit
stopifnot(hpu > 0)
header_weight <- h_header_in/hpu
axis_weight   <- axis_legend_in/hpu

heights <- c(header_weight)   # ADDED: header row first, matches row_of below starting at index 2
row_of  <- integer(vn)    # layout-row index of each predictor's (hist) row
for (i in seq_len(vn)) {
    heights <- c(heights, 1, 1)
    row_of[i] <- length(heights) - 1
    if (i < vn) heights <- c(heights, bg) else heights <- c(heights, axis_weight)
}
stopifnot(abs(sum(heights)*hpu - fig_height_in) < 1e-6)   # weights really fill the page
total_rows <- length(heights)

# Column widths: a predictor-name label column sized to the longest
# predictor name, followed by three equal-width plot columns (zero-order,
# marginal-over-visual, partial) -- measured the same way as lab_load_in in
# the PCA figure (a throwaway pdf() device + strwidth()), rather than
# guessed.
#
# CHANGED: label column moved from the MIDDLE (between the 2 plot columns)
# to the FAR LEFT, and given its full measured width (lab_mid_in) rather
# than the half-width + redistribute-the-rest-to-the-plot-columns trick
# the 2-panel version used. That trick was safe only because the label sat
# between two on-page panels, so text wider than its own (shrunk) column
# overflowed harmlessly (via xpd=NA) into a neighboring panel that was
# still on the page; at the far left edge, the same overflow runs off the
# printable page entirely -- confirmed by test-rendering, which clipped
# the first 1-2 characters of every predictor name until this was reverted
# to the full measured width here.
.tmp <- tempfile(fileext = ".pdf")
pdf(.tmp, width = fig_width_in, height = fig_height_in); plot.new()
lab_mid_in <- max(strwidth(ponames, units = "inches", cex = txt_cex)) + 0.08
dev.off(); unlink(.tmp)
col_w3 <- rep((fig_width_in - lab_mid_in)/3, 3)   # zero-order, marginal-over-visual, partial

# lmat: 4 columns (predictor label | zero-order plot | marginal-over-visual
# plot | partial plot). Cells left at 0 are blank (unused by any panel) --
# layout()'s convention for "no figure here". Panel numbers are assigned in
# the exact order the drawing loop below issues its plotting calls, since
# layout() advances to the next-numbered region on each new high-level plot
# call, not by matrix position.
lmat <- matrix(0L, nrow = total_rows, ncol = 4)
npanel <- 0L
npanel <- npanel + 1L; lmat[1, 2] <- npanel   # header, col 2 (A: zero-order)
npanel <- npanel + 1L; lmat[1, 3] <- npanel   # header, col 3 (B: marginal over visual)
npanel <- npanel + 1L; lmat[1, 4] <- npanel   # header, col 4 (C: partial)
for (i in seq_len(vn)) {
    r <- row_of[i]
    npanel <- npanel + 1L; lmat[r, 1] <- npanel; lmat[r+1, 1] <- npanel   # label, spans both rows
    npanel <- npanel + 1L; lmat[r,   2] <- npanel   # zero-order hist
    npanel <- npanel + 1L; lmat[r+1, 2] <- npanel   # zero-order box
    npanel <- npanel + 1L; lmat[r,   3] <- npanel   # marginal-over-visual hist (or blank)
    npanel <- npanel + 1L; lmat[r+1, 3] <- npanel   # marginal-over-visual box (or blank)
    npanel <- npanel + 1L; lmat[r,   4] <- npanel   # partial hist
    npanel <- npanel + 1L; lmat[r+1, 4] <- npanel   # partial box
    npanel <- npanel + 1L; lmat[r+2, 2] <- npanel   # between-gap / axis, col 2 (zero-order)
    npanel <- npanel + 1L; lmat[r+2, 3] <- npanel   # between-gap / axis, col 3 (marginal)
    npanel <- npanel + 1L; lmat[r+2, 4] <- npanel   # between-gap / axis, col 4 (partial)
    if (i == vn) { npanel <- npanel + 1L; lmat[r+2, 1] <- npanel }   # legend, col 1, same row as the axes
}

pdf("combined_regression.pdf", height = fig_height_in, width = fig_width_in)
layout(lmat, widths = c(lab_mid_in, col_w3[1], col_w3[2], col_w3[3]), heights = heights)
xlim <- c(-.4, .4)   # CHANGED: was c(-1,1) -- this is now only the DISPLAYED window
star_off <- 0.05 * diff(xlim)   # significance-star offset, scaled to the new (narrower) range

# CHANGED: hist() bins must span the full data range or it errors ("some
# 'x' not counted; maybe 'breaks' do not span range of 'x'") -- coefficients
# routinely exceed +-0.5 even though we now only want to DISPLAY the
# +-0.5 window (xlim). bk_lim spans the actual data (zcoef, mcoef, and
# ucoef, symmetric about 0); n_breaks is scaled up from the original
# 100-over-range-2 bin width so bars keep the same visual density inside
# the (narrower) displayed window regardless of how wide the data actually
# go.
bk_lim   <- c(-1, 1) * max(abs(c(zcoef, ucoef, mcoef)), na.rm = TRUE)
n_breaks <- max(20, round(100 * diff(bk_lim) / diff(xlim)))

# ADDED: dedicated header row (A/B/C letter + descriptive column header,
# at the same height) -- replaces the old in-panel draw_panel_letter(),
# which placed a bare letter partway down each panel with no text
# identifying which column was which. Mirrors the Python notebook's
# header+letter convention exactly.
draw_header <- function(letter, label, ml, mr) {
    par(mai = c(0, ml, 0.02, mr))
    plot(0, 0, col = "white", xlim = c(0,1), ylim = c(-1,1), axes = F, xlab = "", ylab = "", main = "")
    text(0.5, 0, label,  cex = txt_cex, font = 2, adj = c(0.5,.5), xpd = NA)
    text(0,   0, letter, cex = txt_cex, font = 2, adj = c(0,.5),   xpd = NA)
}

# CHANGED: hist/box drawing pulled into helpers -- the zero-order and
# partial blocks were near-identical code differing only in data/colors,
# and duplicating that logic (as the old stacked version did) is exactly
# the class of bug this file's predictor_registry refactor was meant to
# stop recurring elsewhere. Each takes ml/mr (left/right margin, inches) so
# the caller can reserve side_pad_in on only the page-edge side of a panel.
draw_hist <- function(x, col_fill, m, psym, ml, mr) {
    par(mai = c(0, ml, top_pad_in, mr), lwd = 0.1)
    hh <- hist(x, 100, axes = F, xlab = "", ylab = "", main = "", border = "white",
               col = col_fill, breaks = seq(bk_lim[1], bk_lim[2], length.out = n_breaks),
               freq = F, xlim = xlim)
    if (m > 0) text(quantile(x, .95, na.rm = TRUE) + star_off, max(hh$density)/2, psym, cex = txt_cex)
    else       text(quantile(x, .05, na.rm = TRUE) - star_off, max(hh$density)/2, psym, cex = txt_cex)
}
draw_box <- function(x, col_fill, border_col, ml, mr) {
    par(mai = c(0, ml, top_pad_in, mr), lwd = 1)
    abline(v = 0, lty = 1, col = rgb(0,0,0,.5))
    boxplot(x, horizontal = T, axes = F, xlab = "", ylab = "", main = "", ylim = xlim,
            notch = T, pch = 20, cex = .1, col = col_fill, border = border_col, lwd=.2)
    abline(v = 0, lty = 1, col = rgb(0,0,0,.5))
}
draw_gap <- function(ml, mr) {
    par(mai = c(0, ml, 0, mr))
    plot(0, 0, col = "white", xlim = xlim, axes = F, xlab = "", ylab = "", main = "")
    abline(v = 0, lty = 1, col = rgb(0,0,0,.5))
}
draw_axis_row <- function(ml, mr) {
    par(mai = c(axis_mai_in, ml, 0, mr), mgp = c(1.8, 0.35, 0), tcl = -0.25)
    plot(0, 0, col = "white", xlim = xlim, axes = F, xlab = "Regression coefficient", cex.lab = txt_cex)
    abline(v = 0, lty = 1, col = rgb(0,0,0,.5))
    axis(1, cex.axis = txt_cex)
}
# ADDED: blank placeholder for a predictor x panel combination with no
# data -- currently only vdist/vecc in the marginal-over-visual column
# (there's no "controlling for visual" version of a visual predictor
# itself). Mirrors draw_gap's zero-reference-line continuation; called
# once with a label for the hist slot and once without for the box slot.
draw_blank <- function(ml, mr, label = NULL) {
    par(mai = c(0, ml, top_pad_in, mr))
    plot(0, 0, col = "white", xlim = xlim, axes = F, xlab = "", ylab = "", main = "")
    abline(v = 0, lty = 1, col = rgb(0,0,0,.5))
    if (!is.null(label))
        text(xlim[2]*0.25, 0, label, cex = txt_cex*0.8, font = 3, col = "grey45", adj = c(0,.5), xpd = NA)
}
# CHANGED: legend drawn in the SAME layout row as the three x-axes -- in
# column 1 (the predictor-label column), since the label moved to the far
# left and the legend moved with it (previously both lived in the middle
# column). Back to 2 (stat-type) rows x 4 (category) columns (was briefly
# 3) -- see the leg_need_h comment above.
#
# Three bugs fixed here, all surfacing together once "Network" (the
# longest category name) landed next to a long predictor name in the row
# above: (1) the row-label anchor was "x_lo/2+.1", which for a right-
# justified ("adj=c(1,...)") string wider than x_lo/2 -- "Zero-order" was
# -- put part of the text to the LEFT of that anchor past x=0, off the
# page; now anchored at x_lo itself (flush against where the swatches
# start), the only place a right-justified string can't overflow left of
# 0 as long as x_lo >= its own width (guaranteed, since x_lo is exactly
# sized to fit the longer of "Zero-order"/"Partial"). (2) the swatches
# were leg_sq tall but 2*leg_sq wide (a left-over sizing artifact, not an
# intentional rectangle) -- now true leg_sq x leg_sq squares. (3) the 4
# category swatches were spread across the ENTIRE remaining width
# (lab_mid_in - x_lo), which (combined with bug 2's doubled width) gave
# "Network" enough horizontal room that its 45-degree diagonal text
# reached up into the predictor-label row above; now packed tightly right
# after the row labels (swatch size + a small explicit leg_cgap) and
# rotated steeper (60 instead of 45 degrees, trading a little extra
# vertical reach for noticeably less horizontal reach) -- all three
# together free enough room that the overlap is gone (confirmed by
# test-rendering with "Network" directly above a long predictor name).
draw_legend <- function() {
    par(mai = c(0,0,0,0))
    plot(0, 0, col = "white", xlim = c(0, lab_mid_in), ylim = c(0, axis_legend_in), xlab = "", ylab = "",
         main = "", axes = FALSE, xaxs = "i", yaxs = "i")
    ncat  <- length(catnames)
    x_lo  <- leg_row_lab_w
    y_top <- axis_legend_in - (axis_legend_in - leg_need_h)/2   # vertically centered in the row
    y_hdr <- y_top - leg_hdr_h
    y_z1 <- y_hdr;           y_z2 <- y_z1 - leg_sq
    y_p1 <- y_z2 - leg_rgap; y_p2 <- y_p1 - leg_sq
    text(x_lo - 0.02, (y_z1+y_z2)/2, "Zero-order", cex = txt_cex, adj = c(1,0.5), xpd = NA)
    text(x_lo - 0.02, (y_p1+y_p2)/2, "Partial",    cex = txt_cex, adj = c(1,0.5), xpd = NA)
    for (cc in seq_len(ncat)) {
        xc <- x_lo + leg_cgap + (cc-1)*(leg_sq+leg_cgap) + leg_sq/2
        text(xc, y_z1+.03, catnames[cc], cex = txt_cex, adj = c(0,0), xpd = NA, srt=60)
        rect(xc-leg_sq/2, y_z2, xc+leg_sq/2, y_z1, col = catcol_dark[cc],  border = NA, xpd = NA)
        rect(xc-leg_sq/2, y_p2, xc+leg_sq/2, y_p1, col = catcol_light[cc], border = NA, xpd = NA)
    }
}

# header row -- drawn first, matching the panel numbers assigned to row 1
# of lmat above (columns 2, 3, 4 in that order)
draw_header("A", "Zero-order",           col2_ml, col2_mr)
draw_header("B", "Marginal over visual", col3_ml, col3_mr)
draw_header("C", "Partial",              col4_ml, col4_mr)

for (i in seq_len(vn)) {
    p <- ord[i]

    # ---- predictor label (panel column 1, spans both rows below) ----
    # CHANGED: right-justified (was centered) and using as little
    # horizontal space as the measured lab_mid_in allows -- flush against
    # column 2's left edge, matching the Python notebook's label column.
    par(mai = c(0,0,0,0))
    plot(0, 0, col = "white", xlim = c(0,1), ylim = c(-1,1), axes = F, xlab = "", ylab = "", main = "")
    text(1, 0, ponames[p], adj = c(1,.5), cex = txt_cex, xpd = NA)

    # ---- zero-order hist + box (panel column 2) ----
    draw_hist(zcoef[,p], col_dark[p], mzc[p], zpsym[p], col2_ml, col2_mr)
    draw_box (zcoef[,p], col_dark[p], "grey40", col2_ml, col2_mr)

    # ---- marginal-over-visual hist + box (panel column 3) ----
    # ADDED: blank placeholder for vdist/vecc (no "controlling for visual"
    # version of a visual predictor itself) -- mirrors the Python
    # notebook's blank "(baseline)" row for 'visual' in its own
    # marginal-over-visual panel. CHANGED: col_mid -> col_light and
    # "grey58" -> "grey75" -- no more intermediate shade (see the colors
    # comment above).
    if (visualidx[p]) {
        draw_blank(col3_ml, col3_mr, label = "(baseline)")
        draw_blank(col3_ml, col3_mr)
    } else {
        draw_hist(mcoef[,p], col_light[p], mmc[p], mpsym[p], col3_ml, col3_mr)
        draw_box (mcoef[,p], col_light[p], "grey75", col3_ml, col3_mr)
    }

    # ---- partial hist + box (panel column 4) ----
    draw_hist(ucoef[,p], col_light[p], muc[p], upsym[p], col4_ml, col4_mr)
    draw_box (ucoef[,p], col_light[p], "grey75", col4_ml, col4_mr)

    # ---- between-predictor gap, or bottom axis + legend on the last row ----
    if (i < vn) {
        draw_gap(col2_ml, col2_mr); draw_gap(col3_ml, col3_mr); draw_gap(col4_ml, col4_mr)
    } else {
        draw_axis_row(col2_ml, col2_mr); draw_axis_row(col3_ml, col3_mr); draw_axis_row(col4_ml, col4_mr)
        draw_legend()
    }
}

dev.off()

