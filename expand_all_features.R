# ===========================================================================
# expand_all_features.R
#
# NOT YET RUN. This consolidates expand_data.R + expand_extra_features.R +
# merge_features.R into the single script that should have produced
# alldata.csv all along, had all three rounds of feature work been planned
# together from the start instead of accreting across sessions. It is meant
# to produce byte-for-byte the same alldata.csv that running expand_data.R,
# then expand_extra_features.R, then merge_features.R in sequence produces
# -- NOT a rewrite or an optimization pass. Every computational line below
# is copied unchanged from one of those three scripts; only the boilerplate
# that all three duplicated (loading net_data.csv/graphs.csv, and the final
# merge/fill/reorder/write step that merge_features.R used to do against
# two files on disk) has been de-duplicated into one copy, run once, here.
#
# Specifically NOT changed, on purpose, even though it would be faster:
# expand_data.R's per-trial loop below re-parses each network's graph JSON
# from scratch every time it's visited (fromJSON(graphs[...]) inside the
# participant/trial loop), even though PHASE 1 below has already parsed
# and cached every network's topology once. Unifying that would be a real
# optimization, but it would also mean this script's two halves are no
# longer simply "the old scripts, verbatim" -- and this script is not being
# run or diffed against the old output before being trusted, so verbatim-
# but-slower was judged the safer choice over faster-but-unverified.
#
# The only other thing worth knowing before running this: PHASE 1 below
# calls set.seed(1) immediately before its own random-restart KK layout
# loop, exactly as expand_extra_features.R did standalone. Nothing above it
# in this script consumes any random draws (expand_data.R's own logic below
# is entirely deterministic), so that seed produces the identical vdist/
# vecc values either way -- but if this script is ever edited to add new
# randomness anywhere before PHASE 1, that guarantee breaks.
# ===========================================================================

if(!require(RJSONIO))     install.packages("RJSONIO");     require(RJSONIO)
if(!require(lmerTest))    install.packages("lmerTest");    require(lmerTest)
if(!require(igraph))      install.packages("igraph");      require(igraph)
if(!require(RColorBrewer)) install.packages("RColorBrewer"); require(RColorBrewer)

# ---- load data (shared by both expansion passes below) ---------------------
dat    <- read.csv("net_data.csv")
graphs <- read.csv("../../graphs.csv")
graphs <- as.character(graphs$graphs)

# ===========================================================================
# PART A -- expand_data.R, verbatim (produces `fdat`, the fulldata.csv-
# equivalent data frame, instead of writing it to disk).
# ===========================================================================

# descriptives ------------------------------------------------------------

hist(dat$age)
summary(dat$age)

table(dat$race)/30
table(dat$ethnicity)/30
table(dat$sex)/30
table(dat$gender)/30
table(dat$education)/30
table(dat$nationality)/30

hist(log10(dat$rt))
summary(dat$rt)

# load network graphs -----------------------------------------------------

expand.graph <- function(g){
  g <- fromJSON(g)
  nodes <- do.call(rbind,g$nodes)
  return(data.frame(nodes,density=g$density,assortativity=g$assortativity,apl=g$`average path length`,modularity=g$modularity))
}

glist <- list()
for (i in 1:length(graphs)){
  g <- expand.graph(graphs[i])
  glist[[i]] <- data.frame(nid=(i-1),nnode=dim(g)[1],g)
  print(i)
}

gdat <- do.call(rbind,glist)

# expand data for analysis ------------------------------------------------

# number of nodes for each trial
nnode <- tapply(gdat$nnode,gdat$nid,mean)
nnode <- nnode[dat$netid+1]
ntrial <- dim(dat)[1]
sel <- unlist(lapply(1:ntrial,function(x) rep(x,nnode[x])))

# expand behavioral data
edat <- dat[sel,]
edat <- data.frame(rid=1:length(sel),edat)

# expand network data
cgd <- do.call(rbind,glist[dat$netid+1])

# combine data
cdat <- cbind(edat,cgd)
rm(edat,cgd)

# dummy code responses
binres <- rep(0,length(sel))
shf <- c(0,cumsum(nnode)[-ntrial])+1
binres[shf+dat$response] <- 1

cdat <- cbind(cdat,binres)
cdat <- cdat[,-15]

# add columns for properties of unode
unodedat<- do.call(rbind,lapply(1:dim(cdat)[1],function(x) glist[[cdat$netid[x]+1]][cdat$unode[x]+1,4:8]))
colnames(unodedat) <- c("udegree","ucloseness","ubetweenness","ueigenvector","ucommunity")
cudat <- cbind(cdat,unodedat)


# remove unode and friends from dataframe rows (not actual options)
# measure dyadic relationship between unode and possible choices (distance, edge betweenness, & same community)
# measure overall network change for each possible choice (modularity, assortativity, & average path length)
# measure specific user change for each possible choice (centralities)
cudat <- cudat[cudat$id!=cudat$unode,]
frm <- rep(0,dim(cudat)[1])
udists <- rep(NA,dim(cudat)[1])
edgbet <- rep(NA,dim(cudat)[1])
chgmod <- rep(NA,dim(cudat)[1])
chgast <- rep(NA,dim(cudat)[1])
chapl <- rep(NA,dim(cudat)[1])
chub <- rep(NA,dim(cudat)[1])
chuc <- rep(NA,dim(cudat)[1])
chue <- rep(NA,dim(cudat)[1])
mutfr <- rep(NA,dim(cudat)[1])
clust <- rep(NA,dim(cudat)[1])
uclust <- rep(NA,dim(cudat)[1])
chuclust <- rep(NA,dim(cudat)[1])
constr <- rep(NA,dim(cudat)[1])
uconstr <- rep(NA,dim(cudat)[1])
chuconstr <- rep(NA,dim(cudat)[1])
usid <- unique(cdat$sid)
sn <- length(usid)
ntr <- 30
start <- Sys.time()
for (i in 1:sn){
  sel1 <- cudat$sid==usid[i]
  trials <- unique(cudat$trial[sel1])
  for (j in trials){
    sel2 <- cudat$trial[sel1]==j
    cunode <- cudat$unode[sel1][sel2][1]
    cg <- fromJSON(graphs[cudat$netid[sel1][sel2][1]+1])
    cg <- graph_from_edgelist(do.call(rbind,cg$links)+1,directed = F)
    friends <- which(as_adjacency_matrix(cg,type = "both")[cunode+1,]==1)-1
    frm[which(((sel1) & (cudat$trial==j) & (cudat$id %in% friends)))] <- 1
    udists[sel1 & cudat$trial==j] <- distances(cg,cudat$id[sel1][sel2]+1,cunode+1)
    mutfr[sel1 & cudat$trial==j] <- sapply(cudat$id[sel1][sel2]+1, function(v) length(intersect(friends+1, as.integer(neighbors(cg,v)))))
    clust[sel1 & cudat$trial==j] <- transitivity(cg,type="local",vids=cudat$id[sel1][sel2]+1)
    constr[sel1 & cudat$trial==j] <- constraint(cg,nodes=cudat$id[sel1][sel2]+1)
    uclust[sel1 & cudat$trial==j] <- transitivity(cg,type="local",vids=cunode+1)
    uconstr[sel1 & cudat$trial==j] <- constraint(cg,nodes=cunode+1)
    cnodes <- cudat$id[sel1][sel2]+1
    cno <- sum(sel1 & cudat$trial==j)
    ceb <- rep(NA,cno)
    ccm <- rep(NA,cno)
    cca <- rep(NA,cno)
    cpl <- rep(NA,cno)
    cud <- rep(NA,cno)
    cub <- rep(NA,cno)
    cuc <- rep(NA,cno)
    cue <- rep(NA,cno)
    ccl <- rep(NA,cno)
    ccc <- rep(NA,cno)
    for (k in 1:cno){
      if ((cnodes[k]-1) %in% friends){
          ceb[k] <- NA
          ccm[k] <- NA
          cca[k] <- NA
          cpl[k] <- NA
          cud[k] <- NA
          cub[k] <- NA
          cuc[k] <- NA
          cue[k] <- NA
          ccl[k] <- NA
          ccc[k] <- NA
      } else {
          ccg <- add_edges(cg,c(cunode+1,cnodes[k]))
          m0 <- ecount(cg)
          ceb[k] <- edge_betweenness(ccg, e = m0 + 1L)
          ccm[k] <- modularity(cluster_optimal(ccg))
          cca[k] <- assortativity_degree(ccg)
          cpl[k] <- mean_distance(ccg)
          cud[k] <- degree(ccg)[cunode+1]
          cub[k] <- betweenness(ccg)[cunode+1]
          cuc[k] <- closeness(ccg)[cunode+1]
          cue[k] <- eigen_centrality(ccg)$vector[cunode+1]
          ccl[k] <- transitivity(ccg,type="local",vids=cunode+1)
          ccc[k] <- constraint(ccg,nodes=cunode+1)
      }
    }
    edgbet[sel1 & cudat$trial==j] <- ceb
    chgmod[sel1 & cudat$trial==j] <- ccm - cudat$modularity[sel1 & cudat$trial==j]
    chgast[sel1 & cudat$trial==j] <- cca - cudat$assortativity[sel1 & cudat$trial==j]
    chapl[sel1 & cudat$trial==j] <- cpl - cudat$apl[sel1 & cudat$trial==j]
    chub[sel1 & cudat$trial==j] <- cub - cudat$ubetweenness[sel1 & cudat$trial==j]
    chuc[sel1 & cudat$trial==j] <- cuc - cudat$ucloseness[sel1 & cudat$trial==j]
    chue[sel1 & cudat$trial==j] <- cue - cudat$ueigenvector[sel1 & cudat$trial==j]
    chuclust[sel1 & cudat$trial==j] <- ccl - uclust[sel1 & cudat$trial==j]
    chuconstr[sel1 & cudat$trial==j] <- ccc - uconstr[sel1 & cudat$trial==j]
  }
  print(i)
  print(Sys.time() -start)
}
samecom <- as.numeric(cudat$community == cudat$ucommunity)
odat <- cbind(edgbet,chgmod,chgast,chapl,udists,chub,chuc,chue,
              mutfr,clust,uclust,chuclust,constr,uconstr,chuconstr)


fdat<- cbind(cudat,samecom,odat)
fdat<-fdat[frm!=1,]
fdat<-fdat[,-c(21,31)]
# CHANGED: no longer written to disk here (was write.csv(fdat,"fulldata.csv",...))
# -- `fdat` is consumed directly by PART C's merge below instead.
cat(sprintf("\nPart A done: fdat has %d rows, %d columns\n", nrow(fdat), ncol(fdat)))

# ===========================================================================
# PART B -- expand_extra_features.R, verbatim (produces `out`, the
# extra_features.csv-equivalent data frame, instead of writing it to disk).
# `dat`/`graphs` are already loaded above; not reloaded here.
# ===========================================================================

# ---- knobs -----------------------------------------------------------------
KK_RESTARTS <- 10     # visual-layout random restarts averaged per network
SR_GAMMA    <- 0.7    # Son, Bhandari & FeldmanHall (PNAS 2023): best-fitting
                       # value for community-like ("Random Walks") structure.
                       # Canonical SR formulation (row-normalized adjacency,
                       # M = (I - gamma*T)^-1) -- PNAS blocked direct access
                       # to their supplement, so double check this matches
                       # their exact transition-matrix definition if it
                       # matters for the paper.
KATZ_BETA_FRAC <- 0.5 # katz beta = KATZ_BETA_FRAC / leading eigenvalue(A),
                       # network-specific so it's always safely convergent
                       # (beta < 1/spectral radius) regardless of that
                       # network's own connectivity.

set.seed(1)

# ---- effective size (Burt 1992, binary/unweighted simplification) ---------
# effective_size(v) = |N(v)| - 2*|E(N(v))| / |N(v)|, i.e. degree minus a
# redundancy term counting ties AMONG v's neighbors (not to v itself).
# Degree 0/1 nodes have no possible redundancy, so effective size = degree.
effective_size <- function(A, deg, v) {
  sapply(v, function(vi) {
    nb <- which(A[vi, ] == 1)
    k  <- length(nb)
    if (k <= 1) return(k)
    m <- sum(A[nb, nb]) / 2
    k - (2 * m) / k
  })
}

# ===========================================================================
# PHASE 1: per-network precompute. One pass over all 10,000 networks; every
# quantity below depends only on network topology.
# ===========================================================================
n_networks <- length(graphs)
netfeat <- vector("list", n_networks)
nnode_vec <- integer(n_networks)

cat("Phase 1: per-network precompute (", n_networks, " networks)\n", sep = "")
start <- Sys.time()

for (i in seq_len(n_networks)) {
  gj <- fromJSON(graphs[i])
  n  <- length(gj$nodes)
  nnode_vec[i] <- n

  g <- graph_from_edgelist(do.call(rbind, gj$links) + 1, directed = FALSE)
  A <- as.matrix(as_adjacency_matrix(g))
  deg <- degree(g)

  # -- graph eccentricity, effective size, community size (per node) --------
  ecc_vec   <- eccentricity(g)
  esize_vec <- effective_size(A, deg, seq_len(n))
  comm_vec  <- sapply(gj$nodes, function(x) x[["community"]])   # fromJSON simplifies
                                                                 # all-numeric node objects
                                                                 # to atomic vectors, so $
                                                                 # fails but [[ works on both
  comsize_lookup <- table(comm_vec)
  comsize_vec <- as.numeric(comsize_lookup[as.character(comm_vec)])

  # -- pairwise: Adamic-Adar, Katz, structural similarity, degree homophily -
  w_aa <- ifelse(deg > 1, 1 / log(deg), 0)          # guards deg<=1 (Inf/undefined);
                                                      # common neighbors always have
                                                      # deg>=2 so this never actually
                                                      # zeroes out a real AA term
  AA <- A %*% diag(w_aa, nrow = n) %*% A
  diag(AA) <- NA

  eig1 <- max(Re(eigen(A, symmetric = TRUE, only.values = TRUE)$values))
  beta <- if (eig1 > 0) KATZ_BETA_FRAC / eig1 else 0
  Katz <- if (beta > 0) solve(diag(n) - beta * A) - diag(n) else matrix(0, n, n)

  normA <- sqrt(rowSums(A^2))
  denomA <- outer(normA, normA)
  cosmat <- (A %*% t(A)) / denomA
  cosmat[denomA == 0] <- 0    # isolated node(s): similarity undefined -> 0
  diag(cosmat) <- NA

  dhom <- outer(deg, deg, function(a, b) abs(a - b))

  # -- successor representation cosine similarity ----------------------------
  Tm <- A / pmax(deg, 1)      # row-normalized; isolated nodes keep an all-
                               # zero row rather than dividing by zero
  M  <- solve(diag(n) - SR_GAMMA * Tm)
  normM <- sqrt(rowSums(M^2))
  denomM <- outer(normM, normM)
  srcos <- (M %*% t(M)) / denomM
  srcos[denomM == 0] <- 0
  diag(srcos) <- NA

  # -- network-level centralization / degree variance (baseline) ------------
  netcent <- centr_degree(g)$centralization
  netvar  <- var(deg)

  # -- visual layout: Kamada-Kawai, averaged over KK_RESTARTS restarts -------
  distsum <- matrix(0, n, n)
  eccsum  <- numeric(n)
  for (r in seq_len(KK_RESTARTS)) {
    xy <- layout_with_kk(g)
    distsum <- distsum + as.matrix(dist(xy))
    centroid <- colMeans(xy)
    eccsum <- eccsum + sqrt(rowSums(sweep(xy, 2, centroid)^2))
  }
  vdist <- distsum / KK_RESTARTS
  vecc  <- eccsum / KK_RESTARTS

  netfeat[[i]] <- list(g = g, A = A, deg = deg, ecc = ecc_vec, esize = esize_vec,
                       comsize = comsize_vec, AA = AA, Katz = Katz, cosmat = cosmat,
                       dhom = dhom, srcos = srcos, netcent = netcent, netvar = netvar,
                       vdist = vdist, vecc = vecc)

  if (i %% 500 == 0) { cat(i, " "); print(Sys.time() - start) }
}
cat("\nPhase 1 done. "); print(Sys.time() - start)

# ===========================================================================
# PHASE 2: expansion bookkeeping (rid/sid/trial/netid/unode/id), mirroring
# expand_data.R's expansion exactly so rid lines up for the later merge.
# ===========================================================================
ntrial <- dim(dat)[1]
nnode  <- nnode_vec[dat$netid + 1]
sel    <- unlist(lapply(1:ntrial, function(x) rep(x, nnode[x])))
id     <- unlist(lapply(1:ntrial, function(x) 0:(nnode[x] - 1)))
rid    <- 1:length(sel)

sid   <- dat$sid[sel]
trial <- dat$trial[sel]
netid <- dat$netid[sel]
unode <- dat$unode[sel]

# drop the self row only (see header note: friend rows are left in and will
# simply not match anything when merged against fulldata.csv's rid set)
keep  <- id != unode
rid <- rid[keep]; sid <- sid[keep]; trial <- trial[keep]
netid <- netid[keep]; unode <- unode[keep]; id <- id[keep]

n_out <- length(rid)
vdist <- ecc <- uecc <- chuecc <- vecc <- uvecc <- dhom <- structsim <- rep(NA, n_out)
esize <- uesize <- chuesize <- aa <- katz <- srcos <- rep(NA, n_out)
comsize <- ucomsize <- netcent <- chgcent <- netvar <- chgvar <- rep(NA, n_out)

# ===========================================================================
# PHASE 3: per-trial lookups (vectorized -- no candidate loop needed) plus
# the per-candidate "change" loop for the four genuinely hypothetical-edge-
# dependent measures.
# ===========================================================================
usid <- unique(sid)
sn <- length(usid)
cat("\nPhase 3: per-trial features (", sn, " participants)\n", sep = "")
start <- Sys.time()

for (i in 1:sn) {
  sel1 <- sid == usid[i]
  trials <- unique(trial[sel1])
  for (j in trials) {
    sel2 <- trial[sel1] == j
    rows <- which(sel1)[sel2]     # indices into the *filtered* vectors

    cunode <- unode[rows[1]]
    cnetid <- netid[rows[1]]
    cids   <- id[rows]
    nf <- netfeat[[cnetid + 1]]

    su1 <- cunode + 1
    ci1 <- cids + 1

    vdist[rows]     <- nf$vdist[su1, ci1]
    ecc[rows]       <- nf$ecc[ci1]
    uecc[rows]      <- nf$ecc[su1]
    vecc[rows]      <- nf$vecc[ci1]
    uvecc[rows]     <- nf$vecc[su1]
    dhom[rows]      <- nf$dhom[su1, ci1]
    structsim[rows] <- nf$cosmat[su1, ci1]
    esize[rows]     <- nf$esize[ci1]
    uesize[rows]    <- nf$esize[su1]
    aa[rows]        <- nf$AA[su1, ci1]
    katz[rows]      <- nf$Katz[su1, ci1]
    srcos[rows]     <- nf$srcos[su1, ci1]
    comsize[rows]   <- nf$comsize[ci1]
    ucomsize[rows]  <- nf$comsize[su1]
    netcent[rows]   <- nf$netcent
    netvar[rows]    <- nf$netvar

    # -- change measures: require the hypothetical graph per candidate ------
    cno <- length(rows)
    cce <- ccc_ <- ccent <- cvar <- rep(NA, cno)
    for (k in 1:cno) {
      ccg  <- add_edges(nf$g, c(su1, ci1[k]))
      Accg <- as.matrix(as_adjacency_matrix(ccg))
      dccg <- degree(ccg)
      cce[k]   <- effective_size(Accg, dccg, su1)
      ccc_[k]  <- eccentricity(ccg, vids = su1)
      ccent[k] <- centr_degree(ccg)$centralization
      cvar[k]  <- var(dccg)
    }
    chuesize[rows] <- cce   - nf$esize[su1]
    chuecc[rows]   <- ccc_  - nf$ecc[su1]
    chgcent[rows]  <- ccent - nf$netcent
    chgvar[rows]   <- cvar  - nf$netvar
  }
  if (i %% 100 == 0) { cat(i, " "); print(Sys.time() - start) }
}
cat("\nPhase 3 done. "); print(Sys.time() - start)

out <- data.frame(rid, sid, trial, netid, id,
                  vdist, ecc, uecc, chuecc, vecc, uvecc,
                  dhom, structsim, esize, uesize, chuesize,
                  aa, katz, srcos, comsize, ucomsize,
                  netcent, chgcent, netvar, chgvar)
# CHANGED: no longer written to disk here (was write.csv(out,"extra_features.csv",...))
# -- `out` is consumed directly by PART C's merge below instead.
cat(sprintf("\nPart B done: out has %d rows, %d columns\n", nrow(out), ncol(out)))

# ===========================================================================
# PART C -- merge_features.R's logic, applied directly to the in-memory
# `fdat`/`out` data frames above instead of re-reading fulldata.csv/
# extra_features.csv from disk. See merge_features.R for the full
# rationale (NaN source, why 0-fill, why this column order); kept in sync
# with that script by hand since there are now two copies of this block --
# if you change one, change the other.
# ===========================================================================
shared_cols <- intersect(names(fdat), names(out))
stopifnot(setequal(shared_cols, c("rid", "sid", "trial", "netid", "id")))
alldata <- merge(fdat, out, by = shared_cols, sort = FALSE)
stopifnot(nrow(alldata) == nrow(fdat))

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
numeric_check_cols <- setdiff(names(alldata), id_cols)
stopifnot(all(sapply(alldata[numeric_check_cols], function(x)
    !is.numeric(x) || sum(is.na(x)) == 0)))

col_order <- c(
    "rid", "sid", "trial", "netid", "unode", "id", "nnode",
    "response", "rt", "binres",
    "age", "race", "ethnicity", "sex", "education", "nationality", "gender",
    "density", "assortativity", "apl", "modularity", "netcent", "netvar",
    "udegree", "ucloseness", "ubetweenness", "ueigenvector",
    "uclust", "uconstr", "uecc", "uesize", "ucomsize",
    "degree", "closeness", "betweenness", "eigenvector",
    "clust", "constr", "ecc", "esize", "comsize",
    "samecom", "edgbet", "udists", "mutfr", "dhom", "structsim", "aa", "katz", "srcos",
    "chgmod", "chgast", "chapl", "chgcent", "chgvar",
    "chub", "chuc", "chue", "chuclust", "chuconstr", "chuecc", "chuesize",
    "vdist", "vecc", "uvecc"
)
stopifnot(setequal(col_order, names(alldata)))
alldata <- alldata[, col_order]

write.csv(alldata, "alldata.csv", row.names = FALSE)
cat(sprintf("\nwrote alldata.csv: %d rows, %d columns\n", nrow(alldata), ncol(alldata)))
