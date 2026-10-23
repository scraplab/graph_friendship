
# ===========================================================================
# expand_extra_features.R
#
# Produces a SEPARATE csv of additional features, one row per (sid, trial,
# candidate id) -- the same grain as fulldata.csv -- meant to be joined to
# fulldata.csv by `rid` (a merge, NOT a positional cbind: rid is assigned
# from the same pre-filter expansion order used in expand_data.R, so a
# merge is robust to any row-count/order differences between this script
# and expand_data.R's friend-filtering, whereas cbind would silently
# misalign on the first discrepancy). sid/trial/id are carried along too,
# as a redundant sanity check for the merge.
#
# Deliberate simplification vs. expand_data.R: this script does NOT
# replicate friend-detection at all. It only drops the self row
# (id == unode); friend rows are left in with (possibly meaningless, e.g.
# duplicate-edge) "change" values. That's fine -- those rows don't exist in
# fulldata.csv's rid set, so an inner merge drops them for free. This
# sidesteps the exact index-juggling that caused both bugs found in
# expand_data.R's friend-detection logic.
#
# Architecture: everything that depends only on network topology (not on
# which trial/participant is using that network) is computed ONCE PER
# NETWORK in the precompute phase below, then just looked up per trial.
# Only the four genuine "change" measures (chuesize, chuecc, chgcent,
# chgvar) need fresh computation per candidate, since those depend on the
# hypothetical graph with that specific candidate's edge added.
# ===========================================================================

if(!require(RJSONIO)) install.packages("RJSONIO"); require(RJSONIO)
if(!require(igraph))  install.packages("igraph");  require(igraph)

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

# ---- load data (same sources as expand_data.R) -----------------------------
dat    <- read.csv("net_data.csv")
graphs <- read.csv("../../graphs.csv")
graphs <- as.character(graphs$graphs)

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

# ===========================================================================
# write output
# ===========================================================================
out <- data.frame(rid, sid, trial, netid, id,
                  vdist, ecc, uecc, chuecc, vecc, uvecc,
                  dhom, structsim, esize, uesize, chuesize,
                  aa, katz, srcos, comsize, ucomsize,
                  netcent, chgcent, netvar, chgvar)
write.csv(out, "extra_features.csv", row.names = FALSE)
cat("\nwrote extra_features.csv:", nrow(out), "rows\n")
