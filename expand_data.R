

# load packages -----------------------------------------------------------

if(!require(RJSONIO)) install.packages("RJSONIO"); require(RJSONIO)
if(!require(lmerTest)) install.packages("lmerTest"); require(lmerTest)
if(!require(igraph)) install.packages("igraph"); require(igraph)
if(!require(RColorBrewer)) install.packages("RColorBrewer"); require(RColorBrewer)

# load data ---------------------------------------------------------------

dat <- read.csv("net_data.csv")
head(dat)

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

graphs <- read.csv("../../graphs.csv")
graphs <- as.character(graphs$graphs)

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
write.csv(fdat,"fulldata.csv",row.names = F)
