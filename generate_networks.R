
# Load packages -----------------------------------------------------------

if(!require(igraph)) install.packages("igraph"); require(igraph)
if(!require(RJSONIO)) install.packages("RJSONIO"); require(RJSONIO)

# simulate graphs ---------------------------------------------------------


# generate initial population of graphs
set.seed(1)
nnet <- 100000
glist <- list()
for (i in 1:nnet){
  nnode <- sample(12:24,1)
  radius <- .3 + (.4-.3)*runif(1)
  glist[[i]]<-list()
  glist[[i]]$nnode <- nnode
  glist[[i]]$radius <- radius
  glist[[i]]$g <- sample_grg(nodes = nnode, radius = radius, torus = F)
  print(i)
}

# graphs must be connected, low-density, with high modularity and degree assortativity
connected <- unlist(lapply(glist,function(x) is.connected(x$g)))
glist <- glist[connected]
dense <- unlist(lapply(glist,function(x) edge_density(x$g)))
glist <- glist[dense < .25]
assort <- unlist(lapply(glist,function(x) assortativity_degree(x$g)))
glist <- glist[assort > 0]
modul <- unlist(lapply(glist,function(x) infomap.community(x$g)$modularity))
glist <- glist[modul > .1]

glist <- glist[1:10000]
#save(glist,file = "glistnetwork.Rdata")
plot(glist[[sample(length(glist),1)]]$g)

ng <- 10000
jsons <- rep(NA,ng)
for (i in 1:ng){
  cg <- glist[[i]]$g
  ccg <- cluster_optimal(cg)
  jlist <- list()
  
  # network level properties
  jlist$id <- i-1
  jlist$density <- edge_density(cg)
  jlist$assortativity <- assortativity.degree(cg,F)
  jlist$"average path length" <- average.path.length(cg,F) 
  jlist$modularity <- ccg$modularity
  
  # node level properties
  nodes <- cbind(id=1:vcount(cg)-1,degree=degree(cg),
                      closeness=closeness(cg),betweenness=betweenness(cg),
                      eigenvector=eigen_centrality(cg)$vector,community=ccg$membership)
  jlist$nodes <- nodes
  
  # edge level properties
  links <- get.edgelist(cg)
  colnames(links) <- c("target","source")
  jlist$links <- links-1
  
  # save as json
  jsons[i] <- toJSON(jlist)
  print(i)
}
out <- data.frame(id=1:ng-1,graphs=jsons)
#write.csv(out,file = "graphs.csv",row.names=F)

# write out multiple files to avoid import timeout
for (i in 1:10){
  sel <- ((i-1)*1000 + 1):(i*1000)
  fname <- paste("graphs_",as.character(i),".csv",sep="")
  write.csv(out[sel,],file=fname,row.names=F)
}


# calculate states for normalization --------------------------------------

dg <- unlist(lapply(glist,function(x) mean(degree(x$g))))
cl <- unlist(lapply(glist,function(x) mean(closeness(x$g))))
bt <- unlist(lapply(glist,function(x) mean(betweenness(x$g))))
ei <- unlist(lapply(glist,function(x) mean(eigen_centrality(x$g)$vector)))

mean(dg)
mean(cl)
mean(bt)
mean(ei)

dg <- unlist(lapply(glist,function(x) sd(degree(x$g))))
cl <- unlist(lapply(glist,function(x) sd(closeness(x$g))))
bt <- unlist(lapply(glist,function(x) sd(betweenness(x$g))))
ei <- unlist(lapply(glist,function(x) sd(eigen_centrality(x$g)$vector)))

sd(dg)
sd(cl)
sd(bt)
sd(ei)
