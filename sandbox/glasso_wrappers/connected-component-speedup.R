load_all("unifips")
library(igraph)

indat  <- read.csv("data/prostate.csv")
dat  <- as.matrix(indat[,1:2000])
S <- cov2cor(cov(dat))

LAM <- .9

amat <- abs(S) > LAM
diag(amat) <- 0

ig <- graph_from_adjacency_matrix(amat)
##plot(ig)
com <- components(ig)


grp <- lapply(1:com$no, function(j){
    unname(which(com$membership == j))
})
grp



cl <- cheap_glasso(dat, lambda=LAM)
#' nodes with no parents
np <- which(sapply(cl, length)==0)
length(np)



#' Singleton connected components
sc <- unlist(grp[sapply(grp, length) == 1])
length(sc)


intersect(np, sc)
