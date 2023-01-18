library(unifips)
library(microbenchmark)
library(gRbase)

##source("helpers.R")
##load_all("unifips")

data(math)
dat <- scale(math)

data(carcass, package="gRbase")
dat <- scale(carcass)

data(carcassall, package="gRbase")
dat <- scale(carcassall[,1:15])

LAMBDA  <- .7
g1 <- glasso_lh(dat, lambda=LAMBDA)      ## Lauritzen-Højsgaard                     
g2 <- glasso_fht(dat, lambda=LAMBDA)     ## Friedman-Hastie-Tibshirani              
g3 <- glasso_ns(dat, lambda=LAMBDA)      ## Neighbourhood selection

par(mfcol=c(2,2))
plot(g1); title("Lau-Hoj")
plot(g2); title("FHT")
plot(g3); title("Neigh-selection")








