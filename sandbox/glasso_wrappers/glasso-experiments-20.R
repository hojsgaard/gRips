library(unifips)
library(gRbase)
library(microbenchmark)
load_all("unifips")

##build("unifips");install("unifips")
## math

data(math)
dat <- scale(math)
LAMBDA  <- .3
g <- glasso_ns(dat, lambda=LAMBDA)
plot(g)

## prostate
USE   <- 1:10
LAMBDA  <- .3

data(prostate)
dat   <- as.data.frame(scale(prostate$x[,USE])) ## FIXME: prostate$x is matrix; not dataframe
S <- cov(dat); attr(S, "n") <- nrow(dat)   ## FIXME: Ugly
## Saturated model

r2 <- glasso_ns(dat, lambda=LAMBDA, fit=F) ## Assumes data to be dataframe
r3 <- glasso_lh(S, edges=r2$emat)          ## model argument missing; input dataframe as well.
r4 <- glasso_lh(S, lambda=LAMBDA)
r5 <- glasso_lh(S, lambda=LAMBDA, penalize.diagonal=F)


## Should be the same (and they are)
r2$emat
r3$emat

## Should be similar to above when lambda is small
r4$emat
r5$emat

emat_compare(r2$emat, r3$emat)
emat_compare(r4$emat, r5$emat)
c(ncol(r2$emat), ncol(r3$emat), ncol(r4$emat), ncol(r5$emat))

load_all("unifips")
gl1 <- glasso_fht(S, lambda=LAMBDA, penalize.diagonal=!TRUE)
un1 <- glasso_lh(S, lambda=LAMBDA, penalize.diagonal=!TRUE)
cg1 <- glasso_ns(dat, lambda=LAMBDA)

K.gl1 <- gl1$K
K.un1 <- un1$K

gl1$emat
un1$emat

e1b = unifips18(esat, S, aux=list(method="lasso", lambda=LAMBDA))
e1d = unifips18(esat, S, check=2, aux=list(method="lasso", lambda=LAMBDA))

load_all("unifips");

ed = unifips18(esat, S, aux=list(method="mtp2"))


LAMBDA=.3
microbenchmark::microbenchmark(
e1a = glasso_fht(S, lambda=LAMBDA, penalize.diagonal=!FALSE),
e1b = unifips18(esat, S, aux=list(method="lasso", lambda=LAMBDA)),
e1d = unifips18(esat, S, check=2, aux=list(method="lasso", lambda=LAMBDA)),
e1c = glasso_ns(dat, lambda=LAMBDA),
ed = unifips18(esat, S, aux=list(method="mtp2")),
times=5)


e1b = unifips18(esat, S, aux=list(method="lasso", lambda=LAMBDA))
e1d = unifips18(esat, S, check=2, aux=list(method="lasso", lambda=LAMBDA))

e1b
e1d


## Caroline data

load("data/lcm.RData")
load("data/ldm.RData")

USE <- 1:500
dat <- ldm[, USE]
S <- cov2cor(lcm)[USE, USE]
esat <- t(names2pairs(1:nrow(S), result="matrix"))
load_all("unifips")

LAMBDA <- .3

system.time({
    gl1 <- glasso::glasso(S, rho=LAMBDA, penalize.diagonal=TRUE) %>% iznogood()})

system.time({
    un1 <- unifips18(esat, S, aux=list(method="lasso", lambda=LAMBDA, penalize.diagonal=TRUE))})


system.time({
    un2 <- unifips18(esat, S, aux=list(method="mtp2"))
})


cg1 <- cheap_glasso(dat, lambda=LAMBDA)

lpen <- function(K, S, LAMBDA){
   log(det(K)) - sum(K*S) - sum(LAMBDA * abs(K))
}


L <- matrix(LAMBDA, nrow=nrow(S), ncol=ncol(S)); diag(L) <- 0

lpen(gl1$K, S, L)
lpen(un1$K, S, L)

d <- gl1$K - un1$K

which(abs(d) > 1e-4, arr.ind=T)


e1  <- gl1$emat
e2  <- un1$emat


n1 <- apply(e1 * c(1, 10000),2, sum)

n2 <- apply(e2 * c(1, 10000),2, sum)

intersect(n1, n2)
