

library(gRbase)
make_initc("unifips")
compileAttributes("unifips")

load_all("unifips"); 

data(math)
S  <- cov(math)
nobs <- nrow(math)
gl <- list(1:3, 3:5)
em <- matrix(c(1,2, 2,3, 1,3, 3,4, 3,5, 4,5), nrow=2)


data(personality)
pers <- personality[,1:10]
S <- cov(pers)
nobs <- nrow(S)
em  <- emat_loop_model(1:nrow(S))

em  <- emat_saturated_model(1:nrow(S))

load_all("unifips"); 
e <- matrix(NA, nrow=2, ncol=0)
g <- list() 


em <- matrix(c(1,2, 2,3, 1,3, 3,4, 3,5, 4,5), nrow=2)
gl <- list(1:3, 3:5)

ed <- em

load_all("unifips"); fit_ggm(S, ed, nobs=nobs, method="ips")             
load_all("unifips"); fit_ggm(S, ed, nobs=nobs, method="gips")            
load_all("unifips"); fit_ggm(S, ed, nobs=nobs, engine="R", method="ips") 
load_all("unifips"); fit_ggm(S, ed, nobs=nobs, engine="R", method="gips")

microbenchmark::microbenchmark(
                    fit_ggm(S, ed, nobs=nobs, method="ips"),             
                    fit_ggm(S, ed, nobs=nobs, method="gips"),
                    fit_ggm(S, ed, nobs=nobs, engine="R", method="ips"), 
                    fit_ggm(S, ed, nobs=nobs, engine="R", method="gips")
)










load_all("unifips"); make_clist_(S, glist)

glist <- list(1:2, 2:3, 3:4, 4:5, c(5,1)) ## Pas på; den modificeres af c-kode
K <- diag(1/diag(S))
.c_ips_ggm_(S, glist, K, nobs=nrow(math))

glist <- list(1:2, 2:3, 3:4, 4:5, c(5,1)) ## Pas på; den modificeres af c-kode
K <- diag(1/diag(S))
.c_gips_ggm_(S, glist, K, nobs=nrow(math))



aa <- (1:nrow(S))


clist <- lapply(glist, function(g) aa[-g])





glist <- lapply(glist, `-`, 1)
clist <- lapply(clist, `-`, 1)

make_initc("unifips")
compileAttributes("unifips")

S <- cov(math)
K <- diag(1/diag(S))
load_all("unifips")


Si <- solve(K)

lapply(glist, function(g){
    g1 <- g + 1
    max(abs(Si[g1, g1] - S[g1, g1]))


})

E <- do.call(cbind, lapply(glist, combn_prim, 2)) + 1

max



























document("unifips")






library(unifips)
library(gRim)



devtools::load_all("unifips");
data(math)



mm1 <- cmod(~.^., data=pers)
mm2 <- cmod2(~.^., data=pers)

vn <- names(pers)
glist <- emat_loop_model(vn)  %>% as.data.frame %>% as.list


microbenchmark::microbenchmark(
    cmod(glist, data=pers),
    cmod2(glist, data=pers)    
)

mm1 <- cmod(glist, data=pers)
mm2 <- cmod2(glist, data=pers)









## Global settings
EPS  <- 1e-6
ITER <- 5000

data(math)
S <- cov(math)

data(personality)
S <- cov(personality)

devtools::load_all("unifips")
m01 <- unifips18(S, nobs=nrow(math), edges=list())
m02 <- unifips18(S, edges=matrix(nrow=2, ncol=0), iter=ITER, eps=EPS)

logLik(m01)
logLik(m02)

AIC(m01)
BIC(m01)

m03 <- unifips18(S, edges=list(1:5), iter=ITER, eps=EPS)
m04 <- unifips18(S, edges=combn(1:5, 2), iter=ITER, eps=EPS)

logLik(m03)
logLik(m04)

m01$K  %>% as_K2graph %>% plot
m02$K  %>% as_K2graph %>% plot
m03$K  %>% as_K2graph %>% plot
m04$K  %>% as_K2graph %>% plot

## Data
S <- .make_breastS()
USE <- 1:20

## Tomme model:
m01 <- unifips18(list(), mm3$S, iter=ITER, eps=EPS)
m02 <- unifips18(matrix(nrow=2, ncol=0), mm3$S, iter=ITER, eps=EPS)
logLik(m01)
logLik(m02)

## Mættet model (på en marginal)
m03 <- unifips18(list(1:5), mm3$S, iter=ITER, eps=EPS)
m04 <- unifips18(combn(1:5, 2), mm3$S, iter=ITER, eps=EPS)
logLik(m03)
logLik(m04)

## NB: logLik beregnes "on the fly" og findes altså ikke model
## objekt. Det kan laves om.

m05 <- unifips18(mm3$sat, mm3$S, iter=ITER, eps=EPS)
m06 <- unifips18(mm3$sat, mm3$S, iter=ITER, eps=EPS, aux=list(method="mtp2"))
logLik(m05)
logLik(m06)

m07 <- unifips18(mm3$loop, mm3$S, iter=ITER, eps=EPS)
m08 <- unifips18(mm3$loop, mm3$S, iter=ITER, eps=EPS, aux=list(method="mtp2"))
logLik(m07)
logLik(m08)

m05$K  %>% K2graph %>% plot
m06$K  %>% K2graph %>% plot
m07$K  %>% K2graph %>% plot
m08$K  %>% K2graph %>% plot


