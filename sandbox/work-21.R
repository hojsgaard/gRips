make_initc("unifips")
compileAttributes("unifips")
document("unifips")

devtools::load_all("unifips")

## load_all("unifips")
## build("unifips");install("unifips")


library(gRbase)
devtools::load_all("unifips")
#library(unifips)
library(microbenchmark)
data(math)
S    <- cov(math)
nobs <- nrow(math)
gl   <- list(1:3, 3:5)
em   <- matrix(c(1,2, 2,3, 1,3, 3,4, 3,5, 4,5), nrow=2)

EPS = 1e-12

devtools::load_all("unifips")
rg1 = gips_fit(S, gl, nobs=nobs, eps=EPS, engine="R")
rg2 = gips_fit(S, em, nobs=nobs, eps=EPS, engine="R")

rg1
rg2

(rg1$K - rg2$K)  %>% zapsmall

## Why so many itersations
devtools::load_all("unifips")
cg1 = gips_fit(S, gl, nobs=nobs, eps=EPS, engine="cpp")
cg2 = gips_fit(S, em, nobs=nobs, eps=EPS, engine="cpp")

cg1
cg2

(cg1$K - cg2$K)  %>% zapsmall


rg1$K - cg1$K


gg1 = fit_ggm(S, gl, nobs=nobs, eps=EPS, engine="cpp")
gg2 = fit_ggm(S, em, nobs=nobs, eps=EPS, engine="cpp")

gg1
gg2

K <- solve(S)
(ci$K - K)  %>% abs %>% max
(cg$K - K)  %>% abs %>% max
(ri$K - K)  %>% abs %>% max
(rg$K - K)  %>% abs %>% max

load_all("unifips")

library(gRbase)
library(magrittr)
edges <- list(c(1,2,3), c(2,3,4,5))


ed <- lapply(edges, combn_prim, 2)  %>% do.call(cbind, .)  %>% t
ed


ed <- lapply(edges, combn_prim, 2)

ff <- function(zz)do.call(cbind, zz)

ed  %>% ff

ed  %>% (function(zz)do.call(cbind, zz))



%>% {function(zz) 



%>% t




compileAttributes("unifips")
make_initc("unifips")
document("unifips")
load_all("unifips")
## install("unifips")

S1  <- cov(iris[,1:4])
S2  <- cov(iris[,1:4])

S <- matrix(rnorm(100*100), nr=100, nc=100)
S <- S + t(S)
S1 <- S
S2 <- S

    
S1 - S2

cc0 <- c(0,1)
H <- 1.01 * diag(1,2)
update_Sigma_(S1, cc0, H)
S1

update_Sigma2_(S2, cc0, H)
S2

S1 - S2

microbenchmark::microbenchmark(
                    update_Sigma_(S1, cc0, H),
                    update_Sigma2_(S1, cc0, H),
                    times=50
                )





library(gRbase)



library(unifips)

data(carcassall)
dat <- carcassall[,1:15]
dat <- dat[1:20,]
S <- cov2cor(cov(dat))
n <- nrow(dat)

f <- ~1:2+2:3+3:4+4:5+5:6+6:7+7:8+8:9+9:10+10:11+11:12+12:13+13:14+14:15+15:1
e <- rhsf2list(f) %>% do.call(rbind, .) 
storage.mode(e) <- "numeric"
e <- t(e)

f <- fit_ggm(S, e, n)
unifips:::mean_abs_diff_on_Emat(f$Sigma, S, e)

g <- fit_ggm(S, e, n, eps=1e-4, method="glasso")
unifips:::mean_abs_diff_on_Emat(g$Sigma, S, e)

f
g





















enquote(do.call(f, list(n=10000)))


do.call(benchmark,
  c(tests, list(replications=10000,
    columns=c('test', 'elapsed', 'replications'),
    order='elapsed')))

ee <- lapply(tests, function(tt) expression(tt))

e <- expression(s)
e2 <- substitute(e, list(s=tests[[1]]))


lapply(tests)

eval_one(tests[[1]], 10000)
eval_one(tests[[2]], 10000)
eval_one(tests[[3]], 10000)
    

eval_one(tests[[i]], 10000)










times <- 1000
dots <- tests
do.call(rbind, lapply(dots, function(fn) eval_one(fn, times=times)))




benchmark(tests)








lapply(tests, eval)


library(microbenchmark)

do.call(microbenchmark,
  c(tests, list(times=5))
)
  
microbenchmark::microbenchmark(
  tests,
  times=5 
)







do.call(benchmark,
  c(tests, )










lapply()

do.call(benchmark,
  c(lapply(a, function(a.) eval(e, list(x=a.))),
    list(replications=10000,
      columns=c('test', 'elapsed', 'replications'),
      order='elapsed'))
)


eval(expression(f(a[[1]])))
  


  # simple test functions used in subsequent examples
     random.array = function(rows, cols, dist=rnorm) 
                       array(dist(rows*cols), c(rows, cols))
     random.replicate = function(rows, cols, dist=rnorm)
                           replicate(cols, dist(rows))



tests = list(rep=expression(random.replicate(100, 100)), 
                  arr=expression(random.array(100, 100)))

do.call(benchmark,
             c(tests, list(replications=100,
                           columns=c('test', 'elapsed', 'replications'),
                           order='elapsed')))





library(microbenchmark)

microbenchmark(
  f(a[[1]]),
  f(a[[2]]),
  f(a[[3]]),
  times=10000)




tic()
f(a[[1]])
f(a[[2]])
f(a[[3]])
toc()


tests = list(f(a[[1]]), f(a[[2]]), f(a[[3]]))








tests = list(rep=expression(random.replicate(100, 100)), 
  arr=expression(random.array(100, 100)))
do.call(benchmark,
  c(tests, list(replications=100,
    columns=c('test', 'elapsed', 'replications'),
    order='elapsed')))



lapply(a, function(a.) do.call(f, list(a.)))

do.call(f, a)






devtools::load_all("unifips")

data(prostate)
names(prostate)

NCOL  <- 8
dat   <- prostate$x[,1:NCOL]

S <- cov(dat)
n.obs <- nrow(dat)

devtools::load_all("unifips")
f <- fit_ggm_calif(S, edges=NULL, n.obs, K=NULL, iter=1000L, eps=1e-6, print=FALSE)

solve(f$S)

S / f$w

solve(S) / f$wi


emat_sat  <- emat_saturated_model(1:NCOL)
emat_sat


emat <- emat_sat[, c(1,2,5,6,7,8)]

get_zero_edges(emat_sat, 8)





(1:nrow(S))
S <- cov(dat)

load_all("unifips")

res <- fit_ggm_all(S, emat, nobs=nrow(dat))

tt <- sapply(res, function(z) unname(z$details["time.elapsed"] ))
kk <- lapply(res, function(z) unname(z$K ))


max(abs(kk[[1]] - kk[[2]]))
max(abs(kk[[1]] - kk[[3]]))
max(abs(kk[[1]] - kk[[4]]))

r_ips = fit_ggm(S, emat, nobs=nrow(dat), engine="R",   method="ips")

do.call(fit_ggm, list(S, emat, nobs=nrow(dat), engine="c", method="ips"))















fit_ggm(S, emat, nobs=nrow(dat), engine="cpp")
fit_ggm(S, emat, nobs=nrow(dat), engine="R")


fit_ggm(S, emat, nobs=nrow(dat), engine="cpp", method="ips")
fit_ggm(S, emat, nobs=nrow(dat), engine="R", method="ips")



    



load_data <- function(...)
{
    if (length(unlist(list(...))) > 1) stop("Only one dataset can be loaded")
    e <- new.env()
    name <- data(..., envir = e)[1]
    e[[name]]
}
