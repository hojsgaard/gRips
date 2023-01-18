USECOL <- 1:4
E    <- emat_saturated_model(USECOL)
E0 <- E - 1

Elist <- matrix2list(E, F)
E0list <- matrix2list(E0, F)

S1 <- S
S2 <- cf$Sigma

devtools::load_all("../unifips")
diff_on_Emat(S1, S2, E)
diff_on_Emat_(S1, S2, E)
diff_on_Emat_(S1, S2, E0, 0)

diff_on_Elist(S1, S2, Elist)
diff_on_Elist_(S1, S2, Elist)
diff_on_Elist_(S1, S2, E0list, 0)

max_abs_diff_on_Emat(S1, S2, E)
max_abs_diff_on_Emat_(S1, S2, E)
max_abs_diff_on_Emat_(S1, S2, E0, 0)

max_abs_diff_on_Elist(S1, S2, Elist)
max_abs_diff_on_Elist_(S1, S2, Elist)
max_abs_diff_on_Elist_(S1, S2, E0list, 0)

mean_abs_diff_on_Emat(S1, S2, E)
mean_abs_diff_on_Emat_(S1, S2, E)
mean_abs_diff_on_Emat_(S1, S2, E0, 0)

#mean_abs_diff_on_Elist(S1, S2, Elist)
mean_abs_diff_on_Elist_(S1, S2, Elist)
mean_abs_diff_on_Elist_(S1, S2, E0list, 0)


#' ## Check matrix norms
#' 

#' del: Change in parameter value at convergence
#'
#' thr: Threshold for convergence. Default value is 1e-4.  Iterations
#'          stop when average absolute parameter change is less than
#'          thr * ave(abs(offdiag(s)))

library(unifips)
library(doBy)
library(microbenchmark)
make_initc("unifips")
compileAttributes("unifips")
document("unifips")
load_all("unifips")
a<-unifips:::doit(c("aa", "bb", "asas"))
a 

a<-unifips:::doit(c("aa", "bb", "asas1"))
a



make_initc("unifips")
load_all("unifips")
unifips:::find_str_("asas", c("aa", "bb", "asas"))
unifips:::find_str_("asas", c("aa", "asas", "bb", "asas"))
unifips:::find_str_("asas2", c("aa", "bb", "asas"))

library(MASS)

data(prostate)
dat <- prostate$x

## data(math)
## dat <- math

USECOL <- 1:32
p <- max(USECOL)
S <- cov(dat[, USECOL])
S <- cov2cor(S)
d <- nrow(S)
n <- nrow(dat)

set.seed(141164)
ed <- emat_random_model(USECOL, p=0.7)
ne <- ncol(ed)
ne
THR <- 1e-4


load_all("unifips")
m1 <- fit_ggm(S, ed, nobs=n, eps=THR, iter=1222)
nobs(m1)
logLik(m1)
AIC(m1)
BIC(m1)


load_all("unifips"); unifips::list_names_(list())
unifips::list_names_(list(a=1, b=list(1,2,3)))
unifips::list_names_(list(a=1, b=list(1,2,3), 9999))


load_all("unifips")
m1 <- fit_ggm(S, ed, nobs=n, eps=THR, iter=1222, method="glasso")
m2 <- fit_ggm(S, ed, nobs=n, eps=THR, iter=1222, convcrit = 1)
m3 <- fit_ggm(S, ed, nobs=n, eps=THR, iter=1222, convcrit = 3, aux=list(conv_ref=logLik(m1)))

m1; m2; m3



base <- restrict(fit_ggm, list(S=S, edges=ed, nobs=n, eps=THR, iter=1222))
base()

gbase <- restrict(base, list(method="glasso"))

g <- gbase()

gbase()

uuu <- restrict(base, list(method="gips", convcrit=3, aux=list(conv_ref=logLik(g))))

uuu()

microbenchmark(
  gbase(),
  uuu(),
  m1 = fit_ggm(S, ed, nobs=n, eps=THR, iter=1222, method="glasso"),
  m3 = fit_ggm(S, ed, nobs=n, eps=THR, iter=1222, convcrit = 3, aux=list(conv_ref=logLik(m1)))
)




load_all("unifips")
m2 <- fit_ggm(S, ed, nobs=n, eps=THR, iter=1222, convcrit = 1)

load_all("unifips")
m2 <- fit_ggm(S, ed, nobs=n, eps=THR, iter=1222, convcrit = 1, aux=list(a=10, b=12, conv_ref=99.0, convref=-111)
)

m2 <- fit_ggm(S, ed, nobs=n, eps=THR, iter=1222, convcrit = 2)




m2 <- fit_ggm(S, ed, nobs=n, eps=THR, iter=m1$details$iter - 1, method="glasso")




del <- m1$details["del"]

S1 <- m1$Sigma
S2 <- m2$Sigma
#S1[1,1]=S1[1,1]+.1
S1
S2

load_all("unifips")
del
nu <- 0.5*sum(abs(S1-S2)) / (p^2-p)
nu

de <- sum(abs(S[lower.tri(S)])) / (p*(p-1)/2)
de

zz <- nu/de
zz


zz/del


1 + 2/(p-1)







max_abs_diff_(S1, S2)








max_abs_(S1)
max_abs(S1)
max(abs(S1))

max_abs_diag_(S1)
max_abs_diag(S1)
max(abs(diag(S1)))

load_all("unifips")
del
max_abs_diff_(S1, S2)


# max_abs_diff(S1, S2) ## Does not exist






max_abs_diff_rel_(S1, S2)
# max_abs_rel_diff(S1, S2) ## Does not exist
max(abs(S1-S2)) / max(abs(S1))

load_all("unifips")
max_abs_diag_diff_(S1, S2)
max_abs_diag_diff(S1, S2)
max(abs(diag(S1 - S2)))

max_diag_diff_(S1, S2)
max_diag_diff(S1, S2)
max(diag(S1 - S2))




library(microbenchmark)
microbenchmark(
  max_abs_diff_(S1, S2),
  # max_abs_diff(S1, S2) ## Does not exist
  max(abs(S1-S2)),
  ##
  max_abs_diag_diff_(S1, S2),
  max_abs_diag_diff(S1, S2),
  max(abs(diag(S1 - S2))),
  ##
  max_diag_diff_(S1, S2),
  max_diag_diff(S1, S2),
  max(diag(S1 - S2))
)


system.time({for (i in 1:10000) })










abs(S1-S2)  %>% sum




mad <- unifips:::max_abs_diff_(S1, S2)



(del/(d-0)) / mad


unifips:::max_abs_diff_on_E_(S1, S2, ed)
max_abs_diff_on_E(S1, S2, ed)

avg_abs_diff_on_E(S1, S2, ed)




unifips:::matrix_locs(S1 - S2, ed-1)

diff_on_E(S1, S2, ed)

diff_on_E_(S1, S2, ed)

crit_value(S1, S2, ed)

c(unifips:::matrix_locs(S1 - S2, ed-1), diag(S1-S2))


S[t(ed)]

make_initc("unifips")
load_all("unifips")






ee <- matrix(c(1,2,2,3,2,4), byrow=T,ncol=2)
S[ee]





