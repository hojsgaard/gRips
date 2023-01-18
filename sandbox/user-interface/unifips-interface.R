library(unifips)
load_all("unifips")

##build("unifips");install("unifips")
## math

data(math)
dat <- scale(math)
## FIXME :gRim:::.extract_cmod_data fails if data is matrix (e.g. after scaling data)
load_all("unifips"); out <- cmod2(~.^., data=math)

out1 <- gRim::cmod(~me:ve:al + ve:al + an, data=math, fit=F)

load_all("unifips"); out2 <- cmod2(~me:ve:al + ve:al + an, data=math)

out1
out2

str(out1)
str(out2)

terms(out1); terms(out2)
formula(out1); formula(out2)

out <- cmod2(~.^., data=math)
load_all("unifips");

#out1 <- gRim::cmod(~.^1, data=math)  ## FAILS
#out2 <- cmod2(~.^1, data=math)  ## FAILS

out$glist




## Switch to indices
glist <- lapply(out$glist, match, out$varNames)

## Genetors consisting of one edge must be removed (otherwise
## combnPrim acts differently)
glist <- glist[sapply(glist, length) > 1]

## Turn into edges (not necessary for standard models, though)
glist <- lapply(glist, combnPrim, 2)
## Gather in 2 x p matrix
glist <- do.call(cbind, glist)
## Remove duplicates
glist <- glist[, !duplicated(glist, MARGIN=2)]
## 2 x p edge matrix
glist


ff <- unifips18(glist, S=out$datainfo$S, aux=list(method="mtp2"))
ff
ff$K

out$datainfo


## Bruger alene de variable, der er specificerede
gRim::cmod(~me:ve:al + ve:al:an, data=math)
