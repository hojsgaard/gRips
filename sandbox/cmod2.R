#' @title Fit graphical Gaussian models
#' 
#' @description Functions for fitting graphical Gaussian models
#'     (ggms); deliberately with a simple interface. These are not
#'     intended to be called directly by the user but by more proper
#'     "modelling functions".
#'
#' @name cmod2
#'
#' @param formula A specification of a graphical Gaussian model.
#' @param data A dataframe
#' @param marginal A subset of the variables.
#' @param fit Should the model be fitted?
#' @param details The amount of details.
#' @param aux List specifying various restrictions on models.
#' 

## 'unified ips' (c++); requires glist
## 'unified ips' (R); requires emat (2 x p)
## Functions here call functions in fit-worker.R

#' @export
cmod2 <- function (formula, data, marginal=NULL, fit=TRUE, details=0, aux=list()) 
{
    #cat(" ... 1) Find covariance matrix of correct dimension, find n.obs\n")

    dd  <- gRim::extract_cmod_data(data)
    vn  <- colnames(dd$S)
    ans <- gRim::parse_gm_formula(formula, vn, marginal)
    vn  <- vn[sort(match(ans$varNames, vn))]

    #str(list(dd=dd, vn=vn))
    
    #cat(" ... 2) Find list of generators\n")

    ## Switch to indices
    glist <- lapply(ans$glist, match, ans$varNames)
    ## Generators consisting of one edge must be removed (otherwise
    ## combn_prim acts differently)
    glist <- glist[sapply(glist, length) > 1]
    ## Turn into edges (not necessary for standard models, though)
    glist <- lapply(glist, combn_prim, 2)
    #print(list(glist=glist))

    #cat(" ... 3) Convert generators to 2xp edge matrix\n")

    ## Gather in 2 x p matrix
    glist <- do.call(cbind, glist)
    ## Remove duplicates
    glist <- glist[, !duplicated(glist, MARGIN=2)]

    ##print(list(glist=glist))

    #cat(" ... 4) Put everything into 'model object'\n")
    
    out <- list(S = dd$S[vn, vn], iter=1000, nobs = dd$n.obs,
                edges=glist, zero=NULL, aux=aux)
    
    ##datainfo <- list(S = dd$S[vn, vn], n.obs = dd$n.obs)
    
    #cat(" ... 5) Call fit if requested.\n")

    #oo <<- out
    if (fit) out <- do.call(gips18, out)

    out$glist <- glist
    out
}
