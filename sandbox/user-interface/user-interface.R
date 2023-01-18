##########################################################
##
## Continuous interaction model (graphical Gaussian model)
##
##########################################################

#' @title Graphical Gaussian Model
#'
#' @description Graphical Gaussian Model based on modern fitting algorithm.
#' 
#' @name ggmu
#'
#' @param formula ...
#' @param data ...
#' @param marginal ...
#' @param fit ...
#' @param details ...
#' 
ggmu <- function(formula, data, marginal=NULL, fit=TRUE, details=0){

    dd <- .extract_cmod_data(data)
    
    vn <- colnames(dd$S)
    ans <- parse_gm_formula(formula, vn, marginal)
    ## Get varNames in the order matching to the data:
    vn <- vn[sort(match(ans$varNames, vn))]
    
    datainfo <- list(S=dd$S[vn, vn],
                     n.obs=dd$n.obs,
                     data=data)
    
    res <- list(glist          = ans$glist,
                varNames       = vn,
                datainfo       = datainfo,
                fitinfo        = NULL,
                isFitted       = FALSE
                )



    
    ## upd   <- .cModel_finalize(ans$glist, vn)  
    ## res[names(upd)] <- upd  
    ## class(res) <- c("cModel", "iModel")
    
    ## if (fit) fit(res) else res

    print(res$glist)
    print(res$datainfo$S)
    vv <- colnames(res$data$S)
    glisti <- lapply(res$glist, match, vv)
    glisti <- unlist(lapply(glisti, function(g) combnPrim(g, 2, FALSE)), recursive=FALSE)
    ff <- ips_modern(glisti, res$datainfo$S)
    ff
}

## FIXME: Fragile: .extract_cmod_data is copied directly from gRim

.extract_cmod_data <- function(data){
    if (inherits(data, "data.frame")){
        data <- cov.wt(data, method="ML")
    } else
        if (inherits(data, "list") && identical(names(data), c("cov", "center", "n.obs"))){
            ## OK
        } else
            stop("Can not proceed...")
    
            
    names(data)[1] <- "S"
    data
}


##ips_modern2_<- function (glist, S, iter=1000, eps=1e-3, method="ips", lambda=0, pos=FALSE, print=FALSE)     

