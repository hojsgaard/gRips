.update_cc_parm <- function(cc, parm, S, eps, method="ips", lambda=0, pos=FALSE){
    ##cat(".update_cc_parm : ", toString(cc), "method : ", method, "\n")
    Scc       <- S[cc, cc]
    Kcc       <- parm$K[cc, cc]
    Sigma.old <- parm$Sigma ## FIXME Careful here
    Sigmacc   <- parm$Sigma[cc, cc]

    if (identical(method, "lasso"))
        lambdacc <- lambda[cc, cc] ##if (is.matrix(lambda)) lambda[cc, cc] else lambda
    else lambdacc <- NULL
    ## list(Kcc=Kcc, Scc=Scc, Sigmacc=Sigmacc)
    
    Kstar <- .getKstar(Sigmacc)   
    Laux  <- .getLaux(Kcc, Kstar)
    Saux  <- .getSaux(Scc, Laux, method=method, lambda=lambdacc, pos=pos)
    
    Ktilde <- .getKtilde(Saux)
    Kupd   <- .getKupd(Ktilde, Laux)
    ##print(list(Kupd=Kupd, Kupdinv=.sol(Kupd)))
    Haux   <- .getHaux(Saux, Kstar)
    ##print(list(cc=cc, Kaux=Kaux, Haux=Haux, det(Haux)))    
    str(list(Kstar=Kstar, Haux=Haux, Scc=Scc))    
    parm$K     <- .updateK(parm$K, Kupd, cc)
    parm$Sigma <- .updateSigma(parm$Sigma, Haux, cc)
    print(parm)
    parm    
}
