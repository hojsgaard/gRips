
## -----------------------------------------------
##
## API like functions below here
##
## -----------------------------------------------

.covips_ggm_update_cc_parm <- function(S, cc, parm, Scc_inv_list, j)
{

    ## sprintf("parm-start:\n") %>% cat(); print(parm)
    
    Scc       <- S[cc, cc] ## FIXED
    Ktilde    <- Scc_inv_list[[j]]
    
    Kcc       <- parm$K[cc, cc]
    Sigmacc   <- parm$Sigma[cc, cc]
    
    Kstar  = .sol(Sigmacc);
    Laux   = Kcc - Kstar;
    
    Kupd   = Ktilde + Laux;
    Haux   = Kstar - Kstar %*% Scc %*% Kstar;
    
    parm$K[cc, cc] = Kupd;
    ## str(list(Kstar=Kstar, Haux=Haux, Scc=Scc))
    
    parm$Sigma     = parm$Sigma - parm$Sigma[, cc] %*% Haux %*% parm$Sigma[cc, ] 

    ## sprintf("parm-end:\n") %>% cat(); print(parm)
    parm
}


.maxKdiff <- function(K_old, K_new){
   max(abs(K_old - K_new)) / max(abs(K_old))
}

.sumKdiff <- function(K_old, K_new){
   sum(abs(K_old - K_new)) / sum(abs(K_old))
}

.getKstar <- function(Sigmacc){
    .sol(Sigmacc)
}

.getLaux <- function(Kcc, Kstar){
    Kcc - Kstar   
}

.getSaux <- function(Scc, Laux, method="ips", ...){
    ##cat("method : ", method, "\n")
    switch(method,
           "ips"    = {.getSaux_ips   (Scc, Laux)},
           "mtp2"   = {.getSaux_mtp2  (Scc, Laux)},
           "lasso"  = {.getSaux_lasso (Scc, Laux, ...)},
           "hybrid" = {.getSaux_hybrid(Scc, Laux, ...)}
           )
}

.getKtilde <- function(Saux){
  .sol(Saux)
}

.getKupd <- function(Ktilde, Laux){
    Ktilde + Laux
}

.getHaux <- function(Saux, Kstar){
    Kstar - Kstar %*% Saux %*% Kstar
}

.updateK <- function(K, Kupd, cc){
    K[cc, cc] <-  Kupd
    K
}

.updateSigma <- function(Sigma, Haux, cc){
    Sigma - Sigma[, cc] %*% Haux %*% Sigma[cc, ] 
    ##print(list(cc=cc, Sigma.old=Sigma.old[, cc], Sigma.cc=out))
}

.updateSigma2 <- function(Sigma, Haux, cc){
    SS <- Sigma[, cc]
    Sigma - SS %*% Haux %*% t.default(SS)
}

## 
## getSaux 'methods' 
## 

.get_suv_adjust <- function(Scc, luv){
    (2 * Scc[1, 1] * Scc[2, 2] * luv)/(sqrt(1 + 4 * luv^2 * Scc[1, 1] * Scc[2, 2]) + 1)
}

.getSaux_ips <- function(Scc, Laux){
    ## cat(".getSaux_i \n"); print(Scc)
    Scc
}

## changed this to new expression not involving ifelse.

.getSaux_lasso <- function(Scc, Laux, lambdacc, pos=FALSE){ 
    ##cat(".getSaux_lasso", "lambda : ", lambda, "\n")
    suv <- Scc[1, 2]
    luv <- Laux[1, 2]
    lambdauv <- lambdacc[1, 2]
    suv_star <- .get_suv_adjust(Scc + lambdacc, luv)

    #cat("suv_star : ", suv_star, "\n"); print(list(Scc))
    
    suv_tilde <- suv_star
    if (suv + lambdauv < suv_star){
        suv_tilde <- suv + lambdauv
    } else if (suv - lambdauv > suv_star){
        suv_tilde <- suv - lambdauv
    }

    Saux <- Scc + lambdacc
    Saux[c(2, 3)] <- suv_tilde
    Saux
}

## changed expression for updating

.getSaux_mtp2 <- function(Scc, Laux){ 

    ##cat(".getSaux_mtp2 \n")  
    suv <- Scc[1, 2]
    luv <- Laux[1, 2]
   
    limuv <- suv / (Scc[1, 1] * Scc[2, 2] - suv^2)

    if (luv <= limuv){
        suv_tilde <- suv
        #cat("case 1 : ", suv_tilde, "\n")
    } else {
        suv_tilde <- .get_suv_adjust(Scc, luv)
        #cat("case 2 : ", suv_tilde, "\n")
    }
    
    Saux <- Scc
    Saux[c(2, 3)] <- suv_tilde
    #print(list(Saux=Saux))
    Saux
}

.getSaux_hybrid <- function(Scc, Laux, lambdacc, pos=FALSE){
    ##cat(".getSaux_hybrid \n")

    ##print(list(Scc=Scc, Laux=Laux, lambdacc=lambdacc, pos=pos))
    suv <- Scc[1, 2]
    luv <- Laux[1, 2]
    lambdauv <- lambdacc[1, 2]
  
    suv_star <- .get_suv_adjust(Scc + lambdacc, luv)
   
    suv_tilde <- suv_star
    if (suv + lambdauv < suv_star){
      suv_tilde <- suv + lambdauv
    } else if (suv - lambdauv > suv_star){
      suv_tilde <- suv - lambdauv
    }
    
    Saux <- Scc + lambdacc
    
    if (pos) {  
        limuv <- suv_tilde / (Saux[1, 1] * Saux[2, 2] - suv_tilde^2)
        if (luv <= limuv){
            suv_tilde <- suv_tilde       ## FIXME : looks strange
        } else {
            suv_tilde <- suv_star
        }
    }
    
    Saux[c(2, 3)] <- suv_tilde
    Saux    
}







