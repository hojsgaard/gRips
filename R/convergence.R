
max_abs <- function(S){
    max(abs(S))
}

max_abs_diag <- function(S){
    max(abs(diag(S)))
}

max_abs_diag_diff <- function(S, Sigma){
    max(abs(diag(S) - diag(Sigma)))
}

max_diag_diff <- function(S, Sigma){
    max(diag(S) - diag(Sigma))
}


diff_on_emat <- function(Sigma, S, E){
  dif <- Sigma - S
  out <- c(diag(dif), dif[t(E)])
  out
}

diff_on_Elist <- function(Sigma, S, E){
    ## E_ <- do.call(cbind, E)
    E <- E[sapply(E, length) > 1]
    E_ <- t.default(unique(t.default(do.call(cbind, lapply(E, combn_prim, 2)))))
    diff_on_emat(Sigma, S, E_)
}

max_abs_diff_on_emat <- function(Sigma, S, E){
    max(abs(diff_on_emat (Sigma, S, E)))
}

max_abs_diff_on_Elist <- function(Sigma, S, E){
    E_ <- do.call(cbind, E)
    max_abs_diff_on_emat(Sigma, S, E_)
}

mean_abs_diff_on_emat <- function(Sigma, S, E){
  mean(abs(diff_on_emat (Sigma, S, E)))
}

max_diff_on_emat <- function(Sigma, S, E){
  out = max_diag_diff(S, Sigma)
  for (j in 1:ncol(E)){
      u = E[1, j]; v = E[2, j]
      d = Sigma[u, v] - S[u, v]
      out = max(out, d)
  }
  out
}

avg_diff_on_emat <- function(Sigma, S, E){
  max(abs(diff_on_emat (Sigma, S, E)))
}

max_abs_diff_on_EK <- function(Sigma, S, E, K){
    out = 0
    for (i in 1:nrow(Sigma))
        out = max(out, abs(Sigma[i, i] - S[i, i]) * K[i, i])

    for (j in 1:ncol(E)){
        u = E[1, j]; v = E[2, j]
        d = abs((Sigma[u, v] - S[u, v]) * K[u, v])
        out = max(out, d)
        ##cat("out: ", out, "\n")
    }
    out
}


## ##############################################
##
## ### does model fit to data
##
## ##############################################

## .does_model_fit_to_data_ips <- function(EE, parm, data, eps){
##     d <- max_abs_diff_on_emat_(data, parm$Sigma, EE)
##     s <- max(abs(diag(data))) 
##     out <- d < s * eps
##     ##cat(sprintf(".does_model_fit_to_data_ips (R) d=%15e, s=%15e, out=%i\n", d, s, out))
##     out
## }

## .does_model_fit_to_data_mtp2 <- function(EE, parm, data, eps){
##     ## cat(".does_model_fit_to_data_mtp2\n")
##     d1 <- max_diff_on_emat_(data, parm$Sigma, EE) 
##     d2 <- max_abs_diff_on_EK_(data, parm$Sigma, EE, parm$K)
##     out <- !((d1 > 0) || (d2 > eps)) 
##     ##cat(sprintf(".does_model_fit_to_data_mtp2 (R) d1=%15e, d2=%15e, out=%i\n", d1, d2, out))
##     out
## }

## .does_model_fit_to_data_lasso <- function(EE, parm, data, eps){
##     cat("## NOT IMPLEMENTED\n")
## }

## .does_model_fit_to_data_hybrid <- function(EE, parm, data, eps){
##     cat("## NOT IMPLEMENTED\n")
## }

## .does_model_fit_to_data <- function(EE, parm, data, eps, method){
##     switch(method,
##            "ips"    ={.does_model_fit_to_data_ips   (EE, parm, data, eps)},
##            "mtp2"   ={.does_model_fit_to_data_mtp2  (EE, parm, data, eps)},
##            "lasso"  ={.does_model_fit_to_data_lasso (EE, parm, data, eps)},
##            "hybrid" ={.does_model_fit_to_data_hybrid(EE, parm, data, eps)})
## }

## ##############################################
##
## ### does edge fit to data
##
## ##############################################


## .does_edge_fit_to_data_ips <- function(cc, parm, data, eps){
##     ##cat(".does_edge_fit_to_data_ips:  edge : ", toString(cc), "\n")
##     d <- max(abs(parm$Sigma[cc, cc] - data[cc, cc]))
##     d < eps
## }

## .does_edge_fit_to_data_mtp2 <- function(cc, parm, data, eps){
##     ##cat(".does_edge_fit_to_data_mtp2:  edge : ", toString(cc), "\n")
##     dif <- data[cc, cc] - parm$Sigma[cc, cc]
##     d1 <- max(dif)
##     d2 <- max(abs(dif) * parm$K[cc, cc])
##     out <- !(d1 > 0 || d2 > eps)
##     ## cat(sprintf(".does_edge_fit_to_data_mtp2 (R) d1=%15e, d2=%15e, out=%i\n",
##     ##             d1, d2, out))
##     out
## }

## .does_edge_fit_to_data_lasso <- function(cc, parm, data, eps){
##     cat("## NOT IMPLEMENTED\n")
## }

## .does_edge_fit_to_data_hybrid <- function(cc, parm, data, eps){
##     cat("## NOT IMPLEMENTED\n")
## }

## .does_edge_fit_to_data <- function(cc, parm, data, eps, method){
##     switch(method,
##            "ips"   ={.does_edge_fit_to_data_ips   (cc, parm, data, eps)},
##            "mtp2"  ={.does_edge_fit_to_data_mtp2  (cc, parm, data, eps)},
##            "lasso" ={.does_edge_fit_to_data_lasso (cc, parm, data, eps)},
##            "hybrid"={.does_edge_fit_to_data_hybrid(cc, parm, data, eps)})  
## }



