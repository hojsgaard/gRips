

get_bd <- function(A, u, amat){
    ## str(list(A=A, u=u, amat=amat))
    amat2 <- amat[A, A, drop=FALSE]
    which(amat2[u, ] > 0)
}


update_ <- function(A, u, b) {
    r <- setdiff(1:nrow(A), c(u, b))
    
    Arb <- A[r, b, drop=FALSE]
    Abb <- A[b, b, drop=FALSE]
    Abu <- A[b, u, drop=FALSE]
    
    Aru = Arb %*% solve(Abb, Abu)
    A[r, u] <- Aru
    A[u, r] <- t.default(Aru)
    return(A)
}




small_last <- function(amat, f) {    
  slo <- smallest_last_ordering(amat)
  nms <- rownames(amat)
  slo_inv <- match(1:nrow(amat), slo)
  amat2 <- amat[slo, slo]
  
  d <- nrow(amat)
  A <- 1:d
  b <- vector("list", d-1)
  b_nms <- vector("list", d-1)
  
  for (i in d:1) {
    bi <- intersect(get_bd(A, i, amat2), A)
    deg <- length(bi)
 ##   if (f <= deg){
 ##     stop("MLE may not exist\n") SLL: redundant naar coreness er checket
 ##   }
    b[[i]] <- bi
    b_nms[[i]] <- nms[slo][bi]
    A <- A[-i]
  }
  names(b) <- nms[slo]
  names(b_nms) <- nms[slo]
  return(list(slo=slo, slo_nms=nms[slo],
              slo_inv=slo_inv,
              slo_inv_nms=nms[slo_inv],
              b_nms=b_nms,
              b=b))
}



smallest_last_ordering <- function(adj_matrix) {
  # Calculate the degrees of each node (number of neighbors)
  degrees <- rowSums(adj_matrix)
  
  ## Initialize an empty vector to store the smallest-last ordering
  ordering <- numeric(length = ncol(adj_matrix))
  passive <- numeric(length = ncol(adj_matrix))
  
  for (i in 1:length(ordering)) {
    ## Find the node with the smallest degree
    degrees
    smallest_degree_node <- which(degrees == min(degrees))[1]
    smallest_degree_node
    
    ## Update the ordering and remove the node and its edges
    ordering[i] <- smallest_degree_node
    adj_matrix[smallest_degree_node,] <- 0
    adj_matrix[,smallest_degree_node] <- 0
    degrees <- rowSums(adj_matrix)
    degrees[smallest_degree_node] <- Inf  # Set the degree of the removed node to Inf
    passive[smallest_degree_node] <- Inf
    degrees <- degrees + passive
    ## str(list(passive=passive, degrees=degrees, ordering=ordering))
    
  }
  
  return(rev(ordering))
}

smallest_first_ordering <- function(adj_matrix) {
  # Calculate the degrees of each node (number of neighbors)
  degrees <- rowSums(adj_matrix)

  ## Initialize an empty vector to store the smallest-last ordering
  ordering <- numeric(length = ncol(adj_matrix))
  passive <- numeric(length = ncol(adj_matrix))

  for (i in 1:length(ordering)) {
    ## Find the node with the smallest degree
    degrees
    smallest_degree_node <- which(degrees == min(degrees))[1]
    smallest_degree_node

    ## Update the ordering and remove the node and its edges
    ordering[i] <- smallest_degree_node
    adj_matrix[smallest_degree_node,] <- 0
    adj_matrix[,smallest_degree_node] <- 0
    degrees <- rowSums(adj_matrix)
    degrees[smallest_degree_node] <- Inf  # Set the degree of the removed node to Inf
    passive[smallest_degree_node] <- Inf
    degrees <- degrees + passive
    ## str(list(passive=passive, degrees=degrees, ordering=ordering))

  }

  return(ordering)
}



small_first <- function(amat) {
  sfo <- smallest_first_ordering(amat)
  nms <- rownames(amat)
  sfo_inv <- match(1:nrow(amat), sfo)

  return(list(sfo=sfo, sfo_nms=nms[sfo],
              sfo_inv=sfo_inv,
              sfo_inv_nms=nms[sfo_inv]))
}


smart_start <- function(S, amat) {
  d <- nrow(amat)
  tmp <- small_first(amat)
  S2 <- S[tmp$sfo, tmp$sfo]
  amat2 <- amat[tmp$sfo, tmp$sfo]
  for (i in 1:d) {
    b <- get_bd(1:d, i, amat2)
    S2 <- update_(S2, i, b)
  }
  return(S2[tmp$sfo_inv, tmp$sfo_inv])
}




symsweep <- function (V, b, eps = 1e-8)
{
  out <- matrix(0, nrow(V), nrow(V))
  dimnames(out) <- dimnames(V)
  Saa <- V[-b,-b, drop = FALSE]
  Sab <- V[-b, b, drop = FALSE]
  Sba <- V[b,-b, drop = FALSE]
  Sbb <- V[b, b, drop = TRUE]
  
  
  if (abs(Sbb) > eps) {
    Sbbi = -1 / Sbb
    
    B <- Sba * Sbbi
    out[-b,-b] <- Saa + Sab %*% B
    
    C <- Sba * abs(Sbbi)
    out[-b, b] <- t(C)
    out[b,-b] <- C
    out[b, b] <- Sbbi
    return(out)
  } else
    return(V)
}



push <- function(Sigma, b, f,eps=1e-8) {
  d <- nrow(Sigma)
  Sigma2 <- Sigma
  for (i in d:(f+2)){
    tmp <- symsweep(Sigma2[1:i, 1:i, drop=FALSE], i,eps=eps)
    Sigma2[1:i, 1:i] <- tmp
    tmp2 <- update_(Sigma2[1:(i-1), 1:(i-1)], (i-1), b[[i-1]]) ## A, u, b
    Sigma2[1:(i-1), 1:(i-1)] <- tmp2
  }
  return(Sigma2)
}


pull <- function(Sigma, f,eps=1e-8) {
  d <- nrow(Sigma)
  Sigma2 <- Sigma    
  for (i in (f+2):d) {
    tmp <- symsweep(Sigma2[1:i, 1:i, drop=FALSE], i,eps=eps)
    Sigma2[1:i, 1:i] <- tmp        
  }
  return(Sigma2)    
}


sweep_start <- function(S, amat, f,eps=1e-8) {
  d <- nrow(amat)
  tmp <- small_last(amat, f=f)
  S2 <- S[tmp$slo, tmp$slo]
  Sigma2 <- update_(S2, d, tmp$b[[d]])
  #Sigma  <- Sigma2[tmp$slo_inv, tmp$slo_inv]
  if (f < d - 1) {
    Sigma2 <- push(Sigma2, tmp$b, f,eps=eps)
    Sigma2 <- pull(Sigma2,  f, eps=eps)
  }
  return(Sigma2[tmp$slo_inv,tmp$slo_inv])
}


auto_start <- function(S, good, amat) {
  d <- nrow(amat)
  Sigma <- S
  
  for (i in good)
  {
    b <- get_bd(1:d, i, amat)
    
    Sigma <- update_(Sigma, i, b)
  }
  return(Sigma)
}



