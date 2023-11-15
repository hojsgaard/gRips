

get_bd <- function(A, u, amat){
    ## str(list(A=A, u=u, amat=amat))
    amat2 <- amat[A, A, drop=FALSE]
    which(amat2[u, ] > 0)
}


update_ <- function(A, u, b) {
    r <- setdiff(1:nrow(A), c(u, b))
   ##  print(list(r=r, u=u, b=b)) 
    Arb <- A[r, b, drop=FALSE]
    Abb <- A[b, b, drop=FALSE]
    Abu <- A[b, u, drop=FALSE]
    ## print(Abb)
    Aru = Arb %*% solve(Abb, Abu)
    A[r, u] <- Aru
    A[u, r] <- t.default(Aru)
    return(A)
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








