

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







reorder <- function(S, amat) {
  tmp <- small_first(amat)
  S2 <- S[tmp$sfo, tmp$sfo]
  amat2 <- amat[tmp$sfo, tmp$sfo]
  sfo_inv <- match(1:nrow(amat), tmp$sfo)
  
  return(list(S2=S2,amat2=amat2, sfo=tmp$sfo, sfo_inv=sfo_inv))
}

