#include "RcppArmadillo.h"
#include <iostream>
#include <vector>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::interfaces(r,cpp)]]

using namespace Rcpp;
using namespace arma;

// https://github.com/petewerner/misc/wiki/RcppArmadillo-cheatsheet

// //[[Rcpp::export]]
inline mat inv22_(mat A){
  double d = A(0,0) * A(1,1) - A(0,1)*A(1,0);
  mat B(2,2);
  B(0) = A(3);
  B(3) = A(0);
  B(1) = -A(1);
  B(2) = -A(2);
  return(B/d);
} 


// FIXME: Need not send A as argument.
// //[[Rcpp::export]]
arma::mat woodsolve2_(mat A, mat Ai, mat D, uvec cc){
  
  uvec cc0 = cc - 1;
  //Rcout << cc0 << std::endl;
  mat Di = inv22_(D);
  //Rcout << "Di" << std::endl << Di << std::endl;
  mat Ais = Ai.submat(cc0, cc0);
  //Rcout << "Ais" << std::endl << Ais << std::endl;
  mat H = inv22_(Ais + Di);
  //Rcout << "H" << std::endl << H << std::endl;
  return( Ai - Ai.cols(cc0) * H * Ai.rows(cc0) );    
}



/*** R

load_all("unifips")
di <- 4
A <- rWishart(1, 40, diag(1, di))[,,1]
A
Ai <- solve(A)
Ai

D <- matrix(c(1,.1, .1, 1), nr=2)
cc <- c(2,4)

solve(D)
Ai[cc, cc]
H <- solve(Ai[cc, cc] + solve(D))
H

Ai - Ai[, cc] %*% H %*% Ai[cc, ]    

r2 = woodsolve2(A, Ai, D, cc)
r2_ = woodsolve2_(A, Ai, D, cc)

r2 - r2_


library(microbenchmark)
di <- 1000
A2 <- rWishart(1, 4000, diag(1, di))[,,1]
A2i <- solve(A2)

microbenchmark(
    ## woodsolve0(A2, A2i, D, cc),
    ## woodsolve1(A2, A2i, D, cc),
    woodsolve2(A2, A2i, D, cc), 
    woodsolve2_(A2, A2i, D, cc), 
times=100)




***/



//    H <- sol(Ai[cc, cc] + sol(D))
//    Ai - Ai[, cc] %*% H %*% Ai[cc, ]    
