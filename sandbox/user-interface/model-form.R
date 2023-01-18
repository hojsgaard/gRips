library(unifips)
data(math)
S <- cov(math)

load_all("unifips")


do_data <- function(data){
    n <- nrow(data)
    data <- cov(data)
    attr(data, "n") <- n
    data
}

mod <- function(form, data){
    UseMethod("mod")
}

mod.list <- function(form, data){
    if (inherits(data, "data.frame")) data <- do_data(data)

    if (length(form) > 0){
        em <- do.call(rbind, lapply(form, function(g)
            if (length(g) > 1) t.default(combnPrim(g, 2))))
        em <- t.default(unique(em))
    } else
        em  <- matrix(nrow=2, ncol=0)

    str(list(em, data))
    mod.matrix(em, data)
}

mod.matrix <- function(form, data){
    if (inherits(data, "data.frame")) data <- do_data(data)

    vn <- unique.default(form)
    n <- attr(data, "n")
    ## Some checking needed
    S <- S[vn, vn]
    attr(S, "n") <- n
    unifips18_r(form, S)
}

gl  <- list(1:3, 2:4)
mod(gl, S)

gl <- matrix(c(1,2,2,3,1,3,2,4,3,4), nr=2)
mod(gl, S)

gl <- matrix(nc=0,nr=2)
mod(gl, S)

gl <- list()
mod(gl, S)






