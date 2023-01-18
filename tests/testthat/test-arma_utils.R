context("gRips - arma utils")

##############################


nr <- 6
nc <- 8
SS <- matrix(rnorm(nr * nc), nrow=nr)

rr <- c(2, 3)
cc <- c(4, 5)

val <- c(999, 9999)

extract_uv <- function(M, u, v){
    if (is.na(u[1])) u <- TRUE
    if (is.na(v[1])) v <- TRUE    
    M[u, v]
}

replace_uv <- function(M, u, v, value){
    if (is.na(u[1])) u <- TRUE
    if (is.na(v[1])) v <- TRUE    
    M[u, v] <- value
    M
}

test_that("arma_utils", {
expect_equal(extract_uv_(SS,  rr,  cc),       extract_uv(SS,  rr,  cc))
expect_equal(extract_uv_(SS,  rr, -cc),       extract_uv(SS,  rr, -cc))
expect_equal(extract_uv_(SS, -rr,  cc),       extract_uv(SS, -rr,  cc))
expect_equal(extract_uv_(SS, -rr, -cc),       extract_uv(SS, -rr, -cc))
expect_equal(extract_uv_(SS,  NA,  NA),       extract_uv(SS,  NA,  NA))
expect_equal(extract_uv_(SS,  rr,  NA),       extract_uv(SS,  rr,  NA))
expect_equal(extract_uv_(SS,  NA,  cc),       extract_uv(SS,  NA,  cc))
expect_equal(extract_uv_(SS, -rr,  NA),       extract_uv(SS, -rr,  NA))
expect_equal(extract_uv_(SS,  NA, -cc),       extract_uv(SS,  NA, -cc))

expect_equal(replace_uv_(SS,  rr,  cc, val),  replace_uv(SS,  rr,  cc, val))
expect_equal(replace_uv_(SS,  rr, -cc, val),  replace_uv(SS,  rr, -cc, val))
expect_equal(replace_uv_(SS, -rr,  cc, val),  replace_uv(SS, -rr,  cc, val))
expect_equal(replace_uv_(SS, -rr, -cc, val),  replace_uv(SS, -rr, -cc, val))
expect_equal(replace_uv_(SS,  NA,  NA, val),  replace_uv(SS,  NA,  NA, val))
expect_equal(replace_uv_(SS,  rr,  NA, val),  replace_uv(SS,  rr,  NA, val))
expect_equal(replace_uv_(SS,  NA,  cc, val),  replace_uv(SS,  NA,  cc, val))
expect_equal(replace_uv_(SS, -rr,  NA, val),  replace_uv(SS, -rr,  NA, val))
expect_equal(replace_uv_(SS,  NA, -cc, val),  replace_uv(SS,  NA, -cc, val))
})


