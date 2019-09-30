#!/usr/bin/env Rscript

library(rjson)

source("../math/dual_simplex.r")

# input data
res <- fromJSON(file = "examples/ex0.json")

rows <- res$rows
cols <- res$cols
A <- matrix(res$A, nr = rows, nc = cols, byrow = TRUE)
b <- res$b
c <- res$c
d_d <- res$d_d
d_u <- res$d_u

result <- dual_simplex(A, b, c, d_d, d_u)

if (result$solved) {
    print("Iterations:")
    print(result$iterations)
    print("Plan:")
    print(result$plan)
    print("Value:")
    print((c %*% result$plan)[1])
} else {
    print("No solutions")
}
