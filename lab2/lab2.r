#!/usr/bin/env Rscript

library(rjson)
library(collections)

source("../math/dual_simplex.r")
source("../math/is_whole.r")

# input data
res <- fromJSON(file = "examples/ex10.json")

rows <- res$rows
cols <- res$cols
A <- matrix(res$A, nr = rows, nc = cols, byrow = TRUE)
b <- res$b
c <- res$c
d_d <- res$d_d
d_u <- res$d_u

t <- 0
nu <- rep(0, cols)
nu0 <- 0
r0 <- -Inf

task0 <- list(lower = d_d, upper = d_u)

q <- Queue$new()
q$push(task0)

while (TRUE) {
    t <- t + 1

    # step 1
    if (q$size() == 0) {
        break
    }

    # step 2
    task <- q$pop()
    result <- dual_simplex(A, b, c, task$lower, task$upper)

    if (!result$solved || (c %*% result$plan)[1] <= r0) {
        next
    }

    # step 3
    check_whole <- is_whole(result$plan, EPSILON)
    if (all(check_whole)) {
        nu <- result$plan
        nu0 <- 1
        r0 <- (c %*% nu)[1]
        next
    }

    # step 4
    float_idx <- which(check_whole == FALSE)[1]

    float_lower <- floor(result$plan[float_idx])
    float_upper <- ceiling(result$plan[float_idx])

    u <- task$upper
    u[float_idx] <- float_lower

    l <- task$lower
    l[float_idx] <- float_upper

    q$push(list(lower = task$lower, upper = u))
    q$push(list(lower = l, upper = task$upper))
}


if (!nu0) {
    print("No solutions")
    stop()
}

print("Iterations:")
print(t)
print("Plan:")
print(nu)
print("Value:")
print((c %*% nu)[1])
