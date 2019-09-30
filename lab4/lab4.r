#!/usr/bin/env Rscript

library(rjson)

# input data
ex_num <- commandArgs(trailingOnly=TRUE)
filename <- paste("examples/ex", ex_num, ".json", sep="")
res <- fromJSON(file = filename)

c <- res$c
n <- res$n
F <- matrix(res$F, nr = n, nc = c + 1, byrow = TRUE)

x <- matrix(0, nr = n, nc = c + 1)
B <- matrix(0, nr = n + 1, nc = c + 1)

# step 1 - forward bellman count
for (k in c(1:n)) {
    f <- F[k,]
    b <- B[k,]

    for (i in c(1:c+1)) {
        s <- f[c(1:i)] + b[c(i:1)]
        B[k+1, i] <- max(s)
        x[k, i] <- which.max(s) - 1
    }
}

# step 2 - backward plan count
plan <- rep(NA, n)
left <- c

profit_idx <- which.max(B[n+1,])
plan[n] <- x[n, profit_idx]
left <- left - plan[n]

for (k in c((n-1):1)) {
    spent <- x[k, left+1]
    left <- left - spent
    plan[k] <- spent
}
print(paste("Plan:", toString(plan)))

profit <- max(b)
print(paste("Profit:", profit))
