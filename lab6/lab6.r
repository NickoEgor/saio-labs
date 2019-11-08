#!/usr/bin/env Rscript

library(rjson)

source("../math/adj_mat.r")

# input data
ex_num <- commandArgs(trailingOnly=TRUE)
# ex_num <- "1"
filename <- paste("examples/ex", ex_num, ".json", sep="")
input <- fromJSON(file = filename)
invisible(list2env(input, .GlobalEnv))

U <- c()
I <- c(s)
B <- rep(Inf, n)
B[s] <- 0
f <- rep(0, n)
f[s] <- s

n_edges <- length(g_from)
m <- adj_mat(n, g_from, g_to, g_weight)

cur <- s

while (TRUE) {
    to <- which(m[cur,] > 0)
    if (length(to) == 0)
        break

    for (i in to) {
        if (B[i] > B[cur] + m[cur, i]) {
            B[i] <- B[cur] + m[cur, i]
            f[i] <- cur
            if (!(i %in% I))
                I <- c(I, i)
        }
    }

    min_w <- which.min(B[to])
    if (min_w == Inf) {
        stop(paste("No paths from", s))
    }

    U <- c(U, cur)
    opts <- setdiff(I, U)
    if (length(opts) == 0) {
        break
    }

    cur <- opts[which.min(B[opts])]
}

print(paste("Bellman:", toString(B)))
print(paste("Expected:", toString(expected)))

# for checking
find_path <- function(t) {
    s <- 1
    path <- c(t)
    while (s != t) {
        t <- f[t]
        path <- c(t, path)
    }
    return(path)
}
