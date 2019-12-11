#!/usr/bin/env Rscript

library(rjson)

# input data
# ex_num <- commandArgs(trailingOnly=TRUE)
ex_num <- "1"
filename <- paste("examples/ex", ex_num, ".json", sep="")
input <- fromJSON(file = filename)
invisible(list2env(input, .GlobalEnv))

solve <- function(mat) {
    n <- nrow(mat)
    for (k in seq(n)) {
        for (i in seq(n)) {
            for (j in seq(n)) {
                mat[i, j] <- min(mat[i, j], mat[i, k] + mat[k, j])
            }
        }
    }
    return(mat)
}

items <- as.numeric(unlist(fromJSON(items)))
for (k in seq(length(items))) {
    if (items[k] == "Inf") {
        items[k] <- Inf
    }
}

n <- sqrt(length(items))
mat <- matrix(nrow=n, ncol=n, data=items)
res <- solve(mat)

print(res)
