#!/usr/bin/env Rscript

library(rjson)

# input data
ex_num <- commandArgs(trailingOnly=TRUE)
# ex_num <- "8"
filename <- paste("examples/ex", ex_num, ".json", sep="")
input <- fromJSON(file = filename)
invisible(list2env(input, .GlobalEnv))

I <- c(s)
B <- rep(0, n)
f <- rep(0, n)

n_edges <- length(g_from)
m_w <- matrix(0, nrow=n, ncol=n)
m_1 <- matrix(0, nrow=n, ncol=n)
for (i in c(1:n_edges)) {
    m_w[g_from[i], g_to[i]] <- g_weight[i]
    m_1[g_from[i], g_to[i]] <- 1
}

usagesAll <- colSums(m_1)

process_new <- function(v) {
    indeces <- intersect(I, which(m_1[,v] == 1))
    costs <- sapply(indeces, function(i) B[i] + m_w[i, v])
    B[v] <<- max(costs)
    f[v] <<- indeces[which.max(costs)]
}

while (!(t %in% I)) {
    usages <- colSums(matrix(m_1[I,], nrow=length(I)))
    filled <- which(usages == usagesAll)
    new <- filled[which(!(filled %in% I))][1]

    I <- c(I, new)
    process_new(new)

    if (t %in% I) { break }
}

print(paste("I:", toString(I)))
print(paste("B:", toString(B)))
print(paste("f:", toString(f)))

path <- c(t)
i <- t
while (i != s) {
    i <- f[i]
    path <- c(i, path)
}

print(paste("path: ", toString(path)))

print(paste("max:", B[t]))
print(paste("expected:", expected))
