#!/usr/bin/env Rscript

library(rjson)

source("../math/adj_mat.r")

# input data
ex_num <- commandArgs(trailingOnly=TRUE)
# ex_num <- "8"
filename <- paste("examples/ex", ex_num, ".json", sep="")
input <- fromJSON(file = filename)
invisible(list2env(input, .GlobalEnv))

improve_flow <- function(s, t, d, x=NULL) {
    Ic <- 1 # iteration counter
    It <- 1 # tag counter
    L <- c(s) # tagged nodes
    g <- rep(-1, n) # first tags
    g[s] <- 0
    p <- rep(-1, n) # second tags
    p[s] <- 1
    i <- s # current

    while (TRUE) {
        to <- which(d[i,] > 0)
        to <- setdiff(to, L)
        for (j in to) {
            if (x[i,j] < d[i,j]) {
                g[j] <- i
                It <- It + 1
                p[j] <- It
                L <- c(L, j)
            }
        }

        from <- which(d[,i] > 0)
        from <- setdiff(from, L)
        for (j in from) {
            if (x[j,i] > 0) {
                g[j] <- -i
                It <- It + 1
                p[j] <- It
                L <- c(L, j)
            }
        }

        if (t %in% L)
            break

        Ic <- Ic + 1
        j0 <- which(p == Ic)
        if (length(j0) == 0) {
            return(NULL)
        }

        i <- j0
    }

    v <- 0 # flow

    i1 <- g[t]
    a <- d[i1,to] - x[i1,to]
    changed <- c(t)
    while (i1 != s) {
        i2 <- g[i1]

        if (i2 > 0) {
            a <- min(a, d[i2,i1] - x[i2,i1])
        } else if (i2 < 0) {
            i2 <- abs(i2)
            a <- min(a, x[i1,i2])
        } else {
            stop("Impossible?")
        }

        changed <- c(changed, i1)
        i1 <- i2
    }

    invisible(lapply(changed, function(to) {
        s <- sign(g[to])
        from <- abs(g[to])
        # print(paste(from,to,s))
        if (s > 0) {
            x[from,to] <<- x[from,to] + a*s
        } else {
            x[to,from] <<- x[to,from] + a*s
        }
    }))

    return(x)
}

next_x <- matrix(0, nrow=n, ncol=n) # flows
x <- NULL
d <- adj_mat(n, g_from, g_to, g_cap)
while (!is.null(next_x)) {
    x <- next_x
    next_x <- improve_flow(s, t, d, x)
}

print("Flow")
print(x)
print(paste("Max flow:", sum(x[,t])))
print(paste("Expected:", expected))
