#!/usr/bin/env Rscript

library(matlib) # inv
library(gtools) # combinations

EPSILON <- 1e-6

dual_simplex <- function(A, b, c, d_d, d_u) {

    rows <- nrow(A)
    cols <- ncol(A)

    # starting basis
    J <- c(1:cols)
    combs <- combinations(cols, rows, J)
    for (i in 1:nrow(combs)) {
        Ab <- A[,combs[i,]]
        if (det(Ab) != 0) {
            Jb <- combs[i,]
            break
        }
    }

    B <- inv(Ab)

    # step 1
    y <- c[Jb] %*% B
    coplan <- y %*% A - c
    Jn <- setdiff(J, Jb)
    Jn_plus <- Jn[which(coplan[Jn] >= 0)]
    Jn_minus <- setdiff(Jn, Jn_plus)

    hasSolutions <- FALSE

    # iteration
    iter <- 0
    while (TRUE) {
        iter <- iter + 1
        # step 2
        nu <- rep(0, cols)
        nu[Jn_plus] <- d_d[Jn_plus]
        nu[Jn_minus] <- d_u[Jn_minus]

        s <- 0
        for (i in Jn) {
            s <- s + A[,i]*nu[i]
        }

        nu[Jb] <- B %*% (b - s)

        # step 3
        if (all(d_d[Jb] - EPSILON < nu[Jb]) && all(d_u[Jb] + EPSILON > nu[Jb])) {
            hasSolutions <- TRUE
            break
        }

        # step 4
        k <- which(nu[Jb] < d_d[Jb] - EPSILON | nu[Jb] > d_u[Jb] + EPSILON)[1]
        jk <- Jb[k]

        # step 5
        if (nu[jk] < d_d[jk])
            u <- 1
        else
            u <- -1

        e <- rep(0, length(Jb))
        e[k] <- 1
        deltaY <- u * e %*% B
        uv <- deltaY %*% A

        # step 6
        sigma <- rep(Inf, length(Jn))
        idx <- union(Jn_plus[which(uv[Jn_plus] < -EPSILON)],
                     Jn_minus[which(uv[Jn_minus] > EPSILON)])

        for (index in 1:length(Jn)) {
            if (Jn[index] %in% idx) {
                sigma[index] <- -coplan[Jn[index]]/uv[Jn[index]]
            }
        }

        sigma0 <- min(sigma)
        if (sigma0 == Inf) {
            break
        }
        j_star <- Jn[which.min(sigma)]

        # step 7
        coplan <- coplan + sigma0*uv

        # step 8
        Jb[k] <- j_star
        Ab <- A[,Jb]
        B <- inv(Ab)

        # step 9
        Jn = setdiff(J, Jb)

        if (j_star %in% Jn_plus) {
            if (u == 1) {
                Jn_plus <- union(setdiff(Jn_plus, j_star), jk)
            }
            else {
                Jn_plus <- setdiff(Jn_plus, j_star)
            }
        } else {
            if (u == 1) {
                Jn_plus <- union(Jn_plus, jk)
            }
        }

        Jn_minus = setdiff(Jn, Jn_plus)
    }

    return(list(hasSolutions = hasSolutions, iterations = iter, plan = nu))
}

