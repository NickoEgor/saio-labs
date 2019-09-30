#!/usr/bin/env Rscript

library(matlib) # inv
library(MASS)   # ginv
library(gtools) # combinations

dual_simplex <- function(A, b, c, d_d, d_u, eps = 1e-6) {

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

    if (det(Ab) < eps) {
        stop("det(Ab) = 0")
    }
    if (nrow(Ab) == 1 && ncol(Ab) == 1) {
        B <- ginv(Ab)
    }
    else {
        print(paste("dims", toString(dim(Ab))))
        B <- inv(Ab)
    }

    # step 1
    y <- c[Jb] %*% B
    coplan <- y %*% A - c
    Jn <- setdiff(J, Jb)
    Jn_plus <- Jn[which(coplan[Jn] >= 0)]
    Jn_minus <- setdiff(Jn, Jn_plus)

    solved <- FALSE

    # iteration
    iter <- 0
    while (TRUE) {
        iter <- iter + 1
        # step 2
        nu <- rep(0, cols)
        nu[Jn_plus] <- d_d[Jn_plus]
        nu[Jn_minus] <- d_u[Jn_minus]
        print(paste("J:", toString(J)))
        print(paste("d_d:", toString(d_d)))
        print(paste("Jn+:", toString(Jn_plus)))
        print(paste("d_u:", toString(d_u)))
        print(paste("Jn-:", toString(Jn_minus)))
        print(paste("nu1:", toString(nu)))

        s <- 0
        for (i in Jn) {
            s <- s + A[,i]*nu[i]
        }

        print(paste("s:", toString(s)))
        nu[Jb] <- B %*% (b - s)

        # step 3
        print(paste("nu2:", toString(nu)))
        print(paste("wtf1", all(d_d[Jb] - eps < nu[Jb])))
        print(paste("wtf2", all(d_u[Jb] + eps > nu[Jb])))
        if (all(d_d[Jb] - eps < nu[Jb]) && all(d_u[Jb] + eps > nu[Jb])) {
            solved <- TRUE
            break
        }

        # step 4
        k <- which(nu[Jb] < d_d[Jb] - eps | nu[Jb] > d_u[Jb] + eps)[1]
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
        idx <- union(Jn_plus[which(uv[Jn_plus] < -eps)],
                     Jn_minus[which(uv[Jn_minus] > eps)])

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
        if (det(Ab) < eps) {
            stop("det(Ab) = 0")
        }
        if (nrow(Ab) == 1 && ncol(Ab) == 1) {
            B <- ginv(Ab)
        }
        else {
            print(paste("dims", toString(dim(Ab))))
            B <- inv(Ab)
        }

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

    return(list(solved = solved, iterations = iter, plan = nu, basis = Jb))
}

