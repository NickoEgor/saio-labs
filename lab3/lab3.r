#!/usr/bin/env Rscript

options(warn=1)

library(lpSolveAPI) # solve
library(matlib)     # inv
library(MASS)       # ginv
library(gtools)     # permutations
library(rjson)

source("../math/is_whole.r")
source("../math/dual_simplex.r")

PRECISION <- 6
EPSILON <- '^'(10,-PRECISION)
MAX_ITERS <- 1500

# input data
res <- fromJSON(file = "examples/ex7.json")

rows <- res$rows
cols <- res$cols
A <- matrix(res$A, nr = rows, nc = cols, byrow = TRUE)
b <- res$b
c <- res$c
J <- c(1:ncol(A))
mode <- res$mode

Jb <- guess_basis(A)
Jn_plus <- NULL

mism <- cols - rows

iter <- 0
while (TRUE) {
    iter <- iter + 1
    if (iter > MAX_ITERS)
        stop("Infinite loop")

    # step 1
    result <- dual_simplex(A, b, c, rep(0, ncol(A)), rep(1e8, ncol(A)), EPSILON, Jb, Jn_plus)
    if (!result$solved) {
        stop("Can't solve task")
    }
    plan <- result$plan

    J_all <- c(1:length(plan))
    Jb <- result$basis
    Jn_plus <- result$nonbasis_plus
    J_art <- setdiff(J_all, J)

    # step 2
    common <- intersect(Jb, J_art)

    while (length(common) > 0) {
        extra <- common[1]
        idx <- which.max(Jb %in% extra)

        shifted_extra <- extra - mism

        row <- A[shifted_extra,]
        row <- row[-extra]
        col <- A[,extra]
        col <- col[-shifted_extra]
        val <- b[shifted_extra]

        coeff <- A[shifted_extra,extra]
        if (abs(coeff) - 1 > EPSILON) {
            stop(paste("wrong coeff:", coeff))
        }

        A <- A[-shifted_extra,] # delete row
        A <- A[,-extra]         # delete column
        b <- b[-shifted_extra]  # delete row
        c <- c[-extra]          # delete column
        plan <- plan[-extra]    # delete column

        J_all <- c(1:length(plan))
        Jb <- Jb[-idx]
        for (i in 1:length(Jb)) {
            if (Jb[i] > extra) {
                Jb[i] <- Jb[i] - 1
            }
        }

        J_art <- setdiff(J_all, J)

        b <- b - col*val

        sub_rows <- matrix(rep(row,length(col)), nr = length(col), byrow = TRUE)
        sub_rows <- sub_rows * col
        A <- A - sub_rows

        common <- intersect(Jb, J_art)
    }

    # step 3
    if (all(is_whole(plan[J], EPSILON))) {
        plan <- plan[J]
        c <- c[J]

        print(paste("Iterations:", toString(iter)))
        print(paste("Plan:", toString(round(plan, PRECISION))))
        print(paste("Value:", toString(round((c %*% plan)[1], PRECISION))))
        break
    }

    # step 4
    Ab <- A[,Jb]
    B <- inv_matrix(Ab,EPSILON)

    fract <- !is_whole(plan[J_all], eps = EPSILON)
    non_art_basis <- J_all %in% intersect(J, Jb)

    js <- fract & non_art_basis
    i0 <- which.max(js)
    k <- which.max(Jb == i0)

    y <- B[k,]
    alpha <- y %*% A
    beta <- y %*% b
    fa <- c(alpha %% 1, -1)
    fa[Jb] <- 0
    fb <- beta %% 1

    Jn <- setdiff(J_all, Jb)
    if (all(fa[Jn] < EPSILON)) {
        stop("No whole plans")
    }

    # step 5
    A <- cbind(A, rep.int(0, nrow(A)))
    A <- rbind(A, -fa)
    b <- c(b, -fb)
    c <- c(c, 0)
    Jb <- c(Jb, ncol(A))
}

