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

# get_basis <- function(plan) {
    # Jb <- which(plan > 0)
    # if (length(Jb) < nrow(A)) {
        # row_diff <- nrow(A) - length(Jb)

        # Jn <- which(plan == 0)

        # combs <- permutations(length(Jn), row_diff, Jn)
        # for (i in 1:nrow(combs)) {
            # indeces <- sort(c(Jb, combs[i,]))
            # Ab <- A[,indeces]
            # if (abs(det(Ab)) > EPSILON) {
                # return(indeces)
            # }
        # }
        # stop("Can't be solved")
    # }
    # return(Jb)
# }

# input data
res <- fromJSON(file = "examples/ex3.json")

rows <- res$rows
cols <- res$cols
A <- matrix(res$A, nr = rows, nc = cols, byrow = TRUE)
b <- res$b
c <- res$c
J <- c(1:ncol(A))
mode <- res$mode

# print("A:")
# print(A)
# print(paste("b:", toString(b)))
# print(paste("c:", toString(c)))

mism <- cols - rows

iter <- 0
while (TRUE) {
    iter <- iter + 1
    print(paste("Iter:", iter))

    # print(paste("A rows:", nrow(A), ", cols:", ncol(A)))
    # print(paste("b len:", length(b)))
    # print(paste("c len:", length(c)))

    # if (!all.equal(dim(A), c(length(b), length(c))) == TRUE) {
        # print(dim(A))
        # print(length(b))
        # print(length(c))
        # stop("The ERROR")
    # }

    # print("A:")
    # print(A)
    # print(paste("b:", toString(b)))
    # print(paste("c:", toString(c)))

    # step 1
    # lprec <- make.lp(0, ncol=ncol(A))
    # set.objfn(lprec, obj = c)
    # for (row_idx in 1:nrow(A)) {
        # add.constraint(lprec, xt = A[row_idx,], type = "=", rhs = b[row_idx])
    # }
    # lp.control(lprec,
               # simplextype = "dual",
               # pivoting = "firstindex",
               # verbose = "important",
               # sense = mode)
    # solve(lprec)
    # plan <- get.variables(lprec)
    # plan <- round(plan, PRECISION)
    # write.lp(lprec, filename = paste("test", iter, ".lp", sep = ""))

    result <- dual_simplex(A, b, c, rep(0, ncol(A)), rep(1e8, ncol(A)), eps = EPSILON)
    if (!result$solved) {
        stop("Can't solve task")
    }
    plan <- result$plan

    # print(paste("basis", toString(get.basis(lprec))))
    # print(paste("non-basic basis", toString(get.basis(lprec))))
    # print(paste("guess", toString(guess.basis(lprec, plan))))
    # print(get.primal.solution(lprec))

    J_all <- c(1:length(plan))
    # Jb <- get_basis(plan)
    Jb <- result$basis
    print(paste("Jb:", toString(Jb)))
    # Jb <- which(plan > EPSILON)
    J_art <- setdiff(J_all, J)

    print(paste("Plan:", toString(plan)))
    # print(paste("Count:", toString(length(plan))))
    # print(paste("Jb:", toString(Jb)))
    # print(paste("J_art:", toString(J_art)))
    # print("A")
    # print(A)

    # step 2
    common <- intersect(Jb, J_art)

    while (length(common) > 0) {
        extra <- common[1]
        idx <- which.min(extra %in% Jb)

        shifted_extra <- extra - mism

        row <- A[shifted_extra,]
        row <- row[-extra]
        col <- A[,extra]
        col <- col[-shifted_extra]
        val <- b[shifted_extra]

        coeff <- A[shifted_extra,extra]
        if (abs(coeff) - 1 > EPSILON) {
            print(paste("wrong coeff:", coeff))
            stop()
        }

        A <- A[-shifted_extra,] # delete row
        A <- A[,-extra]         # delete column
        b <- b[-shifted_extra]  # delete row
        c <- c[-extra]          # delete column
        plan <- plan[-extra]    # delete column

        J_all <- c(1:length(plan))
        # Jb <- get_basis(plan)
        Jb <- Jb[-idx]
        Jb_len <- length(Jb)
        Jb <- Jb[c(idx:Jb_len)] + 1
        print(paste("Jb:", toString(Jb)))
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
    # print("Ab")
    # print(Ab)
    if (nrow(Ab) == 1 && ncol(Ab) == 1)
        B <- ginv(Ab)
    else
        print(paste("dims", toString(dim(Ab))))
        B <- inv(Ab)
    # print("Ab^-1")
    # print(B)

    fract <- is_whole(plan[J_all], eps = EPSILON) == FALSE
    # print(paste("fract:", toString(fract)))
    non_art_basis <- J_all %in% intersect(J, Jb)
    # print(paste("non_art_basis:", toString(non_art_basis)))

    js <- fract & non_art_basis
    # print(paste("js", toString(js)))

    i0 <- which.max(js)
    # print(paste("i0", toString(i0)))

    e0 <- rep(0, nrow(B))
    e0[i0] <- 1
    # print(paste("ed:", toString(e0)))
    y <- e0 %*% B
    # print(paste("y:", toString(y)))
    alpha <- y %*% A
    beta <- y %*% b
    fa <- c(alpha %% 1, -1)
    fa[Jb] <- 0
    fb <- beta %% 1

    # print(paste("Alpha: [", toString(alpha), "]"))
    # print(paste("New cond fa: [", toString(-fa), "]"))
    # print(paste("New cond fb: [", -fb, "]"))

    Jn <- setdiff(J_all, Jb)
    # print(paste("Jn", toString(Jn)))
    if (all(fa[Jn] < EPSILON)) {
        print("No whole plans")
        stop()
    }

    # step 5
    A <- cbind(A, rep.int(0, nrow(A)))
    A <- rbind(A, -fa)
    b <- append(b, -fb)
    c <- append(c, 0)

    # print("_________________________________________________________________________________________________________")
}

