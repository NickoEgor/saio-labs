library(gtools) # combinations

# source("./inv_matrix.r")
source("../math/inv_matrix.r")

guess_basis <- function(A) {
    cols <- ncol(A)
    rows <- nrow(A)

    J <- c(1:cols)
    combs <- combinations(cols, rows, J)
    for (i in 1:nrow(combs)) {
        Ab <- A[,combs[i,]]
        if (det(Ab) != 0) {
            return(combs[i,])
        }
    }

    return(NULL)
}

dual_simplex <- function(A, b, c, d_d, d_u, eps = 1e-6, Jb = NULL, Jn_plus = NULL) {
    rows <- nrow(A)
    cols <- ncol(A)

    # starting basis
    if (is.null(Jb))
        Jb <- guess_basis(A)
    Ab <- A[,Jb]
    B <- inv_matrix(Ab, eps)

    # step 1
    y <- c[Jb] %*% B
    coplan <- y %*% A - c
    J_all <- c(1:cols)
    Jn <- setdiff(J_all, Jb)

    # print(paste("y", toString(y)))
    # print(paste("c", toString(c)))
    # print("A")
    # print(A)
    # print(paste("coplan:", toString(coplan)))

    # Jn_plus <- Jn[which(coplan[Jn] > eps)]
    Jn_plus <- Jn[which(coplan[Jn] >= -eps)]
    Jn_minus <- setdiff(Jn, Jn_plus)

    solved <- FALSE

    # iteration
    iter <- 0
    while (TRUE) {
        iter <- iter + 1

        # print(paste("INT iter:", iter))
        # print(paste("Jb", toString(Jb)))
        # print(paste("Jn", toString(Jn)))
        # print(paste("Jn+", toString(Jn_plus)))
        # print(paste("Jn-", toString(Jn_minus)))
        # print("--------")

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

        deltaY <- u * B[k,]
        # print(paste("delta:", toString(deltaY)))
        # print(paste("B[k,]:", toString(B[k,])))
        # uv <- deltaY %*% A
        uv <- deltaY %*% A

        # print(paste("uv:", toString(uv)))

        # step 6
        sigma <- rep(Inf, length(Jn))
        idx <- union(Jn_plus[which(uv[Jn_plus] < -eps)],
                     Jn_minus[which(uv[Jn_minus] > eps)])
        # print(paste("j's", toString(idx)))

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
        # print(paste("sigma", toString(sigma)))
        # print(paste("j_star", j_star))

        # step 7
        coplan <- coplan + sigma0*uv

        # step 8
        Jb[k] <- j_star
        Ab <- A[,Jb]

        B <- inv_matrix(Ab, eps)

        # step 9
        Jn = setdiff(J_all, Jb)

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


        # print(paste("Jb", toString(Jb)))
        # print(paste("Jn", toString(Jn)))
        # print(paste("Jn+", toString(Jn_plus)))
        # print(paste("Jn-", toString(Jn_minus)))
        # print("____________________________________________")
    }

        # print("++++++++ END ALL ++++++++")
        # print(paste("Jb", toString(Jb)))
        # print(paste("Jn", toString(Jn)))
        # print(paste("Jn+", toString(Jn_plus)))
        # print(paste("Jn-", toString(Jn_minus)))
        # print("+++++++++++++++++++++++++")

    return(list(solved = solved, iterations = iter, plan = nu, basis = Jb, nonbasis_plus = Jn_plus))
}

