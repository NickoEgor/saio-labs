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

dual_simplex <- function(A, b, c, d_d, d_u, eps = 1e-6, Jb = NULL) {

    rows <- nrow(A)
    cols <- ncol(A)

    # starting basis
    if (is.null(Jb))
        Jb <- guess_basis(A)
    Ab <- A[,Jb]

    print(paste("Jb", toString(Jb)))
    print("Ab")
    print(Ab)
    # print(paste("b:", toString(b)))
    # print(paste("c:", toString(c)))

    print("ALIVE")
    B <- inv_matrix(Ab, eps)
    print("DEAD")

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

        s <- 0
        for (i in Jn) {
            s <- s + A[,i]*nu[i]
        }

        nu[Jb] <- B %*% (b - s)

        # step 3
        # print(paste("nu:", toString(nu)))
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
        # print(paste("Jb:", toString(Jb)))
        # print("Ab")
        # print(Ab)
        Jb[k] <- j_star
        Ab <- A[,Jb]

        # print(paste("new Jb:", toString(Jb)))
        # print("new Ab")
        # print(Ab)

        # print("A")
        # print(A)
        # print(paste("b:", toString(b)))
        # print(paste("c:", toString(c)))
        # print("ALIVE")
        B <- inv_matrix(Ab, eps)
        # print("DEAD")
        # print("______________________________________________")

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

