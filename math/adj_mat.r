adj_mat <- function(n, from, to, weights=NULL) {
    n_edges <- length(from)

    if (length(to) != n_edges)
        stop("lengths of arrays not equal")
    if (!is.null(weights) && length(weights) != n_edges)
        stop("lengths of arrays not equal")

    if (is.null(weights)) {
        m <- matrix(FALSE, nrow=n, ncol=n)
        g_weights <- rep(TRUE, n)
    } else {
        m <- matrix(0, nrow=n, ncol=n)
        g_weights <- weights
    }

    for (i in c(1:n_edges))
        m[from[i], to[i]] <- g_weights[i]

    return(m)
}
