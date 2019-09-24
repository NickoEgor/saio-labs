is_whole <- function(x, eps) {
    return(apply(abs(matrix(c(x%%1, x%%1 - 1), nc = 2)), 1, FUN = min) < eps)
}

