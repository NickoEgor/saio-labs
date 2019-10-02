library(matlib) # inv
library(MASS)   # ginv

inv_matrix <- function(A, eps) {
    if (abs(det(A)) < eps) {
        stop("det(Ab) < eps")
    }
    if (nrow(A) == 1 && ncol(A) == 1) {
        return(ginv(A))
    }
    else {
        # print(paste("dims", toString(dim(A))))
        return(inv(A))
    }
}

