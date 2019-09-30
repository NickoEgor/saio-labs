# Three way in R to solve LP task

Given:
+ `rows, cols - integer`
+ `A - matrix[rows][cols]`
+ `b - vector[rows]`
+ `c - vector[cols]`

## library(boot)

```R
result <- simplex(c, b3 = b, A3 = A, maxi = TRUE)
if (result$solved != 1) {
    stop()
}
plan <- result$soln
```

## library(linprog)

```R
result <- solveLP(c, b, A, maximum = TRUE)
if (result$status != 0) {
    stop()
}
plan <- result$solution
```

## library(lpSolveAPI)

```R
lprec <- make.lp(0, ncol=ncol(A))
set.objfn(lprec, obj = c)

for (row_idx in 1:nrow(A)) {
    add.constraint(lprec, xt = A[row_idx,], type = "=", rhs = b[row_idx])
}
lp.control(lprec,
           simplextype = "primal",
           # simplextype = "dual",
           pivoting = "firstindex",
           verbose = "important",
           sense = "max")
solve(lprec)
plan <- get.variables(lprec)
```
