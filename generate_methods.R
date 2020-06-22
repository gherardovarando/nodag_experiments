library(pcalg)
library(igraph)

gen_gmat <- function(p, k, n, r) {
  dag <- gmat::rgraph(p, k / p, dag = TRUE, ordered = TRUE)
  L <- gmat:::mh_u(1, p = p, dag = dag)[, , 1]
  Sigma <- solve(tcrossprod(L))
  true <- as_graphnel(dag)
  X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigma)
  return(list(
    true = true,
    coef = L,
    x = X
  ))
}

gen_pcalg <- function(p, k, n, r, ...) {
  gGtrue <- randomDAG(p, prob = k / p, ...)
  x <- rmvDAG(n, gGtrue)
  true <- gGtrue
  return(list(
    true = true,
    coef = wgtMatrix(gGtrue),
    x = x
  ))
}

gen_pcalg_exp <- function(p, k, n, r) {
  gGtrue <- randomDAG(p, prob = k / p)
  errMat <- matrix(nrow = n,
                   ncol = p,
                   data = rexp(n * p))
  x <- rmvDAG(n, gGtrue, errMat = errMat)
  true <- gGtrue
  return(list(
    true = true,
    coef = wgtMatrix(gGtrue),
    x = x
  ))
}




methods <- list(
  "gmat_mh_u" = gen_gmat,
  "randomDAG_gaus" = gen_pcalg,
  "randomDAG_gaus_2" = function(p,k,n,r) gen_pcalg(p,k,n,r,lB=-1, uB = 1),
  "randomDAG_exp" = gen_pcalg_exp
)