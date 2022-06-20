library(pcalg)
library(igraph)

gen_gmat <- function(p, k, n, r) {
  dag <- gmat::rgraph(p, k / p, dag = TRUE, ordered = TRUE)
  V(dag)$name <- make.names(1:p)
  L <- gmat:::mh_u(1, p = p, dag = dag, h = 500, eps = 0.01)[, , 1]
  Sigma <- solve(tcrossprod(L))
  true <- as_graphnel(dag)
  X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigma)
  colnames(X) <- V(dag)$name
  return(list(
    true = true,
    coef = L,
    x = X
  ))
}

gen_pcalg <- function(p, k, n, r) {
  gGtrue <- randomDAG(p, prob = k / p, V = make.names(1:p), lB = -1, uB = 1)
  x <- rmvDAG(n, gGtrue)
  x <- x
  true <- gGtrue
  return(list(
    true = true,
    coef = wgtMatrix(gGtrue),
    x = x
  ))
}

gen_pcalg_exp <- function(p, k, n, r) {
  gGtrue <- randomDAG(p, prob = k / p, V = make.names(1:p), lB = -1, uB = 1)
  errMat <- matrix(nrow = n,
                   ncol = p,
                   data = rexp(n * p))
  x <- rmvDAG(n, gGtrue, errMat = errMat)
  x <- x
  true <- gGtrue
  return(list(
    true = true,
    coef = wgtMatrix(gGtrue),
    x = x
  ))
}


gen_pcalg_gumb <- function(p, k, n, r) {
  gGtrue <- randomDAG(p, prob = k / p, V = make.names(1:p), lB = -1, uB = 1)
  errMat <- matrix(nrow = n,
                   ncol = p,
                   data = evd::rgumbel(n * p))
  x <- rmvDAG(n, gGtrue, errMat = errMat)
  x <- x
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
  "randomDAG_exp" = gen_pcalg_exp,
  "randomDAG_gumb" = gen_pcalg_gumb
)
