############# estimate func
library(igraph)
library(Matrix)
library(pcalg)
library(bnlearn)
library(reticulate)
reticulate::use_python("/usr/lib/python3")
notears <- import("notears")
dyn.load("nodag.so")

est_pcalg <- function(x, ...) {
  time <- system.time({
    suffStat <- list(C = cor(x), n = nrow(x))
    pcout <-
      pc(suffStat = suffStat,
         indepTest = gaussCItest,
         #p = ncol(x),
         maj.rule = TRUE,
         solve.confl = TRUE,
         labels = colnames(x),
         ...)
  })[3]
  return(list(time = time, graph = pcout@graph))
}


est_nodag <- function(x, ...) {
  arg <- list(...)
  lambda <- arg$lambda
  toll <- ifelse(is.null(arg$toll), 1e-3, arg$toll)
  maxitr <- ifelse(is.null(arg$maxitr), 100, arg$maxitr)
  time <-
    system.time(
      out <- .Fortran(
        "NODAG",
        as.integer(ncol(x)),
        as.double(cor(x)),
        as.double(diag(ncol(x))),
        as.double(lambda),
        as.double(toll),
        as.double(0.5),
        as.integer(maxitr)
      )
    )[3]
  A <- matrix(nrow = ncol(x), out[[3]], dimnames = list(colnames(x), colnames(x)))
  g <- graph_from_adjacency_matrix(sign(abs(A)), diag = FALSE)
  return(list(
    time = time,
    graph = as_graphnel(g),
    coef = A,
    itr = out[[7]]
  ))
}

est_tabu <- function(x, ...) {
  time <- system.time(bn <- tabu(x, maxp = 10, ...))[3]
  return(list(time = time, graph  = as.graphNEL(bn)))
}

est_ges <- function(x, ...) {
  time <- system.time({
    score <- new("GaussL0penObsScore", x)
    ges.fit <- ges(score)
  })[3]
  return(list(
    time = time,
    graph = as(ges.fit$essgraph, "graphNEL")
  ))
}


est_chowliu <- function(x, ...) {
  time <- system.time({
    res <- bnlearn::chow.liu(as.data.frame(x))
  })
  return(list(time = time, graph  = as.graphNEL(res)))
}

est_lingam <- function(x, ...) {
  time <- system.time({
    res <- pcalg::lingam(x)
  })
  g <- as(abs(t(res$Bpruned)), "graphNEL")
  return(list(
    time = time,
    graph  = g,
    coef = res$Bpruned
  ))
}

est_notears <- function(x,  ...) {
  arg <- list(...)
  lambda <- arg$lambda
  w_threshold <- ifelse(is.null(arg$w_threshold), 0.3, arg$w_threshold)
  max_iter <- ifelse(is.null(arg$max_iter), 100, arg$max_iter)
  time <- system.time({
    res <- notears$notears_linear(X = as.matrix(x), lambda1 = as.double(lambda), loss_type = "l2", w_threshold = as.double(w_threshold),
                                  max_iter = as.integer(max_iter), h_tol = 1e-4, rho_max = 10)
  })
  A <- matrix(nrow = ncol(x), res, dimnames = list(colnames(x), colnames(x)))
  g <- graph_from_adjacency_matrix(sign(abs(A)), diag = FALSE)
  return(list(
    time = time,
    graph = as_graphnel(g),
    coef = res
  ))
}