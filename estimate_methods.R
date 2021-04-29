############# estimate func
library(igraph)
library(Matrix)
library(pcalg)
library(bnlearn)
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
