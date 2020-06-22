library(graph)
library(igraph)
library(pcalg)
library(Matrix)
library(graph)
dyn.load("nodag.so")
source("util.R")

set.seed(2020)
p <- 10
n <- 1000
true <- graphNEL(nodes = paste0(1:10), edgeL = list(
  "1" = list(edges = c(5, 6), weights=runif(2, min = 0.5,1)),
  "2" = list(),
  "3" = list(),
  "4" = list(edges = 5, weights= runif(1, 0.5, 1)),
  "5" = list(edges = c(9, 10), weights=runif(2, 0.5, 1)), 
  "6" = list(edges = 10, weights=runif(1, 0.5, 1)),
  "7" = list(edges = 8, weights= runif(1, 0.5, 1)),
  "8" = list(),
  "9" = list(edges = 10, weights=runif(1, 0.5, 1)),
  "10" = list()
), edgemode = "directed")
X <- rmvDAG(n, true)
true_amat <- as(true, "graphAM")@adjMat
plot(true)


res <- list()
for (lambda in c(0.1,0.2,0.3)){
  out <- .Fortran("NODAG", as.integer(p), as.double(cor(X)),
                              as.double(diag(p)), as.double(lambda), 
                              as.double(1e-8), as.double(0.5), 
                              as.integer(1000))
  A <- matrix(nrow =p, out[[3]])
  g <- graph_from_adjacency_matrix(sign(abs(A)), diag = FALSE)
  res[[paste0(lambda)]] <- list(
    A = A,
    lambda = lambda,
    shd = pcalg::shd((as_graphnel(g)), true),
    g = g
  )
}

sapply(res, function(x) x$shd)
plot(res$`0.2`$g)

plot(true)
plot(as_graphnel(res$`0.2`$g))

layout <- matrix(ncol = 2, byrow = TRUE, 
          c(210.0000, 360.0000,
            290.0000, 360.0000,
            380.0000, 360.0000,
             50.0000, 320.0000,
            100.0000, 250.0000,
            250.0000, 170.0000,
            500.0000, 350.0000,
            500.0000, 150.0000,
             10.0000, 170.0000,
            100.0000,   0.0000))

igraph.to.tikz(graph_from_adjacency_matrix(true_amat),layout = layout, file = "examplep10.true.txt")
igraph.to.tikz(res$`0.2`$g,layout = layout, file = "examplep10.0.2.txt")
