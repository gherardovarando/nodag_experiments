library(graph)
library(igraph)
library(pcalg)
dyn.load("../src/nodag.so")
source("util.R")

set.seed(2020)
p <- 10
n <- 1000
## true DAG:
gGtrue <- randomDAG(p, prob = 0.2)
gmG <- list(x = rmvDAG(n, gGtrue), g = gGtrue)
X <- gmG$x
true <- dag2cpdag(gGtrue)
true_amat <- as(true, "graphAM")@adjMat
plot(true)


res <- list()
for (lambda in c(0.1,0.2,0.3)){
  out <- .Fortran("NODAG", as.integer(p), as.double(cor(X)),
                              as.double(diag(p)), as.double(lambda), 
                              as.double(1e-5), as.double(0.5), 
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
plot(res$`0.3`$g)

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

plot(res$`0.2`$g, layout = layout)

igraph.to.tikz(graph_from_adjacency_matrix(true_amat),layout = layout, file = "examplep10.true.txt")
igraph.to.tikz(res$`0.2`$g,layout = layout, file = "examplep10.0.2.txt")

Aest <-   res$`0.2`$A %*% diag(sqrt(diag(cov(X))))
West <- t(diag(p) - diag(1/diag(Aest)) %*% Aest)

W <- wgtMatrix(gGtrue)

plot(W[W!=0], West[W!=0])
abline(0,1)
