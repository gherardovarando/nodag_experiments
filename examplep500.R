library(pcalg)
library(igraph)
dyn.load("../src/nodag.so")

set.seed(2020)
p <- 1000
n <- 1000
## true DAG:
gGtrue <- randomDAG(p, prob =  2/p)
gmG <- list(x = rmvDAG(n, gGtrue), g = gGtrue)
X <- gmG$x
true <- dag2cpdag(gGtrue)
true_amat <- as(true, "graphAM")@adjMat

system.time(out <- .Fortran("NODAG", as.integer(p), as.double(cor(X)),
                as.double(diag(p)), as.double(0.2), 
                as.double(1e-5), as.double(0.5), 
                as.integer(1000)))
out[[7]]
A <- matrix(nrow =p, out[[3]])
g <- graph_from_adjacency_matrix(sign(abs(A)), diag = FALSE)

system.time({
  suffStat <- list(C = cor(X), n = nrow(X))
  pcout <- pc(suffStat = suffStat, indepTest = gaussCItest, 
              p = ncol(X), alpha = 0.001)
}
)

pcalg::shd((as_graphnel(g)), true)
pcalg::shd(pcout, true)
