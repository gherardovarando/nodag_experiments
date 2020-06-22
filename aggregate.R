args <- commandArgs(trailingOnly = TRUE)
out <- "TABLE.rds"
est_methods <- c(
  "pc-0.01" ,
  "pc-0.005" ,
  "pc-0.001" ,
  "nodag-0.3",
  "nodag-0.2" ,
  "nodag-0.1" ,
  "tabu",
  "ges"
)
if (length(args) > 0){
  out <- args[1]
}
if (length(args) > 1){
  est_methods <- args[-1]
}

library(pcalg)
library(Matrix)
library(igraph)

datapath <- "simulations/data"
gtpath <- "simulations/gt"
respath <- "simulations/results"

M <- 20
ks <- c(1,2,3,4)

#### sizes
ps <- c(5, 10, 20, 50, 100)
ns <- c(100, 1000, 10000)

#### generation methods 
gen_methods <- c(
  "randomDAG_gaus",
  "randomDAG_exp"
)

TABLE <- array(dim = c(9,
                       length(gen_methods), 
                       length(est_methods) + 1, 
                       length(ns),
                       length(ps),
                       length(ks),
                       M), 
                       dimnames = list(
                         stat = c("shd-dag", "shd-cpdag", "time", "tpr", "tnr", "fpr", "f1", "fdr", "iter"),
                         gen = gen_methods,
                         meth = c(est_methods, "empty"),
                         n = ns,
                         p = ps,
                         k = ks,
                         rep = 1:M
                       ))
for (n in ns){
for (p in ps){
  ## make the empty graph
  empty <- as_graphnel(make_empty_graph(p, TRUE))
  for (k in ks){
    databasepath <- paste(datapath,n, p, k, sep = "/")
    gtbasepath <- paste(gtpath,n, p, k, sep = "/")
    resbasepath <- paste(respath,n, p, k, sep = "/")
    for (r in 1:M){
      for (gen in gen_methods){
        filename <- paste(gen,n, p, k, r, sep = "_")
        gt <- readRDS(file = paste0(gtbasepath,"/dag_", filename, ".rds"))
        gt_cf <- read.table(file = paste0(gtbasepath,"/coeff_", filename))
        gtamat <- as(gt, "graphAM")@adjMat
        gtsk <- sign(gtamat + t(gtamat))
        for (est in est_methods){
          if (file.exists(paste0(resbasepath, "/", est,"_",filename, ".rds"))){
            message(filename)
            results <- readRDS(file = paste0(resbasepath, "/", est,"_",filename, ".rds"))
            TABLE["shd-dag", gen, est,paste0(n), paste0(p), paste0(k), r] <- pcalg::shd(results$graph, gt)
            TABLE["shd-cpdag", gen, est,paste0(n), paste0(p), paste0(k), r] <- pcalg::shd(results$graph, pcalg::dag2cpdag(gt))
            TABLE["time", gen, est,paste0(n), paste0(p), paste0(k), r] <- ifelse(is.null(results$time), NA, results$time)  
            amat <- as(results$graph, "graphAM")@adjMat
            sk <- sign(amat + t(amat))
            tpr <- sum(sk!=0 & gtsk!=0) / sum(gtsk!=0)
            tnr <- sum(sk==0 & gtsk==0) / sum(gtsk==0)
            fpr <- sum(sk!=0 & gtsk==0) / sum(gtsk==0)
            fdr <- sum(sk!=0 & gtsk==0) / sum(sk == 0)
            pr <-  sum(sk!=0 & gtsk!=0) / sum(sk!=0)
            if (is.nan(pr)) pr <- 1
            if (is.nan(fpr)) fpr <- 0
            if (is.nan(tpr)) tpr <- 1
            f1 <- 2 * (pr * tpr) / (pr + tpr)
            TABLE["tpr", gen, est,paste0(n), paste0(p), paste0(k), r] <- tpr
            TABLE["tnr", gen, est,paste0(n), paste0(p), paste0(k), r] <- tnr
            TABLE["fpr", gen, est,paste0(n), paste0(p), paste0(k), r] <- fpr
            TABLE["fdr", gen, est,paste0(n), paste0(p), paste0(k), r] <- fdr
            TABLE["f1", gen, est,paste0(n), paste0(p), paste0(k), r] <- f1
            TABLE["iter", gen, est,paste0(n), paste0(p), paste0(k), r] <- ifelse(is.null(results$itr), NA, results$itr)
          }
        }
        TABLE["shd-cpdag", gen, "empty",paste0(n), paste0(p), paste0(k), r] <- pcalg::shd(empty, gt)
        TABLE["shd-dag", gen, "empty",paste0(n), paste0(p), paste0(k), r] <- pcalg::shd(empty, gt)
        
      }
    }
  }
}
}

saveRDS(TABLE, file = paste0("simulations/", out))
