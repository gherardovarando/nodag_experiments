args <- commandArgs(trailingOnly = TRUE)
out <- "TABLE.rds"
est_methods <- c(
  #"pc-0.1" ,
  #"pc-0.05" ,
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

# copied from the pcalg code
myshd <- function(m1,m2){
  shd <- 0
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  m1[ind] <- 0
  shd <- shd + length(ind)/2
  ind <- which(ds < 0)
  m1[ind] <- m2[ind]
  shd <- shd + length(ind)/2
  d <- abs(m1 - m2)
  shd + sum((d + t(d)) > 0)/2
}

datapath <- "simulations/data"
gtpath <- "simulations/gt"
respath <- "simulations/results"

M <- 20
ks <- c(1,2,3,4)

#### sizes
ps <- c(5, 10, 20, 50, 100) 
ns <- c(100, 200, 500, 1000)

#### generation methods 
gen_methods <- c(
  "gmat_mh_u",
  "randomDAG_gaus",
  "randomDAG_exp",
  "randomDAG_gumb",
  "randomDAG_gaus_2"
)

TABLE <- array(dim = c(15,
                       length(gen_methods), 
                       length(est_methods) + 1, 
                       length(ns),
                       length(ps),
                       length(ks),
                       M), 
                       dimnames = list(
                         stat = c("sid", "sid-lower", "sid-upper", 
				  "is-dag", "is-cpdag", "is-pdag", 
				  "shd-dag", "shd-cpdag", "time", "tpr", "tnr", "fpr", 
				  "f1", "fdr", "iter"),
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
        gtcpdag <- pcalg::dag2cpdag(gt)
        gtamatcpdag <- as(gtcpdag, "graphAM")@adjMat
        vn <- colnames(gtamat)
        gtsk <- sign(gtamat + t(gtamat))
        for (est in est_methods){
          if (file.exists(paste0(resbasepath, "/", est,"_",filename, ".rds"))){
            message(filename)
            results <- readRDS(file = paste0(resbasepath, "/", est,"_",filename, ".rds"))
            amat <- as(results$graph, "graphAM")@adjMat
            ix <- sapply(vn, function(x) which(x == colnames(amat))) 
            amat <- amat[ix,ix] ### reorder amat as gtamat
            TABLE["time", gen, est, paste0(n), paste0(p), paste0(k), r] <- ifelse(is.null(results$time), NA, results$time)  
            tmp <- SID::structIntervDist(gtamat, amat)
            TABLE["sid", gen, est, paste0(n), paste0(p), paste0(k), r] <- tmp$sid
            TABLE["sid-lower", gen, est, paste0(n), paste0(p), paste0(k), r] <- tmp$sidLowerBound
            TABLE["sid-upper", gen, est, paste0(n), paste0(p), paste0(k), r] <- tmp$sidUpperBound
            TABLE["is-dag", gen, est, paste0(n), paste0(p), paste0(k), r] <- isValidGraph(amat, type = 'dag')
            TABLE["is-cpdag", gen, est, paste0(n), paste0(p), paste0(k), r] <- isValidGraph(amat, type = 'cpdag')
            TABLE["is-pdag", gen, est, paste0(n), paste0(p), paste0(k), r] <- isValidGraph(amat, type = 'pdag')
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
            TABLE["shd-dag", gen, est, paste0(n), paste0(p), paste0(k), r] <- myshd(amat, gtamat)
            TABLE["shd-cpdag", gen, est, paste0(n), paste0(p), paste0(k), r] <- myshd(amat, gtamatcpdag)
            TABLE["iter", gen, est,paste0(n), paste0(p), paste0(k), r] <- ifelse(is.null(results$itr), NA, results$itr)
          }
        }
        TABLE["shd-cpdag", gen, "empty",paste0(n), paste0(p), paste0(k), r] <- myshd(empty, gtamatcpdag)
        TABLE["shd-dag", gen, "empty",paste0(n), paste0(p), paste0(k), r] <- myshd(empty, gtamat)
      }
    }
  }
}
}
message("saving table")
saveRDS(TABLE, file = paste0("simulations/", out))
