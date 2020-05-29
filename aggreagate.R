args <- commandArgs(trailingOnly = TRUE)
out <- "TABLE.rds"
est_methods <- c(
  "pc-0.01" ,  
  "nodag-0.3",
  "nodag-0.2", 
  "nodag-0.1" ,
  "tabu",
  "ges"#,
#  "chowliu"
)
if (length(args) > 0){
  out <- args[1]
}
if (length(args) > 1){
  est_methods <- args[-1]
}

source("util.R")

datapath <- "simulations/data"
gtpath <- "simulations/gt"
respath <- "simulations/results"

M <- 20
ks <- c(1,2,3,4)

#### sizes
ps <- c(5, 10, 20, 50, 100, 200)
ns <- c(100, 1000, 10000)

#### generation methods 
gen_methods <- c(
   "randomDAG_gaus",
   "randomDAG_exp",
  "randomDAG_gaus_2"
)


TABLE <- array(dim = c(9,
                       length(gen_methods), 
                       length(est_methods), 
                       length(ns),
                       length(ps),
                       length(ks),
                       M), 
                       dimnames = list(
                         stat = c("shd simple", "shd", "time", "tpr", "tnr", "fpr", "f1", "fdr", "iter"),
                         gen = gen_methods,
                         meth = est_methods,
                         n = ns,
                         p = ps,
                         k = ks,
                         rep = 1:M
                       ))
for (n in ns){
for (p in ps){
  for (k in ks){
    databasepath <- paste(datapath,n, p, k, sep = "/")
    gtbasepath <- paste(gtpath,n, p, k, sep = "/")
    resbasepath <- paste(respath,n, p, k, sep = "/")
    for (r in 1:M){
      for (gen in gen_methods){
        filename <- paste(gen,n, p, k, r, sep = "_")
        gt <- readRDS(file = paste0(gtbasepath,"/cpdag_", filename, ".rds"))
        gt_cf <- read.table(file = paste0(gtbasepath,"/coeff_", filename))
        gtamat <- as(gt, "graphAM")@adjMat
        gtsk <- sign(gtamat + t(gtamat))
        for (est in est_methods){
          if (file.exists(paste0(resbasepath, "/", est,"_",filename, ".rds"))){
            results <- readRDS(file = paste0(resbasepath, "/", est,"_",filename, ".rds"))
            TABLE["shd", gen, est,paste0(n), paste0(p), paste0(k), r] <- pcalg::shd(results$graph, gt)
            TABLE["shd simple", gen, est,paste0(n), paste0(p), paste0(k), r] <- pcalg::shd(results$graph, gt)
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
      }
    }
  }
}
}

saveRDS(TABLE, file = paste0("simulations/", out))
