args <- commandArgs(trailingOnly = TRUE)
torun <- c(
  #"pc-0.1" ,
  #"pc-0.05" ,
  "pc-0.01" ,
  "pc-0.005" ,
  "pc-0.001" ,
  "nodag-0.3",
  "nodag-0.2" ,
  "nodag-0.1" ,
  "tabu",
  "ges",
  "notears-0",
  "notears-0.1",
  "notears-0.2",
  "notears-0.3"
)
if (length(args) > 0){
  torun <- args
}

source("estimate_methods.R")

datapath <- "simulations/data"
gtpath <- "simulations/gt"
respath <- "simulations/results"
dir.create(respath, showWarnings = FALSE, recursive = TRUE)

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
  "randomDAG_gumb"
)

est_methods <- list(
  "pc-0.1" = function(x) est_pcalg(x, alpha = 0.1),
  "pc-0.05" = function(x) est_pcalg(x, alpha = 0.05),
  "pc-0.01" = function(x) est_pcalg(x, alpha = 0.01),
  "pc-0.005" = function(x) est_pcalg(x, alpha = 0.005),
  "pc-0.001" = function(x) est_pcalg(x, alpha = 0.001),
  "nodag-0.3" = function(x) est_nodag(x, lambda = 0.3),
  "nodag-0.2" = function(x) est_nodag(x, lambda = 0.2),
  "nodag-0.1" = function(x) est_nodag(x, lambda = 0.1),
  "nodag-0.05" = function(x) est_nodag(x, lambda = 0.05),
  "notears-0.3" = function(x) est_notears(x, lambda = 0.3),
  "notears-0.2" = function(x) est_notears(x, lambda = 0.2),
  "notears-0.1" = function(x) est_notears(x, lambda = 0.1),
  "notears-0" = function(x) est_notears(x, lambda = 0),
  "tabu" = function(x) est_tabu(x),
  "ges" = function(x) est_ges(x),
  "chowliu" = function(x) est_chowliu(x),
  "lingam" = function(x) est_lingam(x)
)

for (n in ns){
  for (p in ps){
    for (k in ks){
      databasepath <- paste(datapath,n, p, k, sep = "/")
      gtbasepath <- paste(gtpath,n, p, k, sep = "/")
      resbasepath <- paste(respath,n, p, k, sep = "/")
      dir.create(resbasepath, showWarnings = FALSE, recursive = TRUE)
      for (r in 1:M){
        for (gen in gen_methods){
          filename <- paste(gen,n, p, k, r, sep = "_")
          x <- read.table(file = paste0(databasepath,"/", filename),header = TRUE)  
          for (est in torun){
            alg <- est_methods[[est]] 
            results <- alg(x)
            saveRDS(results, file = paste0(resbasepath, "/", est,"_",filename, ".rds"))
            message(paste(est, filename))
          }
        }
      }
    }
  }
}
