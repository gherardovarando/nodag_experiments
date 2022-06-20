source("generate_methods.R")
## make dir
datapath <- "simulations/data"
gtpath <- "simulations/gt"
dir.create(datapath, recursive = TRUE, showWarnings = FALSE)
dir.create(gtpath, recursive = TRUE, showWarnings = FALSE)


## number of replicates
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

for (n in ns){
  for (p in ps){
    for (k in ks){
      message("n ", n," p ", p, " k ", k)
      databasepath <- paste(datapath, n, p, k, sep = "/")
      gtbasepath <- paste(gtpath,n, p, k, sep = "/")
      dir.create(databasepath, recursive = TRUE, showWarnings = FALSE)
      dir.create(gtbasepath, recursive = TRUE, showWarnings = FALSE)
      for (r in 1:M){
        ### uniform sampling with gmat
        for (gen in gen_methods){
          res <- methods[[gen]](p, k, n, r)
          filename <- paste(gen,n, p, k, r, sep = "_")
          write.table(res$x, file = paste0(databasepath,"/", filename), 
                      row.names = FALSE, col.names = TRUE)  
          ## save the true graph
          saveRDS(object = res$true,  file = paste0(gtbasepath,"/dag_", filename, ".rds"))
          ## save the true coeff matrix
          write.table(res$coef, file = paste0(gtbasepath,"/coeff_", filename), 
                      row.names = FALSE, col.names = FALSE)  
        }
      }
    }
  }
}


