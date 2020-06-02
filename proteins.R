library(igraph)
source("util.R")

filenames <- c(
  "cd3cd28",
  "cd3cd28icam2",
  "cd3cd28+aktinhib",
  "cd3cd28+g0076",
  "cd3cd28+psitect",
  "cd3cd28+u0126",
  "cd3cd28+ly",
  "pma",
  "b2camp"
)
extension <- ".csv"
basepath <- "data/Sachs/Data Files/"
completepaths <- paste0(basepath, 1:9, ". ", filenames, extension)

bpath <- "proteins/all/"
dir.create(bpath, showWarnings = FALSE, recursive = TRUE)
conditions <- 1:9

########################## LOADING DATA 

datalist <- lapply(
  conditions,
  FUN = function(cond) {
    fn <- completepaths[cond]
    tmp <- read.csv(fn)
    colnames(tmp) <- tolower(colnames(tmp))
    return(tmp)
  }
)

p <- ncol(datalist[[1]])


all <- datalist[[1]]
if (length(datalist) > 1){
  for (cond in 2:length(datalist)){
    all <- rbind(all, datalist[[cond]])
  }
}

message("data loaded correctly ", dim(all))

########################


system.time(res <- est_nodag(all, lambda = 0.2))

A <- res$coef 
colnames(A) <- rownames(A) <- colnames(all)
adj <- sign(abs(A)) -diag(p)
g <- graph_from_adjacency_matrix(adjmatrix = adj)
################## 

nicenames <- c("Raf", "Mek$1/2$", "PLC$_\\gamma$", "PIP2", 
               "PIP3", "Erk$1/2$", "Akt", "PKA", "PKC", "p38", "JNK" )
layout <- matrix(nrow = 11, ncol = 2, byrow = TRUE,
                 data = c(10, 6, #RAF 
                          10, 3, #MEK
                          0, 7, #PLCG
                          0, 2, #PIP2
                          2, 4, #PIP3
                          10, 0, #ERK (P44.42)
                          6, 0, #AKT (PAKTS473)
                          7, 6.5, #PKA
                          5, 8, #PKC
                          5, 5, #P38
                          3, 5  #JNK
                 ))

plot(g, layout = layout)
igraph.to.tikz(g,layout = layout, file = "proteins.nodag-0.2.txt")
sum(A!=0) - p

