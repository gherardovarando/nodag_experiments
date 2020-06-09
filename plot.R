args <- commandArgs(trailingOnly = TRUE)
tablefile <- "TABLE.rds"
if (length(args) > 0){
  tablefile <- args[1]
}
library(ggplot2)
library(data.table)
source("util.R")
cols <- c(palette.colors(8, palette = "Classic Tableau"), "black")
types <- c(rep("solid", 3), rep("dashed", 5), "dotted")

methods <- c(
  "nodag-0.3",
  "nodag-0.2" ,
  "nodag-0.1" ,
  "pc-0.01" ,
  "pc-0.005" ,
  "pc-0.001" ,
  "tabu",
  "ges",
  "empty"
)
names(types) <- names(cols) <- methods

TABLE <- readRDS(paste0("simulations/", tablefile))

AVG <- apply(TABLE, c(1,3,4,5), mean, na.rm = TRUE)


data <- as.data.table(AVG)
data$linetype <- 2
data[meth %in% c("nodag-0.1","nodag-0.2", "nodag-0.3"), "linetype"] <- 1
data$linetype <- as.factor(data$linetype)
data$p <- as.numeric(data$p)
data[meth == "nodag-0.1" & stat %in% c("f1", "shd-cpdag", "shd-graph", "fpr", "tpr") & n == "100", "value"] <- NA


plot_select(data, methods, stats = c("f1", "tpr", "fpr") ,
                file = "plot_skeleton.pdf", cols = cols, types = types, height = 4)

plot_select_log(data, methods, stats = c("shd-graph","shd-cpdag") ,
            file = "plot_shd.pdf", cols = cols, types = types, height = 3.5)

plot_select_log(data, methods, stats = c("time") ,
                file = "plot_time.pdf", cols = cols, types = types, height = 2.5)

plot_select(data[p <= 100], c('nodag-0.1', "nodag-0.2", "nodag-0.3"), stats = c("iter") ,
            file = "plot_iter.pdf", cols = cols, types = types, height = 5)


alldata <- as.data.table(TABLE)

plot_density(alldata[p == "50"], methods, stats = c("shd-cpdag", "shd-graph") ,
             file = "plot_denisty.pdf", cols = cols, types = types)
