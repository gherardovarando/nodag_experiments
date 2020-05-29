args <- commandArgs(trailingOnly = TRUE)
tablefile <- "TABLE.rds"
if (length(args) > 0){
  tablefile <- args[1]
}
library(ggplot2)
library(data.table)
source("util.R")
cols <- palette.colors(6, palette = "Classic Tableau")
types <- c(2,1,1,1,2,2)

names(types) <- names(cols) <-   c(
  "pc-0.01",
"nodag-0.3",
"nodag-0.2", 
"nodag-0.1",
"tabu",
"ges" 
)

TABLE <- readRDS(paste0("simulations/", tablefile))

AVG <- apply(TABLE, c(1,3,4,5), mean, na.rm = TRUE)


data <- as.data.table(AVG)
data$linetype <- 2
data[meth %in% c("nodag-0.1","nodag-0.2", "nodag-0.3"), "linetype"] <- 1
data$linetype <- as.factor(data$linetype)
data$p <- as.numeric(data$p)
data[meth == "nodag-0.1" & stat == "fpr" & n == "100", ] <- NA

algsel <- c("pc-0.01", "nodag-0.1", "nodag-0.2", "nodag-0.3", "ges", "tabu")
plot_select(data[p <= 100], algsel, stats = c("f1", "tpr", "fpr") ,
                file = "plot_skeleton.pdf", cols = cols, types = types)

plot_select_log(data[p <= 100], algsel, stats = c("shd", "time") ,
            file = "plot_shd_time.pdf", cols = cols, types = types)

plot_select(data[p <= 100], c('nodag-0.1', "nodag-0.2", "nodag-0.3"), stats = c("iter") ,
            file = "plot_iter.pdf", cols = cols, types = types)


alldata <- as.data.table(TABLE)
hist(alldata[stat == "shd" & meth == "pc-0.01"  & n == "1000"]$value)
hist(alldata[stat == "shd" & meth == "nodag-0.2"  & n == "1000"]$value)

plot_density(alldata[p == "100"], algsel, stats = c("f1") ,
             file = "plot_denisty.pdf", cols = cols, types = types)
