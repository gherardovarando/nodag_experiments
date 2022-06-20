args <- commandArgs(trailingOnly = TRUE)
tablefile <- "TABLE.rds"
if (length(args) > 0){
  tablefile <- args[1]
}
library(ggplot2)
library(RColorBrewer)
library(data.table)
source("util.R")
cols <- c(brewer.pal(n = 8, name = "Paired"), "black")
types <- c(rep("solid", 3), rep("dashed", 5), "dotted")

methods <- c(
  "nodag-0.3",
  "nodag-0.2" ,
  #"nodag-0.1" ,
  "notears-0.3",
  "notears-0.2" , 
  "notears-0.1" , 
  # "pc-0.1" ,
  # "pc-0.05" ,   
  "pc-0.01" ,
  #"pc-0.005" ,
  "pc-0.001" ,
  "tabu",
  "ges",
  "empty"
)
names(types) <- names(cols) <- methods

TABLE <- readRDS(paste0("simulations/", tablefile))

AVG <- apply(TABLE, c(1,2,3,4,5,6), mean, na.rm = TRUE)
data <- as.data.table(AVG)
data$p <- as.numeric(data$p)
data$n <- as.numeric(data$n)

ggplot(data[stat %in% c("shd-dag", "shd-cpdag", "time") & p == 50  & meth %in% methods 
            & gen == "randomDAG_exp"]) +
  geom_line(aes(
    x = n,
    y = value,
    group = meth,
    col = meth
  )) + theme_bw() + 
  #scale_y_log10()+
  facet_grid(rows = vars(stat),
             cols = vars(k),
             scales = "free_y"
  ) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 30),
        legend.box.spacing = unit(0.5, "lines"),
        legend.box.margin = margin(
          t = 0,
          r = 0,
          b = 0,
          l = 0,
          unit = "pt"
        ))  +
  ggsave(
    file = 'plot_1.pdf',
    path = "simulations",
    width = 6,
    height = 3,
    units = "in"
  )
