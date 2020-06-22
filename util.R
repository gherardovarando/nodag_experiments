library(igraph)
library(ggplot2)
library(Matrix)


igraph.to.tikz <- function (graph, layout, file = "") {
  ## Here we get the matrix layout
  if ("function" %in% class(layout))
    layout <- layout(graph)
  
  layout <- layout / max(abs(layout))
  mycat <- function(...)
    cat(..., file = file, append = TRUE)
  ##TikZ initialisation and default options (see pgf/TikZ manual)
  cat("\\tikzset{\n", file = file, append = FALSE)
  mycat(
    "\tnode/.style={circle,inner sep=1mm,minimum size=0.8cm,draw,
      very thick,black,fill=red!20,text=black},\n"
  )
  mycat("\tnondirectional/.style={very thick,black},\n")
  mycat("\tunidirectional/.style={nondirectional,shorten >=2pt,-stealth},\n")
  mycat("\tbidirectional/.style={unidirectional,bend right=10}\n")
  mycat("}\n")
  mycat("\n")
  
  ##Size
  mycat("\\begin{tikzpicture}[scale=5]\n")
  
  for (i in 1:length(V(graph))) {
    vertex <- V(graph)[i]
    label <- V(graph)[vertex]$name
    if (is.null(label))
      label <- ""
    
    ##drawing vertices
    mycat (
      sprintf (
        "\t\\node [node] (v%d) at (%f, %f)\t{%s};\n",
        vertex,
        layout[i, 1],
        layout[i, 2],
        label
      )
    )
  }
  mycat("\n")
  
  adj = igraph::get.adjacency(graph)
  print(adj)
  bidirectional = adj & t(adj)
  
  if (!is.directed(graph))
    ##undirected case
    for (line in 1:nrow(adj)) {
      for (col in line:ncol(adj)) {
        if (adj[line, col] & col > line) {
          mycat(sprintf (
            "\t\\path [nondirectional] (v%d) edge (v%d);\n",
            line,
            col
          )) ##edges drawing
        }
      }
    }
  else
    ##directed case
    for (line in 1:nrow(adj)) {
      for (col in 1:ncol(adj)) {
        if (bidirectional[line, col] & line > col)
          mycat(
            sprintf ("\t\\path [bidirectional] (v%d) edge (v%d);\n", line, col),
            sprintf ("\t\\path [bidirectional] (v%d) edge (v%d);\n", col, line)
          ) ##edges drawing
        else if (!bidirectional[line, col] & adj[line, col])
          mycat(sprintf ("\t\\path [unidirectional] (v%d) edge (v%d);\n", line, col)) ##edges drawing
      }
    }
  
  mycat("\\end{tikzpicture}\n")
}


plot_select <-
  function(data,
           algs,
           stats,
           cols = NULL,
           types = NULL,
           file = "plot.pdf",
           path = "simulations/",
           height = 5) {
    ggplot(data[stat %in% stats & meth %in% algs,],
           aes(
             x = p,
             y = value,
             group = meth,
             color = meth,
             linetype = meth
           )) +
      facet_grid(cols = vars(n),
                 rows = vars(stat),
                 scales = "free_y") +
      geom_line() +
      theme_bw() +
      # scale_y_log10()+
      labs(color = "Algorithm") +
      guides(col = guide_legend(nrow = 2)) +
      scale_color_manual('', values = cols) +
      scale_linetype_manual('', values = types) +
      theme(
        legend.position = "top",
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
        ),
        plot.margin = margin(
          t = 0,
          r = 0,
          b = 0,
          l = 0,
          unit = "pt"
        )
      ) +
      ggsave(
        file = file,
        path = path,
        width = 6,
        height = height,
        units = "in"
      )
  }

plot_density <- function(data,
                            algs,
                            stats,
                            file = "plot_density.pdf",
                            path = "simulations/",
                            cols = NULL,
                            types = NULL,
                            height = 5) {
  ggplot(data[stat %in% stats & meth %in% algs,],
         aes(
           x = value + 1,
           group = meth,
           color = meth,
           linetype = meth
         )) +
    facet_grid(cols = vars(n),
               rows = vars(stat),
               scales = "free_y") +
    geom_density() +
    theme_bw() +
    #ggtitle(title) +
    scale_x_log10()+
    labs(color = "Algorithm") +
    guides(col = guide_legend(nrow = 2)) +
    scale_color_manual("", values = cols) +
    scale_linetype_manual("", values = types) +
    theme(
      legend.position = "top",
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
      ),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0,
        unit = "pt"
      )
    ) +
    ggsave(
      file = file,
      path = path,
      width = 6,
      height = height,
      units = "in"
    )
}


plot_select_log <- function(data,
                            algs,
                            stats,
                            file = "plot.pdf",
                            path = "simulations/",
                            cols = NULL,
                            types = NULL,
                            height = 5) {
  ggplot(data[stat %in% stats & meth %in% algs,],
         aes(
           x = p,
           y = value,
           group = meth,
           color = meth,
           linetype = meth
         )) +
    facet_grid(cols = vars(n),
               rows = vars(stat),
               scales = "free_y") +
    geom_line() +
    theme_bw() +
    #ggtitle(title) +
    scale_y_log10()+
    labs(color = "Algorithm") +
    guides(col = guide_legend(nrow = 2)) +
    scale_color_manual("", values = cols) +
    scale_linetype_manual("", values = types) +
    theme(
      legend.position = "top",
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
      ),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0,
        unit = "pt"
      )
    ) +
    ggsave(
      file = file,
      path = path,
      width = 6,
      height = height,
      units = "in"
    )
}
