library(glasso)
library(igraph)
library(ggplot2)
library(pcalg)
library(bnlearn)
dyn.load("../src/nodag.so")

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
  
  adj = get.adjacency(graph)
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

gen_gmat <- function(p, k, n, r) {
  dag <- rgraph(p, k / p, dag = TRUE, ordered = TRUE)
  L <- gmat:::mh_u(1, p = p, dag = dag)[, , 1]
  Sigma <- solve(tcrossprod(L))
  true <- dag2cpdag(as_graphnel(dag))
  X <- MASS::mvrnorm(n, rep(0, p), Sigma = Sigma)
  name <- paste("gmat_unif", p, k, r, sep = "_")
  return(list(
    true = true,
    coef = L,
    x = X,
    filename = name
  ))
}

gen_pcalg <- function(p, k, n, r, ...) {
  gGtrue <- randomDAG(p, prob = k / p, ...)
  x <- rmvDAG(n, gGtrue)
  true <- dag2cpdag(gGtrue)
  name <- paste("pcdag_randomDAG", p, k, r, sep = "_")
  return(list(
    true = true,
    coef = wgtMatrix(gGtrue),
    x = x,
    filename = name
  ))
}

gen_pcalg_exp <- function(p, k, n, r) {
  gGtrue <- randomDAG(p, prob = k / p)
  errMat <- matrix(nrow = n,
                   ncol = p,
                   data = rexp(n * p))
  x <- rmvDAG(n, gGtrue, errMat = )
  true <- dag2cpdag(gGtrue)
  name <- paste("pcdag_randomDAG", p, k, r, sep = "_")
  return(list(
    true = true,
    coef = wgtMatrix(gGtrue),
    x = x,
    filename = name
  ))
}

############# estimate func

est_pcalg <- function(x, ...) {
  time <- system.time({
    suffStat <- list(C = cor(x), n = nrow(x))
    pcout <-
      pc(suffStat = suffStat,
         indepTest = gaussCItest,
         p = ncol(x),
         maj.rule = TRUE,
         solve.confl = TRUE,
         ...)
  })[3]
  return(list(time = time, graph = pcout@graph))
}

est_fnodag <- function(x, ...) {
  arg <- list(...)
  lambda <- arg$lambda
  toll <- ifelse(is.null(arg$toll), 1e-4, arg$toll)
  maxitr <- ifelse(is.null(arg$maxitr), 1000, arg$maxitr)
  time <-
    system.time(
      out <- .Fortran(
        "FNODAG",
        as.integer(ncol(x)),
        as.double(cor(x)),
        as.double(diag(ncol(x))),
        as.double(lambda),
        as.double(toll),
        as.double(0.5),
        as.integer(maxitr)
      )
    )[3]
  A <- matrix(nrow = ncol(x), out[[3]])
  g <- graph_from_adjacency_matrix(sign(abs(A)), diag = FALSE)
  return(list(
    time = time,
    graph = as_graphnel(g),
    coef = A,
    itr = out[[7]]
  ))
}

est_nodag <- function(x, ...) {
  arg <- list(...)
  lambda <- arg$lambda
  toll <- ifelse(is.null(arg$toll), 1e-4, arg$toll)
  maxitr <- ifelse(is.null(arg$maxitr), 1000, arg$maxitr)
  time <-
    system.time(
      out <- .Fortran(
        "NODAG",
        as.integer(ncol(x)),
        as.double(cor(x)),
        as.double(diag(ncol(x))),
        as.double(lambda),
        as.double(toll),
        as.double(0.5),
        as.integer(maxitr)
      )
    )[3]
  A <- matrix(nrow = ncol(x), out[[3]])
  g <- graph_from_adjacency_matrix(sign(abs(A)), diag = FALSE)
  return(list(
    time = time,
    graph = as_graphnel(g),
    coef = A,
    itr = out[[7]]
  ))
}

est_tabu <- function(x, ...) {
  time <- system.time(bn <- tabu(x, maxp = 10, ...))[3]
  return(list(time = time, graph  = dag2cpdag(as.graphNEL(bn))))
}

est_ges <- function(x, ...) {
  time <- system.time({
    score <- new("GaussL0penObsScore", x)
    ges.fit <- ges(score)
  })[3]
  return(list(
    time = time,
    graph = as(ges.fit$essgraph, "graphNEL")
  ))
}


est_chowliu <- function(x, ...) {
  time <- system.time({
    res <- bnlearn::chow.liu(as.data.frame(x))
  })
  return(list(time = time, graph  = dag2cpdag(as.graphNEL(res))))
}

est_lingam <- function(x, ...) {
  time <- system.time({
    res <- pcalg::lingam(x)
  })
  g <- as(abs(t(res$Bpruned)), "graphNEL")
  return(list(
    time = time,
    graph  = dag2cpdag(g),
    coef = res$Bpruned
  ))
}

plot_select <-
  function(data,
           algs,
           stats,
           cols = NULL,
           types = NULL,
           file = "plot.pdf",
           path = "simulations/") {
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
      guides(col = guide_legend(nrow = 1)) +
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
        height = 3,
        units = "in"
      )
  }

plot_density <- function(data,
                            algs,
                            stats,
                            file = "plot_density.pdf",
                            path = "simulations/",
                            cols = NULL,
                            types = NULL) {
  ggplot(data[stat %in% stats & meth %in% algs,],
         aes(
           x = value,
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
    guides(col = guide_legend(nrow = 1)) +
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
      height = 3,
      units = "in"
    )
}


plot_select_log <- function(data,
                            algs,
                            stats,
                            file = "plot.pdf",
                            path = "simulations/",
                            cols = NULL,
                            types = NULL) {
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
    guides(col = guide_legend(nrow = 1)) +
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
      height = 3,
      units = "in"
    )
}