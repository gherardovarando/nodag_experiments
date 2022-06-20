source("generate_methods.R")
source("estimate_methods.R")

gt <- gen_pcalg_exp(p = 10,k =  2, n = 1000)

gt$coef
colnames(gt$x)
a <- est_notears(as.data.frame(gt$x), lambda = 0)

estpc <- est_pcalg(gt$x, alpha = 0.01)

estges <- est_ges(x = gt$x)

estnodag <- est_nodag(x = gt$x, lambda = 0.2)

esttabu <- est_tabu(as.data.frame(gt$x))

plot(estpc$graph)
plot(estges$graph)
plot(estnodag$graph)
plot(esttabu$graph)
plot(a$graph)
plot(gt$true)

