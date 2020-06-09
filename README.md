# Experiments for *Learning DAGs without imposing acyclicyty* 


This repository contains the code and the data to replicate the experiments in 
Varando G, *Learning DAGs without imposing acyclicyty* (2020) 
([arXiv:2006.03005](https://arxiv.org/abs/2006.03005)). 


## Requirements 

The following packages are needed: 

* `ggplot2`
* `pcalg`
* `bnlearn`
* `igraph`
* `Matrix`   

Moreover `nodag.f` (v0.0.2 from 
[gherardovarando/nodag](https://github.com/gherardovarando/nodag)) 
needs to be compiled with `R CMD SHLIB nodag.f -llapack -lblas`.

## Example 

Run the code in `example10p.R`, the estimated graph is saved in 
tikz format. 

## Simulations 

Run the scripts in the following order:

* `generate_data.R` to generate the data (~6 GB) 
* `estimate.R` to run the methods  
* `aggregate.R` to evaluate and collect results 
* `plot.R` to generate the plots  

## Protein signaling network

Run the code in `proteins.R` the estimated graph is saved in tikz format.  

## License 

`nodag.f` is available from 
[gherardovarando/nodag](https://github.com/gherardovarando/nodag) under  
BSD-3-Clause license. 
