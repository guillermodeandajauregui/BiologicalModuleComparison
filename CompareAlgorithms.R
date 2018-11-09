library(data.table)
library(tidyverse)
library(igraph)
library(plyr)

#read graphs 
g     = read.graph(file = "results/example_01.gml", "gml")

#detect modules 



modulesGlist = list(infomap           = infomap.community(g, nb.trials = 1000, e.weights = NULL),
                    fast.greedy       = cluster_fast_greedy(g, weights = NULL),
                    louvain           = cluster_louvain(g, weights = NULL),
                    walktrap          = cluster_walktrap(g, weights = NULL),
                    spinglass         = cluster_spinglass(g, weights = NULL),
                    edge.betweenness  = cluster_edge_betweenness(g, weights = NULL),
                    leading.eigen     = cluster_leading_eigen(g, weights = NULL),
                    label.prop        = cluster_label_prop(g, weights = NULL)
                    )


algo_comp_matrix =          sapply(X = modulesGlist, FUN = function(i){
                              sapply(X = modulesGlist, FUN = function(j){
                                igraph::compare(i,j, "nmi")
                              })
                            })

save(algo_comp_matrix, modulesGlist, g, file = "results/algorithmComparisons.RData")
