library(data.table)
library(tidyverse)
library(igraph)
library(plyr)

#detect modules
g = read.graph(file = "results/example_01.gml", "gml")
set.seed(1)
g_rewired <- rewire(g, each_edge(prob = 0.5))
write.graph(graph = g_rewired, file = "results/example_04.gml", format = "gml")

?rewire
set.seed(1)
g_rewired2 <- rewire(g, keeping_degseq(loops = FALSE, niter = 100))
write.graph(graph = g_rewired2, file = "results/example_05.gml", format = "gml")

set.seed(1)
g_rewired3 <- rewire(g, each_edge(prob = 0.25))
write.graph(graph = g_rewired3, file = "results/example_06.gml", format = "gml")
