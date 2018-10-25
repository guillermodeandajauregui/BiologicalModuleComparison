source("enrichmentFunction.R")
library(data.table)
library(tidyverse)
library(igraph)
library(plyr)

#read graphs 
g     = read.graph(file = "results/example_01.gml", "gml") #main graph, from merging five known pathways
g_alt = read.graph(file = "results/example_02.gml", "gml") #alternative graph, from five different known pathways
g_mix = read.graph(file = "results/example_03.gml", "gml") #graph from two pathways used in g, and two pathways in g0
g_rew = read.graph(file = "results/example_06.gml", "gml") #rewired g (using each_edge(prob=0.25))

#detect modules with Infomap

modules_g     = infomap.community(g,     nb.trials = 1000)
modules_g_alt = infomap.community(g_alt, nb.trials = 1000)
modules_g_mix = infomap.community(g_mix, nb.trials = 1000)
modules_g_rew = infomap.community(g_rew, nb.trials = 1000)

#save
saveRDS(object = modules_g,     file = "results/modules_g.RDS")
saveRDS(object = modules_g_alt, file = "results/modules_g_alt.RDS")
saveRDS(object = modules_g_mix, file = "results/modules_g_mix.RDS")
saveRDS(object = modules_g_rew, file = "results/modules_g_rew.RDS")

#plotting
plot(modules_g, g, vertex.label="", main = "main")
plot(modules_g_alt, g_alt, vertex.label="", main = "alt")
plot(modules_g_mix, g_mix, vertex.label="", main = "mix")
plot(modules_g_rew, g_rew, vertex.label="", main = "rew")

#enrichment 

#for this work, we will use an ad-hoc list of pathways 
#containing the 10 pathways used for the construction of the example networks 
#plus 2 other pathways 

source("auxiliar_pathwayList.R") #as auxiliar script to keep the main script readable
                                 #generates my_pathways object


enrichment_g = ModuleEnrichment(g = g, 
                                moduleList = communities(modules_g), 
                                pathwayList = my_pathways, 
                                threshold = 0.05)

enrichment_g_alt = ModuleEnrichment(g = g_alt, 
                                moduleList = communities(modules_g_alt), 
                                pathwayList = my_pathways, 
                                threshold = 0.05)

enrichment_g_mix = ModuleEnrichment(g = g_mix, 
                                    moduleList = communities(modules_g_mix), 
                                    pathwayList = my_pathways, 
                                    threshold = 0.05)

enrichment_g_rew = ModuleEnrichment(g = g_rew, 
                                    moduleList = communities(modules_g_rew), 
                                    pathwayList = my_pathways, 
                                    threshold = 0.05)

#write out enrichments

fwrite(enrichment_g,     file = "results/enrichment_g.txt",     quote = FALSE, sep = "\t", col.names = TRUE)
fwrite(enrichment_g_alt, file = "results/enrichment_g_alt.txt", quote = FALSE, sep = "\t", col.names = TRUE)
fwrite(enrichment_g_mix, file = "results/enrichment_g_mix.txt", quote = FALSE, sep = "\t", col.names = TRUE)
fwrite(enrichment_g_rew, file = "results/enrichment_g_rew.txt", quote = FALSE, sep = "\t", col.names = TRUE)
