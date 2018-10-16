library(data.table)
library(tidyverse)
library(igraph)
library(plyr)

#detect modules
g = read.graph(file = "results/example_01.gml", "gml")

#infomap
modules_infomap = infomap.community(g, nb.trials = 1000)
communities(modules_infomap)
plot(modules_infomap, g, vertex.label="")

comms_infomap = communities(modules_infomap)

#enrichment
enrichment_infomap = lapply(X = seq_along(comms_infomap), FUN = function(i){
  #enrichment_list_cases = lapply(X = seq_along(l_comm_cases), FUN = function(i){
  
  nomen = names(comms_infomap)[i]
  
  my_enrichment = HTSanalyzeR::multiHyperGeoTest(collectionOfGeneSets = my_pathway_shortlist, 
                                                 universe = V(g)$name, 
                                                 hits = comms_infomap[[i]], 
                                                 minGeneSetSize = 1, 
                                                 pAdjustMethod = "BH", 
                                                 verbose = TRUE
  )
  print(nrow(my_enrichment))
  my_2 = tibble::rownames_to_column(as.data.frame(my_enrichment))#%>%filter(Adjusted.Pvalue<10)
  my_3 = cbind(comm = nomen, my_2)
  return(my_3)
  
})

enrichment_infomap = ldply(enrichment_infomap, data.frame)
enrichment_infomap$Adjusted.Pvalue2 = p.adjust(p = enrichment_infomap$Pvalue, method = "BH")
filter(enrichment_infomap, Adjusted.Pvalue<0.05)


#louvain
modules_louvain = cluster_louvain(graph = g)
communities(modules_louvain)
plot(modules_louvain, g, vertex.label="")

comms_louvain = communities(modules_louvain)

enrichment_louvain = lapply(X = seq_along(comms_louvain), FUN = function(i){
  #enrichment_list_cases = lapply(X = seq_along(l_comm_cases), FUN = function(i){
  
  nomen = names(comms_louvain)[i]
  
  my_enrichment = HTSanalyzeR::multiHyperGeoTest(collectionOfGeneSets = my_pathway_shortlist, 
                                                 universe = V(g)$name, 
                                                 hits = comms_louvain[[i]], 
                                                 minGeneSetSize = 1, 
                                                 pAdjustMethod = "BH", 
                                                 verbose = TRUE
  )
  print(nrow(my_enrichment))
  my_2 = tibble::rownames_to_column(as.data.frame(my_enrichment))#%>%filter(Adjusted.Pvalue<10)
  my_3 = cbind(comm = nomen, my_2)
  return(my_3)
  
})

enrichment_louvain = ldply(enrichment_louvain, data.frame)
enrichment_louvain$Adjusted.Pvalue2 = p.adjust(p = enrichment_louvain$Pvalue, method = "BH")
filter(enrichment_louvain, Adjusted.Pvalue2<0.05)

#walktrap
modules_walktrap = cluster_walktrap(graph = g)
communities(modules_walktrap)
plot(modules_walktrap, g, vertex.label="")

#label prop
modules_label_prop = cluster_label_prop(graph = g)
communities(modules_label_prop)
plot(modules_label_prop, g, vertex.label="")
