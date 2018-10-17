#explore KEGG database to build an example network of known pathways 

library(data.table)
library(tidyverse)
library(igraph)
library(graphite)

#make a list of KEGG pathways
x = graphite::pathways("hsapiens", "kegg")
#convert ids to symbol
x = convertIdentifiers(x, "symbol")

#make a list of nodes in pathways
#also remove the SYMBOL prefix (not necessary now)
y = lapply(X = x, FUN = function(i){
  tmp = nodes(i)
  tmp = gsub(pattern = "SYMBOL:", 
             replacement = "", 
             x = tmp)
  return(tmp)
})

#make a Jaccard matrix
j_matrix = sapply(X = y, FUN = function(i){
  sapply(X = y, FUN = function(j){
    jei = length(intersect(i,j))/length(union(i,j))
  })
})

#I will be more interested in the size of the intersection regardless of similarity
int_matrix = sapply(X = y, FUN = function(i){
  sapply(X = y, FUN = function(j){
    zei = length(intersect(i,j))
  })
})

#which pathways intersect in exactly one node?
w = which(int_matrix==1, arr.ind = TRUE)

#pick two y[133] y[100]
#hedgehog, NF-kB

#Now I want a non-crosstalking pathway
#want to find some that are actually called pathway 

##grep(pattern = "pathway", x = names(which(int_matrix[133,]==0)), value = TRUE)

#and two pathways that crosstalk with this one

##which(int_matrix["VEGF signaling pathway",]==1)
##which(int_matrix["p53 signaling pathway",]==1)

#Hedgehog signaling pathway 133
#NF-kappa B signaling pathway 100
#VEGF signaling pathway #136
#p53 signaling pathway 109
#Insulin signaling pathway 185
#RIG-I-like receptor signaling pathway 154

#Now that I have my pathways, I would like to get my graphs 

my_ids = c(100, 109, 133, 136, 154)
my_pathways = x[my_ids]

my_pw_graphs = lapply(X = my_pathways, FUN = function(i){
  igraph::igraph.from.graphNEL(graphite::pathwayGraph(pathway = i))
})

my_nw = my_pw_graphs[[1]]
my_nw
for(i in my_pw_graphs[2:5]){
  my_nw = igraph::union(my_nw, i)
}
my_nw = as.undirected(my_nw)
plot(my_nw, vertex.label="")
my_nw = remove.edge.attribute(my_nw, "weight")

V(my_nw)$name = gsub(pattern = "SYMBOL:", replacement = "", x = V(my_nw)$name)
write.graph(graph = my_nw, file = "results/example_01.gml", format = "gml")

#for enrichment purposes, I want all the pathways that are actually called "pathway"

my_pathway_biglist = y[grep(pattern = "pathway", x = names(y))]
saveRDS(object = my_pathway_biglist, file = "results/my_pathway_biglist.RDS")
