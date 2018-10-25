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

pws = y[grep(pattern = "pathway", x = names(y))]



#We will take two pathways from graph 01 and two pathways from graph 02
ThirdSet =     x[c("Pentose phosphate pathway",
                    "Notch signaling pathway",
                    "Hedgehog signaling pathway",
                    "p53 signaling pathway"
                    )]


my_ThirdGraphs = lapply(X = ThirdSet, FUN = function(i){
  igraph::igraph.from.graphNEL(graphite::pathwayGraph(pathway = i))
})

myThird_NW = my_ThirdGraphs[[1]]
myThird_NW
for(i in my_ThirdGraphs[2:5]){
  myThird_NW = igraph::union(myThird_NW, i)
}
myThird_NW = as.undirected(myThird_NW)
plot(myThird_NW, vertex.label="")
myThird_NW = remove.edge.attribute(myThird_NW, "weight")
V(myThird_NW)$name = gsub(pattern = "SYMBOL:", replacement = "", x = V(myThird_NW)$name)
write.graph(graph = myThird_NW, file = "results/example_03.gml", format = "gml")
