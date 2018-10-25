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

pw_j_matrix = sapply(X = pws, FUN = function(i){
  sapply(X = pws, FUN = function(j){
    jei = length(intersect(i,j))/length(union(i,j))
  })
})

#I will be more interested in the size of the intersection regardless of similarity
pw_int_matrix = sapply(X = pws, FUN = function(i){
  sapply(X = pws, FUN = function(j){
    zei = length(intersect(i,j))
  })
})

sort(rowSums(pw_int_matrix))%>%head(10)

pw_int_matrix["Pentose phosphate pathway",c("Hedgehog signaling pathway",
                                            "NF-kappa B signaling pathway",
                                            "VEGF signaling pathway",
                                            "p53 signaling pathway",
                                            "Insulin signaling pathway",
                                            "RIG-I-like receptor signaling pathway"
                                            )]

pw_int_matrix["Notch signaling pathway",c("Hedgehog signaling pathway",
                                          "NF-kappa B signaling pathway",
                                          "VEGF signaling pathway",
                                          "p53 signaling pathway",
                                          "Insulin signaling pathway",
                                          "RIG-I-like receptor signaling pathway",
                                          "Pentose phosphate pathway"
)]

#next one NOPE
pw_int_matrix["Sphingolipid signaling pathway",c("Hedgehog signaling pathway",
                                                 "NF-kappa B signaling pathway",
                                                 "VEGF signaling pathway",
                                                 "p53 signaling pathway",
                                                 "Insulin signaling pathway",
                                                 "RIG-I-like receptor signaling pathway",
                                                 "Pentose phosphate pathway",
                                                 "Notch signaling pathway"
)]

pw_int_matrix["mRNA surveillance pathway",c("Hedgehog signaling pathway",
                                                 "NF-kappa B signaling pathway",
                                                 "VEGF signaling pathway",
                                                 "p53 signaling pathway",
                                                 "Insulin signaling pathway",
                                                 "RIG-I-like receptor signaling pathway",
                                                 "Pentose phosphate pathway",
                                                 "Notch signaling pathway"
)]


pw_int_matrix["TGF-beta signaling pathway",c("Hedgehog signaling pathway",
                                          "NF-kappa B signaling pathway",
                                          "VEGF signaling pathway",
                                          "p53 signaling pathway",
                                          "Insulin signaling pathway",
                                          "RIG-I-like receptor signaling pathway",
                                          "Pentose phosphate pathway",
                                          "Notch signaling pathway",
                                          "mRNA surveillance pathway"
)]


pw_int_matrix["IL-17 signaling pathway",c("Hedgehog signaling pathway",
                                            "NF-kappa B signaling pathway",
                                            "VEGF signaling pathway",
                                            "p53 signaling pathway",
                                            "Insulin signaling pathway",
                                            "RIG-I-like receptor signaling pathway",
                                            "Pentose phosphate pathway",
                                            "Notch signaling pathway",
                                            "mRNA surveillance pathway",
                                            "TGF-beta signaling pathway"
)]



############################################

SecondSet =     x[c("Pentose phosphate pathway",
                "Notch signaling pathway",
                "mRNA surveillance pathway",
                "TGF-beta signaling pathway",
                "IL-17 signaling pathway")]


my_SecondGraphs = lapply(X = SecondSet, FUN = function(i){
  igraph::igraph.from.graphNEL(graphite::pathwayGraph(pathway = i))
})

mySecond_NW = my_SecondGraphs[[1]]
mySecond_NW
for(i in my_SecondGraphs[2:5]){
  mySecond_NW = igraph::union(mySecond_NW, i)
}
mySecond_NW = as.undirected(mySecond_NW)
plot(mySecond_NW, vertex.label="")
mySecond_NW = remove.edge.attribute(mySecond_NW, "weight")
V(mySecond_NW)$name = gsub(pattern = "SYMBOL:", replacement = "", x = V(mySecond_NW)$name)
write.graph(graph = mySecond_NW, file = "results/example_02.gml", format = "gml")
