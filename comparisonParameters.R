library(data.table)
library(tidyverse)
library(igraph)
library(plyr)

#read graphs 
g     = read.graph(file = "results/example_01.gml", "gml") #main graph, from merging five known pathways
g_alt = read.graph(file = "results/example_02.gml", "gml") #alternative graph, from five different known pathways
g_mix = read.graph(file = "results/example_03.gml", "gml") #graph from two pathways used in g, and two pathways in g0
g_rew = read.graph(file = "results/example_06.gml", "gml") #rewired g (using each_edge(prob=0.25))

#read modules detected

modules_g     = readRDS(file = "results/modules_g.RDS")
modules_g_alt = readRDS(file = "results/modules_g_alt.RDS")
modules_g_mix = readRDS(file = "results/modules_g_mix.RDS")
modules_g_rew = readRDS(file = "results/modules_g_rew.RDS")

#read enrichments

enrichment_g     = fread(input = "results/enrichment_g.txt")
enrichment_g_alt = fread(input = "results/enrichment_g_alt.txt")
enrichment_g_mix = fread(input = "results/enrichment_g_mix.txt")
enrichment_g_rew = fread(input = "results/enrichment_g_rew.txt")

#Make List (useful for some analysis)
graphLists = list(main = g,
                  alt  = g_alt,
                  mix  = g_mix,
                  rew  = g_rew)

##########
#Step A ##
##########

#are they topologically similar
ParametersTopo =          rbind(
                                sapply(X = graphLists, FUN = vcount),
                                sapply(X = graphLists, FUN = ecount),
                                sapply(X = graphLists, FUN = no.clusters)
                                )


rownames(ParametersTopo) = c("No.Nodes", "No.Edges", "No.Components")
ParametersTopo = rownames_to_column(as.data.frame(ParametersTopo))
ParametersTopo

#degree distributions look the same? 

degree_df = data.frame(degree = 0:46,
                       main =as.numeric(table(factor(degree(g), levels = 0:46))),
                       alt = as.numeric(table(factor(degree(g_alt), levels = 0:46))),
                       mix = as.numeric(table(factor(degree(g_mix), levels = 0:46))),
                       rew = as.numeric(table(factor(degree(g_rew), levels = 0:46)))
)

p = ggplot(data = degree_df) + 
  geom_point(mapping = aes(degree, main, colour = "main")) + 
  geom_point(mapping = aes(degree, alt, colour = "alternative")) +
  geom_point(mapping = aes(degree, mix, colour = "mixed")) +
  geom_point(mapping = aes(degree, rew, colour = "rewired")) +
  labs(y = "frequency") +  
  guides(colour = guide_legend(title = "")) 
  
p



##########
#Step B ##
##########


#do they have similar nodes?

jaccard_nodes <- function(g1,g2){
  a = sort(vertex.attributes(graph = g1)[["name"]])
  b = sort(vertex.attributes(graph = g2)[["name"]])
  
  deLaCueva = length(intersect(a,b))/length(union(a,b))
  return(deLaCueva)
}

# NodesJaccard = c(
#           jaccard_nodes(g, g_alt),
#           jaccard_nodes(g, g_mix),
#           jaccard_nodes(g, g_rew)
# )

NodesJaccard = sapply(X = graphLists, FUN = jaccard_nodes, g1 = g)


#do they have similar edges?

jaccard_edges <- function(g1, g2){
  return(length(E(igraph::intersection(g1, g2)))/length(E(igraph::union(g1, g2))))
}

# EdgesJaccard = c(
#                   jaccard_edges(g, g_alt),
#                   jaccard_edges(g, g_mix),
#                   jaccard_edge(g, g_rew)
#                   )
EdgesJaccard = sapply(X = graphLists, FUN = jaccard_edges, g1 = g)
EdgesJaccard

##########
#Step C ##
##########

#graphs with equal node set

#main vs rew

possible_algos = c("vi", "nmi", "split.join", "rand", "adjusted.rand")

comparison_methods = sapply(X = possible_algos, FUN = function(i){
  igraph::compare(comm1 = modules_g,
                  comm2 = modules_g_rew,
                  method = i
  )
})




#graphs with unequal node set

# communities(modules_g)
# 
# testing  = sapply(X = communities(modules_g_mix), function(i){#columns
#   sapply(X = communities(modules_g), function(j){#rows
#     length(intersect(i,j))/length(union(i,j))
#   })
# })
# 
# sum(rowSums(testing))/25
# which(rowSums(testing)!=0)
# which(rowSums(testing)==1)

moduleComparer = function(reference_modules, test_modules){
jmat =  sapply(X = communities(test_modules), function(i){#columns
          sapply(X = communities(reference_modules), function(j){
            length(intersect(i,j))/length(union(i,j))
    })
  })

jmat_rowsums = rowSums(x = jmat)
jmat_nonzero = apply(jmat, 1, function(i) length(which(i!=0)))
simScore_mod = ifelse(test = jmat_rowsums==0, 0, jmat_rowsums/jmat_nonzero)

#SimScore = sum(rowSums(jmat))/length(communities(reference_modules))
SimScore = sum(simScore_mod)/length(simScore_mod)
matchingModules = length(which(rowSums(jmat)!=0))
#perfectMatches  = length(which(rowSums(jmat)==1))
perfectMatches  = length(which(jmat_rowsums==1 & jmat_nonzero==1))
resultados = list(similarity_matrix = jmat,
                  Similarity_Score  = SimScore,
                  Matching.Modules  = matchingModules,
                  Perfect.Matches   = perfectMatches
                  )

return(resultados)
}

moduleComparer(reference_modules = modules_g, test_modules = modules_g_alt)$similarity_matrix%>%rowSums
length(communities(modules_g))

moduleComparer(reference_modules = modules_g, test_modules = modules_g_alt)$similarity_matrix%>%as.data.frame%>%xtable(digits = 2)
moduleComparer(reference_modules = modules_g, test_modules = modules_g_mix)$similarity_matrix%>%as.data.frame%>%xtable(digits = 2)
moduleComparer(reference_modules = modules_g, test_modules = modules_g_rew)$similarity_matrix%>%as.data.frame%>%xtable(digits = 2)

temporal1 = moduleComparer(reference_modules = modules_g, test_modules = modules_g_alt)[2:4]%>%as.data.frame
temporal2 = moduleComparer(reference_modules = modules_g, test_modules = modules_g_mix)[2:4]%>%as.data.frame
temporal3 = rbind(temporal1, temporal2)
temporal3
xtable(temporal3, digits = 3)
rm(temporal1, temporal2, temporal3)





moduleComparer_edges = function(reference, 
                                reference_modules, 
                                test,
                                test_modules){
  
  graphs_reference = lapply(X = communities(reference_modules), function(i){
    return(induced.subgraph(graph = reference, vids = i))
  })
  
  graphs_test = lapply(X = communities(test_modules), function(i){
    induced.subgraph(graph = test, vids = i)
  })
  
  jmat =  sapply(X = graphs_test, function(i){#columns
    sapply(X = graphs_reference, function(j){
      rrr = length(E(igraph::intersection(i, j)))/length(E(igraph::union(i, j)))
    })
  })
  
  jmat_rowsums = rowSums(x = jmat)
  jmat_nonzero = apply(jmat, 1, function(i) length(which(i!=0)))
  simScore_mod = ifelse(test = jmat_rowsums==0, 0, jmat_rowsums/jmat_nonzero)
  
  #SimScore = sum(rowSums(jmat))/length(communities(reference_modules))
  SimScore = sum(simScore_mod)/length(simScore_mod)
  matchingModules = length(which(rowSums(jmat)!=0))
  #perfectMatches  = length(which(rowSums(jmat)==1))
  perfectMatches  = length(which(jmat_rowsums==1 & jmat_nonzero==1))
  resultados = list(similarity_matrix = jmat,
                    Similarity_Score  = SimScore,
                    Matching.Modules  = matchingModules,
                    Perfect.Matches   = perfectMatches
  )
  
  return(resultados)
}

moduleComparer_edges(reference = g, 
                     reference_modules = modules_g, 
                     test = g_alt, 
                     test_modules = modules_g_alt)

moduleComparer_edges(reference = g, 
                     reference_modules = modules_g, 
                     test = g_mix, 
                     test_modules = modules_g_mix)

temporal1 = moduleComparer_edges(reference = g, reference_modules = modules_g, test = g_alt, test_modules = modules_g_alt)[2:4]%>%as.data.frame
temporal2 = moduleComparer_edges(reference = g, reference_modules = modules_g, test = g_mix, test_modules = modules_g_mix)[2:4]%>%as.data.frame
temporal3 = rbind(temporal1, temporal2)
temporal3
xtable(temporal3, digits = 3)
rm(temporal1, temporal2, temporal3)
