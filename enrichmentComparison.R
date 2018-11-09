library(data.table)
library(tidyverse)
library(igraph)
library(plyr)

#Functions to use

jaccard_simplex = function(a,b){
  length(intersect(a,b))/length(union(a,b))
}


EnrichmentEuclidianDistance = function(enrichment_1, 
                                       enrichment_2){
  
  tally_1 = enrichment_1 %>% group_by(rowname) %>% tally
  tally_2 = enrichment_2 %>% group_by(rowname) %>% tally
  join.df = full_join(x = tally_1, y = tally_2, "rowname")
  join.df = mutate(join.df, n.x = ifelse(is.na(n.x), 0, n.x))
  join.df = mutate(join.df, n.y = ifelse(is.na(n.y), 0, n.y))
  join.df = mutate(join.df, diff.squared = (n.x - n.y)**2)
  #join.df = mutate(join.df, diff = (n.x - n.y))
  resulta = sum(join.df$diff.squared)
  #resulta = sum(join.df$diff)
  resulta = sqrt(resulta)
  #print(nrow(join.df))
  return(resulta)
  
}

EnrichmentManhattanDistance = function(enrichment_1, 
                                       enrichment_2){
  
  tally_1 = enrichment_1 %>% group_by(rowname) %>% tally
  tally_2 = enrichment_2 %>% group_by(rowname) %>% tally
  join.df = full_join(x = tally_1, y = tally_2, "rowname")
  join.df = mutate(join.df, n.x = ifelse(is.na(n.x), 0, n.x))
  join.df = mutate(join.df, n.y = ifelse(is.na(n.y), 0, n.y))
  join.df = mutate(join.df, abs.diff = abs(n.x - n.y))
  resulta = sum(join.df$abs.diff)
  return(resulta)
  
}

ModuleFunctionalSimilarity = function(enrichment_1, enrichment_2){
  list_one = split(enrichment_1, enrichment_1$comm)
  list_one = lapply(list_one, FUN=function(i) i$rowname)
  
  list_two = split(enrichment_2, enrichment_2$comm)
  list_two = lapply(list_two, FUN=function(i) i$rowname)
  
  rr =sapply(X = list_one, function(i){#columns
    sapply(X = list_two, function(j){#rows
      delaCueva = jaccard_simplex(i,j)
      delaCueva = ifelse(is.na(delaCueva), 0, delaCueva)
    })
  })
  
  return(rr)
}

#read enrichments

enrichment_g     = fread(input = "results/enrichment_g.txt")
enrichment_g_alt = fread(input = "results/enrichment_g_alt.txt")
enrichment_g_mix = fread(input = "results/enrichment_g_mix.txt")
enrichment_g_rew = fread(input = "results/enrichment_g_rew.txt")

enrichment_List = list(enrichment_g = enrichment_g,
                       enrichment_g_alt = enrichment_g_alt,
                       enrichment_g_mix = enrichment_g_mix,
                       enrichment_g_rew = enrichment_g_rew )

#Similarity of enriched sets 

EnrichedProcessJ = sapply(X = enrichment_List, FUN = function(i){
                   jaccard_simplex(a = unique(enrichment_g[,"rowname"]),
                                  b = unique(i[,"rowname"]))
                })

#number of communities by process

EnrichmentEuclidianDistances = sapply(X = enrichment_List, FUN = function(i){
  EnrichmentEuclidianDistance(enrichment_1 = enrichment_g, enrichment_2 = i)
})

EnrichmentManhattanDistances = sapply(X = enrichment_List, FUN = function(i){
  EnrichmentManhattanDistance(enrichment_1 = enrichment_g, enrichment_2 = i)
})

EnrichmentManhattanDistances
EnrichmentEuclidianDistances

#simmilarity of communities process sets 

ModuleFunctionalSimilarity(enrichment_g, enrichment_g)
ModuleFunctionalSimilarity(enrichment_g, enrichment_g)%>%apply(MARGIN = 2, FUN = function(i){length(which(i>0.5))})


ModuleFunctionalSimilarity(enrichment_g, enrichment_g_alt)
ModuleFunctionalSimilarity(enrichment_g, enrichment_g_alt)%>%apply(MARGIN = 2, FUN = function(i){length(which(i>0.5))})

ModuleFunctionalSimilarity(enrichment_g, enrichment_g_mix)
ModuleFunctionalSimilarity(enrichment_g, enrichment_g_mix)%>%apply(MARGIN = 2, FUN = function(i){length(which(i>0.5))})

ModuleFunctionalSimilarity(enrichment_g, enrichment_g_rew)
ModuleFunctionalSimilarity(enrichment_g, enrichment_g_rew)%>%apply(MARGIN = 2, FUN = function(i){length(which(i>0.5))})
