graphLists = list(main = g,
                  alt  = g_alt,
                  mix  = g_mix,
                  rew  = g_rew)

#are they topologically similar
rbind(
sapply(X = graphLists, FUN = vcount),
sapply(X = graphLists, FUN = ecount),
sapply(X = graphLists, FUN = no.clusters)
)


#do they have similar nodes?

jaccard_nodes <- function(g1,g2){
  a = sort(vertex.attributes(graph = g1)[["name"]])
  b = sort(vertex.attributes(graph = g2)[["name"]])
  
  deLaCueva = length(intersect(a,b))/length(union(a,b))
  return(deLaCueva)
}

jaccard_nodes(g, g_alt)
jaccard_nodes(g, g_mix)
jaccard_nodes(g, g_rew)

#do they have similar edges?

jaccard_edges <- function(g1, g2){
  return(length(E(igraph::intersection(g1, g2)))/length(E(igraph::union(g1, g2))))
}

jaccard_edges(g, g_alt)
jaccard_edges(g, g_mix)
jaccard_edge(g, g_rew)


table(degree(g))
table(degree(g_alt))
table(degree(g_mix))

hist(degree(g), breaks = 50, xlim = c(0,50))

plot(x = names(table(degree(g))), y = (table(degree(g))))
plot(x = names(table(degree(g_alt))), y = (table(degree(g_alt))))


data.frame(degree = 0:46,
           main =as.numeric(table(factor(degree(g), levels = 0:46))),
           alt = as.numeric(table(factor(degree(g_alt), levels = 0:46))),
           mix = as.numeric(table(factor(degree(g_mix), levels = 0:46))),
           rew = as.numeric(table(factor(degree(g_rew), levels = 0:46)))
           )
