#modular centrality
#as seen in https://export.arxiv.org/abs/1810.05101
library(igraph)

#input: graph, centrality measure, modules
g
modules_infomap
#1) Choose a standard centrality measure
##degree

#2) remove all inter-community edges to obtain local network Gl

igraph::crossing(modules_infomap, g)

Gl = delete.edges(g, edges = E(g)[igraph::crossing(modules_infomap, g)])

#3) compute local measure Bl 
degree_local = degree(Gl)

#4) remove intra community edges
Gg = delete.edges(g, edges = E(g)[!(igraph::crossing(modules_infomap, g))])

#5) compute global measure Bg

degree_global = degree(Gg)

#
tan_fi = degree_global/degree_local
tan_fi

degree_modular_modulus = sqrt((degree_local**2) + (degree_global**2))
sort(degree_modular_modulus)

degree_local[1:5]
degree_global[1:5]

V(g)$dmm = degree_modular_modulus
V(g)$infomap = membership(modules_infomap)

?arrange
by_info = igraph::get.data.frame(g, "vertices")%>%
  dplyr::group_by(infomap)

by_info%>%
  filter(dmm==max(dmm))%>%
  arrange(infomap)

by_info%>%
  top_n(1, wt = dmm)%>%
      arrange(infomap, by_group = TRUE)

by_info = igraph::get.data.frame(g, "vertices")%>%
  arrange(infomap)%>%
  group_by(infomap)%>%
  top_n(1, wt = dmm)
by_info

by_info = igraph::get.data.frame(g, "vertices")%>%
  mutate(row_number()==1L)%>%
  arrange(infomap)%>%
  group_by(infomap)%>%
  top_n(1, wt = dmm)%>%
  top_n(-1, wt = id)

paste0("COMM_", by_info$name)

V(g)$pageRank = page.rank(g)$vector

by_pagerank = igraph::get.data.frame(g, "vertices")%>%
  mutate(row_number()==1L)%>%
  arrange(infomap)%>%
  group_by(infomap)%>%
  top_n(1, wt = pageRank)%>%
  top_n(-1, wt = id)

paste0("COMM_", by_pagerank$name)

modular_degree = function(g, commz){
  
  #calculate degree local and global
  GL = delete.edges(g, edges = E(g)[igraph::crossing(commz, g)])
  degree_local = degree(GL)
  GG = delete.edges(g, edges = E(g)[!(igraph::crossing(commz, g))])
  degree_global = degree(GG)
  
  #make output dataframe
  r = data.frame(degree_local  = degree_local,
                 degree_global = degree_global,
                 degree_modular_modulus = sqrt((degree_local**2) + (degree_global**2)),
                 degree_modular_tan     = degree_global/degree_local
                 
  )
  
  return(r)
}

modular_degree(g, modules_infomap)


modular_betweenness = function(g, commz){
  
  #calculate betweenness local and global
  GL = delete.edges(g, edges = E(g)[igraph::crossing(commz, g)])
  betweenness_local = betweenness(GL)
  GG = delete.edges(g, edges = E(g)[!(igraph::crossing(commz, g))])
  betweenness_global = betweenness(GG)
  
  #make output dataframe
  r = data.frame(betweenness_local  = betweenness_local,
                 betweenness_global = betweenness_global,
                 betweenness_modular_modulus = sqrt((betweenness_local**2) + (betweenness_global**2)),
                 betweenness_modular_tan     = betweenness_global/betweenness_local
                 
  )
  #r[is.na(r)] <- 0
  return(r)
}

modular_betweenness(g, modules_infomap)
