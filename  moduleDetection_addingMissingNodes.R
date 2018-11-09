############
#Unless this is justified by the biological experimental setting
#adding missing nodes to a network to make them comparable 
#is not necessarily a good idea (tm)
############

#read graphs 
g     = read.graph(file = "results/example_01.gml", "gml") #main graph, from merging five known pathways
g_alt = read.graph(file = "results/example_02.gml", "gml") #alternative graph, from five different known pathways
g_mix = read.graph(file = "results/example_03.gml", "gml") #graph from two pathways used in g, and two pathways in g0
g_rew = read.graph(file = "results/example_06.gml", "gml") #rewired g (using each_edge(prob=0.25))

#get a list of all nodes to add them to networks 


all_nodes =               unique(c(vertex.attributes(g)[["name"]],
                               vertex.attributes(g_alt)[["name"]],
                               vertex.attributes(g_mix)[["name"]],
                               vertex.attributes(g_rew)[["name"]]
                               )
                             )
g_ori_plus = add.vertices(g, 
                          nv = length(setdiff(all_nodes, 
                                              vertex.attributes(g)[["name"]])),
                          attr = list(name = setdiff(all_nodes, 
                                                     vertex.attributes(g)[["name"]]))
                            )

g_alt_plus = add.vertices(g_alt, 
                          nv = length(setdiff(all_nodes, 
                                              vertex.attributes(g_alt)[["name"]])),
                          attr = list(name = setdiff(all_nodes, 
                                                     vertex.attributes(g_alt)[["name"]]))
)

g_mix_plus = add.vertices(g_mix, 
                          nv = length(setdiff(all_nodes, 
                                              vertex.attributes(g_mix)[["name"]])),
                          attr = list(name = setdiff(all_nodes, 
                                                     vertex.attributes(g_mix)[["name"]]))
)

g_rew_plus = add.vertices(g_rew, 
                          nv = length(setdiff(all_nodes, 
                                              vertex.attributes(g_rew)[["name"]])),
                          attr = list(name = setdiff(all_nodes, 
                                                     vertex.attributes(g_rew)[["name"]]))
)

g_ori_plus
g_alt_plus
g_mix_plus
g_rew_plus

#detect modules with Infomap

modules_g_ori_plus = infomap.community(g_ori_plus,     nb.trials = 1000)
modules_g_alt_plus = infomap.community(g_alt_plus,     nb.trials = 1000)
modules_g_mix_plus = infomap.community(g_mix_plus,     nb.trials = 1000)
modules_g_rew_plus = infomap.community(g_rew_plus,     nb.trials = 1000)


igraph::compare(comm1 = modules_g_ori_plus, 
                comm2 = modules_g_alt_plus, 
                method = "nmi")

igraph::compare(comm1 = modules_g_ori_plus, 
                comm2 = modules_g_mix_plus, 
                method = "nmi")

igraph::compare(comm1 = modules_g_ori_plus, 
                comm2 = modules_g_rew_plus, 
                method = "nmi")

igraph::compare(comm1 = modules_g_ori_plus, 
                comm2 = modules_g_rew_plus, 
                method = "nmi")
