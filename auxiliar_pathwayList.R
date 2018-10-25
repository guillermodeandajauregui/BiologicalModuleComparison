library(graphite)

#make a list of KEGG pathways
x = graphite::pathways("hsapiens", "kegg")

#pathways of interest

pathway_names =         c("Hedgehog signaling pathway",#main
                          "NF-kappa B signaling pathway",#main
                          "VEGF signaling pathway",#main
                          "p53 signaling pathway",#main
                          "RIG-I-like receptor signaling pathway",#main
                          "Pentose phosphate pathway",#alt
                          "Notch signaling pathway",#alt
                          "mRNA surveillance pathway",#alt
                          "TGF-beta signaling pathway",#alt
                          "IL-17 signaling pathway",#alt
                          "Sphingolipid signaling pathway",#extra
                          "Insulin signaling pathway"#extra
)

x = x[pathway_names]

#convert ids to symbol
x = convertIdentifiers(x, "symbol")

#make a list of nodes in pathways
#also remove the SYMBOL prefix (not necessary now)
my_pathways = lapply(X = x, FUN = function(i){
  tmp = nodes(i)
  tmp = gsub(pattern = "SYMBOL:", 
             replacement = "", 
             x = tmp)
  return(tmp)
})

my_pathways
rm(x)
