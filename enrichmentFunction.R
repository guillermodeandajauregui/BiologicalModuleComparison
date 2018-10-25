library(HTSanalyzeR)
library(tidyverse)


ModuleEnrichment = function(g, moduleList, pathwayList, threshold=0.05){
  
  enrichmentResults = lapply(X = seq_along(moduleList), FUN = function(i){
    
    
    nomen = names(moduleList)[i]
    
    my_enrichment = HTSanalyzeR::multiHyperGeoTest(collectionOfGeneSets = pathwayList, 
                                                   universe = V(g)$name, 
                                                   hits = moduleList[[i]], 
                                                   minGeneSetSize = 1, 
                                                   pAdjustMethod = "BH", 
                                                   verbose = TRUE
    )
    print(nrow(my_enrichment))
    my_2 = tibble::rownames_to_column(as.data.frame(my_enrichment))#%>%filter(Adjusted.Pvalue<10)
    my_3 = cbind(comm = nomen, my_2)
    return(my_3) 
  })
  
  enrichmentResults = ldply(enrichmentResults, data.frame)
  enrichmentResults$Adjusted.Pvalue = p.adjust(p = enrichmentResults$Pvalue, method = "BH")
  enrichmentResults = filter(enrichmentResults, Adjusted.Pvalue<threshold)
  
  return(enrichmentResults)
}



