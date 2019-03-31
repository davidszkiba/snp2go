##This function returns all offspring terms of its first parameter 'goterm'
## that are found in the list of goterms provided by the second parameter 
## 'significantGOTerms' as a comma-seperated list.
## The function is used to indicate for each significant GO term the significant
## offspring terms.

AddOffspringInformation <-
function(goterm,significantGOTerms) {
  
  offspring <- switch(Ontology(goterm), 
    MF=AnnotationDbi::get(goterm, GOMFOFFSPRING),
    BP=AnnotationDbi::get(goterm, GOBPOFFSPRING),
    CC=AnnotationDbi::get(goterm, GOCCOFFSPRING))
    
  offspring<-offspring[offspring %in% significantGOTerms]
  
  ## Return string
  paste(offspring,collapse=",")
}
