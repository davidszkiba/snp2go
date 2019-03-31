## This function returns for a given GO term all the
## offspring terms that are found in the input file

GetCandidateTerms <-
function(goterm, allGOTerms) {
  if(is.na(Ontology(goterm))) {	#no Ontology for this term is found. This can happen if the term is obsolete or has a synonymous term.
    return(NA)
  }
  else{
    offspringTerms <- switch(Ontology(goterm), 
      MF=AnnotationDbi::get(goterm, GOMFOFFSPRING), 
      BP=AnnotationDbi::get(goterm, GOBPOFFSPRING),
      CC=AnnotationDbi::get(goterm, GOCCOFFSPRING))
    if(is.na(offspringTerms[1])) {
      return(goterm)
    } else {
      offspringTerms <- offspringTerms[offspringTerms %in% allGOTerms]
      return(c(goterm, offspringTerms))
    }
  }
}
