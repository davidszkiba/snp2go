## This function returns the level (i.e. longest path to root node) of a
## GO term 'goterm' in the ontology 'ontology'.
## It iterates through the list from bottom to top to get the longest path.

FindLevel <-
function(goterm, ontology){
  for(i in length(ontology):1){
    if(goterm %in% ontology[[i]]){
      return(i)
    }
  }	
}

