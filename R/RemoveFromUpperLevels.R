## This functions modifies lists that contain GO terms of the ontologies
## ordered by level. 
## It returns those lists with each GO term associated to its deepest possible
## level. E.g. if a GO term can be found at levels 2, 3 and 4, the returned 
## lists contains the GO term at level 4.

RemoveFromUpperLevels <-
function(levels){
  for(i in 1:(length(levels)-1)){
    for(j in (i+1):length(levels)){
      levels[[i]] <- 
        levels[[i]][(levels[[i]] %in% levels[[j]]) == FALSE]
    }
  }	
  ## Return the levels
  levels
}

