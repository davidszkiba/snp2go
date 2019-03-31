## This function stores the GO structure in a list variable. The level of a GO
## term is calculated as the longest path to the root node. In the returned 
## list, all GO terms associated with level i are stored at index i.
## Note: A GO term is also associated with level i, if one of its ancestors
## is located at level i.
ComputeLevels <-
function(goterms) {
	
  count <- 1	
  GO_Levels <- list()
  termlevel <- list()
  
  ## Generate a list for each of the three ontologies:
  ## MF - Molecular Function
  ## BP - Biological Process
  ## CC - Cellular Component
  for(onto in c("MF", "BP", "CC")) {
    Levels <- list() 
    
    repeat {
      levelTerms <- getGOLevel(onto,count)
      if(length(levelTerms) == 0)  {
        count <- 1
        break
      }
      Levels[[onto]]<-c(Levels[[onto]], list(levelTerms))
      count <- count + 1	
    }


    ## process Levels, so that each GO-Term appears at the most specific Level
    ## E.g. if a GO term was previously associated with levels 2, 3 and 4, it 
    ## is assigned to level 4.
    Levels[[onto]] <- RemoveFromUpperLevels(Levels[[onto]])


    ## For each Level, remove GO terms that are not found in the GFF- 
    ## or gene-association file 
    for(i in 1:length(Levels[[onto]])){
      Levels[[onto]][[i]] <- Levels[[onto]][[i]][Levels[[onto]][[i]] %in% goterms]
    }


  #level 1 should contain GO-Terms from level 1 up to the last level
  #level 2 should contain GO-Terms from level 2 up to the last level
  #and so on
  oo <- paste(onto,".new",sep="")
  Levels[[oo]]<-list()
  for(i in 1:length(Levels[[onto]])) {
    Levels[[oo]][[i]]<-character()
      for(j in i:length(Levels[[onto]])) {
        Levels[[oo]][[i]] <- c(Levels[[oo]][[i]], Levels[[onto]][[j]])	
      }
    Levels[[oo]][[i]] <- unique(Levels[[oo]][[i]])

    ## Set the level of the GO terms in this level.
    termlevel[ Levels[[oo]][[i]] ]  <- i
  }

  GO_Levels[[onto]] <- Levels[[oo]]
  }
  
  return(
    list(
      "MF"=GO_Levels[["MF"]], 
      "CC"=GO_Levels[["CC"]], 
      "BP"=GO_Levels[["BP"]],
      "termlevel"=termlevel
    )
  )
}
