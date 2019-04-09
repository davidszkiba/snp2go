nestedResults <- function(snp2goResult) {
  dd <- snp2goResult$enriched
  
  if(length(dd$GO) == 0) {
    warning("Stopping, because the result does not contain any significant GO Terms.")
    return(NULL)
  }

  ## Add GO Term ontology and level and order data frame by Ontology and level.
  dd$Level <- as.character(snp2goResult$termlevel[as.character(dd$GO)])
  dd$Ontology <- Ontology(as.character(dd$GO))  
  dd <- dd[with(dd, order(dd$Ontology, as.numeric(dd$Level))), ]

  sorted_rows <- c()
  processed <- list()

  ## Reorder the data frame, so that each GO term is followed by its significant
  ## child terms.
  for(i in 1:length(dd$GO)) {
    if(dd$GO[i] %in% names(processed)) {
      next;
    }
    child_terms <- strsplit(as.character(dd[i, c("child.GOs")]), ",")[[1]]
    child_term_rows <- which(sapply(dd$GO, function(goterm) (goterm %in% child_terms)))
    processed[child_terms] <- TRUE
    sorted_rows <- c(sorted_rows, i, child_term_rows)
  }
  return(dd[sorted_rows, c(12, 13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)])
}

