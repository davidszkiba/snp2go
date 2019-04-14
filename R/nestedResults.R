#' Creates a data.frame that sorts enriched GO terms by domain and level.
#'
#' For each enriched GO term, the domain and level inside of the domain is added.
#' Additionally, GO terms are sorted by their level inside the GO DAG and nested.
#'
#' Nesting of GO terms means, that all child terms of a GO term at level L are
#' displayed, before the next significant GO term of level L is shown.
#' This is best illustrated by an example:\cr
#' \code{
#' GO:1 level:1 childterms: GO:2, GO:3, GO:4 \cr
#' GO:2 level:2 childterms: GO:4 \cr
#' GO:3 level:2 childterms: GO:4 \cr
#' GO:4 level:3 childterms: --- \cr
#' GO:5: level:1 childterms: GO:6, GO:7 \cr
#' }
#' GO term GO:1 has 3 significant child terms (GO:2, GO:3 and GO:4). The next
#' GO term at the same level as GO:1, GO:5, is displayed after all child terms
#' of GO:1 were displayed.
#' Additionally, GO terms are sorted by their level inside the GO DAG and
#' nested.
#'
#' @param snp2goResult The list returned by the snp2go function.
#' @return Returns a data.frame.
#' @examples
#'   c <- GenomicRanges::GRanges(
#'       seqnames=rep(1, 100),
#'       ranges=IRanges(
#'           rep(1, 100),
#'           runif(min=2, max=1000, n=100)
#'       )
#'   )
#'   nc <- GenomicRanges::GRanges(
#'       seqnames=rep(1, 1000),
#'       ranges=IRanges(
#'           rep(1, 100),
#'           runif(min=2, max=1000, n=1000)
#'       )
#'   )
#'   \dontrun{x <- snp2go(gff="dmel.gff", candidateSNPs=c, noncandidateSNPs=nc)}
#'   \dontrun{nested <- nestedResults(x)}
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

