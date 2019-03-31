## This function is used for the inclusive analysis to update the 
## mapping of a GO term to its regions.
## It assigns to each GO term all the regions and SNPs associated
## with any of this GO term's offspring terms.

GetAllRegionsAndSNPs <- function(goterm, allGOTerms, go2region){
  
  ## candidateTerms are all offspring terms that are found in the input file
  ## and the given GO term itself
  candidateTerms <- GetCandidateTerms(goterm, allGOTerms)
  if(is.na(candidateTerms[1])) { ## Term is obsolete or has a synonym
    return(
      list(
        "allCandidateRegions"=NA, 
        "allCandidateSNPs"=NA, 
        "allNoncandidateSNPs"=NA
      )
    )
  }
  else{
    ## store the indices of regions/genes, candidate SNPs and non-candidate SNPs
    ## associated with the given GO term. These indices can then be used to
    ## access the regions (ranges[indices[["regions"]]]), candidate SNPs 
    ## c.ranges[indices[["candidates"]]], or non-candidate SNPs 
    ## nc.ranges[indices[["noncandidates"]]] associated with the given GO term. 
    ## -> 'indices' is a list with three named elements:
    ##    list(regions, candidates, noncandidates)
    ## We pre-allocate a vector of size 1024 and increase the vector size if 
    ## needed. See chapter 2 of Patrick Burns' "The R Inferno". 
    ## Additionally, 'upper' points to the number of indices already 
    ## stored in the vector. This may seem a bit complicated, but prevents
    ## memory fragmentation as it would be caused by appending the vector to 
    ## itself (e.g.: regions <- c(regions, newRegions) ) 
 
    indices <- list() ## will contain the three vectors of indices
    vlen <- 1024 ## starting vector length for regions vector
    
    for(type in c("regions", "candidates", "noncandidates")) {
      ## store the indices in the 'regions' variable
      regions <- vector("integer", vlen)
      upper <- 0 ## number of indices stored
      for(ct in candidateTerms){
        ind <- go2region[[ct]][[type]]
        if(length(ind) > 0) {
          lower <- upper + 1
          upper <- lower + length(ind) - 1
          if(upper > length(regions)) {
            tmp <- regions  ## temporarily store old regions
            regions <- vector("integer", upper * 2) 
            regions[1:length(tmp)] <- tmp
          }
          regions[lower:upper] <- ind
        }
      }
     
      if(upper > 0) {  ## indices were assigned to regions
        indices[[type]] <- unique(regions[1:upper])
      } else {  ## if no indices were assigned, return empty vector
        indices[[type]] <- vector("integer")
      }
    }

    ## Return a list		
    list(
      "allCandidateRegions"=indices[["regions"]], 
      "allCandidateSNPs"=indices[["candidates"]],
      "allNoncandidateSNPs"=indices[["noncandidates"]]
    )
  } # end else
}
