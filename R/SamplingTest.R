## This function samples SNPs to produce a null distribution of the
## number of regions containing c SNPs, where c is the number of
## candidate SNPs associated with the respective GO term. The SNPs are
## sampled from all SNPs (candidates + non-candidates) associated with
## the respective GO term.
## The number of regions, g, containing c candidate SNPs is calculated
## and the proportion of iterations containing the c sampled SNPs
## in less (p.L) or more (p.G) regions than g is returned.

SamplingTest <-
function(goterm, iterations, ranges, c.ranges, nc.ranges, regions2candidates, 
  regions2noncandidates, go2allRegions, go2allCandidateSNPs, 
  go2allNoncandidateSNPs) {
	
  ## get the number of candidate SNPs, c
  c <- length(
    unique(
      go2allCandidateSNPs[[goterm]]
    )
  )
  
  ## get the number of genes, g, containing the candidate SNPs
  g <- length(
    reduce(ranges[queryHits(
      regions2candidates[
        subjectHits(regions2candidates) %in% go2allCandidateSNPs[[goterm]]
      ]
    )])
  )
  
  ## get the set of SNPs (allSNPs) from which c SNPs are sampled
  if(c == 0){
    allSNPs <- nc.ranges[go2allNoncandidateSNPs[[goterm]]]
  }
  else if(length(unique(go2allNoncandidateSNPs[[goterm]]))==0) {
    allSNPs <- c.ranges[go2allCandidateSNPs[[goterm]]]
  }
  else{
    allSNPs <- c(c.ranges[go2allCandidateSNPs[[goterm]]], 
      nc.ranges[go2allNoncandidateSNPs[[goterm]]])
  }

  ## get overlaps of SNPs and genes
  overlappingSNPs <- subjectHits(findOverlaps(allSNPs, 
    ranges[go2allRegions[[goterm]]]))
  ## preallocate a vector to store the number of regions found in each iteration
  numRegions <- vector("integer",iterations)

  ## In each iteration, sample c SNPs from all SNPs and count
  ## the number (length) of distinct (unique) regions those SNPs are
  ## overlapping with.
  numSNPs <- length(allSNPs)
  numRegions <- sapply(
    1:iterations,  ## apply the next lines 'iteration' times 
    function(x) 
      length(unique(overlappingSNPs[sample(1:numSNPs, c)]))
  )
 
  ## If the sampled SNPs were always found in the same number of genes
  ## set p.L and p.G to NA, 
  ## else, get the proportion of iterations in which the sampled SNPs
  ## were found in <= (p.L) or >= (p.G) g regions.
  if(length(unique(numRegions)) == 1) {
    lower <- NA;
    bigger <- NA;
  }
  else{
    lower <- sum(numRegions <= g) / iterations
    bigger <- sum(numRegions >= g) / iterations
  }

  ## Return list...
  list(lower, bigger, g)
}
