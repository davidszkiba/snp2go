## The snp2go function implements the SNP2GO method.
## This is the only function exported by the SNP2GO package,
## all other functions are private.
   
snp2go <-
function(gtf, gff, goFile, FDR=0.05, runs=100000, candidateSNPs, 
  noncandidateSNPs, extension=0, verbose=TRUE, min.regions=10) {

  
  #############################################################################
  ## Extract gene regions and mappings of genes to GO terms from input files
  
  ## Check which kind of input file is used (GFF or GTF + Mapping File)
  ## and import Gene <-> GOterm mappings using the
  ## "ProcessGFF" or "ProcessMartExport" function accordingly.
  if(!missing("gff")) {
    if(!missing("goFile")){
      warning("Parameter 'goFile' is ignored in GFF mode")
    }
    if(!missing("gtf")){
      warning("Paramter 'gtf' is ignored in GFF mode")
    }
    ## See below for a description of inData
    inData <- ProcessGFF(gff.file=gff, c.ranges=candidateSNPs, 
      nc.ranges=noncandidateSNPs, updownregion=extension, 
      verbose=verbose)
  } else if(!missing("gtf")) {
    if(!missing("gff")){
      warning("Parameter 'gff' is ignored in GTF mode")
    }
    ## See below for a description of inData
    inData <- ProcessMartExport(gtf.file=gtf, mart.file=goFile, 
      c.ranges=candidateSNPs, nc.ranges=noncandidateSNPs, 
      region.type="gene", updownregion=extension, verbose=verbose)
  } else {
    stop("No GTF or GFF file specified")
}

  ## The 'inData' variable is used to store data extracted/calculated from input
  ## files (GO terms, location of GO term regions (i.e. genes or gene + given 
  ## number of nucleotides upstream/downstream associated with a GO 
  ## term), overlaps between SNPs and GO term regions, genes/SNPs associated
  ## with GO terms)
  ## 'inData' is a named list variable: 
  ##   list(goterms, ranges, regions2candidates, regions2noncandidates, go2region) 
  ## 
  ## + The 'goterms' elemet stores the total set of GO terms used
  ## 
  ## + The 'regions2candidates' and 'regions2noncandidates' elements store
  ##   mappings of GO term regions to candidate and non-candidate SNPs, 
  ##   respectively. They are "Hits" Objects as returned by the findOverlaps 
  ##     method of the GenomicRanges package.
  ##
  ## + The ranges variable contains a GRanges object of the GO term regions 
  ##   (i.e. genes annotated with at least one GO term).
  ##
  ## + go2region stores for each  GO terms a list variable: 
  ##   list(regions, candidates, noncandidates).
  ##   "regions", "candidates", and "noncandidates" are integer vectors 
  ##   containing the indices of the genomic regions, candidate SNPs and
  ##   non-candidate SNPs associated with the respective GO term.
  ##   E.g. x <- go2region[["GO:0000001"]][["candidates"]] returns the indices 
  ##   of all candidate SNPs associated with GO:0000001, and candidateSNPs[x] 
  ##   returns those SNPs as a GenomicRanges object.
  ##   noncandidateSNPs[go2region[["GO:0000001"]][["noncandidates"]]] returns all
  ##   non-candidate SNPs associated with GO:0000001.
  ##   ranges[go2region[["GO:0000001"]][["regions"]]] returns all regions
  ##   associated with GO:0000001.
  
  ## We store the number of candidate- and non-candidate SNPs,
  ## as well as the number of informative SNPs (i.e. SNPs in GO term
  ## regions) just to provide this information.
  numSNPs <- list()
  numInformativeSNPs <- list()
  numSNPs[["candidates"]] <- length(candidateSNPs)
  numInformativeSNPs[["candidates"]] <- 
    length(
      unique(subjectHits(inData$regions2candidates))
    )
  numSNPs[["noncandidates"]] <- length(noncandidateSNPs)
  numInformativeSNPs[["noncandidates"]] <- 
    length(
      unique(subjectHits(inData$regions2noncandidates))
    )

  
  #############################################################################
  ## Initialize data structures representing the directed acyclic graphs (DAGs) 
  ## as lists
  
  ## For each of the three ontologies "Molecular Function", 
  ## "Biological Process", and "Cellular Component", we construct a list that
  ## contains the GO terms associated with the different levels. For the 
  ##inclusive analysis, those lists do not only contain the GO term associated 
  ## with each level, but also the offspring terms of the associated GO terms.
  ## The lists of the 3 ontologies are gotermsByLevel$MF, gotermsByLevel$BP 
  ## and gotermsByLevel$CC.
  ## E.g., gotermsByLevel$MF[[1]] contains all GO terms found at any 
  ## level >= 1 of the MF ontology.
  ## gotermsByLevel$MF[[2]] contains all GO terms found at any level >= 2 of the
  ## MF ontology, and so on.
  if(verbose) { message("Organizing GO categories by level") }
  gotermsByLevel <- ComputeLevels(inData$goterms)
  termlevel <- gotermsByLevel[["termlevel"]]
  ## To speed up the calculation of the number of SNPs associated with 
  ## background GO terms, we create lists that contain the _number_ of SNPs 
  ## associated with GO terms from the different levels of the previously 
  ## created gotermsByLevel lists.
  ## E.g. snpsByLevel$MF$LengthCandidatesByLevel[[2]] returns the
  ## number of SNPs associated with gotermsByLevel$MF[[2]]
  if(verbose) {
    message("Calculating number of candidate SNPs and non-candidate SNPs\
associated with GO terms from different ontology levels:")
  }
  snpsByLevel <- list()
  for(onto in c("MF", "BP", "CC")) {
  snpsByLevel[[onto]] <- 
    ComputeSNPsByLevel(
      gotermsByLevel=gotermsByLevel[[onto]], 
      gotermMappings=inData$go2region, 
      numInformativeSNPs=numInformativeSNPs,
      ontoName=onto,
      verbose=verbose
    )
  }


  #############################################################################
  ## Update mappings of GO terms for the inclusive analysis

  ## Due to the inclusive analysis, we need to update the mappings of GO terms 
  ## to regions, candidate SNPs and non-candidate SNPs as found in go2region, so
  ## that each GO term contains mappings to all regions and SNPs associated with
  ## any of its offspring terms.

  goterm2ranges <- list()
  goterm2ranges[["regions"]] <- hash::hash()
  goterm2ranges[["candidates"]] <- hash::hash()
  goterm2ranges[["noncandidates"]] <- hash::hash()

  if(verbose) { 
    message(
      "Updating mappings of GO terms to SNPs and genomic regions\
for inclusive analysis"
    ) 
  }
  progEnvir <- new.env()
  progIncr <-
    GetProgressBarIncrementer(
      verbose=verbose,
      environment=progEnvir,
      min=0,
      max=length(inData$goterms),
      style=3
    )


  for(g in inData$goterms){
    progIncr(progEnvir)
    x <- GetAllRegionsAndSNPs(g, inData$goterms, inData$go2region)
    goterm2ranges[["regions"]][[g]] <- x[["allCandidateRegions"]]
    goterm2ranges[["candidates"]][[g]] <- x[["allCandidateSNPs"]]
    goterm2ranges[["noncandidates"]][[g]] <- x[["allNoncandidateSNPs"]]
  }
  rm(progEnvir)
  
  #############################################################################
  ## Test GO terms for over-representation of candidate SNPs

  ## The progEnvir environment is needed for the progress bar
  if(verbose) { 
    message("Testing GO terms for over-representation of candidate SNPs") 
  }
  progEnvir <- new.env()
  progIncr <- 
    GetProgressBarIncrementer(
      verbose=verbose, 
      environment=progEnvir, 
      min=0, 
      max=length(inData$goterms),
      style=3
    )
  
  ## In the lapply function, the GetStatistic function is invoked for each GO
  ## term and the results are stored in the enrichmentPvals variable 
  ## which is a list.
  
  enrichmentPvals <- lapply(inData$goterms, function(g) {
    progIncr(progEnvir)  ## Increment the progress bar
    ## Return NA if the GO term is deprecated or has a synonymous term
    if(is.na(goterm2ranges[["regions"]][[g]][1])) { 
      return(NA)
    }
    ## If GO term is not deprecated, return the result of GetStatistic
    GetStatistic(
      goterm=g, 
      gotermsByLevel=gotermsByLevel, 
      snpsByLevel=snpsByLevel, 
      cc=length(goterm2ranges[["candidates"]][[g]]), 
      cn=length(goterm2ranges[["noncandidates"]][[g]])
    )
  }) #End of lapply
  rm(progEnvir)

  ## Keep only GO terms for which we got a meaningful result (i.e. 
  ## non-deprecated GO terms)
  enrichmentPvals <- enrichmentPvals[!is.na(enrichmentPvals)]

  ## Calculate FDR adjusted p-values according to Benjamini-Hochberg method
  ## enrichmentPvals is a list(goterm, pvalue)
  ## 
  enrichmentPvalsAdjusted <- 
    p.adjust(unlist(sapply(enrichmentPvals, "[", 2)), method="BH")
  isSignificant <- enrichmentPvalsAdjusted < FDR
  signEnrichmentPvals <- enrichmentPvals[isSignificant]
  signGOTerms <- unlist(sapply(signEnrichmentPvals, "[", 1))
  signEnrichmentPvalsAdjusted <- enrichmentPvalsAdjusted[isSignificant]

  #############################################################################
  ## Test for local effects (e.g. spatial clustering of candidate SNPs)


  ##Perform sampling procedure only for GO terms having at least minRegions
  ## hasMinRegions is a boolean vector with the same length as signGOTerms:
  ## TRUE means that the corresponding GO term has more than 'minRegions'
  ## _non-overlapping_ (i.e. reduced) regions/genes.

  #TODO: change parameter name from 'min.regions' to minRegions
  # and update documentation
  minRegions <- min.regions

  hasMinRegions <- unlist(sapply(signGOTerms,
    function(g){
      length(reduce(inData$ranges[ goterm2ranges[["regions"]][[g]] ])) >= minRegions
    }))
  signGOTerms <- signGOTerms[hasMinRegions]
  signEnrichmentPvals <- signEnrichmentPvals[hasMinRegions]
  signEnrichmentPvalsAdjusted <- signEnrichmentPvalsAdjusted[hasMinRegions]

  if(length(signGOTerms) == 0) {
    if(verbose){ message("No significant GO terms found") }
    result <- data.frame(
      "GO"=character(0), 
      "P"=numeric(0), 
      "FDR"=numeric(0), 
      "p.L"=numeric(0), 
      "p.G"=numeric(0), 
      "g"=numeric(0), 
      "G"=numeric(0), 
      "nc"=numeric(0), 
      "mc"=numeric(0), 
      "GO.def"=character(0), 
      "child.GOs"=character(0)
    )
  } 
  else {
    ## In the lapply function, the SamplingTest function is invoked for each GO
    ## term and the results are stored in the samplingValues variable.
    ## ProgEnvir is again needed for the progress bar.
    if(verbose){ 
      message(paste("Found", length(signGOTerms), "significant GO terms."))
      message("Testing significant GO terms for local effects")
    }
      
      progEnvir <- new.env()
      progIncr <- 
      GetProgressBarIncrementer(
        verbose=verbose, 
        environment=progEnvir, 
        min=0, 
        max=length(signGOTerms), 
        style=3
      )

    samplingValues <-lapply(signGOTerms,
      function(sig){
        progIncr(progEnvir)
        SamplingTest(
          goterm=sig,   
          iterations=runs, 
          ranges=inData$ranges, 
          c.ranges=candidateSNPs, 
          nc.ranges=noncandidateSNPs, 
          regions2candidates=inData$regions2candidates, 
          regions2noncandidates=inData$regions2noncandidates, 
          go2allRegions=goterm2ranges[["regions"]], 
          go2allCandidateSNPs=goterm2ranges[["candidates"]], 
          go2allNoncandidateSNPs=goterm2ranges[["noncandidates"]]
        )
      }
    ) # End of lapply
    #SamplingTest returns list(pvalLower,pvalGreater, number of candidate regions)
    rm(progEnvir)

  
    #############################################################################
    ## Return the results as a list
    ## A data.frame containing information about significant GO terms is returned
    ## as part of this list
    ## Other variables of the list contain parameters used for the analysis.
  
    ## Construct a data.frame:
    ## Columns: goterm, pval, pval_adjusted, pval-lower, pval-greater, 
    ## numCandidateRegions, numTotalRegions, numCandidates, numNonCandidates,  
    ## definition, offspring
  
    ## get the pval-lower, pval-greater and numCandidateRegions data
    samplingDataFrame <- 
      data.frame(
        matrix(unlist(samplingValues), byrow=TRUE, nrow=length(samplingValues))
      )
    colnames(samplingDataFrame) <- c("p.L", "p.G", "g")
  
    ## get numTotalRegions, numCandidates, numNoncandidates, definition and 
    ## offspring data
    signData <- 
      lapply(signGOTerms,
        function(g) {
          list(
            "numTotalRegions"=
              length(reduce(inData$ranges[goterm2ranges[["regions"]][[g]]])),
            "numCandidateSNPs"=length(unique(goterm2ranges[["candidates"]][[g]])),
            "numNoncandidateSNPs"=length(unique(goterm2ranges[["noncandidates"]][[g]])),
            "definition"=Definition(g),
            "significantOffspringTerms"=
              AddOffspringInformation(g, significantGOTerms=signGOTerms)
          )
        }
      )

    signDataFrame <- 
      matrix(unlist(signData), byrow=TRUE, nrow=length(signGOTerms))
    colnames(signDataFrame) <- c("G", "nc", "mc", "GO.def", "child.GOs")
   
    ## create the data.frame by sticking together samplingDataFrame, signDataFrame,
    ## GO term names ("GO"), p-values ("P") and FDR ("FDR")
    result <- 
      data.frame(
        "GO"=signGOTerms, 
        "P"=unlist(sapply(signEnrichmentPvals, "[", 2)), 
        "FDR"=unlist(signEnrichmentPvalsAdjusted), 
        samplingDataFrame, 
        signDataFrame)

    ## The columns "G", "nc" and "mc" are stored as factors but I think it is
    ## better to store them as numeric values, so I convert them:
    for(col in c("G", "nc", "mc")) {
      result[[col]] <- as.numeric( levels(result[[col]])) [result[[col]]]
    }
  } # end else

  ## Return the list
  ## The data frame is stored in the "enriched" element
  return(
    list(
      "FDR"=FDR,
      "runs"=runs,
      "extension"=extension,
      "candidateSNPs"=numSNPs[["candidates"]],
      "noncandidateSNPs"=numSNPs[["noncandidates"]],
      "informative.candidateSNPs"=numInformativeSNPs[["candidates"]],
      "informative.noncandidateSNPs"=numInformativeSNPs[["noncandidates"]],
      "enriched"=result,
      "goterms"=inData$goterms,
      "regions"=inData$ranges,
      "go2ranges"=goterm2ranges,
      "termlevel"=termlevel
    )
  )
}
