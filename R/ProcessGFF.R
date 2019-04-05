## Extract mapping of genes to GO terms from a GFF file

ProcessGFF <-
function(gff.file, c.ranges, nc.ranges, updownregion, verbose){

  if(!(class(c.ranges) == "GRanges" && class(nc.ranges) == "GRanges")) stop ("SNPs must be stored as GRanges objects.")
  ############################################################################# 
  ##Read the GFF file and extract genes associated to GO terms

  if(verbose) { message(paste("Reading",gff.file)) }
  
  ## Using the readLines function is faster than read.delim
  regions <- readLines(gff.file)
  
  ## Keep only gene regions
  regions <- regions[grep("\tgene\t", regions, perl=T)] 
  
  ## Split lines by the tab-character "\t"
  splitRegions <- strsplit(regions, "\t")
  
  ## store everything between Ontology_term and the next ";" of the 9th column
  ## in the 'GO' variable
  GO <- gsub(".*Ontology_term=([^;]+).*", "\\1", sapply(splitRegions, "[", 9))
  
  ## remove lines that do not contain any GO term
  ##GO <- GO[grep("GO:[[:alnum:]]{7}", GO, perl=TRUE)]
 
  ## remove GO terms starting with "!" and SO-Terms
  GO <- gsub("(!GO:[[:alnum:]]{7})|(SO:[[:alnum:]]{7})", "", GO)
 
  ## Keep only lines associated to at least one GO term
  splitRegions<-splitRegions[grep("GO:[[:alnum:]]{7}", GO, perl=T)]
  rm(regions)
  
  ## Store those lines in a GRanges object with associated GO terms as metadata
  ranges_GO <- GRanges(
    seqnames=sapply(splitRegions, "[", 1),
    ranges=IRanges(
      as.numeric(sapply(splitRegions, "[", 4)),
      as.numeric(sapply(splitRegions, "[", 5))
    ),
    GO=GO[grep("GO:[[:alnum:]]{7}", GO, perl=TRUE)]
  )

  ## extend ranges according to updownregion parameter
  if(updownregion > 0){
    start(ranges_GO) <- start(ranges_GO) - updownregion
    end(ranges_GO) <- end(ranges_GO) + updownregion
  }

  ## Get mappings of SNPs and GO term regions (i.e. genes).
  ## The findOverlaps method is provided by the GenomicRanges package
  ## and returns a Hits object.
  regions2candidates <- findOverlaps(ranges_GO, c.ranges)
  regions2noncandidates <- findOverlaps(ranges_GO, nc.ranges)

  ## Store all GOTerms found in the GFF-File in the 'goterms' variable
  goterms <- unique(unlist(strsplit(GO, ",")))
  goterms <- goterms[goterms != ""]
  goterms <- goterms[grep("GO:[[:alnum:]]{7}", goterms, perl=TRUE)]
  
  ## Remove contribtes_to... GOTerms
  goterms <- goterms[grep("contributes", goterms, perl=TRUE, invert=TRUE)]
  rm(GO)
  
  if(length(goterms) == 0) stop("No GO terms found. Please check the format of your GFF file.")
  
  #############################################################################
  ## Store mappings of GO terms to genes and SNPs in a hash table

  go2region <- hash::hash()
  if(verbose) { message("Mapping GO terms to genomic regions\n", appendLF=FALSE) }
  progEnvir <- new.env()
  progIncr <- GetProgressBarIncrementer(verbose=verbose, environment=progEnvir, 
    min=0, max=length(goterms), style=3)
  
  for(g in goterms) {
    progIncr(progEnvir)
    goterm_lines <- grep(g,elementMetadata(ranges_GO)$GO)
    go2region[[g]] <- 
      list(
        "regions"=goterm_lines, 
        "candidates"=subjectHits(
          regions2candidates[queryHits(regions2candidates) %in% goterm_lines]),
        "noncandidates"=subjectHits(
          regions2noncandidates[queryHits(regions2noncandidates) %in% goterm_lines])
      )
  }
  rm(progEnvir)

  ## Return a list
  list(
    "goterms"=goterms, 
    "regions2candidates"=regions2candidates,
    "regions2noncandidates"=regions2noncandidates, 
    "ranges"=ranges_GO, 
    "go2region"=go2region
  )
}
