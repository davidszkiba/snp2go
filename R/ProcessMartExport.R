ProcessMartExport <-
function(gtf.file, mart.file, region.type="gene", updownregion=0, c.ranges, 
  nc.ranges, verbose) {
  
  if(!(class(c.ranges) == "GRanges" && class(nc.ranges) == "GRanges")) stop ("SNPs must be stored as GRanges objects.")

  #############################################################################
  ## Read the GTF file and extract gene regions
  
  ## Keep only lines corresponding to exons
  if(verbose) { message(paste("Processing", gtf.file)) }	
  gtf <- read.delim(gtf.file, comment.char="#", header=F)
  gtf <- gtf[gtf$V3=="exon", c(1, 4, 5, 9)]

  ##Get the gene-ID: the attribute column (9th column) should contain only 
  ## the gene-ID
  gtf$V9 <- as.character(gtf$V9)
  ## Remove everything but the gene ID from the 9th column
  gtf$V9 <- gsub(".*gene_id\\s([^;]+);.*$", "\\1", gtf$V9)

  ## Extract the regions according to the specified region.type (which is GENE 
  ## by default)
  if(toupper(region.type) == "EXON"){
    ranges <- 
      GRanges(
        seqnames=gtf$V1, 
        ranges=IRanges(gtf$V4, gtf$V5), 
        ID=gtf$V9
      )
  }	
  else if(toupper(region.type) == "GENE"){
    a <- split(gtf, gtf$V9)
    chromosomes <- sapply(a, "[", 1)
    gene_starts <- sapply(a, function(x) min(x$V4))
    gene_ends <- sapply(a, function(x) max(x$V5))
    IDs <- sapply(a, "[", 4) 
    ranges <- 
      GRanges(
        seqnames=sapply(chromosomes, "[", 1), 
        ranges=IRanges(gene_starts, gene_ends), 
        ID=sapply(IDs, "[", 1)
      )
  }
  else {
  stop(paste("Unknown region type:", region.type))
  }
  rm(gtf)

  ## Extend the GO term regions by the specified nucleotide number
  if(updownregion > 0) {
  start(ranges) <- start(ranges) - updownregion
  end(ranges) <- end(ranges) + updownregion
  }

  #############################################################################
  ## Read the mapping file provided by biomart, which contains gene IDs in the
  ## first column and GO terms in the second column
  
  if(verbose) { message(paste("Processing", mart.file)) }
  mart.file <- read.delim(mart.file, header=T)
  mart.file[, 1] <- as.vector(mart.file[, 1])
  ## Keep only lines containing a GO term in the second column:
  mart.file <- mart.file[mart.file[, 2] != "", ]
  ## Keep only lines with gene IDs that are also found in the GTF file:
  mart.file <- mart.file[mart.file[, 1] %in% unlist(IDs), ]
  if (dim(mart.file)[1] == 0) stop("No GO mappings for the gene IDs of the GTF file found. Please check the formats of the GTF and mapping files.");
  ## This creates a mapping of GO terms to corresponding gene IDs
  go2id <- split(mart.file[, 1], mart.file[, 2])

  ## Keep only ranges from the GTF file that are mapped to a GO term:
  ranges <- ranges[which(elementMetadata(ranges)$ID %in% mart.file[, 1])]
  ## Store non-empty goterms in a variable
  goterms <- unique(mart.file[, 2])
  goterms <- goterms[goterms != ""]

  ## Find overlaps between GO term regions (i.e. genes) and SNPs
  ## The findOverlaps function is provided by the GenomicRanges package and
  ## returns a 'Hits' object.
  regions2candidates <- findOverlaps(ranges, c.ranges)
  regions2noncandidates <- findOverlaps(ranges, nc.ranges)
  
  #############################################################################
  ## Construct a hash table that contains mappings of GO terms to genes, 
  ## candidate SNPs and non-candidate SNPs.
  if(verbose) { message("Mapping GO terms to genomic regions\n", appendLF=FALSE) }
  
  ## progEnvir and progIncr are used to show the progress with a progress bar
  progEnvir <- new.env()
  progIncr <- 
    GetProgressBarIncrementer(
      verbose=verbose, 
      environment=progEnvir, 
      min=0, 
      max=length(goterms), 
      style=3
    )
  ## go2regions contains the mappings of GO terms to genes and SNPs
  go2region <- hash()
  for(g in goterms) {	
    progIncr(progEnvir)
    g_lines <- which((elementMetadata(ranges)$ID %in% go2id[[g]] == TRUE) == TRUE)
    go2region[[g]] <- list("regions"=g_lines, 
      "candidates"=subjectHits(
        regions2candidates[queryHits(regions2candidates) %in% g_lines]),
      "noncandidates"=subjectHits(
        regions2noncandidates[queryHits(regions2noncandidates) %in% g_lines]))
  }
  rm(progEnvir)
	
  ## Return list...
  list(
    "goterms"=as.character(goterms), 
    "regions2candidates"=regions2candidates,
    "regions2noncandidates"=regions2noncandidates,
    "ranges"=ranges, 
    "go2region"=go2region
  )
}
