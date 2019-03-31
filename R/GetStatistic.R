## This function uses Fisher's exact test to detect GO terms having an
## over-representation of candidate SNPs.
## Parameters: 
## - goterm is the term to be tested,
## - gotermsByLevel contains GO terms ordered by level
## - snpsByLevel contains _number_ of SNPs assigned to different levels
## - cc is the number of candidate SNPs associated with goterm
## - cn is the number of non-candidate SNPs associated with goterm

GetStatistic <-
function(goterm, gotermsByLevel, snpsByLevel, cc, cn) {
    
  ## cc - candidate term and candidate SNP
  ## cn - candidate term and non-candidate SNP

  o <- Ontology(goterm)
  level <- FindLevel(goterm, gotermsByLevel[[o]])
  allCandidates <- snpsByLevel[[o]][["LengthCandidatesByLevel"]][level]
  allNoncandidates <- snpsByLevel[[o]][["LengthNoncandidatesByLevel"]][level]
	
  ## nc - non-candidate terms (background GO terms) and candidate SNP
  ## nn - non-candidate term and non-candidate SNP
  nc <- allCandidates - cc
  nn <- allNoncandidates - cn
  
  ## Return list
  list(
    goterm,
    fisher.test(matrix(c(cc,cn,nc,nn),byrow=T,nrow=2),alternative="g")$p.value
  )
}
