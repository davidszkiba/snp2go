test.snp2go <- function() {
  library("RUnit")
  library("SNP2GO")
  
  fileDir = file.path(path.package("SNP2GO"),"inst","extdata")
  c <- read.delim(file.path(fileDir,"testcandidates.tsv"), header=F)
  nc <- read.delim(file.path(fileDir, "testnoncandidates.tsv"),header=F)
  cr <- GRanges(seqnames=c[, 1], ranges=IRanges(c[, 2], c[, 2]))
  ncr <- GRanges(seqnames=nc[, 1], ranges=IRanges(nc[, 2], nc[, 2]))
  
  s <- 
    snp2go(
      gff=file.path(fileDir, "test.gff"), 
      candidateSNPs=cr, 
      noncandidateSNPs=ncr, 
      FDR=1.1, 
      min.regions=1, 
      runs=1000
    )
  print("check g:")
  checkEqualsNumeric(s$enriched[s$enriched$GO=="GO:0043089", "g"], 2)
  print("check G:")
  checkEqualsNumeric(s$enriched[s$enriched$GO=="GO:0043089", "G"], 2)
  print("check nc:")
  checkEqualsNumeric(s$enriched[s$enriched$GO=="GO:0043089", "nc"], 2)
  print("check mc:")
  checkEqualsNumeric(s$enriched[s$enriched$GO=="GO:0043089", "mc"], 4)
  print("check P:")
  checkEqualsNumeric(s$enriched[s$enriched$GO=="GO:0043089", "P"], 
    fisher.test(matrix(c(2,4,3,6), nrow=2, byrow=T), alternative="g")$p.value)
  
}
