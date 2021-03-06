\name{snp2go}
\alias{snp2go}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Functional Analysis of Genome Wide Association Studies
}
\description{
	The snp2go function detects GO categories that have a significant overrepresentation of candidate SNPs and infers whether local effects, such as linkage disequilibrium, are contributing to the significance.
}
\usage{
snp2go(gtf, gff, goFile, FDR = 0.05, runs = 100000, candidateSNPs, noncandidateSNPs, extension=0, verbose=TRUE, min.regions=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gtf}{
	Path to the GTF file.
}
  \item{gff}{
	 Path to the GFF file. The GFF files must contain the attribute "Ontology_term".
}
  \item{goFile}{
	 Path to a file containing the gene ID -> GO ID associations. The gene IDs must be the same as the gene IDs used in the GTF file. This parameter is needed only when the "gtf" parameter is used.
}
  \item{FDR}{
	Significance level. SNP2GO will report GO categories with FDR values lower than that set by the user. Default: 0.05.
}
  \item{runs}{
	 Number of times the program will repeat the random sampling. Default: 100,000.
}
  \item{candidateSNPs}{
	A GenomicRanges object containing the coordinates of the candidate SNPs.
}
  \item{noncandidateSNPs}{
	 A GenomicRanges object containing the coordinates of the non-candidate SNPs.
}
  \item{extension}{
	 Number of nucleotides up and down the gene that will be added for the extended definition of the genome region. Default: 0 (i.e., gene region = gene)
}
  \item{verbose}{Print progress report of the analysis. Default = 'TRUE'}

  \item{min.regions}{
	  Number of regions associated to a GO category. SNP2GO will report GO categories having at least min.regions. Default: 10.
}
}
\details{
	
SNP2GO performs an inclusive analysis of GO categories, according to which GO categories form the same level are tested together and genes associated to categories that are descendant of the tested categories take the annotation from the parent. SNP2GO calls the function getGOLevel from the Bioconductor package goProfiles to extract the GO structure. The analysis is done from top to bottom and genes exclusively associated to levels above the tested GO category are discarded.

After the inclusive analysis, SNP2GO defines the genome regions associated to the GO categories and links candidate and non-candidate SNPs accordingly. This association is performed by calling the Bioconductor package GenomicRanges. For each GO level, the program computes the total number t = n + m of candidate SNPs, n, and non-candidate SNPs, m. Later, for each GO category, Ci, belonging to that level SNP2GO computes the number of candidate (nc) and non-candidate (mc) SNPs associated to it and calculate the probability to obtain at least nc candidate SNPs in that category. Correction for multiple testing is performed applying false discovery rate (Benjamini and Hochberg, 1995).

Finally, a local test is applied to significant GO categories. To this end, SNP2GO determines the number, g, of genomic regions in Ci with at least one of the nc candidate SNPs. Then it samples nc SNPs from the nc + mc SNPs without replacement and counts the number of genomic regions where the SNPs belong to. After r runs, the resulting distribution is used as null distribution to test if g is bigger or smaller than expected. Non-significance indicates the absence of local effects. If g is smaller then LD may cause the significance; if g is larger, then the candidate SNP distribution is overdispersed. 
}

\value{
  A named list containing the following elements:

    \item{candidate.snps}{Number of candidate SNPs provided by the user.}
    \item{noncandidate.snps}{Number of non-candidate SNPs provided by the user.}
    \item{informative.candidate.snps}{Number of candidate SNPs that can be associated with a GO category.}
    \item{informative.noncandidate.snps}{Number of non-candidate SNPs that can be associated with an annotated GO category.}
    \item{FDR}{FDR value set by the user.}
    \item{runs}{Number of times SNP2GO repeated the hypergeometric sampling.}
    \item{extension}{Number of nucleotides up and down the gene used to extend the definition of the genome region.}
    \item{regions}{A GRanges object containing all regions of the input file that are associated to at least one GO category. The regions are extended by the number of nucleotides specified by the \code{extension} parameter.}
    \item{goterms}{The list of GO categories used in the analysis.}
    \item{go2ranges}{A list of three hashes that contains the mappings of GO categories to GO regions, candidate SNPs and noncandidate SNPs: \code{go2ranges[["regions"]][["g"]]} contains the indices of all regions, \code{go2ranges[["candidates"]][["g"]]} contains the indices of all candidate SNPs and \code{go2ranges[["noncandidates"]][["g"]]} contains the indices of all noncandidate SNPs associated with GO category \code{g}. The indices of the GO regions refer to the \code{regions} element returned by the snp2go function, the indices of the candidate- and noncandidate SNPs refer to the candidate- and noncandidate SNPs passed to the snp2go function.} 
\item{enriched}{A dataframe that reports the significant GO categories.}

	The columns of the \code{enriched} dataframe provide the next information:
        \item{GO: }{ID of the significant GO category.}
        \item{p.L: }{Proportion of iterations in which the hypergeometric sampling found less or equal candidate regions than observed.}
        \item{p.G: }{Proportion of iterations in which the hypergeometric sampling found more or equal candidate regions than observed.}
        \item{g: }{Number of regions associated with the significant GO category having at least one candidate SNP.}
        \item{G: }{Total number of genome regions associated with the significant GO category.}
        \item{nc: }{Number of candidate SNPs located in G.}
        \item{mc: }{Number of non-candidate SNPs located in G.}
        \item{P: }{Hypergeometric probability to obtain nc or more candidate SNPs when you randomly sample t SNPs from T. See http://www.cibiv.at/software/snp2go for more details.}
        \item{FDR: }{Adjusted P-Value after applying the Benjamini-Hochberg method.}
        \item{GO.definition: }{Definition of the significant GO category.}
        \item{Child.GOs: }{Child categories of the significant GO category that are also significantly enriched with candidate SNPs.}
}

\author{
David Szkiba
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

    ## Randomly create 1000 SNPs
    library(GenomicRanges)
    seqnames <- c("dmel_mitochondrion_genome","2RHet","2L","X","3L","XHet","Uextra",
    "4","YHet","U","2LHet","2R","3LHet","3R","3RHet")
    rel.ends <- c(19519,3288763,23011546,22422829,24543559,204113,29004658,
    1351859,347040,10049039,368874,21146710,2555493,27905055,2517509)
    tot.ends <- cumsum(rel.ends)
    snp.pos <- sample(max(tot.ends),1000)
    snp.iv <- findInterval(snp.pos,tot.ends,rightmost.closed=TRUE)+1
    snp.pos[snp.iv>1]<-snp.pos[snp.iv>1]-tot.ends[snp.iv[snp.iv>1]-1]
    snp.chr <- seqnames[snp.iv]
    snps <- unique(GRanges(seqnames=snp.chr,ranges=IRanges(snp.pos,snp.pos)))
    candidates <- snps[1:200]
    noncandidates <- snps[201:length(snps)]
    
    ## Run SNP2GO with a modified/shortened GFF file taken from www.flybase.org
    x <- 
      snp2go(gff=system.file(package="SNP2GO", "extdata/dmel.gff"),
        candidateSNPs=candidates,
        noncandidateSNPs=noncandidates,
        FDR=1.1,
        runs=1000)

    ## Show enriched GO categories
    x$enriched$GO
    
    ## Store the enriched GO categories in a tab-seperated file
    \dontrun{write.table(x$enriched,file="result.tsv",sep="\t",row.names=FALSE)}
}

