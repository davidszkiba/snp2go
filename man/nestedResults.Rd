\name{nestedResults}
\alias{nestedResults}
\title{
    Creates a data.frame that sorts enriched GO terms by domain and level.
}
\description{
    For each enriched GO term, the domain and level inside of the domain is added.
    Additionally, GO terms are sorted by their level inside the GO DAG and nested.

    
    ...
}
\usage{
    nestedResults(snp2goResult)
}
\arguments{
    \item{snp2goResult}{
        The list returned by the snp2go function.
    }
}
\details{
    Nesting of GO terms means, that all child terms of a GO term at level L are 
    displayed, before the next significant GO term of level L is shown.
    This is best illustrated by an example: 
    \code{
    GO:1 level:1 childterms: GO:2, GO:3, GO:4
    GO:2 level:2 childterms: GO:4
    GO:3 level:2 childterms: GO:4
    GO:4 level:3 childterms: 
    GO:5: level:1 childterms: GO:6, GO:7
    }
    GO term GO:1 has 3 significant child terms (GO:2, GO:3 and GO:4). The next 
    GO term at the same level as GO:1, GO:5, is displayed after all child terms 
    of GO:1 were displayed.
}
\value{
    Returns a data.frame.
}
\author{
    David Szkiba
}
\seealso{
   \code{\link{snp2go}}
}
\examples{
   c <- GenomicRanges::GRanges(
       seqnames=rep(1, 100),
       ranges=IRanges(
           rep(1, 100),
           runif(min=2, max=1000, n=100)
       )
   )
   nc <- GenomicRanges::GRanges(
       seqnames=rep(1, 1000),
       ranges=IRanges(
           rep(1, 100),
           runif(min=2, max=1000, n=1000)
       )
   )
\dontrun{x <- snp2go(gff="dmel.gff", candidateSNPs=c, noncandidateSNPs=nc)}
   \dontrun{nested <- nestedResults(x)}
}
\keyword{snp2go, SNP2GO}
