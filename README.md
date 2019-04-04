# SNP2GO

SNP2GO is an [R](https://www.r-project.org/) package that tests [Gene Ontology](http://geneontology.org/) terms for enrichment of candidate SNPs and infers whether enrichment is influenced by local effects.

# Requirements

SNP2GO depends on the *goProfiles*, *hash* and *GenomicRanges* packages.
You can install them typing the following commands in your R console:

```R
install.packages("hash")

# Install Bioconductor.
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()

# Install the bioconductor packages:
BiocManager::install(c("GenomicRanges", "goProfiles"))
```

# Installing SNP2GO

You can install SNP2GO by either downloading the package manually or by using the *devtools*
package, which automatically downloads it from github (recommended).

## Using the *devtools* package.

```R
install.packages("devtools")
library("devtools")
install.github("davidszkiba/snp2go")
```

## Manual download
You can download the program [here](http://www.cibiv.at/software/snp2go/download.php) and install it using R's install.packages function:

```R
install.packages(pkgs="/home/me/downloads/SNP2GO_1.0.2.tar.gz", type="source")
library(SNP2GO) #load the package
?SNP2GO #view the help-page of the package
?snp2go #view the help-page of the snp2go function
```

# Parameters

The `snp2go` function takes the following arguments:

- candidateSNPs: A GenomicRanges object containing the coordinates of the candidate SNPs.
- noncandidateSNPs: A GenomicRanges object containing the coordinates of the non-candidate SNPs.
- gff: Path to the GFF file. The GFF files must contain the attribute "Ontology_term"
- gtf: Path to the GTF file
- goFile: Path to a file containing the gene ID --> GO ID associations. The gene IDs must be the same as the gene IDs used in the GTF file. This parameter is needed only when the "gtf" parameter is used.
- runs: Number of times the program will repeat the random sampling. Default: 100,000.
- FDR: Significance of the Fisher's exact test after false discovery rate. SNP2GO will report GO terms with FDR values lower than the value set by the user. Default: 0.05
- extension: Number of nucleotides up and down the gene that will be included for the extended definition of the genome region. Default: 0 (i.e., gene region = gene)
- min.regions: Minimum number of regions associated to GO terms. SNP2GO will report GO terms having at least min.regions.
- verbose: Progress report. Default: TRUE.

# Output
SNP2GO returns a named list containing the following elements:

- candidate.snps: Number of candidate SNPs provided by the user.
- noncandidate.snps: Number of non-candidate SNPs provided by the user.
- informative.candidate.snps: Number of candidate SNPs that can be associated with a GO term.
- informative.noncandidate.snps: Number of non-candidate SNPs that can be associated with an annotated GO term.
- FDR: FDR value set by the user.
- runs: Number of times SNP2GO repeated the hypergeometric sampling.
- extension: Number of nucleotides up and down the gene used to extend the definition of the genome region.
- goterms: The list of all GO terms analysed by SNP2GO.
- regions: A GenomicRanges object that contains all regions of the input file that are associated to at least one GO term. These regions are extended by the number of nucleotides specified by the extension parameter.
- go2ranges: A list that contains three hash tables that store the mappings of GO terms to (1) GO regions (i.e. genes), (2) candidate SNPs and (3) noncandidate SNPs:
- go2ranges[["regions"]][[g]] contains the indices of all GO regions, go2ranges[["candidates"]][[g]] contains the indices of all candidate SNPs and go2ranges[["noncandidates"]][[g]] contains the indices of all noncandidate SNPs associated with GO term g.
- The indices of the GO regions refer to the regions element returned by the snp2go function, the indices of the candidate- and noncandidate SNPs refer to the candidate- and noncandidate SNPs passed as GenomicRanges objects to the snp2go function.
- enriched: A dataframe that reports the significant GO terms.

The enriched dataframe reporting the significant GO terms provides the next information:
- GO: ID of the significant GO term.
- p.L: Proportion of iterations in which the hypergeometric sampling found less or equal candidate regions than observed.
- p.G: Proportion of iterations in which the hypergeometric sampling found more or equal candidate regions than observed.
- g: Number of genomic regions associated with the significant GO term having at least one candidate SNP.
- G: Total number of genomic regions associated with the significant GO term.
- nc: Number of candidate SNPs located in G.
- mc: Number of non-candidate SNPs located in G.
- P: P-value of the Fisher's exact test.
- FDR: Adjusted P-Value after applying the Benjamini-Hochberg method.
- GO.def: Definition of the significant GO term.
- Child.GOs: Child terms that are also significantly enriched with candidate SNPs.

# Examples
In the next example we use data from *Drosophila melanogaster*, but the same work-flow can be applied for any other organism. SNP2GO needs a GFF file (or GTF) and a VCF file(s) containing the SNP coordinates.

A GFF file can be found at the [Flybase](http://www.flybase.org/) website. This GFF file contains the gene coordinates and the GO terms associated with them. Alternatively, a GTF file can be download from [Ensembl](http://www.ensembl.org/). The GTF format is slightly different to GFF's and does not provide GO annotations. Therefore, if you are using a GTF file, you will also need to provide SNP2GO with a gene association file. To get a gene association file, go to Ensembl's [MartView](http://www.ensembl.org/biomart/martview) website and select the "Ensembl Genes" database and the "Drosophila melanogaster genes" dataset. The only attributes you need are the "Ensembl Gene ID" and "GO Term Accession". Download the data and request Ensembl to save it in a tsv-file. This protocol is very useful since many annotated genomes are provided in GTF format.

Finally, you need to provide SNP2GO with a list of SNPs. In our example we will use a single VCF file containing the SNPs. The file can be downloaded [here](http://datadryad.org/handle/10255/dryad.39713), and consists of a list of SNPs sorted by a P-value. The top 2,000 SNPs are described in the original study as candidate SNPs and the rest as non-candidate SNPs. You can try yourself to import them into R and run SNP2GO with the *D. melanogaster* GFF file.

You can copy and paste the next work-flow in your R console:

    # R console:
    # load the SNP2GO package
    library(SNP2GO)

    # Read the VCF file and construct a GenomicRanges object:
    snps <- read.delim("BF15.vcf",header=FALSE,comment.char="#")
        snps[,2] <- as.numeric(snps[,2])
        snps <- GRanges(seqnames=snps[,1],ranges=IRanges(snps[,2],snps[,2]))

    # Use the first 2000 SNPs as candidate SNPs and the rest as non-candidate SNPs:
        cand <- snps[1:2000]
        noncand <- snps[2001:length(snps)]

    # Case 1: Using a GFF file
        x <- snp2go(gff="dmel-all-no-analysis-r5.49.gff.gz",
             candidateSNPs=cand,
             noncandidateSNPs=noncand)

    # Case 2: Using a GTF file + gene association file
        y <- snp2go(gtf="Drosophila_melanogaster.BDGP5.70.gtf.gz",
             goFile="mart_export_dmel.txt",
             candidateSNPs=cand,
             noncandidateSNPs=noncand,
             FDR=0.05,
             runs=10000,
             extension=50)

    # Get all enriched GO terms of GFF analysis:
    gff.significant.terms <- x$enriched$GO

        # Get the first of the enriched GO terms:
        first.term <- gff.significant.terms[1] # this is "GO:0051726"

        # Print all regions associated with at least one GO term:
        print(x$regions)

        # Print the regions associated with the first of the enriched GO terms:
        # There are two possibilities to do so:
        # version 1:
        print(x$regions[ x$go2ranges[[ "regions" ]][[ "GO:0051726" ]] ])
        # version 2:
        print(x$regions[unlist(as.list(x$go2ranges[["regions"]][first.term]))])
        # Although version 2 seems more complicated, it allows to get the regions
        # associated with more than one term. In the following example, all regions associated
        # with the first ten enriched GO terms (gff.significant.terms[1:10]) are printed:
        print(x$regions[unlist(as.list(x$go2ranges[["regions"]][gff.significant.terms[1:10]]))])

        # Print the candidate SNPs associated with the first of the enriched GO terms:
        # Like for the GO regions, there are also two possibilities to do that:
        # version 1:
        print(cand[ x$go2ranges[["candidates"]][["GO:0051726"]] ])
        # version 2:
        print(cand[unlist(as.list(x$go2ranges[["candidates"]][first.term]))])


        # Print the noncandidate SNPs associated with the first of the enriched GO terms:
        # version 1:
        print(noncand[ x$go2ranges[["noncandidates"]][["GO:0051726"]] ])
        # version 2:
        print(noncand[unlist(as.list(x$go2ranges[["noncandidates"]][first.term]))])

    # Get number of informative candidates of the GTF analysis:
    z <- y$informative.candidate.snps

        # Print the list of all GO terms associated to at least one gene in the
        # GFF analysis:
        print(x$goterms)


        # Store the results in tab-seperated files
    write.table(file="snp2go_gff.tsv",x$enriched,sep="\t",row.names=F)
    write.table(file="snp2go_gtf.tsv",y$enriched,sep="\t",row.names=F)

*Note:* SNP2GO will by default print the significant GO terms having at least 10 genome regions associated to them.This behaviour can be changed by the min.regions parameter.

# Reference
D. Szkiba, M. Kapun, A. von Haeseler, and M. Gallach. SNP2GO: functional analysis of genome-wide association studies. (2014) Genetics. 197: 285-289 (https://dx.doi.org/10.1534/genetics.113.160341)

[CIBIV SNP2GO Website](http://www.cibiv.at/software/snp2go/index.shtml)
