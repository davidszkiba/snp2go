# SNP2GO

SNP2GO is an R package that tests Gene Ontology terms for enrichment of candidate SNPs and infers whether enrichment is influenced by local effects.

# Installing SNP2GO

You can download the program here and install it using R's install.packages function:

```R
install.packages(pkgs="/home/me/downloads/SNP2GO_1.0.2.tar.gz", type="source")
library(SNP2GO)	#load the package
?SNP2GO #view the help-page of the package
?snp2go #view the help-page of the snp2go function
```

Note: The package depends on the packages goProfiles, hash and GenomicRanges.
Please make sure that those packages are installed.
If not, you can install them typing the following commands in your R console:



```R
install.packages("hash")
source("http://bioconductor.org/biocLite.R")
biocLite("goProfiles")
biocLite("GenomicRanges")
```

# Reference
D. Szkiba, M. Kapun, A. von Haeseler, and M. Gallach. SNP2GO: functional analysis of genome-wide association studies. (2014) Genetics. 197: 285-289 (https://dx.doi.org/10.1534/genetics.113.160341)

[CIBIV SNP2GO Website](http://www.cibiv.at/software/snp2go/index.shtml)
