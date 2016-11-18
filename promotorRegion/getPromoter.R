## Falk library
## Caleb, 04/03/2016
## wong10
## for single cell methylaiton


rm(list=ls())
options(stringsAsFactors=FALSE)

WD <- '/home/xiaoguang/Database/promoter'
system(paste('mkdir -p', WD))
setwd(WD)

if(F){
	source("http://bioconductor.org/biocLite.R")
	biocLite("org.Mm.eg.db")
	biocLite("GenomicFeatures")
	biocLite("BSgenome")
	biocLite("BSgenome.Mmusculus.UCSC.mm9")
	biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene")	
}

library(org.Mm.eg.db)
library(GenomicFeatures)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(foreach)
library(doParallel)

genes <- transcriptsBy(TxDb.Mmusculus.UCSC.mm9.knownGene, "gene")
promoter = promoters(genes,upstream=2000, downstream=0)

entrez <- names(promoter)
chr <- as.character(sapply(runValue(seqnames(promoter)),function(x) x[[1]]))
start <- sapply(start(ranges(promoter)),function(x) x[[1]])
end <- sapply(end(ranges(promoter)),function(x) x[[1]])
strand <- as.character(sapply(runValue(strand(promoter)),function(x) x[[1]]))


x <- org.Mm.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
# Get the SYMBOL for the first five genes
xx[1:5]
# Get the first one
xx[[1]]
}



mm_genes <- unlist(xx[entrez])

nonnaIndex <- match(names(mm_genes),entrez)
all(entrez[nonnaIndex] == names(mm_genes))

promoterList <- data.frame(entrez=entrez[nonnaIndex], genes = mm_genes, chr=chr[nonnaIndex], start=start[nonnaIndex], end=end[nonnaIndex], strand=strand[nonnaIndex])

write.csv(promoterList,file='promoterList_mm9_up2000down0.csv',row.names=FALSE)


