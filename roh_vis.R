#Visualize ROH findings

require(dplyr)
require(ggplot2)

roh <- read.table("C:/Users/Mary/PepLab/data/blasto/ROH/ROH_50SNPs_3allowance_50kb.hom.overlap", header=T)

#convert df to tbl_df
roh <- tbl_df(roh)
glimpse(roh)
roh
class(roh)

chr7 <- filter(roh, CHR == 7, FID != "UNION")
sum <- filter(roh, CHR == 7, FID == "UNION" | FID == "CON")
pchr7 <- filter(chr7, BP1 >= 22618338 & BP1 <= 22919561 & BP2 >= 22618338 & BP2 <= 22919561 ) 

ggplot(sum) + 
  geom_segment(aes(x=BP1, y=FID, xend=BP2, yend=FID))

ggplot(pchr7) + 
  geom_segment(aes(x=BP1, y=IID, xend=BP2, yend=IID)) +
  geom_vline(xintercept = 22766246, colour="red") +
  geom_vline(xintercept = 22768219, colour="red") + 
  geom_vline(xintercept = 22769249, colour="red") + 
  geom_vline(xintercept = 22799526, colour="red") +
  ylab("") + 
  xlab("Genomic Positon") + 
  


require(vcfR)
vcf_file <- "C:/Users/Mary/PepLab/data/blasto/ROH/chr7_22618338-22919561.pilot.vcf"
dna_file <- "C:/Users/Mary/PepLab/data/blasto/ROH/chr7_22618338-22919561.fasta"
gff_file <- "C:/Users/Mary/PepLab/data/blasto/ROH/chr7_22618338-22919561.gtf"

bdh001 <- read.vcfR("C:/Users/Mary/PepLab/data/blasto/ROH/BDH001_chr7_22618338-22919561.pilot.vcf")

vcf <- read.vcfR(vcf_file, verbose=FALSE) 
dna <- ape::read.dna(dna_file, format="fasta")
gff <- read.table(gff_file, sep="\t", quote="")

chrom <- create.chromR(vcf=bdh001, name="Chr7", seq=dna, ann=gff)
plot(chrom)
chromoqc(chrom, dp.alpha=20, xlim=c(22618338, 22919561))

library(VariantAnnotation)
vcf <- readVcf(vcf_file, "hg19")

vr <- as(vcf[, 1:3], "VRanges")
#vr <- renameSeqlevels(vr, value = c("7", "chr7"))

gr7 <- GRanges("chr7", IRanges(22618338, 22919561))
p.vr <- autoplot(vr, which = wh)


autoplot(vr, which = wh, geom="rect", arrow = FALSE)
##################3
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomeGraphs")

require(GenomeGraphs)

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

plusStrand <- makeGeneRegion(chromosome = 19, start = 12050000, end = 12230000, strand = "+", biomart = mart)
minStrand <- makeGeneRegion( chromosome = 19, start = 12050000, end = 12230000, strand = "-", biomart = mart)
ideogram <- makeIdeogram(chromosome = 19)
genomeAxis <- makeGenomeAxis(add53=TRUE, add35=TRUE)
gdPlot(list(ideogram, plusStrand, genomeAxis, minStrand))



minbase <- 22618338
maxbase <- 22919561

genesplus <- makeGeneRegion(start = minbase, end = maxbase, 
                            strand = "+", chromosome = "7", biomart=mart)
genesmin <- makeGeneRegion(start = minbase, end = maxbase, 
                           strand = "-", chromosome = "7", biomart=mart)
ideogram <- makeIdeogram(chromosome = 7)
genomeAxis <- makeGenomeAxis(add53=TRUE, add35=TRUE)
gdPlot(list(ideogram, genesplus, genesmin, genomeAxis))

seg <- makeSegmentation(segStart[[1]], segEnd[[1]], segments[[1]], 
                        dp = DisplayPars(color = "black", lwd=2,lty = "solid"))



################ggbio
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggbio")

library(ggbio)
p.ideo <- Ideogram(genome = "hg19")
p.ideo

library(GenomicRanges)
p.ideo + xlim(GRanges("chr7", IRanges(22618388, 22919561)))

#source("https://bioconductor.org/biocLite.R")
#biocLite("Homo.sapiens")
library(Homo.sapiens)
class(Homo.sapiens)

data(genesymbol, package = "biovizBase")
#wh <- genesymbol[c("IL6")]
wh <- GRanges("chr7", IRanges(22618388, 22919561))
wh <- range(wh, ignore.strand = TRUE)

p.txdb <- autoplot(Homo.sapiens, which = wh)
p.txdb

autoplot(Homo.sapiens, which = wh, label.color = "black", color = "brown", fill = "brown")


#source("https://bioconductor.org/biocLite.R")
#biocLite("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
ensdb <- EnsDb.Hsapiens.v75
autoplot(ensdb, GenenameFilter("IL6"))

gr <- GRanges(seqnames=16, IRanges(30768000, 30770000), strand="*")
autoplot(ensdb, GRangesFilter(gr, "overlapping"), names.expr="gene_name")

#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
bg <- BSgenome.Hsapiens.UCSC.hg19
p.bg <- autoplot(bg, which = wh)
## no geom
p.bg
## segment
p.bg + zoom(1/100)
## rectangle
p.bg + zoom(1/1000)
## text
p.bg + zoom(1/2500)



###karyogram

data(ideoCyto, package = "biovizBase")
autoplot(seqinfo(ideoCyto$hg19), layout = "karyogram")

biovizBase::isIdeogram(ideoCyto$hg19)
autoplot(ideoCyto$hg19, layout = "karyogram", cytoband = TRUE)
