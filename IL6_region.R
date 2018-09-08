#IL6 figure


library(ggbio)
#chr 1 is automatically drawn by default (subchr="chr1")
p.ideo <- Ideogram(genome = "hg19")
p.ideo
#Highlights a region on "chr7"
library(GenomicRanges)
p.chrom <- p.ideo + xlim(GRanges("chr7", IRanges(22618338, 22919561)))

library(Homo.sapiens)
#load gene symbol : GRanges, one gene/row
data(genesymbol, package = "biovizBase")
#retrieve information of the gene of interest
gr7 <- GRanges("chr7", IRanges(22618338, 22919561))
wh <- range(gr7, ignore.strand = TRUE)
#Plot the different transcripts  for our gene of interest
p.txdb <- autoplot(Homo.sapiens, which = wh)
p.txdb
#Change inton geometry, use gap.geom
autoplot(Homo.sapiens, which = wh, gap.geom = "chevron")




### Grab SNPs in region from 1k.R
snps <- filter(same, CHROM == "chr7", POS <= 22919561 & POS >= 22618338)
snp <- snps[,c('CHROM', 'POS', 'EUR_AF', 'ID', 'REF', 'ALT')] 
snp.min <- snps[,c('CHROM', 'POS')]
#make into a grange object
var.gr <- makeGRangesFromDataFrame(snp, seqnames.field="CHROM", start.field="POS", end.field="POS", keep.extra.columns=TRUE, ignore.strand=TRUE)
var.gr

var.min <- makeGRangesFromDataFrame(snp.min, seqnames.field="CHROM", start.field="POS", end.field="POS", keep.extra.columns=FALSE, ignore.strand=TRUE)


### Grab ROH in region
#Read in ROH data from PLINK

#roh <- read.table("C:/Users/Mary/PepLab/data/blasto/ROH/ROH_50SNPs_3allowance_50kb.hom.overlap", header=T)

roh <- read.table("C:/Users/Mary/PepLab/data/blasto/postdoc/ROH_IL6region.overlap", header=F, sep="\t")
names(roh) <- c("POOL", "FID", "IID", "PHE", "CHR", "SNP1", "SNP2", "BP1", "BP2", "KB", "NSNP", "NSIM", "GRP")


gar <- read.table("C:/Users/Mary/PepLab/data/blasto/postdoc/auto_test_garlic_il6region.txt", header=F, sep="\t")
names(gar) <- c("IID", "CHR", "BP1", "BP2", "Size_Class", "ROH_length", "placeholder", "rgb_track_color")
gar$method <- "garlic"

roh.il6$method <- "plink"
#convert df to tbl_df
roh <- tbl_df(roh)
glimpse(roh)
roh
class(roh)
roh <- mutate(roh, CHR = paste("chr",as.character(CHR), sep=""))

chr7 <- filter(roh, CHR == "chr7", FID != "UNION" & FID != "CON")
roh.il6 <- chr7[,c(3:12)]
roh.il6 <- distinct(roh.il6)

pchr7 <- filter(roh.il6, BP1 >= 22618338 & BP1 <= 22919561 & BP2 >= 22618338 & BP2 <= 22919561 ) 

roh.gr <- makeGRangesFromDataFrame(roh.il6, seqnames.field="IID", start.field="BP1", end.field="BP2", keep.extra.columns=TRUE, ignore.strand=TRUE)

#Easier just to make the plot with ggplot
p.roh <- ggplot(roh.il6) + 
  geom_segment(roh.il6, aes(x=BP1, y=IID, xend=BP2, yend=IID)) +
  #scale_x_continuous(limits=c(22618338, 22919561)) +
  theme_bw()

p.roh.gar <- ggplot(gar) + 
  geom_segment(gar, aes(x=BP1, y=IID, xend=BP2, yend=IID)) +
  #scale_x_continuous(limits=c(22618338, 22919561)) +
  theme_bw()



gar$id <- paste(gar$IID, gar$method, sep="-") 
roh.il6$id <- paste(roh.il6$IID, roh.il6$method, sep="-")

comb <- full_join(roh.il6, gar, by="id")

rohP + xlim(GRanges("chr7", IRanges(22618338, 22919561)))


### Test adding plots together
tks <- tracks(p.chrom, gene = p.txdb, roh = p.roh)
heights = c(1, 2, 3) + theme_tracks_sunset()
tks



####
#source("https://bioconductor.org/biocLite.R")
#biocLite("Gviz")
library(Gviz)
library(GenomicFeatures)


### Gene Region Track
gff_file <- "C:/Users/Mary/PepLab/data/blasto/ROH/chr7_22618338-22919561.gtf"
gff <- read.table(gff_file, sep="\t", quote="")
tt <- makeTxDbFromGFF(gff_file)
txTr <- GeneRegionTrack(tt, chromosome="chr7", start=22618338, end=22919561)

gtrack <- GenomeAxisTrack()

grtrack <- GeneRegionTrack(tt, genome = "hg19",
                           chromosome = "chr7", name = "Gene Model", 
                           transcriptAnnotation = "symbol",
                           #background.panel = "#FFFEDB",
                           #background.title = "darkblue",
                           start=22618338, end=22919561)


itrack <- IdeogramTrack(genome = "hg19", chromosome = "chr7")
plotTracks(itrack, from = 22618338, to = 22919561, showId = FALSE)

plotTracks(list(itrack, gtrack, grtrack), from = 22618338, to = 22919561)

atrack <- AnnotationTrack(var.min, name = "SNPs",
                          background.panel = "#FFFEDB",
                          background.title = "darkblue",
                          start=22618338, end=22919561)
plotTracks(atrack)




ind.roh <- rbind(c("bdh001", "chr7", 22683279, 22813060),
c("bdh002", "chr7", 22708276, 22771738),
c("bdh002", "chr7", 22795077, 22919561),
c("bdh003", "chr7", 22755012, 22811537),
c("bdh004", "chr7", 22710432, 22813225),
c("bdh005", "chr7", 22754768, 22812324),
c("bdh006", "chr7", 22618338, 22803560),
c("bdh006", "chr7", 22807884, 22878104))
ind.roh <- data.frame(ind.roh)
names(ind.roh) <- c("IID", "CHR", "BP1", "BP2")
ind.roh$BP1 <- as.numeric(as.character(ind.roh$BP1))
ind.roh$BP2 <- as.numeric(as.character(ind.roh$BP2))

roh.gr <- makeGRangesFromDataFrame(ind.roh, seqnames.field="CHR", start.field="BP1", end.field="BP2", keep.extra.columns=TRUE, ignore.strand=TRUE)

roh.gr <- makeGRangesFromDataFrame(roh.il6, seqnames.field="CHR", start.field="BP1", end.field="BP2", keep.extra.columns=TRUE, ignore.strand=TRUE)


rohtrack <- AnnotationTrack(roh.gr, name = "ROH")

plotTracks(list(itrack, gtrack, grtrack, rohtrack, atrack), from = 22618338, to = 22919561)

plotTracks(list(itrack, gtrack, grtrack, rohtrack), from = 22618338, to = 22919561)

snptrack <- DataTrack(var.gr, name="EUR AF of SNPs")
plotTracks(snptrack)

plotTracks(list(itrack, gtrack, grtrack, rohtrack, snptrack), from = 22618338, to = 22919561)



#highlight
ht <- HighlightTrack(trackList = list(grtrack, rohtrack),
                     start = c(22755012, 22795077, 22807884), end = c(22771738, 22803560, 22811537), chromosome = 7)

ht2 <- HighlightTrack(trackList = list(grtrack, snptrack),
                     start = c(22755012, 22795077, 22807884), end = c(22771738, 22803560, 22811537), chromosome = 7)

plotTracks(list(itrack, gtrack, ht), from = 22618338, to = 22919561)


plotTracks(list(itrack, gtrack, ht), from = 22755012, to = 22811537)
plotTracks(list(itrack, gtrack, ht2), from = 22755012, to = 22811537)


#setEPS()

cairo_ps("cairo_IL6_zoom.eps", width = 8, heigh = 4)
plotTracks(ht2, from = 22755012, to = 22811537)
dev.off()

cairo_ps("cairo_IL6_zoom3.eps", width = 8, heigh = 2)
plotTracks(ht2, from = 22755012, to = 22811537)
dev.off()

#cairo_ps("cairo_IL6.eps", width = 8, height = )
#plotTracks(list(itrack, gtrack, ht), from = 22618338, to = 22919561)
#dev.off()




###

ggplot() + geom_point(aes(x=snps$POS, y=snps$EUR_AF, color=snps$REG))
ggplot() + geom_point(aes(x=snps$POS, y=snps$EUR_AF, color=snps$GTEx))


snps[is.na(snps$GTEx),]
