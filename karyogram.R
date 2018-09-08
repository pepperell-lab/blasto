

require(dplyr)
require(ggplot2)


#Read in ROH data from PLINK

roh <- read.table("C:/Users/Mary/PepLab/data/blasto/ROH/ROH_50SNPs_3allowance_50kb.hom.overlap", header=T)
pilot <- read.table("C:/Users/Mary/PepLab/data/blasto/postdoc/pilot_CON_fiveOrMore.overlap", header=F)
names(pilot) <- names(roh)
blasto <- read.table("C:/Users/Mary/PepLab/data/blasto/postdoc/blasto_11_CON_nineCasesOrMore.overlap", header=F)
names(blasto) <- names(roh)

#convert df to tbl_df
roh <- tbl_df(roh)
glimpse(roh)
roh
class(roh)
roh <- mutate(roh, CHR = paste("chr",as.character(CHR), sep=""))
pilot <- mutate(pilot, CHR = paste("chr",as.character(CHR), sep=""))
blasto <- mutate(blasto, CHR = paste("chr",as.character(CHR), sep=""))


sum <- filter(roh, FID == "CON", PHE == "5:1" | PHE == "5:0") #| FID == "CON")
#sum <- arrange(sum, desc(FID))

bdh001 <- filter(roh, IID == "bdh001")

ggplot(sum) + 
  geom_segment(aes(x=BP1, y=FID, xend=BP2, yend=FID))


ggplot() + geom_boxplot(aes(x="KB", y=sum$KB))


sum %>%
  group_by(CHR) %>%
  summarise(min(KB), max(KB), mean(KB), sum(KB))

pilot %>%
  group_by(CHR) %>%
  summarise(min(KB), max(KB), mean(KB), sum(KB))

blasto %>%
  group_by(CHR) %>%
  summarise(min(KB), max(KB), mean(KB), sum(KB))



combROH <- rbind(data.frame(sample="pilot", sum[,c("CHR", "KB")]),
                            data.frame(sample="pilot_BROAD", pilot[,c("CHR", "KB")]),
                            data.frame(sample="all_BROAD", blasto[,c("CHR", "KB")]))

combROH$sampCHR <- paste(combROH$CHR, combROH$sample, sep="_")
  
ggplot(combROH) + 
  geom_boxplot(aes(x=sampCHR, y=KB, colour=sample)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("") +
  ylab("Consensus ROH Segment Length (KB)")


avg <- combROH %>%
  group_by(sample, CHR) %>%
  summarise(min(KB), max(KB), mean(KB), sum(KB))

avg$CHRsamp <- paste(avg$CHR, avg$sample, sep="_")

ggplot(avg, aes(x=CHRsamp, y=`sum(KB)`, fill=sample)) +
  geom_bar(stat="identity") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("") +
  ylab("Sum of Consensus ROH (KB)")

#Create Karyogram
library(ggbio)
library(GenomicRanges)
library(Homo.sapiens)

data(ideoCyto, package = "biovizBase")
autoplot(seqinfo(ideoCyto$hg19), layout = "karyogram")

biovizBase::isIdeogram(ideoCyto$hg19)
autoplot(ideoCyto$hg19, layout = "karyogram", cytoband = TRUE)

data(darned_hg19_subset500, package = "biovizBase")
dn <- darned_hg19_subset500
library(GenomicRanges)
seqlengths(dn)
## add seqlengths
## we have seqlegnths information in another data set
seqlengths(dn) <- seqlengths(ideoCyto$hg19)[names(seqlengths(dn))]
## then we change order
dn <- keepSeqlevels(dn, paste0("chr", c(1:22, "X")))
seqlengths(dn)
autoplot(dn, layout = "karyogram")

test <- makeGRangesFromDataFrame(sum, seqnames.field="CHR", start.field="BP1", end.field="BP2", keep.extra.columns=TRUE, ignore.strand=TRUE)
pilotGR <- makeGRangesFromDataFrame(pilot, seqnames.field="CHR", start.field="BP1", end.field="BP2", keep.extra.columns=TRUE, ignore.strand=TRUE)
blastoGR <- makeGRangesFromDataFrame(blasto, seqnames.field="CHR", start.field="BP1", end.field="BP2", keep.extra.columns=TRUE, ignore.strand=TRUE)

#rohP <- autoplot(test, layout = "karyogram", aes(color = FID, fill = FID), alpha=0.5) + scale_color_manual(values=c("black", "#66CCCC")) + scale_fill_manual(values=c("black", "#66CCCC"))

rohP <- autoplot(test, layout = "karyogram") + scale_color_manual(values=c("black")) + scale_fill_manual(values=c("black")) 

ggbio::ggsave(filename="C:/Users/Mary/PepLab/data/blasto/figs_for_ASHG/karyorgam_roh.eps", rohP, width = 6, height = 6, units="in", dpi = 300 )




bdh1 <- makeGRangesFromDataFrame(bdh001, seqnames.field="CHR", start.field="BP1", end.field="BP2", keep.extra.columns=FALSE, ignore.strand=TRUE)
autoplot(bdh1, layout = "karyogram")




