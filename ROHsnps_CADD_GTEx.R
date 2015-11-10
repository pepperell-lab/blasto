#manipulating ROH tables
rm(list=ls())

filename = "C:/Users/Mary/Desktop/CLEANUP/50SNPs_3allowance_50kb.ann.gz"
hist <- read.delim(filename, strip.white = TRUE, fill = TRUE, header = FALSE, sep = "\t")
hist <- hist[-1,]
COLS <- c("Chrom","Pos","Ref","Anc","Alt","Type","Length","isTv","isDewhrived","AnnoType","Consequence","ConsScore","ConsDetail","GC","CpG","mapAbility20bp",  "mapAbility35bp",  "scoreSegDup","priPhCons","mamPhCons","verPhCons","priPhyloP","mamPhyloP","verPhyloP","GerpN","GerpS", "GerpRS",  "GerpRSpval","bStatistic","EncExp","EncH3K27Ac","EncH3K4Me1","EncH3K4Me3","EncNucleo","EncOCC", "EncOCCombPVal", "EncOCDNasePVal","EncOCFairePVal", "EncOCpolIIPVal","EncOCctcfPVal", "EncOCmycPVal","EncOCDNaseSig", "EncOCFaireSig",   "EncOCpolIISig",   "EncOCctcfSig","EncOCmycSig","Segway",  "tOverlapMotifs",  "motifDist","motifECount","motifEName","motifEHIPos","motifEScoreChng", "TFBS","TFBSPeaks","TFBSPeaksMax","isKnownVariant", "ESP_AF", "ESP_AFR", "ESP_EUR", "TG_AF", "TG_ASN", "TG_AMR",  "TG_AFR", "TG_EUR",  "minDistTSS","minDistTSE","GeneID", "FeatureID","CCDS","GeneName","cDNApos", "relcDNApos","CDSpos", "relCDSpos","protPos", "relProtPos","Dst2Splice","Dst2SplType","Exon","Intron","oAA","nAA","Grantham","PolyPhenCat","PolyPhenVal","SIFTcat","SIFTval", "RawScore","Cscore","Genotypes","snpID","CandidateGene","CONstart","CONstop","SegID","NumInd")
colnames(hist) <- COLS
#novel <- subset(hist, isKnownVariant == "FALSE" & (NumInd == "5:0" | NumInd == "5:1"))
#tempGeno <- hist$Genotypes
which(COLS=="GeneName")
#unique positions
hist_nodups <- hist[!duplicated(hist[,'Pos']),]

library(stringr)
newGeno <- str_split_fixed(hist$Genotypes, " ", 12)
colnames(newGeno) <- paste(c("G"), 1:12, sep="")
hist2 <- cbind(hist, newGeno)
sameA <- subset(hist2, Ref != 'A' & G1 == 'A' & G2 == 'A' & G3 == 'A' & G4 == 'A' & G5 == 'A' & G6 == 'A' & G7 == 'A' & G8 == 'A' & G9 == 'A' & G10 == 'A')
sameC <- subset(hist2, Ref != 'C' & G1 == 'C' & G2 == 'C' & G3 == 'C' & G4 == 'C' & G5 == 'C' & G6 == 'C' & G7 == 'C' & G8 == 'C' & G9 == 'C' & G10 == 'C')
sameT <- subset(hist2, Ref != 'T' & G1 == 'T' & G2 == 'T' & G3 == 'T' & G4 == 'T' & G5 == 'T' & G6 == 'T' & G7 == 'T' & G8 == 'T' & G9 == 'T' & G10 == 'T')
sameG <- subset(hist2, Ref != 'G' & G1 == 'G' & G2 == 'G' & G3 == 'G' & G4 == 'G' & G5 == 'G' & G6 == 'G' & G7 == 'G' & G8 == 'G' & G9 == 'G' & G10 == 'G')
same <- rbind(sameA, sameC, sameT, sameG)
same_nodups <- same[!duplicated(same[,'Pos']),]

levels(same$G1)
same_novel <- subset(same, isKnownVariant == "FALSE")
same_novel_nodups <- same_novel[!duplicated(same_novel[,'Pos']),]
same_novel_sub <- same_novel[,c(1:6,10,11,57,61,65,68,69,71,89,90,92,93,98:109)]
same_novel_sub_CG <- subset(same_novel_sub, CandidateGene != ".")
same_novel_sub_CG_nodups <- same_novel_sub_CG[!duplicated(same_novel_sub_CG[,'Pos']),]
same_novel_sub_CG_nodups

#same_novel_sub
same_novel_sub_reg <- subset(same_novel_sub, AnnoType == 'RegulatoryFeature')
same_novel_sub_reg

same_novel_sub_reg <- subset(same_novel_sub, AnnoType == 'RegulatoryFeature')
same_novel_sub_reg
same_novel_sub_reg_nodups <- same_novel_sub_reg[!duplicated(same_novel_sub_reg[,'Pos']),]
same_novel_sub_reg_nodups

levels(same$AnnoType)
same_novel_sub_code <- subset(same_novel_sub, AnnoType == 'CodingTranscript')

same_rare <- subset(same, TG_AF < 0.2)
same_rare_sub <- same_rare[,c(1:6,10,11,57,61,65,68,69,89,90,92,93,98:109)]
same_rare_sub_CG <- subset(same_rare_sub, CandidateGene != ".")

same_rare_sub_reg <- subset(same_rare_sub, AnnoType == 'RegulatoryFeature')
same_rare_sub_reg

before = data.frame(attr = c(1,30), type=c('roll_and_bar','roll_and_foo'))
after <- str_split_fixed(before$type, "_and_", 2); after
colnames(after) <- c("attr1", "type2")
new <- cbind(before, after); new
new[!duplicated(new[,'attr1']),]

duplicated(new)
subset(new, attr1 == 'foo' | type2 == 'foo')



----
  #Mary modifications
  same_candGene <- subset(same, CandidateGene != ".")
same_candGene_sub <- same_candGene[,c(1:6,10,11,57,61:65,68,69,89,90,92,93,98:109)]
same_cadGene_sub_nodups <- same_candGene_sub[!duplicated(same_candGene_sub[,'Pos']),]


rare <- same[is.na(same$TG_EUR) | same$TG_EUR <0.2,]
rare_reg <- subset(rare, AnnoType == 'RegulatoryFeature')




rare_CG <- subset(rare, CandidateGene != ".")
rare_CG_short <- rare_CG[,c(1:6,10,11,57,61:65,68,69,89,90,92,93,98:109)]





same_candGene_sub_nodups_com <- subset(same_cadGene_sub_nodups, (TG_EUR <0.2 | TG_EUR > 0.8))

same_short <- same[,c(1:6,10,11,57,61:65,68,69,89,90,92,93,98:109)]

subset(same_short, snpID == "rs4838590")
subset(same_short, Pos == 112034062)


#make change to test git