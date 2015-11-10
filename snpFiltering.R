#Explore the newly added eQTL annotations to SNPs found in ROH during the pilot study

filename = "151110_ROHsnps_CADD_GTEx.txt"
dat <- read.table(filename, header=T, sep='\t', na.strings="NA")
#rownames(dat) <- paste(dat$Chrom, dat$Pos, sep="_") #cannot do because it is not unique, i.e. 
#multiple rows per SNP

which(colnames(dat)=="GeneName") #Way to find the column that corresponds

sub_dat <- dat[,c(1:5, 10, 79, 81:90, 96, 115:120, 125)] #get rid of unnecessary columns
sub_dat$Genotypes <- as.character(sub_dat$Genotypes) #convert genotypes to characters
sub_dat$snpID <- paste(sub_dat$Chrom, sub_dat$Pos, sep="_")

sameA <- subset(sub_dat, substr(Genotypes, 1, 10) == "AAAAAAAAAA")
sameT <- subset(sub_dat, substr(Genotypes, 1, 10) == "TTTTTTTTTT")
sameC <- subset(sub_dat, substr(Genotypes, 1, 10) == "CCCCCCCCCC")
sameG <- subset(sub_dat, substr(Genotypes, 1, 10) == "GGGGGGGGGG")
same <- rbind(sameA, sameC, sameT, sameG)

same_nodups <- same[!duplicated(same[,'snpID']),] #changed from Pos to snpID so check chrom too

#to see what is removed 
require(sqldf)
diff <- sqldf('SELECT * FROM blood_rare EXCEPT SELECT * FROM blood_rare_CG')
subset(same, Pos == 109144544)
subset(same_nodups, Pos ==109144544)
#identified bug and fixed - was removing 'duplicated' snps that weren't truly duplicates

rare <- same[is.na(same$TG_EUR) | same$TG_EUR <0.2,]
rare_reg <- subset(rare, AnnoType == 'RegulatoryFeature')

eQTL <- same[!is.na(same$eQTLs),]
lung <- eQTL[grep("lung", eQTL$eQTLs), ]
blood <- eQTL[grep("WholeBlood", eQTL$eQTLs), ]

lung_CG <- lung[!is.na(lung$candReg) | !is.na(lung$candGene), ]
blood_CG <- blood[!is.na(blood$candReg) | !is.na(blood$candGene), ]

lung_rare <- lung[(is.na(lung$TG_EUR) & lung$Alt == substr(lung$Genotypes, 1, 1)) | (lung$TG_EUR <=0.1 & lung$Alt == substr(lung$Genotypes, 1, 1)) | (lung$TG_EUR >= 0.9 & lung$Ref == substr(lung$Genotypes, 1, 1)), ]
lung_rare_CG <- lung_rare[!is.na(lung_rare$candReg) | !is.na(lung_rare$candGene), ]


blood_rare <- blood[(is.na(blood$TG_EUR) & blood$Alt == substr(blood$Genotypes, 1, 1)) | (blood$TG_EUR <=0.1 & blood$Alt == substr(blood$Genotypes, 1, 1)) | (blood$TG_EUR >= 0.9 & blood$Ref == substr(blood$Genotypes, 1, 1)), ]
blood_rare_CG <- blood_rare[!is.na(blood_rare$candReg) | !is.na(blood_rare$candGene), ]



