#Explore the newly added eQTL annotations to SNPs found in ROH during the pilot study

setwd("C:/Users/Mary/PepLab/data/blasto")
filename = "151110_ROHsnps_CADD_GTEx.txt"
dat <- read.table(filename, header=T, sep='\t', na.strings="NA")
rownames(dat) <- paste(dat$Chrom, dat$Pos, sep="_") #cannot do because it is not unique, i.e. 
#multiple rows per SNP
#non-unique values when setting 'row.names': ‘1_195783178’, ‘1_33175285’, ‘1_62982891’, ‘1_73453918’, ‘11_54706240’, ‘14_59835972’, ‘17_16240639’, ‘19_29867469’, ‘2_150422200’, ‘3_25926142’, ‘6_123442126’, ‘6_123824253’, ‘6_124515590’, ‘6_127868211’, ‘6_131826865’, ‘6_166372059’, ‘6_79094468’, ‘7_133500313’, ‘8_27777387’, ‘9_131386135’
#20 snps are causing the problem

rownames(dat) <- paste(dat$Chrom, dat$Pos, dat$isKnownVariant, sep="_") #try adding a third identifier
#non-unique values when setting 'row.names': ‘1_33175285_FALSE’, ‘1_33175285_TRUE’, ‘1_73453918_FALSE’, ‘1_73453918_TRUE’, ‘11_54706240_FALSE’, ‘6_166372059_FALSE’

subset(sub_dat, Chrom == 1 & Pos == 33175285)


which(colnames(dat)=="GeneName") #Way to find the column that corresponds

sub_dat <- dat[,c(1:5, 10, 79, 81:90, 96, 115:120, 125)] #get rid of unnecessary columns
sub_dat$Genotypes <- as.character(sub_dat$Genotypes) #convert genotypes to characters
sub_dat$snpID <- paste(sub_dat$Chrom, sub_dat$Pos, sep="_")

###Alternative method that will include hets
sameGeno <- subset(sub_dat, 
                   (substr(Genotypes, 1, 2) == substr(Genotypes, 3,4) & 
                    substr(Genotypes, 1, 2) == substr(Genotypes, 5,6) &
                    substr(Genotypes, 1, 2) == substr(Genotypes, 7,8) &
                    substr(Genotypes, 1, 2) == substr(Genotypes, 9,10)
                   ))
              

sameGeno2 <- subset(sub_dat, 
                     ((substr(Genotypes, 1, 2) == substr(Genotypes, 3,4) | substr(Genotypes, 1, 2) == substr(Genotypes, 4,3)) & 
                      (substr(Genotypes, 1, 2) == substr(Genotypes, 5,6) | substr(Genotypes, 1, 2) == substr(Genotypes, 6,5)) &
                      (substr(Genotypes, 1, 2) == substr(Genotypes, 7,8) | substr(Genotypes, 1, 2) == substr(Genotypes, 8,7)) &
                      (substr(Genotypes, 1, 2) == substr(Genotypes, 9,10) | substr(Genotypes, 1, 2) == substr(Genotypes, 10,9))
                   ))
                   

###Allow for one missing genotype
sameAm <- subset(sub_dat, (
    substr(Genotypes, 1, 10) == "AAAAAAAAAA" |
    substr(Genotypes, 1, 10) == "00AAAAAAAA" |
    substr(Genotypes, 1, 10) == "AA00AAAAAA" |
    substr(Genotypes, 1, 10) == "AAAA00AAAA" |
    substr(Genotypes, 1, 10) == "AAAAAA00AA" |
    substr(Genotypes, 1, 10) == "AAAAAAAA00"
  ))
  
sameTm <- subset(sub_dat, (
  substr(Genotypes, 1, 10) == "TTTTTTTTTT" |
    substr(Genotypes, 1, 10) == "00TTTTTTTT" |
    substr(Genotypes, 1, 10) == "TT00TTTTTT" |
    substr(Genotypes, 1, 10) == "TTTT00TTTT" |
    substr(Genotypes, 1, 10) == "TTTTTT00TT" |
    substr(Genotypes, 1, 10) == "TTTTTTTT00"
))

sameCm <- subset(sub_dat, (
  substr(Genotypes, 1, 10) == "CCCCCCCCCC" |
    substr(Genotypes, 1, 10) == "00CCCCCCCC" |
    substr(Genotypes, 1, 10) == "CC00CCCCCC" |
    substr(Genotypes, 1, 10) == "CCCC00CCCC" |
    substr(Genotypes, 1, 10) == "CCCCCC00CC" |
    substr(Genotypes, 1, 10) == "CCCCCCCC00"
))

sameGm <- subset(sub_dat, (
  substr(Genotypes, 1, 10) == "GGGGGGGGGG" |
    substr(Genotypes, 1, 10) == "00GGGGGGGG" |
    substr(Genotypes, 1, 10) == "GG00GGGGGG" |
    substr(Genotypes, 1, 10) == "GGGG00GGGG" |
    substr(Genotypes, 1, 10) == "GGGGGG00GG" |
    substr(Genotypes, 1, 10) == "GGGGGGGG00"
))

sameM <- rbind(sameAm, sameCm, sameTm, sameGm)


#####
sameA <- subset(sub_dat, substr(Genotypes, 1, 10) == "AAAAAAAAAA")
sameT <- subset(sub_dat, substr(Genotypes, 1, 10) == "TTTTTTTTTT")
sameC <- subset(sub_dat, substr(Genotypes, 1, 10) == "CCCCCCCCCC")
sameG <- subset(sub_dat, substr(Genotypes, 1, 10) == "GGGGGGGGGG")
same <- rbind(sameA, sameC, sameT, sameG)

same_nodups <- same[!duplicated(same[,'snpID']),] #changed from Pos to snpID so check chrom too

#to see what is removed 
require(sqldf)
diff <- sqldf('SELECT * FROM sameGeno2 EXCEPT SELECT * FROM same')
subset(same, Pos == 109144544)
subset(same_nodups, Pos ==109144544)
#identified bug and fixed - was removing 'duplicated' snps that weren't truly duplicates


same_rare <- same[(is.na(same$TG_EUR) & same$Alt == substr(same$Genotypes, 1, 1)) | (same$TG_EUR <=0.1 & same$Alt == substr(same$Genotypes, 1, 1)) | (same$TG_EUR >= 0.9 & same$Ref == substr(same$Genotypes, 1, 1)), ]
same_rare_reg <- subset(same_rare, AnnoType == 'RegulatoryFeature')


eQTL <- same[!is.na(same$eQTLs),]
lung <- eQTL[grep("lung", eQTL$eQTLs), ]
blood <- eQTL[grep("WholeBlood", eQTL$eQTLs), ]

lung_CG <- lung[!is.na(lung$candReg) | !is.na(lung$candGene), ]
blood_CG <- blood[!is.na(blood$candReg) | !is.na(blood$candGene), ]

lung_rare <- lung[(is.na(lung$TG_EUR) & lung$Alt == substr(lung$Genotypes, 1, 1)) | (lung$TG_EUR <=0.1 & lung$Alt == substr(lung$Genotypes, 1, 1)) | (lung$TG_EUR >= 0.9 & lung$Ref == substr(lung$Genotypes, 1, 1)), ]
lung_rare_CG <- lung_rare[!is.na(lung_rare$candReg) | !is.na(lung_rare$candGene), ]


blood_rare <- blood[(is.na(blood$TG_EUR) & blood$Alt == substr(blood$Genotypes, 1, 1)) | (blood$TG_EUR <=0.1 & blood$Alt == substr(blood$Genotypes, 1, 1)) | (blood$TG_EUR >= 0.9 & blood$Ref == substr(blood$Genotypes, 1, 1)), ]
blood_rare_CG <- blood_rare[!is.na(blood_rare$candReg) | !is.na(blood_rare$candGene), ]


------
  
#Barreiro eQTLs & reQTLs

barreiro <- read.table("C://Users/Mary/PepLab/data/blasto/BarreiroSig_rs.txt", header=F)  
barreiro2 <- as.character(barreiro$V1)

same_barreiro <- subset(same, rs %in% barreiro2)
same_barreiro_rare <- subset(same_rare, rs %in% barreiro2)

same_barreiro2 <- subset(sameM, rs %in% barreiro2)



------
#Novel  
  
same_novel <- subset(same, isKnownVariant == "FALSE")
#same_novel_nodups <- same_novel[!duplicated(same_novel[,'Pos']),]
#same_novel_sub <- same_novel[,c(1:6,10,11,57,61,65,68,69,71,89,90,92,93,98:109)]
same_novel_CG <- subset(same_novel, candGene != "." | candReg != ".")

  
