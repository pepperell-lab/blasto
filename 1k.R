#all my code for the 1K.Rmd (I have since deleted things that didn't work, or where I found a prettier way)
suppressMessages(library(dplyr))

filepath <- "C://Users/Mary/PepLab/data/blasto/1K_genomes/SNPs_overlap_1k_combined_GTEx.txt"

TG <- read.table(filepath, header = F, sep = '\t', stringsAsFactors = FALSE, na.strings = c("NA", "."))
names(TG) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "BDH001", "BDH002", "BDH003", "BDH004", "BDH005", "BDH006", "GTEx")

TG.l <- tbl_df(TG)
glimpse(TG.l)

# A simple function to extract the AF from various super populations
extractAF <- function(pop, vec) {
  info <- unlist((strsplit(vec, ";", fixed=TRUE)))
  AF <- as.numeric(unlist(strsplit((info[grep(pop, (unlist((strsplit(vec, ";", fixed=TRUE)))))]), "=", fixed=TRUE))[2])
  return(AF)
}

# An alternative way of exctracting the AF's
all_pop <- c("AMR_AF", "AFR_AF", "EUR_AF", "SAS_AF", "EAS_AF")
pat <- paste0(".*\\b", all_pop, "=(\\d+(\\.\\d+)?)\\b.*")

out <- as.data.frame(sapply(pat, gsub, replacement="\\1", x=TG.l$INFO))
newdf <- setNames(as.data.frame(out), all_pop)

AF <- cbind.data.frame(TG.l, newdf)

AF$AMR_AF <- as.numeric(as.character(AF$AMR_AF))
AF$AFR_AF <- as.numeric(as.character(AF$AFR_AF))
AF$EUR_AF <- as.numeric(as.character(AF$EUR_AF))
AF$SAS_AF <- as.numeric(as.character(AF$SAS_AF))
AF$EAS_AF <- as.numeric(as.character(AF$EAS_AF))


### A little code to look at the distribution of FunSEQ scores

funseq <- as.data.frame(sapply(c(".*\\bAMR_AF=(\\d+(\\.\\d+)?)\\b.*"), gsub, replacement="\\1", x=TG.l$INFO))
names(funseq) <- "FunSEQ"
funseq <- as.numeric(as.character(funseq$FunSEQ))

hist(funseq)

extractGT <- function(row) {
  gt <- gsub("/", "", paste(unlist(strsplit(AF[row,'BDH001'], ":", fixed=TRUE))[1], unlist(strsplit(AF[row,'BDH002'], ":", fixed=TRUE))[1], unlist(strsplit(AF[row,'BDH003'], ":", fixed=TRUE))[1], unlist(strsplit(AF[row,'BDH004'], ":", fixed=TRUE))[1], unlist(strsplit(AF[row,'BDH005'], ":", fixed=TRUE))[1], sep=""))
  return(gt)
}


extractGT <- function(x) {
  gt <- gsub("/", "", paste(unlist(strsplit(x[,'BDH001'], ":", fixed=TRUE))[1], unlist(strsplit(x[,'BDH002'], ":", fixed=TRUE))[1], unlist(strsplit(x[,'BDH003'], ":", fixed=TRUE))[1], unlist(strsplit(x[,'BDH004'], ":", fixed=TRUE))[1], unlist(strsplit(x[,'BDH005'], ":", fixed=TRUE))[1], sep=""))
  return(gt)
}



#gt <- gsub("/", "", paste(unlist(strsplit(AF[1,'BDH001'], ":", fixed=TRUE))[1], unlist(strsplit(AF[1,'BDH002'], ":", fixed=TRUE))[1], unlist(strsplit(AF[1,'BDH003'], ":", fixed=TRUE))[1], unlist(strsplit(AF[1,'BDH004'], ":", fixed=TRUE))[1], unlist(strsplit(AF[1,'BDH005'], ":", fixed=TRUE))[1], sep=""))

#sampAF <- AF[1:10,]


#apply(sampAF, 1, function(x) gsub("/", "", paste(unlist(strsplit(x[,'BDH001'], ":", fixed=TRUE))[1], unlist(strsplit(x[,"BDH002"], ":", fixed=TRUE))[1], unlist(strsplit(x[,"BDH003"], ":", fixed=TRUE))[1], unlist(strsplit(x["BDH004"], ":", fixed=TRUE))[1], unlist(strsplit(x[,"BDH005"], ":", fixed=TRUE))[1], sep="")))




AFgt <- AF %>%
  rowwise() %>%
  mutate(gt = gsub("/", "", paste(unlist(strsplit(BDH001, ":", fixed=TRUE))[1], unlist(strsplit(BDH002, ":", fixed=TRUE))[1], unlist(strsplit(BDH003, ":", fixed=TRUE))[1], unlist(strsplit(BDH004, ":", fixed=TRUE))[1], unlist(strsplit(BDH005, ":", fixed=TRUE))[1], sep="")))

# Filter for loci where all 5 cases have the same genotype (including het)
sameGT <- subset(AFgt, 
                 ((substr(gt, 1, 2) == substr(gt, 3,4) | substr(gt, 1, 2) == substr(gt, 4,3)) & 
                    (substr(gt, 1, 2) == substr(gt, 5,6) | substr(gt, 1, 2) == substr(gt, 6,5)) &
                    (substr(gt, 1, 2) == substr(gt, 7,8) | substr(gt, 1, 2) == substr(gt, 8,7)) &
                    (substr(gt, 1, 2) == substr(gt, 9,10) | substr(gt, 1, 2) == substr(gt, 10,9))
                 ))

# Lets try to achieve the same thing with dplyr
same <- filter(AFgt, (substr(gt, 1, 2) == substr(gt, 3,4) | substr(gt, 1, 2) == substr(gt, 4,3)) & 
                 (substr(gt, 1, 2) == substr(gt, 5,6) | substr(gt, 1, 2) == substr(gt, 6,5)) &
                 (substr(gt, 1, 2) == substr(gt, 7,8) | substr(gt, 1, 2) == substr(gt, 8,7)) &
                 (substr(gt, 1, 2) == substr(gt, 9,10) | substr(gt, 1, 2) == substr(gt, 10,9)))



eQTL <- same[!is.na(same$GTEx),]
lung <- eQTL[grep("Lung", eQTL$GTEx), ]
blood <- eQTL[grep("WholeBlood", eQTL$GTEx), ]

# Rare (0.1) & demonstrated eQTL from GTEx
eQTL_rare <- same_rare[!is.na(same_rare$GTEx),]
lung_rare <- eQTL_rare[grep("Lung", eQTL_rare$GTEx), ]
blood_rare <- eQTL_rare[grep("WholeBlood", eQTL_rare$GTEx), ]

# Rare (o.2) & demonstrated eQTL from GTEx
eQTL_rare_20 <- same_rare_20[!is.na(same_rare_20$GTEx),]
lung_rare_20 <- eQTL_rare_20[grep("Lung", eQTL_rare_20$GTEx), ]
blood_rare_20 <- eQTL_rare_20[grep("WholeBlood", eQTL_rare_20$GTEx), ]





#IMMVAR

DCs <- read.table("C://Users/Mary/PepLab/data/blasto/1K_genomes/S4_Lee_collated.txt", stringsAsFactors = FALSE)

dcs <- DCs$V1

# test vector with 3 SNPs I know are in the data frame
testimm <- c("rs4648959", "rs662949", "rs1358035")

# Cross ref same_rare with immvar QTL's

data.frame(filter(same, ID %in% imm))




MOs <- read.table("C://Users/Mary/Downloads/tableS11_meta_monocytes_cis_fdr05.tsv", stringsAsFactors = FALSE)

mos <- MOs$V1

mos_rare <- data.frame(filter(same_rare, ID %in% mos))

mos_rare[,c(1:5, 17:22)]



CD4T <- read.table("C://Users/Mary/Downloads/tableS12_meta_cd4T_cis_fdr05.tsv", stringsAsFactors = FALSE)

cd4t <- CD4T$V1

cd4t_rare <- data.frame(filter(same_rare, ID %in% cd4t))

cd4t_rare[,c(1:5, 17:22)]



tcell <- read.table("C://Users/Mary/Downloads/S7_tcell.xlsx", stringsAsFactors = FALSE)

tc <- tcell$V1

tc_rare <- data.frame(filter(same_rare, ID %in% tc))
####

clec6a <- filter(same_rare, POS >= 8605000 & POS <= 8630926)
