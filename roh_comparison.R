roh <- read.table("C:/Users/Mary/PepLab/data/blasto/ROH/ROH_50SNPs_3allowance_50kb.hom.overlap", header=T)

require(dplyr)
chr7 <- filter(roh, CHR == 7, FID != "UNION" & FID != "CON")
il6  <- filter(chr7, BP1 >= 22618338 & BP1 <= 22919561)
il6$method = "pilot_plink"

plink <- read.table("C:/Users/Mary/PepLab/data/blasto/postdoc/ROH_IL6region.overlap", header=F, sep="\t")
names(plink) <- c("POOL", "FID", "IID", "PHE", "CHR", "SNP1", "SNP2", "BP1", "BP2", "KB", "NSNP", "NSIM", "GRP")

garlic <- read.table("C:/Users/Mary/PepLab/data/blasto/postdoc/auto_test_garlic_il6region.txt", header=F, sep="\t")
names(garlic) <- c("IID", "CHR", "BP1", "BP2", "Size_Class", "ROH_length", "placeholder", "rgb_track_color")

bcf <- read.table("C:/Users/Mary/PepLab/data/blasto/postdoc/bcftool_il6region.txt", header=F, sep="\t")
names(bcf) <- c("IID", "BP1", "BP2")
bcf$IID <- gsub('-', '', bcf$IID)

gar$method <- "garlic"
plink$method <- "plink"
bcf$method <- "bcftools"


plink <- filter(plink, CHR == "7", FID != "UNION" & FID != "CON")
plink <- plink[,c(3:12)]
plink <- distinct(plink)

garlic$method <- "garlic"
plink$method <- "plink"
bcf$method <- "bcftools"

comp <- rbind(il6[,c("IID", "BP1", "BP2", "method")],
              garlic[,c("IID", "BP1", "BP2", "method")], 
              plink[,c("IID", "BP1", "BP2", "method")],
              bcf[,c("IID", "BP1", "BP2", "method")]
              )
comp$id <- paste(comp$IID, comp$method, sep="-")
comp$id <- as.factor(comp$id)

require(ggplot2)

comp <- filter(comp, BP2 < 23000000)
roh.p <- ggplot(comp) +
  geom_segment(aes(x=BP1, y=id, xend=BP2, yend=id, colour=method)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept= 22766246, linetype="dotted", color="black") +
  geom_vline(xintercept= 22768219, linetype="dotted", color="black") +
  geom_vline(xintercept= 22768249, linetype="dotted", color="black") +
  geom_vline(xintercept= 22799526, linetype="dotted", color="black") +
  #scale_x_continuous(limits=c(22618338,22919561)) +
  xlab("Chromosome 7") + 
  ylab("")
