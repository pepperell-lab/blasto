#Working with Geuvadis data

dat <- read.table("C:/Users/Mary/PepLab/data/blasto/Geuvadis/IL6_lncRNA_RPKM.txt", header=T, sep='\t')
str(dat)

require(reshape2)

dat.m <- melt(dat, id.vars=c("TargetID", "Gene_Symbol", "Chr", "Coord"))

require(ggplot2)
ggplot(dat.m) + geom_boxplot(aes(x=Gene_Symbol, y=value))


gts <- read.table("C:/Users/Mary/PepLab/data/blasto/Geuvadis/snps_of_interest.ped", header=F, sep="\t", na.strings=NA)

dat.m$GT <- gts[match(dat.m$variable, gts$V1), 7:14]



dat.m$rs1800796 <- (paste(dat.m$GT$V7, dat.m$GT$V8, sep=""))
dat.m$rs1800796[dat.m$rs1800796 == "CG"] <- "GC"
dat.m$rs1800796[dat.m$rs1800796 == "NANA"] <- NA
dat.m$rs1800796 <- as.factor(dat.m$rs1800796)
ggplot(dat.m) + geom_boxplot(aes(x=rs1800796, y=value)) + facet_wrap(~Gene_Symbol)

dat.m$rs7802308 <- paste(dat.m$GT$V9, dat.m$GT$V10, sep="")
dat.m$rs7802308[dat.m$rs7802308 == "AT"] <- "TA"
dat.m$rs7802308[dat.m$rs7802308 == "NANA"] <- NA
dat.m$rs7802308 <- as.factor(dat.m$rs7802308)
ggplot(dat.m) + geom_boxplot(aes(x=rs7802308, y=value)) + facet_wrap(~Gene_Symbol)

dat.m$rs1524107 <- paste(dat.m$GT$V11, dat.m$GT$V12, sep="")
dat.m$rs1524107[dat.m$rs1524107 == "TC"] <- "CT"
dat.m$rs1524107[dat.m$rs1524107 == "NANA"] <- NA
dat.m$rs1524107 <- as.factor(dat.m$rs1524107)
ggplot(dat.m) + geom_boxplot(aes(x=rs1524107, y=value)) + facet_wrap(~Gene_Symbol)

dat.m$rs2066992 <- paste(dat.m$GT$V13, dat.m$GT$V14, sep="")
dat.m$rs2066992[dat.m$rs2066992 == "TG"] <- "GT"
dat.m$rs2066992[dat.m$rs2066992 == "NANA"] <- NA
dat.m$rs2066992 <- as.factor(dat.m$rs2066992)
ggplot(dat.m) + geom_boxplot(aes(x=rs2066992, y=value)) + facet_wrap(~Gene_Symbol)


###
dat.t.p <- dat[,c(2, 5:466)]
dat.t <- as.data.frame(t(dat.t.p))
dat.t <- dat.t[-1,]


ggplot(dat.t) + geom_point(aes(x=V1, y=V2))


##import the population information


samps <- names(dat[5:length(dat)])
info <- read.table("C:/Users/Mary/PepLab/data/blasto/Geuvadis/integrated_call_samples_v3.20130502.ALL.panel", header=F, sep='\t')
names(info) <- c("sample", "pop", "super_pop", "gender")

df <- data.frame(samps)
df$pop <- info[match(df$samps, info$sample), 'pop']
df$super_pop <- info[match(df$samps, info$sample), 'super_pop']
df$gender <- info[match(df$samps, info$sample), 'gender']

keep <- dat.m[, c('Gene_Symbol', 'variable', 'value', 'rs1800796', 'rs7802308', 'rs1524107', 'rs2066992')]
keep <- na.omit(keep)

keep$Gene_Symbol <- as.character(keep$Gene_Symbol)
keep$Gene_Symbol[keep$Gene_Symbol == "ENSG00000136244.7"] <- "IL6"
keep$Gene_Symbol[keep$Gene_Symbol == "ENSG00000179428.2"] <- "AC073072.5"
keep$Gene_Symbol <- as.factor(keep$Gene_Symbol)

keep$pop <- df[match(keep$variable, df$samps), 'pop']
keep$super_pop <- df[match(keep$variable, df$samps), 'super_pop']
keep$gender <- df[match(keep$variable, df$samps), 'gender']


ggplot(keep) + geom_boxplot(aes(x=rs1800796, y=value)) + facet_wrap(~Gene_Symbol)


IL6 <- keep[as.character(keep$Gene_Symbol) == "IL6",]
gg <- IL6[as.character(IL6$rs1800796) == "GG", 'value']
gc <- IL6[as.character(IL6$rs1800796) == "GC", 'value']
  
AC073072.5 <- keep[as.character(keep$Gene_Symbol) == "AC073072.5",]
hom <- AC073072.5[as.character(AC073072.5$rs1800796) == "GG", 'value']
het <- AC073072.5[as.character(AC073072.5$rs1800796) == "GC", 'value']
