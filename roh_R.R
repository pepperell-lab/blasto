roh <- read.table("overlapping_ROH.txt", header=F, sep='\t')
cases <- roh[roh$V5 == 9,]

require(ggplot2)
require(dplyr)

chr.sum <- cases %>%
  group_by(V1) %>%
  summarise(mean = sum(V4))

ggplot(chr.sum) + geom_point(aes(x=V1, y=mean))


#read individual roh identifications

dat <- read.table("bcftools_garlic_combined.roh", header=F, sep="\t", na.strings="NA")
names(dat) <- c("chrom", "start", "stop", "ind", "method", "length", "class", "no_snps", "qual", "D1", "D2", "D3")



total <- dad%
  group_by(ind) %>%
  summarise(total = sum(length))

p <- ggplot(dat, aes(ind, length))
p + geom_violin()
