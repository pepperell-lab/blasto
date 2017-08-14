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

#Using cutoffs from garlic to assign lengths to bcftools segments
dat$bins <- cut(dat$length, breaks=c(0,97801.9,470624,25000000), labels=c("A","B","C"))
dat$clust5 <- cut(dat$length, breaks=c(0,49277.2,122474,314478,1.02379e+06,25000000), labels=c("A","B","C","D","E"))


tot <- dat %>%
  group_by(ind,method) %>%
  summarise(sum = sum(as.numeric(length)))

tot2 <- dat %>%
  group_by(ind,bins,method) %>%
  summarise(sum = sum(as.numeric(length)))

tot3 <- dat %>%
  group_by(ind,clust5,method) %>%
  summarise(sum = sum(as.numeric(length)))

total = rbind(data.frame(ind = tot$ind, bins="ALL", method=tot$method, sum=tot$sum),
              data.frame(tot2))
total2 = rbind(data.frame(ind = tot$ind, clust5="ALL", method=tot$method, sum=tot$sum),
              data.frame(tot3))

p <- ggplot(total, aes(bins, sum))
p + facet_wrap(~method) + geom_violin() + theme_bw()

q <- ggplot(total2, aes(clust5, sum))
q + facet_wrap(~method) + geom_violin() + theme_bw()




tot.c <- dat %>%
  group_by(ind,method) %>%
  summarise(n = length(length))

tot2.c <- dat %>%
  group_by(ind,bins,method) %>%
  summarise(n = length(length))

tot3.c <- dat %>%
  group_by(ind,clust5,method) %>%
  summarise(n = length(length))


total.c = rbind(data.frame(ind = tot.c$ind, bins="ALL", method=tot.c$method, n=tot.c$n),
              data.frame(tot2.c))
total2.c = rbind(data.frame(ind = tot.c$ind, clust5="ALL", method=tot.c$method, n=tot.c$n),
               data.frame(tot3.c))


r <- ggplot(total.c, aes(bins, n))
r + facet_wrap(~method) + geom_violin() + theme_bw()

s <- ggplot(total2.c, aes(clust5, n))
s + facet_wrap(~method) + geom_violin() + theme_bw()
