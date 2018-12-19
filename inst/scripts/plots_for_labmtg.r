library(ggplot2)
library(AMOR)
library(reshape2)

setwd("~/micropopgen/exp/2018/today4")

Map <- read.table("~/micropopgen/data/test_data/midas/map.txt", header = TRUE, sep = "\t")
head(Map)


p1 <- ggplot(Map, aes(x = Group, col = Group, fill = Group)) +
  geom_bar() + 
  AMOR::theme_blackbox +
  theme(axis.text.x = element_text(angle = 90))
p1
ggsave('sample_num_barplot.svg', p1, width = 5, height = 5)

rm(list=ls())


######
Cov <- read.table("~/micropopgen/data/test_data/midas/coverage.txt", header = TRUE, sep = "\t", row.names = 1)
head(Cov)
Map <- read.table("~/micropopgen/data/test_data/midas/map.txt", header = TRUE, sep = "\t")
row.names(Map) <- as.character(Map$ID)
head(Map)

Map <- Map[colnames(Cov), ]
Map$Cov_0 <- colSums(Cov > 0)
Map$Cov_0.1 <- colSums(Cov > 0.1)
Map$Cov_0.5 <- colSums(Cov > 0.5)
Map$Cov_1 <- colSums(Cov > 1)
Map$Cov_3 <- colSums(Cov > 3)
head(Map)


dat <- melt(Map, id.vars = c("ID", "Group"), variable.name = 'Coverage', value.name = "number.species")
head(dat)

p1 <- ggplot(dat,aes(x = Coverage, y = number.species, color = Group, group = ID)) +
  facet_wrap(~ Group) + 
  geom_line(alpha = 0.1) +
  geom_hline(yintercept = 100) +
  AMOR::theme_blackbox +
  theme(axis.text.x = element_text(angle = 90))
p1
ggsave('species_per_coverage.svg', p1, width = 6, height = 6)


Cov <- Cov[ rowSums(Cov >= 3) > 0, ]
Cov <- Cov[ ,colSums(Cov >= 3) > 0 ]
Cov <- Cov >= 3
Map <- Map[colnames(Cov),]

dat <- cbind(Map[,c("ID", "Group")], t(Cov))
head(dat)
dat[1:5,1:5]

dat <- melt(dat, id.vars = c("ID", "Group"),
            variable.name = 'Strain',
            value.name = 'present')
head(dat)
s <- aggregate(present ~ Strain, data = dat, FUN = sum)
dat$Strain <- factor(dat$Strain, levels = as.character(s[order(s$present, decreasing = TRUE),"Strain"]))
dat$Strain <- factor(dat$Strain, levels = sort(levels(dat$Strain)))

p1 <- ggplot(dat,aes(x = ID, y = Strain)) +
  facet_grid(~ Group, scales = "free_x", space = 'free_x') + 
  geom_tile(aes(fill=present), col = NA) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p1
ggsave('heatmap_present.png', p1, width = 12, height = 8, dpi = 200)

rm(list = ls())

####
MK <- read.table("../today3/all_mktest", sep = "\t", header = TRUE)
head(MK)

dir <- "../today3/all_mktests/signalP/"
files <- list.files(dir)

for(f in files){
  # f <- files[1]
  
  comp <- strsplit(sub(pattern = ".txt$", replacement = "", x = f), split = "_")[[1]]
  
  infile <- paste(dir, "/", f, sep = "")
  SigP <- read.table(infile,
                     sep = "", skip = 1, comment.char = "", header = TRUE)
  
  mk <- subset(MK, A == comp[1] & B == comp[2])
  
  
  row.names(SigP) <- as.character(SigP$X.name)
  sum(table(as.character(SigP$X.name)) > 1)
  sum(table(as.character(mk$gene)) > 1)
  sum(is.na(SigP[ as.character(mk$gene), "X." ]))
  mk$SigP <- SigP[ as.character(mk$gene), "X." ]
  # head(mk)
  
  
  p1 <- ggplot(mk, aes(x = SigP, y = log10(ratio))) +
    # facet_wrap(~ Species) +
    geom_violin() +
    AMOR::theme_blackbox +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 16))
  # p1
  cat(f, "\n")
  print(summary(lm(log10(ratio) ~ SigP, data = mk)))
  outfile <- paste("signalp_mk/", comp[1], "_", comp[2], ".svg", sep = "")
  cat(outfile, "\n")
  ggsave(outfile, p1, width = 6, height = 6)
  
  
}

rm(list = ls())


###
Tab <- read.table("merged.species/relative_abundance.txt", sep = "\t",
                  header = TRUE, row.names = 1)
Tab[1:5, 1:5]

Map <- read.table("~/micropopgen/data/test_data/midas/map.txt", header = TRUE, sep = "\t")
row.names(Map) <- as.character(Map$ID)
Map <- Map[colnames(Tab),]
head(Map)

# Cov <- read.table("~/micropopgen/data/test_data/midas/coverage.txt", header = TRUE, sep = "\t", row.names = 1)
# Cov[row.names(Tab),colnames(Tab)]
# Cov[1:5,1:5]

Tab[ Cov < 3 ] <- 0

Dat <- create_dataset(Tab, Map)
Dat <- clean(Dat)
Dat


Dat.pca <- PCA(Dat, cor = FALSE)
p1 <- plotgg(Dat.pca, col = "Group")
p1
ggsave("PCA_RA.svg", p1, width = 6, height = 5)

rm(list=ls())

#####

Tab <- read.table("../today3/all_mktest", header = TRUE, sep = "\t")
head(Tab)


p1 <- ggplot(Tab,aes(x = log2(ratio), y = DoS)) + 
  geom_point(aes(col = -log10(ratio.pval)), size = 3) +
  AMOR::theme_blackbox +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12))
# p1
ggsave("logratio_DoS.png", p1, width = 15, height = 15)
# ggsave("logratio_DoS.svg", p1, width = 6, height = 5)

p1 <- ggplot(Tab,aes(x = log2(ratio), y = DoS)) + 
  geom_point(data = subset(Tab, ratio.qval >= 0.1), size = 5, col = "pink") +
  geom_point(data = subset(Tab, ratio.qval < 0.1), size = 5, col = "blue") +
  AMOR::theme_blackbox +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 26))
# p1
ggsave("logratio_DoS_sig.png", p1, width = 15, height = 15)


p1 <- ggplot(Tab,aes(x = log((Dn * Ps) / (Ds * Pn)), y = DoS)) + 
  geom_point(aes(col = -log10(ratio.pval))) +
  AMOR::theme_blackbox
# p1

rm(list = ls())


###

Res <- read.table("../today3/signifincant_mktest", header = TRUE, sep =)
head(Res)


Res$fisher.pval <- apply(Res[,5:8],1,function(x){
  mat <- matrix(x, ncol = 2)
  t <- fisher.test(mat, alternative = 'two.sided')
  return(t$p.value)
})
head(Res)


head(Res[order(Res$fisher.pval, decreasing = FALSE),])
