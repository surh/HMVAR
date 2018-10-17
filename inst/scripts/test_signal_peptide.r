


#########################################


pos <- read.table(file = "gff_grampos.txt", sep = "\t", stringsAsFactors = FALSE)
head(pos)
neg <- read.table(file = "gff_gramneg.txt", sep = "\t", stringsAsFactors = FALSE)
head(neg)

signalp <- c(neg$V1, pos$V1)
signalp <- unique(signalp)

Tab <- read.table(file = "signifincant_mktest", header = TRUE)
head(Tab)
ftable(A ~ B, Tab)

Tab$signalP <- "No"
Tab$signalP[ Tab$gene %in% signalp ] <- "Yes"
table(Tab$signalP)

p1 <- ggplot(Tab, aes(x = signalP, y = log2(ratio))) +
  # geom_boxplot() +
  geom_violin()
p1

# summary(lm(log2(ratio) ~ signalP, data = Tab))
