library(AMOR)

# cov <- read.am("~/micropopgen/exp/2018/2018-02-20.summarymk_cov/coverage.txt")
Freq <- read.table("~/micropopgen/data/test_data/midas/merged.snps/Veillonella_parvula_57794/snps_freq.txt",
                  header = TRUE, row.names = 1)
Depth <- read.table("~/micropopgen/data/test_data/midas/merged.snps/Veillonella_parvula_57794/snps_depth.txt",
                    header = TRUE, row.names = 1)

# Freq[1:5, 1:5]
# Depth[1:5, 1:5]

Info <- read.table("~/micropopgen/data/test_data/midas/merged.snps/Veillonella_parvula_57794/snps_info.txt",
                   header = TRUE, sep = "\t")
row.names(Info) <- as.character(Info$site_id)
head(Info)

Map <- read.table("~/micropopgen/data/test_data/midas/map.txt", sep = "\t", header = TRUE)
row.names(Map) <- as.character(Map$ID)
head(Map)
Map <- Map[ colnames(Freq), ]

Dat <- create_dataset(Freq, Map, Info)
Dat2 <- create_dataset(Depth, Map, Info)
Dat.pca <- PCA(Dat)
# p1 <- plotgg(Dat.pca, col = "Group")
# p1
# ggsave("Veillonella_parvula_57794.snps_pca.svg", p1, width = 5, height = 5)

Res <- NULL
for(gene in levels(Info$gene_id)){
  # gene <- levels(Info$gene_id)[1]

  # Res <- rbind(Res, data.frame(gene = gene, D = 0, P = 0))
  snps <- taxa(Dat)[ Dat$Tax$gene_id == gene & !is.na(Dat$Tax$gene_id == gene) ]
  to_remove <- setdiff(taxa(Dat), snps)
  
  Dat.sub <- remove_taxons(Dat, taxons = to_remove)
  Dat2.sub <- remove_taxons(Dat2, taxons = to_remove)
  Dat2.sub <- clean(Dat2.sub)
  Dat.sub <- remove_samples(Dat.sub, samples = setdiff(samples(Dat.sub), samples(Dat2.sub)))
  
  Dat.sub <- subset(Dat.sub, Group %in% c("Supragingival plaque", "Buccal mucosa"), drop = TRUE)
  Dat2.sub <- subset(Dat2.sub, Group %in% c("Supragingival plaque", "Buccal mucosa"), drop = TRUE)
  
  
  Dat.pca <- PCA(Dat.sub)
  p1 <- plotgg(Dat.pca, col = "Group")
  # p1
  # barplot(rowSums(Dat2.sub$Tab > 0))
  # 
  
  for(i in 1:length(taxa(Dat.sub))){
    depth <- Dat2.sub$Tab[i,]
    freq <- Dat.sub$Tab[i,]
    groups <- Dat.sub$Map$Group
    
    freq <- freq[ depth > 0 ]
    groups <- groups[ depth > 0 ]
    
    # all(names(depth) == names(freq))
    # (depth > 0) * (groups == "Buccal mucosa")
    
    tab <- ftable(groups ~ factor(freq > 0.5))
    
    if(sum(diag(tab)) == 0 | sum(diag(tab )) == length(groups)){
      # Res$P[ Res$gene == gene ] <- Res$P[ Res$gene == gene ] + 1
      res <- data.frame(Site = taxa(Dat.sub)[i], nsamples = length(groups), DN = "D",
                        major_allele = Dat.sub$Tax[i, "major_allele"],
                        minor_allele = Dat.sub$Tax[i, "minor_allele"],
                        aa = Dat.sub$Tax[i, "amino_acids"],
                        gene_id = gene)
    }else{
      # Res$D[ Res$gene == gene ] <- Res$D[ Res$gene == gene ] + 1
      res <- data.frame(Site = taxa(Dat.sub)[i], nsamples = length(groups), DN = "P",
                        major_allele = Dat.sub$Tax[i, "major_allele"],
                        minor_allele = Dat.sub$Tax[i, "minor_allele"],
                        aa = Dat.sub$Tax[i, "amino_acids"],
                        gene_id = gene)
    }
    Res <- rbind(Res,res)
  }

  
  
  
}

