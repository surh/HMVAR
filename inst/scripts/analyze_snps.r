# (C) Copyright 2019 Sur Herrera Paredes
# 
# This file is part of HMVAR.
# 
# HMVAR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# HMVAR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with HMVAR.  If not, see <http://www.gnu.org/licenses/>.


library(AMOR)
# setwd("~/micropopgen/exp/2018/today4/")

indir <- "~/micropopgen/data/test_data/midas/merged.snps/"
map_file <- "~/micropopgen/data/test_data/midas/map.txt"

indir <- opts[1]
map_file <- opts[2]

# get midas merge output dirs
dirs <- list.dirs(indir, recursive = FALSE)
dirs

Res <- NULL
for(dir in dirs){
  # dir <- dirs[2]
  
  snp_freq_file <- paste(dir, "snps_freq.txt", sep = "/")
  snp_depth_file <- paste(dir, "snps_depth.txt", sep = "/")
  snp_info_file <- paste(dir, "snps_info.txt", sep = "/")
  
  # Read data
  Freq <- read.table(snp_freq_file,
                     header = TRUE, row.names = 1)
  Depth <- read.table(snp_depth_file,
                      header = TRUE, row.names = 1)
  Info <- read.table(snp_info_file,
                     header = TRUE, sep = "\t")
  row.names(Info) <- as.character(Info$site_id)
  
  # Read map
  Map <- read.table(map_file, sep = "\t", header = TRUE)
  row.names(Map) <- as.character(Map$ID)
  Map <- droplevels(Map[ colnames(Freq), ])
  
  # Create datasets
  Freq <- create_dataset(Freq, Map, Info)
  Depth <- create_dataset(Depth, Map, Info)
  
  # Get combinations
  if(length(levels(Map$Group)) > 1){
    combinations <- combn(levels(Map$Group),m = 2)
    
    for(i in 1:ncol(combinations)){
      sites <- combinations[,i]
      
      # Select comparison
      Freq.s <- subset(Freq, Group %in% sites,
                       drop = TRUE)
      Depth.s <- subset(Depth, Group %in% sites,
                        drop = TRUE)
      
      # Clean
      Depth.s <- clean(Depth.s)
      Freq.s <- remove_samples(Freq.s, samples = setdiff(samples(Freq.s), samples(Depth.s)))
      Freq.s <- remove_taxons(Freq.s, taxons = setdiff(taxa(Freq.s), taxa(Depth.s)))
      
      if(all(table(Freq.s$Map$Group) > 2)){
        # PCA
        Dat.pca <- PCA(Freq.s)
        p1 <- plotgg(Dat.pca, col = "Group")
        # p1
        filename <- paste(basename(dir),
                          "_",
                          gsub(pattern = " ", replacement = ".",
                               paste(combinations[,i], collapse = "_")),
                          "_snps.pca.svg", sep = "")
        ggsave(filename, p1, width = 5, height = 5)
        
        m1 <- lm(PC1 ~ Group, data = p1$data)
        m1 <- summary(m1)
        m2 <- lm(PC2 ~ Group, data = p1$data)
        m2 <- summary(m2)
        m3 <- lm(PC3 ~ Group, data = p1$data)
        m3 <- summary(m3)
        
        res <- data.frame(Genome = basename(dir),
                          Group1 = sites[1],
                          Group2 = sites[2],
                          n1 = table(Freq.s$Map$Group)[sites[1]],
                          n2 = table(Freq.s$Map$Group)[sites[2]],
                          PC1.pval = m1$coefficients[2,4],
                          PC2.pval = m2$coefficients[2,4],
                          PC3.pval = m3$coefficients[2,4],
                          row.names = NULL)
      }else{
        res <- data.frame(Genome = basename(dir),
                          Group1 = sites[1],
                          Group2 = sites[2],
                          n1 = table(Freq.s$Map$Group)[sites[1]],
                          n2 = table(Freq.s$Map$Group)[sites[2]],
                          PC1.pval = NA,
                          PC2.pval = NA,
                          PC3.pval = NA,
                          row.names = NULL)
      }
      Res <- rbind(Res, res)
    }
  }
}
write.table(Res, file = "differentiated_strains.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# Res <- NULL
# for(gene in levels(Info$gene_id)){
#   # gene <- levels(Info$gene_id)[1]
# 
#   # Res <- rbind(Res, data.frame(gene = gene, D = 0, P = 0))
#   snps <- taxa(Dat)[ Dat$Tax$gene_id == gene & !is.na(Dat$Tax$gene_id == gene) ]
#   to_remove <- setdiff(taxa(Dat), snps)
#   
#   Dat.sub <- remove_taxons(Dat, taxons = to_remove)
#   Dat2.sub <- remove_taxons(Dat2, taxons = to_remove)
#   Dat2.sub <- clean(Dat2.sub)
#   
#   Dat.sub <- remove_samples(Dat.sub, samples = setdiff(samples(Dat.sub), samples(Dat2.sub)))
#   
#   Dat.sub <- subset(Dat.sub, Group %in% c("Supragingival plaque", "Buccal mucosa"), drop = TRUE)
#   Dat2.sub <- subset(Dat2.sub, Group %in% c("Supragingival plaque", "Buccal mucosa"), drop = TRUE)
#   
#   
#   Dat.pca <- PCA(Dat.sub)
#   p1 <- plotgg(Dat.pca, col = "Group")
#   # p1
#   # barplot(rowSums(Dat2.sub$Tab > 0))
#   # 
#   
#   for(i in 1:length(taxa(Dat.sub))){
#     depth <- Dat2.sub$Tab[i,]
#     freq <- Dat.sub$Tab[i,]
#     groups <- Dat.sub$Map$Group
#     
#     freq <- freq[ depth > 0 ]
#     groups <- groups[ depth > 0 ]
#     
#     # all(names(depth) == names(freq))
#     # (depth > 0) * (groups == "Buccal mucosa")
#     
#     tab <- ftable(groups ~ factor(freq > 0.5))
#     
#     if(sum(diag(tab)) == 0 | sum(diag(tab )) == length(groups)){
#       # Res$P[ Res$gene == gene ] <- Res$P[ Res$gene == gene ] + 1
#       res <- data.frame(Site = taxa(Dat.sub)[i], nsamples = length(groups), DN = "D",
#                         major_allele = Dat.sub$Tax[i, "major_allele"],
#                         minor_allele = Dat.sub$Tax[i, "minor_allele"],
#                         aa = Dat.sub$Tax[i, "amino_acids"],
#                         gene_id = gene)
#     }else{
#       # Res$D[ Res$gene == gene ] <- Res$D[ Res$gene == gene ] + 1
#       res <- data.frame(Site = taxa(Dat.sub)[i], nsamples = length(groups), DN = "P",
#                         major_allele = Dat.sub$Tax[i, "major_allele"],
#                         minor_allele = Dat.sub$Tax[i, "minor_allele"],
#                         aa = Dat.sub$Tax[i, "amino_acids"],
#                         gene_id = gene)
#     }
#     Res <- rbind(Res,res)
#   }
# 
#   
#   
#   
# }

