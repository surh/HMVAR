library(AMOR)
# setwd("~/micropopgen/exp/2018/today4/")

indir <- "~/micropopgen/data/test_data/midas/merged.snps/"
map_file <- "~/micropopgen/data/test_data/midas/map.txt"

# indir <- opts[1]
# map_file <- opts[2]

# get midas merge output dirs
dirs <- list.dirs(indir, recursive = FALSE)
dirs

Res <- NULL
for(dir in dirs){
  dir <- dirs[2]
  
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
      sites <- combinations[,2]
      sites <- combinations[,i]
      
      # Select comparison
      Freq.s <- subset(Freq, Group %in% sites,
                       drop = TRUE)
      Depth.s <- subset(Depth, Group %in% sites,
                        drop = TRUE)
      
      for(contig in levels(Freq.s$Tax$ref_id)){
        contig <- levels(Freq.s$Tax$ref_id)[1]
        
        contig_end <- max(Freq.s$Tax$ref_pos[ Freq.s$Tax$ref_id == contig ])
        windows <- seq(from = 1, to = contig_end, by = 100)
        for(w in windows){
          w <- windows[1]
          
          Freq.c <- remove_taxons(Freq.s,
                                  as.character(subset(Freq.s$Tax, 
                                                      !(ref_id == contig &
                                                          ref_pos >= w &
                                                          ref_pos < w + 1000))$site_id))
          Depth.c <- remove_taxons(Depth.s,
                                  as.character(subset(Depth.s$Tax, 
                                                      !(ref_id == contig &
                                                          ref_pos >= w &
                                                          ref_pos < w + 1000))$site_id))
          
          # Clean
          Depth.c <- clean(Depth.c)
          Freq.c <- remove_samples(Freq.c, samples = setdiff(samples(Freq.c), samples(Depth.c)))
          Freq.c <- remove_taxons(Freq.c, taxons = setdiff(taxa(Freq.c), taxa(Depth.c)))
          
                                  
          dd <- dist(t(Freq.c$Tab), 'manhattan')
          dd <- reshape2::melt(as.matrix(dd), varnames = c('row', 'col'))
          dd <- dd[ as.numeric(dd$row) > as.numeric(dd$col), ]
          dd$row.group <- Freq.c$Map[ as.character(dd$row), "Group"]
          dd$col.group <- Freq.c$Map[ as.character(dd$col), "Group"]
          dd$Intra <- dd$row.group == dd$col.group
          
          pi.b <- mean(dd$value[ !dd$Intra ])
          pi.w <- mean(dd$value[ dd$Intra ])
          
          Fst <- (pi.b - pi.w) / pi.b
        }
        
        
        
      }
      
      
      
      
      
      
      
      
      if(all(table(Freq.s$Map$Group) > 2)){
        # PCA
        Dat.pca <- PCA(Freq.s)
        p1 <- plotgg(Dat.pca, col = "Group")
        # p1
        filename <- paste(basename(dir),
                          "_",
                          gsub(pattern = " ", replacement = ".",
                               paste(combinations[,1], collapse = "_")),
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
