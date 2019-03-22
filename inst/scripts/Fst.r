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

# indir <- opts[1]
# map_file <- opts[2]

# get midas merge output dirs
dirs <- list.dirs(indir, recursive = FALSE)
dirs

window_size <- 10000
step_size <- 10000

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
    
    Res <- NULL
    for(i in 1:ncol(combinations)){
      # sites <- combinations[,2]
      sites <- combinations[,i]
      
      # Select comparison
      Freq.s <- subset(Freq, Group %in% sites,
                       drop = TRUE)
      Depth.s <- subset(Depth, Group %in% sites,
                        drop = TRUE)
      
      for(contig in levels(Freq.s$Tax$ref_id)){
        # contig <- levels(Freq.s$Tax$ref_id)[1]
        
        contig_end <- max(Freq.s$Tax$ref_pos[ Freq.s$Tax$ref_id == contig ])
        windows <- seq(from = 1, to = contig_end, by = step_size)
        for(w in windows){
          # w <- windows[1]
          
          Freq.c <- remove_taxons(Freq.s,
                                  as.character(subset(Freq.s$Tax, 
                                                      !(ref_id == contig &
                                                          ref_pos >= w &
                                                          ref_pos < w + window_size))$site_id))
          Depth.c <- remove_taxons(Depth.s,
                                  as.character(subset(Depth.s$Tax, 
                                                      !(ref_id == contig &
                                                          ref_pos >= w &
                                                          ref_pos < w + window_size))$site_id))
          
          # Clean
          Depth.c <- clean(Depth.c)
          Freq.c <- remove_samples(Freq.c, samples = setdiff(samples(Freq.c), samples(Depth.c)))
          if(length(setdiff(taxa(Freq.c), taxa(Depth.c))) > 0)
            Freq.c <- remove_taxons(Freq.c, taxons = setdiff(taxa(Freq.c), taxa(Depth.c)))
          
          if(length(table(Freq.c$Map$Group)) == 2){
            dd <- dist(t(Freq.c$Tab), 'manhattan')
            dd <- reshape2::melt(as.matrix(dd), varnames = c('row', 'col'))
            dd <- dd[ as.numeric(dd$row) > as.numeric(dd$col), ]
            dd$row.group <- Freq.c$Map[ as.character(dd$row), "Group"]
            dd$col.group <- Freq.c$Map[ as.character(dd$col), "Group"]
            dd$Intra <- dd$row.group == dd$col.group
            
            pi.b <- mean(dd$value[ !dd$Intra ])
            pi.w <- mean(dd$value[ dd$Intra ])
            
            Fst <- (pi.b - pi.w) / pi.b
            
            res <- data.frame(genome = basename(dir),
                              contig = contig,
                              start = w,
                              end = min(w + window_size, contig_end),
                              g1 = sites[1],
                              g2 = sites[2],
                              n1 = table(Freq.c$Map$Group)[sites[1]],
                              n2 = table(Freq.c$Map$Group)[sites[2]],
                              nsites = length(taxa(Freq.c)),
                              pi.b = pi.b,
                              pi.w = pi.w,
                              Fst = Fst,
                              row.names = NULL)
          }else{
            res <- data.frame(genome = basename(dir),
                              contig = contig,
                              start = w,
                              end = min(w + window_size, contig_end),
                              g1 = sites[1],
                              g2 = sites[2],
                              n1 = table(Freq.c$Map$Group)[sites[1]],
                              n2 = table(Freq.c$Map$Group)[sites[2]],
                              nsites = length(taxa(Freq.c)),
                              pi.b = NA,
                              pi.w = NA,
                              Fst = NA,
                              row.names = NULL)
          }
          Res <- rbind(Res, res)
        }
      }
    }
    
    filename <- paste(basename(dir),
                      "_Fst.txt", sep = "")
    write.table(Res, filename, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}
