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


library(HMAR)
library(ggtree)
library(phytools)
library(reshape2)

setwd("~/micropopgen/exp/2018/today/")

# treefile <- "~/micropopgen/exp/2018/today2/RAxML_bestTree.result_reroot.tre"
treefile <- "~/micropopgen/exp/2018/today2/RAxML_bestTree.result"
# treefile <- "~/micropopgen/exp/2018/today2/RAxML_bestTree.result_reroot.nex"
# treefile <- "~/micropopgen/exp/2018/2018-04-15.tested_species_tree/fastree_result_reroot.tre"
tree <- read.tree(treefile)
tree <- ggtree::reroot(tree, node = 110)
tree$tip.label <- sub(pattern = ".CDS_0$", replacement = "", x = tree$tip.label)


p1 <- ggtree::ggtree(tree, ladderize = FALSE) +
  # ggtree::geom_text2(aes(subset=!isTip, label=node)) +
  ggtree::geom_tiplab(align = TRUE, 
                      size = 2, linetype = 'dotted',
                      linesize = 0.2)
p1
ggsave('phylogeny_raxml.svg', p1, width = 10, height = 10)

##### Get neutralitty indice ####
dir <- "~/micropopgen/exp/2018/2018-04-13.analyze_mktest_fraserv/"
files <- list.files(dir)
files <- files[ grep(pattern = "^NI_", x = files) ]

Tab <- NULL
for(f in files){
  # f <- files[1]
  infile <- paste(dir, "/", f, sep = "")
  cat(infile, "\n")
  comp <- sub(pattern = ".txt$", replacement = "", x = sub(pattern = "^NI_", replacement = "", f))
  tab <- read.table(infile, header = TRUE, sep = "\t")
  tab$Comparison <- comp
  Tab <- rbind(Tab, tab)
}
Tab <- acast(Species ~ Comparison, data = Tab, value.var = "NI")
Tab <- as.data.frame(Tab)
head(Tab)
extra <- setdiff(tree$tip.label, row.names(Tab))
Tab <- rbind(Tab, matrix(NA, nrow = length(extra), ncol = ncol(Tab),
                         dimnames = list(extra, colnames(Tab))))
Tab <- Tab[ tree$tip.label, ]
head(Tab)

Tab$Species <- row.names(Tab)
# row.names(Tab) <- NULL
Tab$Genus <- sapply(strsplit(x = Tab$Species, split = "_"), function(x) x[1])
head(Tab)
# Tab$Phylum <- NA
# Tab$Phylum[ Tab$Genus %in% c()]

p1 <- ggtree::ggtree(tree, ladderize = FALSE) +
  geom_tiplab(align = TRUE, size = 1, linetype = 'dotted',
              linesize = 0.2)
p1 <- p1 %>% gheatmap(Tab[,1:6],
                      colnames = TRUE,
                      colnames_position = 'top',
                      colnames_offset_y = 40,
                      colnames_angle = 90,
                      color = NA,
                      offset = 0.1) +
  scale_fill_gradient2(low = "#d01c8b", mid = 'white',
                       high = '#4dac26', midpoint = 1,
                       na.value = 'grey')
p1
ggsave("phylo_NI.svg", p1, width = 10, height = 6)


####### Get Pi0 ####

dir <- "~/micropopgen/exp/2018/2018-04-13.analyze_mktest_fraserv/"
dirs <- list.dirs(dir, recursive = FALSE)
dirs <- basename(dirs)
dirs <- dirs[ grep(pattern = '^pvals_', x = dirs) ]

Tab <- NULL
for(d in dirs){
  # d <- dirs[1]
  cat(d, "\n")
  comp <- paste(strsplit(d, split = "_")[[1]][-1], collapse = "_")
  
  infile <- paste(dir, "/", d, '/Pi0_summmary.txt', sep ="")
  tab <- read.table(infile, sep = "\t", header = TRUE)
  tab <- tab[, c('Species', 'ratio')]
  tab$Comparison <- comp
  Tab <- rbind(Tab, tab)
}
Tab <- acast(Species ~ Comparison, data = Tab, value.var = "ratio")
Tab <- as.data.frame(Tab)
extra <- setdiff(tree$tip.label, row.names(Tab))
Tab <- rbind(Tab, matrix(NA, nrow = length(extra), ncol = ncol(Tab),
                         dimnames = list(extra, colnames(Tab))))
Tab <- Tab[ tree$tip.label, ]
head(Tab)

Tab$Species <- row.names(Tab)
Tab$Genus <- sapply(strsplit(x = Tab$Species, split = "_"), function(x) x[1])
head(Tab)


p1 <- ggtree::ggtree(tree, ladderize = FALSE) +
  geom_tiplab(align = TRUE, size = 1, linetype = 'dotted',
              linesize = 0.2)
p1 <- p1 %>% gheatmap(Tab[,1:6],
                      colnames = TRUE,
                      colnames_position = 'top',
                      colnames_offset_y = 40,
                      colnames_angle = 90,
                      color = NA,
                      offset = 0.1) +
  scale_fill_gradient2(low = "#d01c8b", mid = 'white',
                       high = '#4dac26', midpoint = 1,
                       na.value = 'grey')
p1
ggsave('Pi0_phylo.svg', p1, width = 10, height = 6)
