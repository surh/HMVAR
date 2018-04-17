library(HMAR)
library(ggtree)
library(phytools)

setwd("~/micropopgen/exp/2018/today/")

# treefile <- "~/micropopgen/exp/2018/today2/RAxML_bestTree.result_reroot.tre"
treefile <- "~/micropopgen/exp/2018/today2/RAxML_bestTree.result"
# treefile <- "~/micropopgen/exp/2018/today2/RAxML_bestTree.result_reroot.nex"
# treefile <- "~/micropopgen/exp/2018/2018-04-15.tested_species_tree/fastree_result_reroot.tre"
tree <- read.tree(treefile)
tree <- ggtree::reroot(tree, node = 110)
tree$tip.label <- sub(pattern = ".CDS_0$", replacement = "", x = tree$tip.label)


p1 <- ggtree::ggtree(tree, ladderize = FALSE) +
  ggtree::geom_text2(aes(subset=!isTip, label=node)) +
  ggtree::geom_tiplab()
p1



## Get neutralitty indice
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




