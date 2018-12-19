library(tidyverse)
library(HMVAR)


# 
# mktest_file <- "../2018-12-14.hmp_mktest/results/Granulicatella_adiacens_61980_mktest.txt"
# Genes <- read_tsv(mktest_file)
# gene <- Genes$gene[2]
genes <- c("638301.3.peg.1",
           "638301.3.peg.444",
           "638301.3.peg.884",
           "638301.3.peg.944",
           "638301.3.peg.955",
           "638301.3.peg.1091",
           "638301.3.peg.1273",
           "638301.3.peg.1346",
           "638301.3.peg.1361")

# Parameters
midas_dir <- "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/Granulicatella_adiacens_61980/"
depth_thres <- 1
map_file <- "../2018-12-14.hmp_mktest/hmp_SPvsTD_map.txt"


mktest <- midas_mktest(midas_dir = midas_dir,
                       map_file = map_file,
                       genes = genes,
                       depth_thres = depth_thres)

vmwa <- read_tsv("hmp.vmwa.pvals.genes/Granulicatella_adiacens_61980_vmwa.genes.txt")
vmwa <- vmwa %>% filter(gene_id %in% genes)
vmwa %>% arrange(q.value)
