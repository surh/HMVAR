#!/usr/bin/env Rscript

# (C) Copyright 2018 Sur Herrera Paredes
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

library(HMVAR)
library(tidyverse)
library(argparser)


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
