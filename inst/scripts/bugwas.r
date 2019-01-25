library(HMVAR)
library(tidyverse)
library(bugwas)

Sys.setenv(LD_LIBRARY_PATH="/opt/modules/pkgs/eqtlbma/git/lib/")


args <- list(midas_dir = "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/Actinomyces_odontolyticus_57475/",
             map_file = "hmp_SPvsTD_map.txt",
             outdir = "testout/")


map <- read_tsv(args$map_file, col_types = 'cc')
map <- map %>% select(sample = ID, Group = Group)
Dat <- read_midas_data(midas_dir = args$midas_dir,
                       map = map,
                       genes = NULL)

