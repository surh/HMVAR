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


args <- list(midas_dir = "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/Leptotrichia_shahii_58912/",
             genes = "/godot/users/sur/exp/fraserv/2018/today2/selected_genes.txt",
             outdir = "subset/",
             samples = "/godot/users/sur/exp/fraserv/2018/today2/selected_samples.txt",
             cds_only = TRUE)

# Read genes
if(is.na(args$genes)){
  genes <- NULL
}else{
  cat("Reading list of genes...\n")
  genes <- read_tsv(args$genes, col_names = FALSE, col_types = 'c')
  genes <- genes$X1
}

# Read samples
if(is.na(args$samples)){
  samples <- NULL
}else{
  cat("Reading list of samples...\n")
  samples <- read_tsv(args$samples, col_names = FALSE, col_types = 'c')
  samples <- samples$X1
}

# Read MIDAS data
cat("Reading MIDAS output...\n")
info <- read_tsv(paste0(args$midas_dir, "/snps_info.txt"),
                 col_types = 'ccncccnnnnnccccc',
                 na = 'NA')
depth <- read_midas_abun(paste0(args$midas_dir, "/snps_depth.txt"))
freq <- read_midas_abun(paste0(args$midas_dir, "/snps_freq.txt"))

# Process data
if (!is.null(samples)){
  cat("Selecting samples from freq and depth files...\n")
  # Clean depth and freq
  depth <- depth %>%
    select(site_id, intersect(samples, colnames(depth)) )
  freq <- freq %>%
    select(site_id, intersect(samples, colnames(freq)) )
}

# Select gene/cds data
if(args$cds_only){
  cat("Selecting CDS only...\n")
  info <- info %>% 
    filter(!is.na(gene_id))
}
if(!is.null(genes)){
  cat("Selecting chosen genes...\n")
  info <- info %>%
    filter(gene_id %in% genes)
}
freq <- freq %>% 
  filter(site_id %in% info$site_id)
depth <- depth %>% 
  filter(site_id %in% info$site_id)

# Write output
# Prepare output dir
if(!dir.exists(args$outdir)){
  cat("Preparing output directory...\n")
  dir.create(args$outdir)
}
cat("Writing info file...\n")
filename <- paste0(args$outdir, "/snps_info.txt")
write_tsv(info, filename)
cat("Writing depth file...\n")
filename <- paste0(args$outdir, "/snps_depth.txt")
write_tsv(depth, filename)
cat("Writing freq file...\n")
filename <- paste0(args$outdir, "/snps_freq.txt")
write_tsv(freq, filename)
