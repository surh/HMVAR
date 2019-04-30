#!/usr/bin/env Rscript

# (C) Copyright 2018-2019 Sur Herrera Paredes
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

################ FUNCTIONS ################
process_arguments <- function(){
  p <- arg_parser(paste0("Calculates and plots the number of variable sites ",
                         "per sample within a genome (or a subset of genes in that ",
                         "genome."))
  
  p <- add_argument(p, "midas_dir",
                    help = paste0("Directory with output from merging MIDAS SNPs. ",
                                  "Must contain 'snps_info.txt', 'snps_depth.txt' ",
                                  "and 'snps_freq.txt'."),
                    type = "character")
  p <- add_argument(p, "map_file",
                    help = paste0("Mapping file associating samples with groups. ",
                                  "Must have 'Group' and 'ID' columns."),
                    type = "character")
  p <- add_argument(p, "--genes",
                    help = paste0("File with subset of genes to analyze. ",
                                  "It must be one gene ID per line without headers."),
                    type = "character",
                    default = NULL)
  p <- add_argument(p, "--cds_only",
                    help = paste0("Flag indicating whether only sites at CDSs should ",
                                  "be considered."),
                    flag = TRUE)
  p <- add_argument(p, "--depth_thres",
                    help = paste0("Minumum number of reads at a given site ",
                                  "at a given sample, for that site in that ",
                                  "sample to be included in calculations"),
                    default = 1,
                    type = "integer")
  p <- add_argument(p, "--freq_thres",
                    help = paste0("This value is the threshold that will be ",
                                  "used to determine which allele is present in ",
                                  "a sample. The value from ",
                                  "min(freq_thres, 1 - freq_thres) represents the ",
                                  "distance from 0 or 1, for a site to be assigned ",
                                  "to the major or minor allele respectively. ",
                                  "It must be a value in [0,1]."),
                    default = 0.5,
                    type = "double")
  p <- add_argument(p, "--outdir",
                    help = paste0("Directory to write (and create if needed) to write ",
                                  "the output."),
                    default = "results/",
                    type = "character")
  p <- add_argument(p, "--prefix",
                    help = paste0("Prefix for output filenames."),
                    default = "out",
                    type = "character")

  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(args$freq_thres < 0 || args$freq_thres > 1)
    stop("ERROR: freq_thres must be a value in [0,1].")
  if(args$depth_thres < 0)
    stop("ERROR: depth_thres must be a non-negative integer.")
  if(!dir.exists(args$midas_dir))
    stop("ERROR: midas_dir does not exist.")
  if(!file.exists(args$map_file))
    stop("ERROR: map_file does not exist.")
  
  # Read genes
  if(is.na(args$genes)){
    args$genes <- NULL
  }else{
    cat("Reading list of genes...\n")
    args$genes <- read_tsv(args$genes, col_names = FALSE, col_types = 'c')
    args$genes <- genes$X1
  }
  
  if(is.na(args$prefix)){
    args$prefix <- NULL
  }
  
  return(args)
}

# Parameters
args <- list(midas_dir = "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/genomes/Granulicatella_adiacens_61980/",
             map_file = "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/hmp_SPvsTD_map.txt",
             genes = NA,
             depth_thres = 1,
             freq_thres = 0.5,
             cds_only = FALSE,
             prefix = "Granulicatella_adiacens_61980",
             outdir = "results/")
args <- process_arguments()

# Read map
cat("Reading map...\n")
map <- read_tsv(args$map_file,
                col_types = cols(ID = col_character(),
                                 Group = col_character())) %>%
  select(sample = ID,
         everything())

# Read and process midas sites
midas_data <- read_midas_data(midas_dir = args$midas_dir,
                              map = map, genes = args$genes,
                              cds_only = FALSE)

# Varsites pipeline
Res <- varsites_pipeline(freq = midas_data$freq,
                         depth = midas_data$depth,
                         info = midas_data$info,
                         map = map,
                         depth_thres = args$depth_thres,
                         freq_thres = args$freq_thres,
                         plot = TRUE)

# Write output
# Prepare output dir
if(!dir.exists(args$outdir)){
  cat("Preparing output directory...\n")
  dir.create(args$outdir)
}

cat("Writing output...\n")
filename <- paste0(c(args$prefix, "varsites.txt"), collapse = ".")
filename <- file.path(args$outdir, filename)
write_tsv(Res$varsites, filename)

filename <- paste0(c(args$prefix, "varsites.pos.txt"), collapse = ".")
filename <- file.path(args$outdir, filename)
write_tsv(Res$varsites.pos, filename)

filename <- paste0(c(args$prefix, "dnds.png"), collapse = ".")
filename <- file.path(args$outdir, filename)
ggsave(filename, Res$dnds.plot, width = 6 , height = 4, units = 'in', dpi = 200)

filename <- paste0(c(args$prefix, "variability.png"), collapse = ".")
filename <- file.path(args$outdir, filename)
ggsave(filename, Res$variability.plot, width = 6 , height = 4, units = 'in', dpi = 200)

filename <- paste0(c(args$prefix, "pos.png"), collapse = ".")
filename <- file.path(args$outdir, filename)
ggsave(filename, Res$pos.plot, width = 12 , height = 4, units = 'in', dpi = 200)

