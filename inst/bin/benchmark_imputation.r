#!/usr/bin/env Rscript

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

library(HMVAR)
library(tidyverse)
library(argparser)

################### Functions ###################
process_arguments <- function(){
  p <- arg_parser(paste0("Hides a proportion of the data and compares ",
                         "the results of mice imputation with original ",
                         "observations"))
  
  # Positional arguments
  p <- add_argument(p, "midas_dir",
                    help = paste0("Directory containing ",
                                  "SNPs merged by MIDAS for a ",
                                  "single genome. Must contain files: ",
                                  "'snps_info.txt', 'snps_depth.txt', ",
                                  "and 'snps_freq.txt'."),
                    type = "character")
  
  
  # Optional arguments
  p <- add_argument(p, "--map_file",
                    help = paste("Mapping file for samples. Must contain,",
                                 "columns ID and Group."),
                    type = "character",
                    default = "map.txt")
  p <- add_argument(p, "--outdir",
                    help = paste("Directory to save output."),
                    default = "benchmark_imputation/",
                    type = "character")
  p <- add_argument(p, "--hidden_proportion",
                    help = paste("Proportion of the data to hide"),
                    type = "numeric",
                    default = 0.1)
  p <- add_argument(p, "--m",
                    help = paste("Number of imputations to perform."),
                    default = 5,
                    type = "integer")
  p <- add_argument(p, "--seed",
                    help = paste("Seed for selecting data to hide, ",
                                 "and imputation."),
                    default = 76543,
                    type = "integer")
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(args$hidden_proportion < 0 || args$hidden_proportion > 1){
    stop("ERROR: hidden_proportion must be a value in [0,1].")
  }
  
  return(args)
}
######################3

# indir <- commandArgs(trailingOnly = TRUE)[1]
# spec <- commandArgs(trailingOnly = TRUE)[2]
# indir <- "./"
# spec <- "midas_output_small/"
# 
# # Eventually replace this with argparse
# args <- list(midas_dir = file.path(indir, spec),
#              map_file = "map.txt",
#              outdir = "benchmark_imputation",
#              hidden_proportion = 0.1,
#              seed = 76543,
#              m = 5)
args <- process_arguments()

# Main output directory
cat("Creating output directory...\n")
dir.create(args$outdir)

# Create list for filenames
Files <- list(Dirs = list(),
              Files = list())

# Read map
cat("Reading map...\n")
map <- read_tsv(args$map_file, col_types = 'cc')
map <- map %>% select(sample = ID, Group = Group)

# Convert to bimbam
cat("Converting MIDAS data to BIMBAM...\n")
Files$Dirs$bimbam_dir <- file.path(args$outdir, "original")
midas_bimbam <- midas_to_bimbam(midas_dir = args$midas_dir,
                                map = map,
                                outdir = Files$Dirs$bimbam_dir,
                                focal_group = "",
                                prefix = NULL)
Files$Files$midas_geno_file <- midas_bimbam$filenames$geno_file
Files$Files$pheno_file <- midas_bimbam$filenames$pheno_file
Files$Files$snp_file <- midas_bimbam$filenames$snp_file
rm(map)

### Benchmark imputation
cat("Benchmarking imputation...\n")
Files$Dirs$data_hidden_dir <- file.path(args$outdir, "data_hidden_geno_files/")
Res <- benchmark_imputation(geno = midas_bimbam$Dat$geno,
                            snp = midas_bimbam$Dat$snp,
                            outdir = Files$Dirs$data_hidden_dir,
                            p = args$hidden_proportion,
                            m = args$m,
                            verbose = FALSE,
                            seed = args$seed)
Files$Files$imputed_geno_file <- Res$imputed_geno_file

# Write results
Files$Files$summary_stats_file <- file.path(args$outdir, "summary_stats.txt")
write_tsv(tibble(path = args$midas_dir,
                 r = Res$r,
                 p.imputed = Res$p.imputed,
                 t.user = Res$t[1],
                 t.system = Res$t[2],
                 t.elapsed = Res$t[3]),
          Files$Files$summary_stats_file)
Res$res$path <- args$midas_dir
Files$Files$imputation_results_file <- file.path(args$outdir, "imputation_results.txt")
write_tsv(Res$res, Files$Files$imputation_results_file)
