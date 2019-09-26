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

# setwd("~/micropopgen/exp/2019/today")
# devtools::document(pkg = "~/micropopgen/src/HMVAR/")

library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Calculate Fst/Fst_pool for all sites in MIDAS output"))
  
  # Positional arguments
  p <- add_argument(p, "midas_dir",
                    help = paste("Directory with midas merge output. Must have",
                                 "snps_info.txt, snps_freq.txt and snps_depth.txt.",
                                 "Altenatively, a directory containing subdirectories",
                                 "with midas merge output."),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--type",
                    help = paste("Indicates whether midas_dir is a single species ('single') or",
                                 "a directory containing multiple species directories ('multi')"),
                    type = "character",
                    default = "single")
  p <- add_argument(p, "--map_file",
                    help = paste("Mapping file. Must be tab-delimited and have ",
                                 "ID and Group columns. Only samples in the mapping file",
                                 "will be considered."),
                    type = "character",
                    default = "map.txt")
  p <-add_argument(p, "--method",
                   help = paste("Either 'Fst' or 'Fstpool'"),
                   default = "Fstpool",
                   type = "character")
  p <- add_argument(p, "--outdir",
                    help = paste("Directory for output"),
                    type = "character",
                    default = "output/")
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  # Check type
  if(!(args$type %in% c('sinlge', 'multi')))
    stop("ERROR: --type must be one of 'single' or 'multi'")
  
  # Check method
  if(!(args$method %in% c('Fst', 'Fstpool'))){
    stop("ERROR: --method must be one of 'Fst' or 'Fstpool'")
  }else if(args$method == "Fst"){
    args$method <- "Weir-Cockerham"
  }else if(args$method == "Fstpool"){
    args$method <- 'Fstpool'
  }
  
  return(args)
}

# Read arguments
args <- process_arguments()

# Dependencies for script
library(tidyverse)
library(HMVAR)

# Read map
map <- read_tsv("midas/map.txt") %>% select(sample=ID, Group)

# Get list of midas_dirs
if(args$type == 'multi'){
  args$midas_dir <- list.dirs(args$midas_dir, recursive = FALSE)
}

if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

# Calculate for each midas dir
args$midas_dir %>%
  purrr::map(function(midas_dir, method,
                      support_thres = 1,
                      map, sorted, outdir){
    spec <- basename(midas_dir)
    cat(spec, "\n")
    
    # Read data
    Dat <- read_midas_data(midas_dir = midas_dir, map = map, cds_only = FALSE)
    Dat <- list(info = Dat$info[1:1000, ],
                depth = Dat$depth[1:1000,],
                freq = Dat$freq[1:1000,])
    
    # map <- map %>% filter(sample %in% colnames(Dat$freq))
    # map$Group %>% table
    
    # Calculate Fst
    fst <- calculate_fst(Dat = Dat,
                         method = method,
                         support_thres = support_thres,
                         map = map %>% filter(sample %in% colnames(Dat$freq)),
                         sorted = sorted,
                         verbose = TRUE)
    
    filename <- file.path(outdir, paste0(spec, ".fst.txt"))
    readr::write_tsv(fst$fst, filename)
    
  }, method = args$method,
  support_thres = 1,
  map = map,
  sorted = args$sorted,
  outdir = args$outdir)