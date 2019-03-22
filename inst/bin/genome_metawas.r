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

process_arguments <- function(){
  p <- arg_parser(paste("Perform metawas on a genome"))
  
  # Positional arguments
  p <- add_argument(p, "midas_dir",
                    help = paste("Directory with midas merge output. Must have",
                                 "snps_info.txt, snps_freq.txt and snps_depth.txt"),
                    type = "character")
  p <- add_argument(p, "focal_group",
                    help = paste("Group that is going to be compared against the rest",
                                 "in the LMM. It must correspond to one of the levels",
                                 "in the Group column of the map_file"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--map_file",
                    help = paste("Mapping file. Must be tab-delimited and have ",
                                 "ID and Group columns."),
                    type = "character",
                    default = "map.txt")
  p <- add_argument(p, "--outdir",
                    help = paste("Directory for output"),
                    type = "character",
                    default = "metawas")
  p <- add_argument(p, "--gemma",
                    help = paste("GEMMA executable."),
                    type = "character",
                    default = "gemma")
  p <- add_argument(p, "--impute",
                    help = paste("If passed, imputation will be attempted",
                                 "before LMM."),
                    flag = TRUE)
  p <- add_argument(p, "--pcs",
                    help = paste("File with community abundance  principal",
                                 "components. Must be tab-delimited, the first",
                                 "column must be named 'ID' and correspond to",
                                 "the sample IDs, and there should be one column",
                                 "per principal component."),
                    type = "character",
                    default = NULL)
  p <- add_argument(p, "--pval_thres",
                    help = paste("p-value threshold for considering significance"),
                    type = "numeric",
                    default = 1e-6)
  
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  # Check gemma
  o <- HMVAR:::run_command(paste(args$gemma, '-v'))
  if(o != 1){
    stop("ERROR: GEMMA executable not found", call. = TRUE)
  }
  
  return(args)
}


# Read arguments
# indir <- "./"
# spec <- "midas_output_small/"

# Eventually replace this with argparse
# args <- list(midas_dir = file.path(indir, spec),
#              focal_group = "Supragingival.plaque",
#              map_file = "map.txt",
#              outdir = "script_pcs_noimpute/",
#              gemma = "~/bin/gemma-0.98.1-linux-static",
#              impute = FALSE,
#              pcs = "pcs.txt",
#              pval_thres = 1e-6)
# rm(indir, spec)
args <- process_arguments()
print(args)


# The steps will be performed
# Steps
# 1. Convert midas data to BIMBAM
# 2. Impute genotypes with mice (optional)
# 3. Get kinship matrix with gemma
# 4. Run lmm
# 5. Run lmm with PCs (optional)

# Required in fraserv for the bugwas modified gemma
# Trying to move to newer GEMMA
# Sys.setenv(LD_LIBRARY_PATH="/opt/modules/pkgs/eqtlbma/git/lib/")

# Main output directory
cat("Creating output directory...\n")
dir.create(args$outdir)

# Create list for filenames
Files <- list(Dirs = list(),
              Files = list())

# Read map
cat("Processing map...\n")
map <- read_tsv(args$map_file, col_types = 'cc')
map <- map %>% select(sample = ID, Group = Group)

# Convert to bimbam
cat("Converting to BIMBAM format...\n")
Files$Dirs$bimbam_dir <- file.path(args$outdir, "bimbam")
midas_bimbam <- midas_to_bimbam(midas_dir = args$midas_dir,
                                map = map,
                                outdir = Files$Dirs$bimbam_dir,
                                focal_group = args$focal_group,
                                prefix = NULL)
Files$Files$midas_geno_file <- midas_bimbam$filenames$geno_file
Files$Files$pheno_file <- midas_bimbam$filenames$pheno_file
Files$Files$snp_file <- midas_bimbam$filenames$snp_file
rm(map)


if(args$impute){
  # Impute
  Files$Dirs$imputed_dir <- file.path(args$outdir, "imputed")
  Files$Files$imputed_geno_file <- mice_impute(geno = midas_bimbam$Dat$geno,
                                               snp = midas_bimbam$Dat$snp,
                                               outdir = Files$Dirs$imputed_dir,
                                               m = 5,
                                               verbose = FALSE,
                                               prefix = "imputed",
                                               return_table = FALSE,
                                               seed = 76543)
  Files$Files$midas_geno_file <- Files$Files$imputed_geno_file
}

# Get kinship matrix
# Works with both gemma v0.93b & v0.98.1
# I am ingoring patterns since genotypes are not fixed but frequencies instead
Files$Dirs$kinship_dir <- file.path(args$outdir, "kinship")
Files$Files$kinship_file <- gemma_kinship(geno_file = Files$Files$midas_geno_file,
                                          pheno_file = Files$Files$pheno_file,
                                          snp_file = Files$Files$snp_file,
                                          gemma = args$gemma,
                                          outdir = Files$Dirs$kinship_dir,
                                          prefix = 'kinship')


# Run standard lmm
Files$Dirs$lmm_dir <- file.path(args$outdir, "lmm")
res <- gemma_lmm(geno_file = Files$Files$midas_geno_file,
                 pheno_file = Files$Files$pheno_file,
                 snp_file = Files$Files$snp_file,
                 kinship_file = Files$Files$kinship_file,
                 cov_file = NULL,
                 gemma = args$gemma,
                 outdir = Files$Dirs$lmm_dir,
                 maf = 0,
                 prefix = "lmm")
Files$Files$lmm_log_file <- res[1]
Files$Files$lmm_assoc_file <- res[2]
rm(res)

if(!is.na(args$pcs)){
  # Run lmm with pcs
  # Format PCs as covariates
  pcs <- read_tsv(args$pcs, col_types = cols(ID = col_character(),
                                             .default = col_double()))
  pcs <- pcs %>% slice(match(midas_bimbam$Dat$pheno$id, ID))
  pcs$ID <- 1
  Files$Files$pc_covariates <- file.path(Files$Dirs$bimbam_dir, "pcs.bimbam")
  write_tsv(pcs, Files$Files$pc_covariates, col_names = FALSE)

  # Run lmmpcs
  Files$Dirs$lmmpcs_dir <- file.path(args$outdir, "lmmpcs")
  res <- gemma_lmm(geno_file = Files$Files$midas_geno_file,
                   pheno_file = Files$Files$pheno_file,
                   snp_file = Files$Files$snp_file,
                   kinship_file = Files$Files$kinship_file,
                   cov_file = Files$Files$pc_covariates,
                   gemma = args$gemma,
                   outdir = Files$Dirs$lmmpcs_dir,
                   maf = 0,
                   prefix = "lmmpcs")
  Files$Files$lmmpcs_log_file <- res[1]
  Files$Files$lmmpcs_assoc_file <- res[2]
  rm(res)
  
  # Combine results
  lmm <- read_tsv(Files$Files$lmm_assoc_file,
                  col_types = cols(chr = col_character(),
                                   rs = col_character(),
                                   ps = col_integer(),
                                   n_miss = col_integer(),
                                   allele1 = col_character(),
                                   allele0 = col_character(),
                                   af = col_number(),
                                   logl_H1 = col_number(),
                                   l_mle = col_number(),
                                   p_lrt = col_number()))
  lmmpcs <- read_tsv(Files$Files$lmmpcs_assoc_file,
                     col_types = cols(chr = col_character(),
                                      rs = col_character(),
                                      ps = col_integer(),
                                      n_miss = col_integer(),
                                      allele1 = col_character(),
                                      allele0 = col_character(),
                                      af = col_number(),
                                      logl_H1 = col_number(),
                                      l_mle = col_number(),
                                      p_lrt = col_number())) %>%
    select(-n_miss, -allele1, -allele0, - af)
  lmm <- lmm %>% full_join(lmmpcs, by = c("chr", "rs", "ps"),
                           suffix = c(".lmm", ".lmmpcs"))
  
  # Select interpretation
  res <- rep('none', nrow(lmm))
  res[ lmm$p_lrt.lmm < args$pval_thres & lmm$p_lrt.lmmpcs >= args$pval_thres ] <- "int"
  res[ lmm$p_lrt.lmm >= args$pval_thres & lmm$p_lrt.lmmpcs < args$pval_thres ] <- "env"
  res[ lmm$p_lrt.lmm < args$pval_thres & lmm$p_lrt.lmmpcs < args$pval_thres ] <- "both"
  lmm <- lmm %>% bind_cols(type = res)
  # lmm$type %>% table
  
  # Write results
  Files$Files$results <- file.path(args$outdir, "lmm.results.txt")
  write_tsv(lmm, Files$Files$results)
}