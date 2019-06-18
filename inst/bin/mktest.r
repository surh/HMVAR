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

################### Functions ###################
process_arguments <- function(){
  p <- arg_parser(paste0("Performs McDonald-Kreitman test on a subset ",
                         "of genes from the output of MIDAS."))
  
  # Positional arguments
  p <- add_argument(p, "midas_dir",
                    help = paste0("Directory containing ",
                                  "SNPs merged by MIDAS for a ",
                                  "single genome. Must contain files: ",
                                  "'snps_info.txt', 'snps_depth.txt', ",
                                  "and 'snps_freq.txt'."),
                    type = "character")
  p <- add_argument(p, "map_file",
                    help = paste0("Mapping file that associates samples ",
                                  "with groups. It must be a tab-delimited ",
                                  "file with headers 'ID' and 'Group'"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--outfile",
                    help = paste0("Name of file to write output"),
                    default = "mktest.txt",
                    type = "character")
  p <- add_argument(p, "--genes",
                    help = paste0("File with list of gene IDs to test. ",
                                  "It must contain one gene id per line, ",
                                  "without header, and the IDs must ",
                                  "correspond to the gene_id column from ",
                                  "the snps_info.txt file. If NULL, all ",
                                  "genes will be tested."),
                    default = NULL,
                    type = "character")
  p <- add_argument(p, '--merge',
                    help = paste0("If passed, all the McDonald-Kreitman values ",
                                  "will be combined into one. Useful to test if ",
                                  "a set of genes are under adaptive evolution."),
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
                                  " distance from 0 or 1, for a site to be assigned ",
                                  "to the major or minor allele respectively. ",
                                  "It must be a value in [0,1]."),
                    default = 0.5,
                    type = "double")
  p <- add_argument(p, "--focal_group",
                    help = paste("Group to use as a focal group to compare against all other.",
                                 "If passed, the script will convert all values in the Group",
                                 "column of map into 'non.<focal_group>', effectiveley",
                                 "ensuring a dichotomous grouping. If not passed, the script",
                                 "will expect only two levels in the Group column of map and it",
                                 "will fail if more than one group is found."),
                    default = NULL,
                    type = "character")
  p <- add_argument(p, '--test',
                    help = paste0("If passed, the Fisher's exact test will be performed ",
                                  "on the McDonald-Kreitman contingency table."),
                    flag = TRUE)
  
  # Read arguments
  args <- parse_args(p)
  
  # Process arguments
  if(args$freq_thres < 0 || args$freq_thres > 1){
    stop("ERROR: freq_thres must be a value in [0,1].")
  }
  if(args$depth_thres < 0){
    stop("ERROR: depth_thres must be a non-negative integer.")
  }
  
  return(args)
}
#################################################

# Set up
args <- process_arguments()
if(is.na(args$genes)){
  genes <- NULL
}else{
  genes <- read_tsv(args$genes, col_names = FALSE, col_types = 'c')
  genes <- genes$X1
}

# MK contingency tables
cat("Calculating contingency tables...\n")
mktest <- midas_mktest(midas_dir = args$midas_dir,
                       map_file = args$map_file,
                       genes = genes,
                       depth_thres = args$depth_thres,
                       freq_thres = args$freq_thres,
                       focal_group = args$focal_group)

# Merge
if(args$merge){
  cat("Merging...\n")
  mktest <- mktest %>% summarise(gene_id = 'merged.set',
                                 Dn = sum(Dn),
                                 Ds = sum(Ds),
                                 Pn = sum(Pn),
                                 Ps = sum(Ps))
}

# Test
if(args$test){
  cat("Testing...\n")
  mktest <- mktest %>%
    pmap_dfr(function(gene_id, Dn, Ds, Pn, Ps){
      mat <- matrix(c(Dn,Ds,Pn,Ps), ncol = 2)
      test <- fisher.test(mat)
      return(tibble(gene_id = gene_id,
                    Dn = Dn,
                    Ds = Ds,
                    Pn = Pn,
                    Ps = Ps,
                    OR = test$estimate,
                    p.value = test$p.value))
    }) %>%
    bind_rows(tibble(gene_id = character(0),
                     Dn = integer(0),
                     Ds = integer(0),
                     Pn = integer(0),
                     Ps = integer(0),
                     OR = numeric(0),
                     p.value = numeric(0)))
}

# Write output
cat("Writing output...\n")
mktest %>% print(n = 20)
write_tsv(mktest, args$outfile)
