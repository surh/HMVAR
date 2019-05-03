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

library(tidyverse)
library(argparser)
library(HMVAR)

process_arguments <- function(){
  p <- arg_parser(paste("Calculate DoS statistics on a single genome or a
                        set of pre-processed genomes."))
  
  # Positional arguments
  p <- add_argument(p, "input",
                    help = paste("Input can be either a single file or a directory.",
                                 "If a single file is passed, it hould be a tab-delimited",
                                 "file where each row corresponds to a gene, and it must",
                                 "include columns Dn, Ds, Pn, and Ps corresponding to the",
                                 "values in the McDonald-Kreitman contingency table.\n",
                                 "If a directory is passed, it should either contain a",
                                 "set of tab-delimited files as above, or a set of",
                                 "sub-directories, each one of which is the ouptut from",
                                 "<midas_merge.py snps>. The option --dir_type allows",
                                 "you to specify which type of directory it is, and",
                                 "the option --suffix lets you specify the suffix of",
                                 "the tab-delimited files."),
                    type = "character")
  p <- add_argument(p, "outdir",
                    help = paste("This should be the outpu directory for the output"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--dir_type",
                    help = paste("Either 'tabs' or 'midas' to indicate whether",
                                 "tab-delimited files or midas_merge.py directories",
                                 "within the input directory",
                                 "are to be processed. If <input> is a file, then",
                                 "this option will be ignored."),
                    type = "character",
                    default = "files")
  p <- add_argument(p, "--map_file",
                    help = paste("Mapping file that associates samples",
                                 "with groups. It must be a tab-delimited",
                                 "file with headers 'ID' and 'Group'. Only",
                                 "required if --dir_type=midas"),
                    type = "character")
  p <- add_argument(p, "--genes",
                    help = paste0("File with list of gene IDs to test. ",
                                  "It must contain one gene id per line, ",
                                  "without header, and the IDs must ",
                                  "correspond to the gene_id column from ",
                                  "the snps_info.txt file if --dir_type=midas", 
                                  "or to the same column in the tab-delimited files",
                                  "if --dir_type=tabs or if <input> was a file.",
                                  "If NULL, all ",
                                  "genes will be tested."),
                    default = NULL,
                    type = "character")
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
  p <- add_argument(p, "--suffix",
                    help = paste("The suffix of the tab-delimited files. It will be used",
                                 "to identify which files to process, and also, to",
                                 "determine the output files. Basically if an input",
                                 "file is of the form /path/<prefix><suffix>",
                                 "the ouptut will be of the form <outdir>/<prefix>*"),
                    type = "character",
                    default = ".txt")
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(args$freq_thres < 0 || args$freq_thres > 1){
    stop("ERROR: freq_thres must be a value in [0,1].")
  }
  if(args$depth_thres < 0){
    stop("ERROR: depth_thres must be a non-negative integer.")
  }
  args$suffix <- paste0(args$suffix, "$")
  
  # Read genes
  if(is.na(args$genes)){
    args$genes <- NULL
  }else{
    cat("Reading list of genes...\n")
    args$genes <- read_tsv(args$genes, col_names = FALSE, col_types = 'c')
    args$genes <- genes$X1
  }
  
  return(args)
}

write_outputs <- function(dat, infile, suffix, outdir, pval_col){
  # infile <- args$input
  # outdir <- args$outdir
  
  prefix <- str_replace(string = basename(infile), pattern = suffix, replacement = "")
  
  # Write table
  filename <- paste0(c(prefix, "DoS.table.txt"), collapse = ".")
  filename <- file.path(outdir, filename)
  write_tsv(dos, filename)
  
  # Produce plots
  if(sum(!is.na(dat$DoS)) > 1){
    p1 <- ggplot(dat, aes(x = DoS)) +
      geom_histogram(bins = 15) +
      ggtitle(label = prefix) +
      AMOR::theme_blackbox()
    filename <- paste0(c(prefix, "DoS.histogram.png"), collapse = ".")
    filename <- file.path(outdir, filename)
    ggsave(filename, p1, width = 6, height = 5, dpi = 200)
    
    p1 <- ggplot(dat, aes_string(x = pval_col)) +
      geom_histogram(bins = 20) +
      ggtitle(label = prefix) +
      AMOR::theme_blackbox()
    filename <- paste0(c(prefix, "DoS.pval.histogram.png"), collapse = ".")
    filename <- file.path(outdir, filename)
    ggsave(filename, p1, width = 6, height = 5, dpi = 200)
    
    filename <- paste0(c(prefix, "DoS.pval.qqplot.png"), collapse = ".")
    filename <- file.path(outdir, filename)
    png(filename = filename, width = 6, height = 5, units = "in", res = 200)
    HMVAR::pval_qqplot(dat[ , pval_col, drop = TRUE ])
    dev.off()
  }
}


args <- process_arguments()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

if(dir.exists(args$input)){
  # Input is dir
  cat("Input is a dir\n")
  if(args$dir_type == 'tabs'){
    cat("Inputs are tab files\n")
    inputs <- list.files(args$input)
    inputs <- str_subset(string = inputs, pattern = args$suffix)
    cat("Found ", length(inputs), " files\n")
  }else if(args$dir_type == "midas"){
    cat("Inputs are dirs\n")
    if(!file.exists(args$map_file)){
      stop("ERROR: --map_file must be provided in <input> is dir and --dir_type=midas.")
    }
    inputs <- list.dirs(args$input, recursive = FALSE)
    # inputs <- str_subset(string = inputs, pattern = args$suffix)
    cat("Found ", length(inputs), " dirs\n")
    pval_col <- 'p.value'
  }else{
    stop("ERROR: --dir_type must be tabs or midas.")
  }
  
  Res <- NULL
  for(input in inputs){
    cat("Processing ", input, "...\n")
    if(args$dir_type == 'tabs'){
      # Each input is a file
      # Input is file
      dat <- read_tsv(file.path(args$input, input),
                      col_types = cols(.default = col_character(),
                                       Dn = col_integer(),
                                       Ds = col_integer(),
                                       Pn = col_integer(),
                                       Ps = col_integer()))
      
      if(!is.null(args$genes)){
        dat <- dat %>%
          filter(gene_id %in% args$genes)
      }
      
      if('p.value' %in% colnames(dat)){
        pval_col <- 'DoS.p.value' 
      }else{
        pval_col <- 'p.value'
      }
      
      dos <- calculate_dos(dat = dat, test = TRUE, clean = FALSE)
      write_outputs(dat = dos, infile = input, suffix = args$suffix, outdir = args$outdir, pval_col = pval_col)
    }else if(args$dir_type == 'midas'){
      # Each input is a midas merge dir
      dos <- midas_dos(midas_dir = input,
                       map = args$map_file,
                       genes = args$genes,
                       depth_thres = args$depth_thres,
                       freq_thres = args$freq_thres,
                       focal_group = args$focal_group)
      
      write_outputs(dat = dos, infile = input, suffix = "$", outdir = args$outdir, pval_col = 'p.value')
    }else{
      stop("ERROR: --dir_type must be tabs or midas.")
    }
    dos$input <- input
    Res <- Res %>% bind_rows(dos)
  }
  # Write combined
  write_outputs(dat = Res, infile = "summary", suffix = '$', outdir = args$outdir, pval_col = pval_col)
}else if(file.exists(args$input)){
    # Input is file
    dat <- read_tsv(args$input,
                    col_types = cols(.default = col_character(),
                                     Dn = col_integer(),
                                     Ds = col_integer(),
                                     Pn = col_integer(),
                                     Ps = col_integer()))
    
    if(!is.null(args$genes)){
      dat <- dat %>%
        filter(gene_id %in% args$genes)
    }
    
    if('p.value' %in% colnames(dat)){
      pval_col <- 'DoS.p.value' 
    }else{
      pval_col <- 'p.value'
    }
     
    dos <- calculate_dos(dat = dat, test = TRUE, clean = FALSE)
    write_outputs(dat = dos, infile = args$input, suffix = args$suffix, outdir = args$outdir, pval_col = pval_col)
}else{
  stop("ERROR: input doesn't exist")
}