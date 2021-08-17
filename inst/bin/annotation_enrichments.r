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

gw_test_enrichments <- function(input, annotations, closest,
                                dist_thres = 500, score_column = 'p.value',
                                gene_score = 'min', alternative = 'less',
                                annot_column = 'GO_terms', method = 'gsea',
                                min_size = 3){
  # Read test
  cat("\tReading genome-wide test...\n")
  col_specs <- rlang::list2(chr = col_character(),
                            ps = col_integer(),
                            rs = col_character(),
                            gene_id = col_character(),
                            !!score_column := col_double())
  gw.test <- read_tsv(input,
                      col_types = do.call(cols, col_specs))
  
  # Reading closest
  if(!is.na(closest)){
    
    if(!all(c('chr', 'ps') %in% colnames(gw.test))){
      stop("ERROR: if closest is provided, input must have 'chr' and 'ps' columns.", call. = TRUE)
    }
    if(!(gene_score %in% c('min', 'max'))){
      stop("ERROR: gene_score must be 'min' or 'max'", call. = TRUE)
    }
    
    cat("\tReading closest...\n")
    closest <- read_tsv(closest,
                        col_names = c("chr", "ps", "ps2", "chr2",
                                      "start", "end", "gene_id", "dist"),
                        col_types = cols(.default = col_number(),
                                         chr = col_character(),
                                         chr2 = col_character(),
                                         gene_id = col_character())) %>%
      select(-ps2, -chr2) %>%
      filter(gene_id != '.')
    # closest
    
    # Match genes and tests
    # Remove too distant features
    closest <- closest %>%
      filter(abs(dist) <= dist_thres) %>%
      select(gene_id, chr, ps)
    # Match features to genes
    gw.test <- gw.test %>%
      left_join(closest, by = c('chr', 'ps'))
    # Get lowest score per gene
    
    score_sel_fun <- match.fun(gene_score)
    gw.test <- gw.test %>%
      split(.$gene_id) %>%
      map_dfr(~ tibble(!!score_column := score_sel_fun(.x[,score_column, drop = TRUE])),
              score_column = score_column, .id = 'gene_id')
    
    cat("\tNumber of genes associated: ", nrow(gw.test), "\n")
  }else if(!('gene_id' %in% colnames(gw.test))){
    stop("ERROR: if closest not provided, input must have 'gene_id' column.", call. = TRUE)
  }
  
  # Read annotations
  cat("\tReading annotations...\n")
  annots <- read_eggnog(annotations) %>%
    select(gene_id = query_name, everything()) %>%
    select(gene_id, terms = annot_column)
  #### CUSTOM FOR TESTS ####
  annots <- annots %>%
    mutate(gene_id = str_replace(string = gene_id,
                                 pattern = "\\([+-]\\)_[0-9]",
                                 replacement = ""))
  ##########
  # Match genes with annotations
  gw.test <- annots %>%
    right_join(gw.test, by = "gene_id") %>%
    select(gene_id, terms, score = score_column) %>%
    filter(!is.na(score))
  cat("\tGenes with annotation: ", sum(!is.na(gw.test$terms)), "\n")
  
  # test
  if(method == 'gsea'){
    res <- terms_enrichment(dat = gw.test, method = method, test = 'wilcoxon',
                            alternative = alternative, min_size = min_size)
  }else if(method == 'sign_test'){
    res <- terms_enrichment(dat = gw.test, method = 'sign_test',
                            alternative = alternative, min_size = min_size)
  }else if(method == 'test_go'){
    res <- terms_enrichment(dat = gw.test, method = method,
                            description = '', algorithm = 'weight01',
                            statatistic = 'ks', node_size = min_size)
  }else{
    stop("ERROR: method must be 'gsea' or 'test_go'", call. = TRUE)
  }
  
  return(list(data = gw.test, res = res))
}

process_arguments <- function(){
  p <- arg_parser(paste("Calculate DoS statistics on a single genome or a
                        set of pre-processed genomes."))
  
  # Positional arguments
  p <- add_argument(p, "input",
                    help = paste("Input can be either a single file or a directory.",
                                 "If a single file is passed, it hould be a tab-delimited",
                                 "file where each row corresponds to a tested feature.",
                                 "If a directory is passed, it should contain a",
                                 "set of tab-delimited files as above. In the latter case",
                                 "the option --suffix lets you specify the suffix of",
                                 "the tab-delimited files. In any case, the tab-delimited",
                                 "files must contain columns named 'chr' and 'ps' which",
                                 "indicate the feature's position."),
                    type = "character")
  p <- add_argument(p, "outdir",
                    help = paste("This should be the directory for the output"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--suffix",
                    help = paste("The suffix of the tab-delimited files. It will be used",
                                 "to identify which files to process, and also to",
                                 "determine the output files. Basically if an input",
                                 "file is of the form /path/<prefix><suffix>",
                                 "the ouptut will be of the form <outdir>/<prefix>*"),
                    type = "character",
                    default = ".txt")
  p <- add_argument(p, "--closest", help = paste("This can be either a single file or a",
                                                 "directory. It must match the same type",
                                                 "as <input>. If a single file it should",
                                                 "be a tab-delimited file produced by bedtools'",
                                                 "closestBed.",
                                                 "Followng the BEDtools convention, a dot (.) in the",
                                                 "seventh column (gene_id) is interpreted as no gene",
                                                 "being close to the given SNP and all rows with a dot",
                                                 "on the seventh column are filtered out.",
                                                 "If a directory is passed, then",
                                                 "for each input file in the <input> directory,",
                                                 "a file that has the same prefix will be expected.",
                                                 "E.g. if an inout file is of the form",
                                                 "<input>/<prefix><suffix>, then there should be a",
                                                 "file of the form <closest>/<prefix>*. If nothing",
                                                 "is passed, then each tab-delimited file must have",
                                                 "a 'gene_id' column."),
                    type = "character",
                    default = NA)
  p <- add_argument(p, "--annotations", help = paste("This can be either a single file or a",
                                                     "directory. It must match the same type",
                                                     "as <input>. If a single file it should",
                                                     "be a tab-delimited file produced by eggNOG",
                                                     "mapper. If a directory is passed, then",
                                                     "for each input file in the <input> directory,",
                                                     "a file that has the same prefix will be expected.",
                                                     "E.g. if an inout file is of the form",
                                                     "<input>/<prefix><suffix>, then there should be a",
                                                     "file of the form <annotations>/<prefix>*."),
                    type = "character",
                    default = NA)
  p <- add_argument(p, "--dist_thres",
                    help = paste("If --closest is passed, this is the maximum distance threshold",
                                 "to associate a gene to a feature."),
                    type = "numeric",
                    default = 500)
  p <- add_argument(p, "--min_size", help = paste("This is the minimum number of genes with",
                                                     "a given annotation, for that annotation to",
                                                     "be tested"),
                    type = "numeric",
                    default = 3)
  p <- add_argument(p, "--annot_column",
                    help = paste("Name of the column with the annotation to test in --annotations files"),
                    type = "character",
                    default = 'GO_terms')
  p <- add_argument(p, "--score_column",
                    help = paste("Name of the column with the score to test in <input> files."),
                    type = "character",
                    default = 'p.value')
  p <- add_argument(p, "--gene_score",
                    help = paste("If --closest is passed, this indicates how to select the per-gene score.",
                                 "Either the maximum (max) or minimum (min) score associated with that gene."),
                    default = "min",
                    type = "character")
  p <- add_argument(p, "--alternative",
                    help = paste("Alternative hypothesis to test. See documentation in",
                                 "test_go and gsea functions from HMVAR."),
                    type = "character",
                    default = 'less')
  
  p <- add_argument(p, "--method",
                    help = paste("Which HMVAR function to use. Either 'test_go', 'sign_test', or 'gsea'"),
                    type = "character",
                    default = 'gsea')

  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  args$suffix <- paste0(args$suffix, "$")
  if(!(args$method %in% c('test_go', 'gsea', 'sign_test'))){
    stop("ERROR: --method must be 'test_go', 'sign_test' or 'gsea'", call. = TRUE)
  }
  if(!(args$gene_score %in% c('min', 'max'))){
    stop("ERROR: --gene_score must be 'min' or 'max'")
  }
  
  return(args)
}

args <- process_arguments()
# args <- list(input = "../lmmres/",
#              outdir = "test/",
#              suffix = "_lmm.results.txt",
#              closest = "../closest/",
#              annotations = "../annots/",
#              dist_thres = 500,
#              min_size = 3,
#              annot_column = "GO_terms",
#              score_column = "p_lrt.lmmpcs",
#              gene_score = 'min',
#              alternative = 'less',
#              method = 'gsea')

if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

if(dir.exists(args$input)){
  cat("Input is a dir\n")
  if(!dir.exists(args$closest) || !dir.exists(args$annotations)){
    stop("ERROR: If <input> is a dir, --closest and --annotations must be dirs too.", call. = TRUE)
  }
  
  # Find input files
  # args$input <- "~/micropopgen/exp/2019/2019-03-29.hmp_metawas_data/Supragingival.plaque/metawas/lmmpcs/"
  inputs <- list.files(args$input)
  inputs <- str_subset(string = inputs, pattern = args$suffix)
  cat("Found ", length(inputs), " files\n")
  prefixes <- str_replace(inputs, args$suffix, "")
  prefixes <- paste0("^", prefixes)
  
  # Get closest and annotation files
  closest_files <- list.files(args$closest)
  cat("Total ", length(closest_files), " closest files\n")
  # cat(head(closest_files), "\n")
  closest_files <- prefixes %>% map(~str_subset(closest_files, .)) %>% unlist %>% unique
  cat("Found ", length(closest_files), " closest files\n")
  annot_files <- list.files(args$annotations)
  cat("Total ", length(annot_files), " annotation files\n")
  # cat(head(annot_files), "\n")
  annot_files <- annot_files <- prefixes %>% map(~str_subset(annot_files, .)) %>% unlist %>% unique
  cat("Found ", length(annot_files), " annotation files\n")
  if(length(closest_files) < length(inputs) || length(annot_files) < length(inputs)){
    stop("ERROR: number of closest ot annotation files is less than number of input files.")
  }
  
  Dat <- tibble()
  for(i in 1:length(inputs)){
    input <- inputs[i]
    prefix <- prefixes[i]
    closest <- str_subset(closest_files, pattern = prefix)
    annotations <- str_subset(annot_files, pattern = prefix)
    
    if(length(closest) != 1 || length(annotations) != 1){
      cat("Input:", input, "\n")
      stop("ERROR: missing or extra closest or annotation file")
    }
    cat(input, "\n")
    Res <- gw_test_enrichments(input = file.path(args$input, input),
                               annotations = file.path(args$annotations, annotations),
                               closest = file.path(args$closest, closest),
                               dist_thres = args$dist_thres,
                               score_column = args$score_column,
                               gene_score = args$gene_score,
                               alternative = args$alternative,
                               annot_column = args$annot_column,
                               method = args$method, min_size = args$min_size)
    
    filename <- paste0(str_replace(prefix, '^\\^', ''), ".enrichments.txt")
    filename <- file.path(args$outdir, filename)
    write_tsv(Res$res, path = filename)
    
    Dat <- Dat %>% bind_rows(Res$data)
  }
  
  
  # Test overall
  if(args$method == 'gsea'){
    res <- terms_enrichment(dat = Dat, method = args$method, test = 'wilcoxon',
                            alternative = args$alternative, min_size = args$min_size)
  }else if(args$method == 'sign_test'){
    res <- terms_enrichment(dat = Dat, method = args$method,
                            alternative = args$alternative,
                            min_size = args$min_size)
  }else if(args$method == 'test_go'){
    res <- terms_enrichment(dat = Dat, method = args$method,
                            description = '', algorithm = 'weight01',
                            statatistic = 'ks', node_size = args$min_size)
  }else{
    stop("ERROR: --method must be 'gsea', 'sign_test' or 'test_go'", call. = TRUE)
  }
  filename <- file.path(args$outdir, 'overall.enrichments.txt')
  write_tsv(res, path = filename)
  
  Dat <- Dat %>% bind_rows(Res$data)
  
}else if(file.exists(args$input)){
  cat("Input is a file\n")
  Res <- gw_test_enrichments(input = args$input, annotations = args$annotations,
                             closest = args$closest, dist_thres = args$dist_thres,
                             score_column = args$score_column,
                             gene_score = args$gene_score,
                             alternative = args$alternative,
                             annot_column = args$annot_column,
                             method = args$method, min_size = args$min_size)
  
  cat("============data=================\n")
  print(Res$data)
  cat("=============res=================\n")
  print(Res$res)
  
  cat("Writing output...\n")
  prefix <- str_replace(basename(args$input), args$suffix, "")
  filename <- paste0(prefix, ".enrichments.txt")
  filename <- file.path(args$outdir, filename)
  write_tsv(Res$res, path = filename)
  
}else{
  stop("ERROR: input doesn't exist")
}

warnings()