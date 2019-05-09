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
  col_specs <- rlang::list2(chr = col_character(),
                            ps = col_integer(),
                            rs = col_character(),
                            gene_id = col_character(),
                            !!score_column := col_double())
  gw.test <- read_tsv(args$input,
                      col_types = do.call(cols, col_specs))
  
  # Reading closest
  if(!is.na(closest)){
    
    if(!all(c('chr', 'ps') %in% colnames(gw.test))){
      stop("ERROR: if closest is provided, input must have 'chr' and 'ps' columns.", call. = TRUE)
    }
    if(!(gene_score %in% c('min', 'max'))){
      stop("ERROR: gene_score must be 'min' or 'max'", call. = TRUE)
    }
    
    closest <- read_tsv(closest,
                        col_names = c("chr", "ps", "ps2", "chr2",
                                      "start", "end", "gene_id", "dist"),
                        col_types = cols(.default = col_number(),
                                         chr = col_character(),
                                         chr2 = col_character(),
                                         gene_id = col_character())) %>%
      select(-ps2, -chr2)
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
  }else if(!('gene_id' %in% colnames(gw.test))){
    stop("ERROR: if closest not provided, input must have 'gene_id' column.", call. = TRUE)
  }
  
  # Read annotations
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
    select(gene_id, terms, score = score_column)
  
  # Test
  if(method == 'gsea'){
    res <- gsea(dat = gw.test, test = 'wilcoxon', alternative = alternative, min_size = min_size)
    
    # If terms are GO match with annotations
    if(all(str_detect(res$term, "^GO:[0-9]{7}"))){
      res <- res %>%
        pmap_dfr(function(term, size, statistic, p.value){
          t <- GO.db::GOTERM[[term]]
          if(!is.null(t)){
            ontology <- t@Ontology
            annotation <- t@Term
          }else{
            ontology <- NA
            annotation <- NA
          }
          tibble(term = term,
                 size = size,
                 statistic = statistic,
                 p.value = p.value,
                 ontology = ontology,
                 annotation = annotation)})
    }
    
  }else if(method == 'test_go'){
    genes <- gw.test$score
    names(genes) <- gw.test$gene_id
    bp.res <- test_go(genes = genes,
                      annots = annots,
                      ontology = 'BP',
                      algorithm = 'weight01',
                      statistic = 'ks',
                      node_size = min_size,
                      score_threshold = 1e-5)
    cc.res <- test_go(genes = genes,
                      annots = annots,
                      ontology = 'CC',
                      algorithm = 'weight01',
                      statistic = 'ks',
                      node_size = min_size,
                      score_threshold = 1e-5)
    mf.res <- test_go(genes = genes,
                      annots = annots,
                      ontology = 'MF',
                      algorithm = 'weight01',
                      statistic = 'ks',
                      node_size = min_size,
                      score_threshold = 1e-5)
    
    res <- topGO::GenTable(bp.res$topgo_data,
                           p.value = bp.res$topgo_res,
                           topNodes = length(bp.res$topgo_res@score)) %>%
      bind_cols(ontology = rep('BP', length(bp.res$topgo_res@score))) %>%
      bind_rows(topGO::GenTable(cc.res$topgo_data,
                                p.value = cc.res$topgo_res,
                                topNodes = length(cc.res$topgo_res@score)) %>%
                  bind_cols(ontology = rep('CC', length(cc.res$topgo_res@score)))) %>%
      bind_rows(topGO::GenTable(mf.res$topgo_data,
                                p.value = mf.res$topgo_res,
                                topNodes = length(mf.res$topgo_res@score)) %>%
                  bind_cols(ontology = rep('MF', length(mf.res$topgo_res@score)))) %>%
      as_tibble() %>%
      select(term = GO.ID, size = Annotated, p.value, ontology, annotation = Term) %>%
      mutate(p.value = as.numeric(p.value)) %>%
      arrange(p.value)
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
                                                 "closestBed. If a directory is passed, then",
                                                 "for each input file in the <input> directory,",
                                                 "a file that has the same prefix will be expected.",
                                                 "E.g. if an inout file is of the form",
                                                 "<input>/<prefix><suffix>, then there should be a",
                                                 "file of the form <closest>/<prefix>*. If nothing",
                                                 "is passed, then each tab-delimited file must have",
                                                 "a 'gene_id' column."),
                    class = "character",
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
                    help = paste("Which HMVAR function to use. Either 'test_go' or 'gsea'"),
                    type = "character",
                    default = 'gsea')

  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  args$suffix <- paste0(args$suffix, "$")
  if(!(args$method %in% c('test_go', 'gsea'))){
    stop("ERROR: --method must be 'test_go' or 'gsea'", call. = TRUE)
  }
  if(!(args$gene_score %in% c('min', 'max'))){
    stop("ERROR: --gene_score must be 'min' or 'max'")
  }
  
  return(args)
}

count_vars <- function(d){
  n.variants <- nrow(d)
  n.sig <- sum(d$type %in% c("int", "both"))
  n.outside <- sum(d$dist > 0)
  tibble(chr = unique(d$chr),
         start = unique(d$start),
         end = unique(d$end),
         n.variants = n.variants,
         n.sig = n.sig,
         n.outside = n.outside)
}

metawas_gene_counts <- function(metawas, closest, annot, outdir = "./", prefix = NULL){
  # Match everything per gene
  all <- metawas %>% full_join(closest, by = c("chr", "ps"))
  
  cat("\tCleaning annotation...\n")
  annot <- annot %>%
    select(gene_id, predicted_gene_name, eggNOG_annot)
  
  
  Genes <- all %>%
    split(.$gene_id) %>%
    map_dfr(~count_vars(.), .id = "gene_id") %>%
    arrange(chr, start) %>%
    left_join(annot, by = "gene_id")
  
  filename <- paste0(c(prefix, "gene_metawas_counts.txt"), collapse = ".")
  filename <- file.path(outdir, filename)
  write_tsv(Genes, filename)
  
  return(filename)
}

#' Sums named vectors by name
#'
#' @param a A named vector.
#' @param b A named vector.
#'
#' @return A named vector.
#' @export
#'
#' @examples
sum_vecs <- function(a, b){
  if(is.null(a)){
    return(b)
  }
  if(is.null(b)){
    return(a)
  }
  
  new <- tibble(term = names(a), count = a) %>%
    full_join(tibble(term = names(b), count = b), by = "term")
  new <- new %>%
    bind_cols(count = new %>% select(-term) %>% rowSums(na.rm = TRUE)) %>%
    select(term, count)
  vec <- new$count
  names(vec) <- new$term
  
  return(vec)
}

# args <- list(input = "~/micropopgen/exp/2019/2019-03-29.hmp_metawas_data/Supragingival.plaque/metawas/lmmpcs/Porphyromonas_sp_57899_lmm.results.txt",
#              closest = "~/micropopgen/exp/2019/2019-03-29.hmp_metawas_data/Supragingival.plaque/closest/Porphyromonas_sp_57899.closest",
#              annotations = "~/micropopgen/exp/2019/2019-04-01.hmp_subsite_annotations/hmp.subsite_annotations/Porphyromonas_sp_57899.emapper.annotations",
#              dist_thres = 500,
#              min_size = 3,
#              outdir = "metawas_enrichments",
#              suffix = ".txt$",
#              score_column = 'p_lrt.lmmpcs',
#              annot_column = 'GO_terms',
#              alternative = 'less',
#              method = 'gsea',
#              gene_score = 'min')
# args <- process_arguments()

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
  closest_files <- str_subset(closest_files, pattern = prefixes)
  annot_files <- list.files(args$annotations)
  annot_files <- str_subset(annot_files, pattern = prefixes)
  if(length(closest_files) != length(inputs) || length(annot_files) != length(inputs)){
    stop("ERROR: number of closest, annotation or input files does not match.")
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
    
    Res <- gw_test_enrichments(input = input, annotations = annotations,
                               closest = closest, dist_thres = args$dist_thres,
                               score_column = args$score_column,
                               gene_score = args$gene_score,
                               alternative = args$alternative,
                               annot_column = args$annot_column,
                               method = args$method, min_size = args$min_size)
    
    filename <- paste0(prefix, ".enrichments.txt")
    filename <- file.path(args$outdir, filename)
    write_tsv(Res$res, path = filename)
    
    Dat <- Dat %>% bind_rows(Res$res)
  }
}else if(file.exists(args$input)){
  Res <- gw_test_enrichments(input = args$input, annotations = args$annotations,
                             closest = args$closest, dist_thres = args$dist_thres,
                             score_column = args$score_column,
                             gene_score = args$gene_score,
                             alternative = args$alternative,
                             annot_column = args$annot_column,
                             method = args$method, min_size = args$min_size)
  
  prefix <- str_replace(basename(args$input), args$suffix, "")
  filename <- paste0(prefix, ".enrichments.txt")
  filename <- file.path(args$outdir, filename)
  write_tsv(Res$res, path = filename)
  
}else{
  stop("ERROR: input doesn't exist")
}

# 
# # Metawas gene counts
# cat("Counting genes...\n")
# metawas_gene_counts(metawas = metawas,
#                     closest = closest,
#                     annot = annot,
#                     outdir = args$outdir,
#                     prefix = args$prefix)
# 
# # Test GO
# if(nrow(annot) > 0){
#   cat("Test GO...\n")
#   process_annotation(annot = annot,
#                      sig_genes = sig_genes,
#                      annotation = "GO_terms",
#                      outdir = args$outdir,
#                      prefix = args$prefix,
#                      test = TRUE,
#                      count_thres = args$count_thres,
#                      match_go = TRUE)
#   # Test KO
#   cat("Test KO...\n")
#   process_annotation(annot = annot,
#                      sig_genes = sig_genes,
#                      annotation = "KEGG_KOs",
#                      outdir = args$outdir,
#                      prefix = args$prefix,
#                      test = TRUE,
#                      count_thres = args$count_thres,
#                      match_go = FALSE)
#   # Test eggNOG
#   cat("Test OG...\n")
#   process_annotation(annot = annot,
#                      sig_genes = sig_genes,
#                      annotation = "OGs",
#                      outdir = args$outdir,
#                      prefix = args$prefix,
#                      test = TRUE,
#                      count_thres = args$count_thres,
#                      match_go = FALSE)
# }
