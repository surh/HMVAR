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


# Copyright (C) 2019 Sur Herrera Paredes
# Performs functional enrichment on metawas hits from eggnog annotations
# Based on script from /home/sur/micropopgen/exp/2019/2019-01-29.metawas_enrichment/

library(tidyverse)
library(argparser)

expand_annot <- function(gene_id, annot){
  annot <- str_split(string = annot, pattern = ",") %>% unlist
  tibble(gene_id = gene_id,
         term = annot)
}

process_annotation <- function(annot,
                               sig_genes,
                               annotation = "GO_terms",
                               outdir = "./",
                               prefix = NULL,
                               test = TRUE,
                               count_thres= 3,
                               match_go = TRUE){
  
  # Create background
  BG <- annot %>%
    select(gene_id, annot = annotation) %>%
    pmap_dfr(expand_annot) %>%
    mutate(sig_gene = gene_id %in% sig_genes) %>%
    filter(!is.na(term))
    
  
  if(test){
    res <- test_annotation(BG = BG,
                           count_thres = count_thres,
                           match_go = match_go)
    filename <- paste0(c(prefix, annotation, "enrichments.txt"), collapse = ".")
    filename <- file.path(outdir, filename)
    write_tsv(res, filename)
  }
  
  # Write background
  filename <- paste0(c(prefix, annotation, "BG.txt"), collapse = ".")
  filename <- file.path(outdir, filename)
  write_tsv(BG, filename)
  
  
  return(filename)
}

test_annotation <- function(BG, count_thres = 3, match_go = TRUE){
  # Get list to test
  to_test <- BG %>%
    filter(sig_gene) %>%
    count(term) %>%
    filter(n >= count_thres) %>%
    select(term, n.sig = n)
  
  # Get BG counts
  bg_counts <- BG %>%
    filter(term %in% to_test$term) %>%
    count(term) %>%
    select(term, n.bg = n)
  
  selection.size <- BG %>%
    filter(sig_gene) %>%
    select(gene_id) %>%
    unique() %>%
    nrow()
  
  bg.size <- BG$gene_id %>% unique %>% length
  
  # Test
  test_res <- to_test %>%
    left_join(bg_counts, by = "term") 
  
  if(nrow(test_res) == 0){
    return(test_res)
  }
  
  test_res <- test_res %>%
    pmap_dfr(function(term, n.sig, n.bg, selection.size, bg.size){
      mat <- matrix(c(n.sig, n.bg, selection.size, bg.size), ncol = 2)
      res <- fisher.test(mat)
      tibble(term = term,
             n.sig = n.sig,
             n.bg = n.bg,
             OR = res$estimate,
             p.value = res$p.value)
    }, selection.size = selection.size, bg.size = bg.size) %>%
    mutate(q.value = p.adjust(p.value, 'fdr')) %>%
    arrange(q.value)
  
  if(match_go){
    # Match metadata
    test_res <- test_res %>%
      bind_cols(test_res %>%
                  select(term) %>%
                  unlist %>%
                  AnnotationDbi::select(GO.db::GO.db, .,
                                        columns = c("ONTOLOGY","TERM","DEFINITION")))
  }
  
  return(test_res)
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

process_arguments <- function(){
  p <- arg_parser(paste(""))
  
  # Positional arguments
  p <- add_argument(p, "lmmres",
                    help = paste(""),
                    type = "character")
  p <- add_argument(p, "closest", help = "")
  p <- add_argument(p, "annotations", help = "")
  
  # Optional arguments
  p <- add_argument(p, "--dist_thres",
                     help = paste(""),
                     type = "numeric",
                     default = 500)
  p <- add_argument(p, "--count_thres", help = "",
                    default = 3)
  p <- add_argument(p, "--outdir", help = "",
                    default = "output")
  p <- add_argument(p, "--prefix", help = "",
                     default = NULL,
                     type = "character")
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(is.na(args$prefix)){
    args$prefix <- NULL
  }
  
  return(args)
}

# args <- list(lmmres = "../2019-03-29.hmp_metawas_data/Supragingival.plaque/metawas/lmmpcs/Porphyromonas_sp_57899_lmm.results.txt",
#              closest = "../2019-03-29.hmp_metawas_data/Supragingival.plaque/closest/Porphyromonas_sp_57899.closest",
#              annotations = "../2019-04-01.hmp_subsite_annotations/hmp.subsite_annotations/Porphyromonas_sp_57899.emapper.annotations",
#              dist_thres = 500,
#              count_thres = 3,
#              outdir = "metawas_enrichments",
#              prefix = NULL)
# args <- list(lmmres = "Streptococcus_sp_60488_lmm.results.txt",
#              closest = "Streptococcus_sp_60488.closest",
#              annotations = "Streptococcus_sp_60488.emapper.annotations",
#              dist_thres = 500,
#              count_thres = 3,
#              outdir = "metawas_enrichments",
#              prefix = NULL)
args <- process_arguments()

dir.create(args$outdir)

# manhat_dir <- paste0(outdir, "/manhattan")
# dir.create(manhat_dir)
# genes_dir <- paste0(outdir, "/genes/")
# dir.create(genes_dir)
# enrich_dir <- paste0(outdir, "/enrich/")
# dir.create(enrich_dir)
# 
# genomes <- read_tsv(genomes_file, col_names = FALSE, col_types = 'c')
# genomes <- genomes$X1
# 
# GENES <- NULL
# OG_bg <- NULL
# GO_bg <- NULL
# KO_bg <- NULL

# Read data
cat("Reading lmm...\n")
metawas <- read_tsv(args$lmmres,
                    col_types = cols(.default = col_double(),
                                     chr = col_character(),
                                     rs = col_character(),
                                     allele1 = col_character(),
                                     allele0 = col_character(),
                                     type = col_character())) %>%
  select(-starts_with("logl_H1"), -starts_with("l_mle"))
# metawas
cat("Reading closest...\n")
closest <- read_tsv(args$closest,
                    col_names = c("chr", "ps", "ps2", "chr2", "start", "end", "gene_id", "dist"),
                    col_types = cols(.default = col_number(),
                                     chr = col_character(),
                                     chr2 = col_character(),
                                     gene_id = col_character())) %>%
  select(-ps2, -chr2)
# closest

# Get genes
cat("Getting gene lists...\n")
genes_tested <- closest %>%
  filter(abs(dist) <= args$dist_thres) %>%
  select(gene_id) %>% unique()
sig_genes <- closest %>%
  left_join(metawas %>%
              filter(type %in% c("int", "both")) %>%
              select(chr, rs, ps),
            by = c("chr", "ps")) %>%
  filter(!is.na(rs)) %>%
  filter(abs(dist) <= args$dist_thres) %>%
  select(gene_id) %>%
  unique() %>%
  unlist()

cat("Reading annot file...\n")
# annot <- read_tsv(annot_file, comment = "#", col_names = FALSE)
annot <- read_tsv(args$annotations,
                  comment = "#",
                  col_names = c("query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue",
                                "seed_ortholog_score",	"predicted_gene_name", "GO_terms",
                                "KEGG_KOs", "BiGG_reactions", "Annotation_tax_scope", "OGs",
                                "bestOG|evalue|score", "COG_cat", "eggNOG_annot"),
                  col_types = cols(.default = col_character(),
                                   seed_ortholog_score = col_double(),
                                   seed_ortholog_evalue = col_double()))
# annot
# Reformat gene name
cat("Reformatting annotation...\n")
annot <- annot %>%
  mutate(gene_id = str_replace(string = query_name,
                               pattern = "\\([+-]\\)_[0-9]",
                               replacement = "")) %>%
  select(gene_id, everything(), -query_name, -seed_ortholog_evalue, -seed_ortholog_score,
         -Annotation_tax_scope)
# annot

# Select onl tested genes
# Only these will be considered in the universe background
cat("Selecting annotations from tested genes...\n")
annot <- annot %>% filter(gene_id %in% genes_tested$gene_id)


# Metawas gene counts
cat("Counting genes...\n")
metawas_gene_counts(metawas = metawas,
                    closest = closest,
                    annot = annot,
                    outdir = args$outdir,
                    prefix = args$prefix)

# Test GO
if(nrow(annot) > 0){
  cat("Test GO...\n")
  process_annotation(annot = annot,
                     sig_genes = sig_genes,
                     annotation = "GO_terms",
                     outdir = args$outdir,
                     prefix = args$prefix,
                     test = TRUE,
                     count_thres = args$count_thres,
                     match_go = TRUE)
  # Test KO
  cat("Test KO...\n")
  process_annotation(annot = annot,
                     sig_genes = sig_genes,
                     annotation = "KEGG_KOs",
                     outdir = args$outdir,
                     prefix = args$prefix,
                     test = TRUE,
                     count_thres = args$count_thres,
                     match_go = FALSE)
  # Test eggNOG
  cat("Test OG...\n")
  process_annotation(annot = annot,
                     sig_genes = sig_genes,
                     annotation = "OGs",
                     outdir = args$outdir,
                     prefix = args$prefix,
                     test = TRUE,
                     count_thres = args$count_thres,
                     match_go = FALSE)
}
