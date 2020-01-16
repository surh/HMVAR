#!/usr/bin/env Rscript

# (C) Copyright 2020 Sur Herrera Paredes
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

library(argparser)
library(tidyverse)
library(HMVAR)

process_arguments <- function(){
  p <- arg_parser(paste("Test functional enrichments from RER permutations test"))
  
  # Positional arguments
  p <- add_argument(p, "input",
                    help = paste("Input with RER permutations results. Directory or file."),
                    type = "character")
  p <- add_argument(p, "annots",
                    help = paste("eggnog annotations. Directory or file."),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--plot",
                    help = paste("Plot rerperms?"),
                    flag = TRUE)
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(dir.exists(args$input)){
    if(!dir.exists(args$annots)){
      stop("ERROR: annots must be a direcotry if input is a directory", call. = TRUE)
    }
    args$input <- list.files(args$input, full.names = TRUE)
    spec <- stringr::str_remove(basename(args$input), pattern = "\\.rer\\.fdr\\.txt$")
    args$annots <- file.path(args$annots,paste0(spec, ".emapper.annotations"))
    if(!all(file.exists(args$annots))){
      cat(args$annots[ !file.exists(args$annots) ], "not found\n")
      stop("ERROR: There must be annotation files for all genomes tested", call. = TRUE)
    }
  }else if(file.exists(args$input)){
    if(dir.exists(args$annots)){
      spec <- stringr::str_remove(basename(args$input), pattern = "\\.rer\\.fdr\\.txt$")
      args$annots <- file.path(args$annots, paste0(spec, ".emapper.annotations"))
    }
    if(!file.exists(args$annots)){
      stop("ERROR: annots must be an existing file or a dir containing <spec>.emapper.annotations", call. = TRUE)
    }
  }else{
    stop("ERROR: input must be either an existing file or directory.", call. = TRUE)
  }
  
  return(args)
}


# setwd("/cashew/users/sur/exp/fraserv/2020/today/gut")
# args <- list(input = "rerperms/Bacteroides_ovatus_58035.rer.fdr.txt",
#              annots = "annots/Bacteroides_ovatus_58035.emapper.annotations",
#              outdir = "output/",
#              plot = TRUE)
args <- list(input = "rerperms/Bacteroides_ovatus_58035.rer.fdr.txt",
             annots = "annots/",
             outdir = "output/",
             plot = TRUE)
# args <- list(input = "rerperms/",
#              annots = "annots/",
#              outdir = "output/",
#              plot = TRUE)
# args <- list(input = "rerperms/",
#              annots = "annots/Bacteroides_coprocola_61586.emapper.annotations",
#              outdir = "output/",
#              plot = TRUE)
if(dir.exists(args$input)){
  if(!dir.exists(args$annots)){
    stop("ERROR: annots must be a direcotry if input is a directory", call. = TRUE)
  }
  args$input <- list.files(args$input, full.names = TRUE)
  spec <- stringr::str_remove(basename(args$input), pattern = "\\.rer\\.fdr\\.txt$")
  args$annots <- file.path(args$annots,paste0(spec, ".emapper.annotations"))
  if(!all(file.exists(args$annots))){
    cat(args$annots[ !file.exists(args$annots) ], "not found\n")
    stop("ERROR: There must be annotation files for all genomes tested", call. = TRUE)
  }
}else if(file.exists(args$input)){
  if(dir.exists(args$annots)){
    spec <- stringr::str_remove(basename(args$input), pattern = "\\.rer\\.fdr\\.txt$")
    args$annots <- file.path(args$annots, paste0(spec, ".emapper.annotations"))
  }
  if(!file.exists(args$annots)){
    stop("ERROR: annots must be an existing file or a dir containing <spec>.emapper.annotations", call. = TRUE)
  }
}else{
  stop("ERROR: input must be either an existing file or directory.", call. = TRUE)
}

args <- process_arguments()

args

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}


GO <- KO <- OG <- COG <- BIGG <- NULL
for(i in 1:length(args$input)){
  rertest_file <- args$input[i]
  annots_file <- args$annots[i]
  
  # Get species name
  spec <- stringr::str_remove(basename(rertest_file),
                              pattern = "\\.rer\\.fdr\\.txt$")
  
  # Read data
  RER <- read_tsv(rertest_file,
                  col_types = cols(gene_id = col_character(),
                                   Rho = col_number(),
                                   P = col_number(),
                                   p.adj = col_number(),
                                   FDR = col_number()))
  annots <- read_eggnog(annots_file) %>%
    mutate(query_name = str_replace(string = query_name,
                                    pattern = "\\([+-]\\)_[0-9]+$",
                                    replacement = "")) %>%
    select(gene_id = query_name, predicted_gene_name,
           GO_terms, KEGG_KOs, BiGG_reactions,
           Annotation_tax_scope, OGs,
           'bestOG|evalue|score', COG_cat, eggNOG_annot)
  
  if(args$plot){
    # p1 <- ggplot(RER, aes(x = -log10(P), y = -log10(FDR))) +
    #   geom_point() +
    #   geom_abline(intercept = 0, slope = 1) +
    #   AMOR::theme_blackbox()
    # 
    # pval_qqplot(RER$P)
    # pval_qqplot(RER$FDR)
    
    p1 <- ggplot(RER, aes(x = Rho, y = -log10(FDR))) +
      geom_point(aes(size = N), shape = 19, alpha = 0.1) +
      scale_radius(trans="sqrt") +
      AMOR::theme_blackbox()
    filename <- file.path(args$outdir, paste0(spec, ".volcano.png"))
    ggsave(filename, p1, width = 5, height = 5, dpi = 150)
  }
  
  # Match genes and annots
  RER <- RER %>%
    filter(!is.na(FDR)) %>%
    left_join(annots, by = "gene_id")
  # RER %>%
  #   select(-p.adj,
  #          -BiGG_reactions,
  #          -Annotation_tax_scope,
  #          -"bestOG|evalue|score") %>%
  #   print(n = 30)
  # RER %>%
  #   filter(!is.na(BiGG_reactions))
  
  # Enrichments
  terms <- RER %>% select(gene_id, terms = GO_terms, score = FDR) %>% filter(!is.na(terms))
  res  <- terms_enrichment(dat = terms,
                           method = "gsea",
                           test = "ks",
                           alternative = "greater",
                           min_size = 10)
  filename <-  file.path(args$outdir, paste0(spec, ".go.enrichments.png"))
  write_tsv(res, filename)
  GO <- GO %>%
    bind_rows(terms %>% mutate(spec = spec))
  
  terms <- RER %>% select(gene_id, terms = KEGG_KOs, score = FDR) %>% filter(!is.na(terms))
  res  <- terms_enrichment(dat = terms,
                           method = "gsea",
                           test = "ks",
                           alternative = "greater",
                           min_size = 10)
  filename <-  file.path(args$outdir, paste0(spec, ".ko.enrichments.png"))
  write_tsv(res, filename)
  KO <- KO %>%
    bind_rows(terms %>% mutate(spec = spec))
  
  terms <- RER %>% select(gene_id, terms = OGs, score = FDR) %>% filter(!is.na(terms))
  res  <- terms_enrichment(dat = terms,
                           method = "gsea",
                           test = "ks",
                           alternative = "greater",
                           min_size = 10)
  filename <-  file.path(args$outdir, paste0(spec, ".og.enrichments.png"))
  write_tsv(res, filename)
  OG <- OG %>%
    bind_rows(terms %>% mutate(spec = spec))
  
  terms <- RER %>% select(gene_id, terms = BiGG_reactions, score = FDR) %>% filter(!is.na(terms))
  res  <- terms_enrichment(dat = terms,
                           method = "gsea",
                           test = "ks",
                           alternative = "greater",
                           min_size = 10)
  filename <-  file.path(args$outdir, paste0(spec, ".bigg.enrichments.png"))
  write_tsv(res, filename)
  BIGG <- BIGG %>%
    bind_rows(terms %>% mutate(spec = spec))
  
  terms <- RER %>% select(gene_id, terms = COG_cat, score = FDR) %>% filter(!is.na(terms))
  res  <- terms_enrichment(dat = terms,
                           method = "gsea",
                           test = "ks",
                           alternative = "greater",
                           min_size = 10)
  filename <-  file.path(args$outdir, paste0(spec, ".cog.enrichments.png"))
  write_tsv(res, filename)
  COG <- COG %>%
    bind_rows(terms %>% mutate(spec = spec))
}




  



# Try Rho GSEA
RER %>%
  filter(FDR <= 0.01) %>%
  select(gene_id, terms = BiGG_reactions, score = Rho) %>%
  filter(!is.na(terms))
terms_enrichment(dat = RER %>%
                   filter(FDR <= 0.01) %>%
                   select(gene_id, terms = GO_terms, score = Rho),
                 method = "gsea",
                 test = "ks",
                 alternative = "two.sided",
                 min_size = 2)
terms_enrichment(dat = RER %>%
                   filter(FDR <= 0.01) %>%
                   select(gene_id, terms = KEGG_KOs, score = Rho),
                 method = "gsea",
                 test = "ks",
                 alternative = "two.sided",
                 min_size = 2)
terms_enrichment(dat = RER %>%
                   filter(FDR <= 0.01) %>%
                   select(gene_id, terms = OGs, score = Rho),
                 method = "gsea",
                 test = "ks",
                 alternative = "two.sided",
                 min_size = 2)
terms_enrichment(dat = RER %>%
                   filter(FDR <= 0.01) %>%
                   select(gene_id, terms = BiGG_reactions, score = Rho),
                 method = "gsea",
                 test = "ks",
                 alternative = "two.sided",
                 min_size = 2)
terms_enrichment(dat = RER %>%
                   filter(FDR <= 0.01) %>%
                   select(gene_id, terms = COG_cat, score = Rho),
                 method = "gsea",
                 test = "ks",
                 alternative = "two.sided",
                 min_size = 2)


