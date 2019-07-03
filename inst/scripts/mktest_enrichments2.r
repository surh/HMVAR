library(tidyverse)
# library(argparser)
library(HMVAR)


setwd("/godot/users/sur/exp/fraserv/2019/today2/")
mkres_dir <- "mkres/"
annot_dir <- "annots/"
suffix <- "_mktest.txt$"
filter_fixed <- FALSE
annot_column <- "GO_terms"
outdir <- "enrichments"



mkres_files <- list.files(mkres_dir)
mkres_files <- str_subset(mkres_files, pattern = suffix)
prefixes <- str_replace(mkres_files, suffix, "")
prefixes <- paste0("^", prefixes)
dir.create(outdir)

annot_files <- list.files(annot_dir)
cat("Total ", length(annot_files), " annotation files\n")
annot_files <- prefixes %>% map(~str_subset(annot_files, .)) %>% unlist %>% unique
cat("Found ", length(annot_files), " annotation files\n")
if(length(annot_files) < length(mkres_files)){
  stop("ERROR: number of closest ot annotation files is less than number of input files.")
}



Dat <- tibble()
for(i in 1:length(mkres_files)){
  # i <- 3
  input <- mkres_files[i]
  prefix <- prefixes[i]
  annotations <- str_subset(annot_files, pattern = prefix)
  input
  prefix
  annotations
  
  if(length(annotations) != 1){
    cat("Input:", input, "\n")
    stop("ERROR: missing or extra annotation file")
  }
  cat(input, "\n")
  # Res <- gw_test_enrichments(input = file.path(args$input, input),
  #                            annotations = file.path(args$annotations, annotations),
  #                            closest = file.path(args$closest, closest),
  #                            dist_thres = args$dist_thres,
  #                            score_column = args$score_column,
  #                            gene_score = args$gene_score,
  #                            alternative = args$alternative,
  #                            annot_column = args$annot_column,
  #                            method = args$method, min_size = args$min_size)
  
  # Read test
  cat("\tReading mkres...\n")
  mkres <- read_tsv(file.path(mkres_dir, input),
                    col_types = cols(gene_id = col_character(),
                                     .default = col_double())) %>%
    filter(!is.na(p.value))
  
  if(filter_fixed){
    mkres <- mkres %>%
      filter(Dn > 0 & Ds > 0)
  }
  mkres
  if(nrow(mkres) == 0) next
  
  
  cat("\tReading annotations...\n")
  annots <- read_eggnog(file.path(annot_dir, annotations)) %>%
    select(gene_id = query_name, everything()) %>%
    select(gene_id, terms = annot_column)
  #### CUSTOM FOR TESTS ####
  annots <- annots %>%
    mutate(gene_id = str_replace(string = gene_id,
                                 pattern = "\\([+-]\\)_[0-9]",
                                 replacement = ""))
  ##########
  annots
  
  # Calculate score
  cat("Calculate signed -log10(p-value)...\n")
  mkres <- mkres %>%
    mutate(score = sign(log(OR)) * (-log10(p.value)))
  
  # Match genes with annotations
  cat("Match genes with annotations \n")
  gw.test <- annots %>%
    right_join(mkres, by = "gene_id") %>%
    select(gene_id, terms, score) %>%
    filter(!is.na(score))
  gw.test
  cat("\tGenes with annotation: ", sum(!is.na(gw.test$terms)), "\n")
  if(nrow(gw.test %>% filter(!is.na(terms)) %>% filter(score != 0)) > 0){
    cat("Testing\n")
    res <- gsea(gw.test, test = 'wilcoxon', alternative = 'two.sided', min_size = 3)
    filename <- paste0(str_replace(prefix, '^\\^', ''), ".enrichments.txt")
    filename <- file.path(outdir, filename)
    filename
    write_tsv(res, path = filename)
  }
  
  Dat <- Dat %>% bind_rows(gw.test)
}



term_gsea <- function(genes, scores, test = "wilcoxon",
                      alternative = "greater", min_size = 3){
  
  if(!is.character(genes)){
    stop("ERROR: genes must be a character vector", call. = TRUE)
  }
  if(!is.numeric(scores) || is.null(attr(scores, "names"))){
    stop("ERROR: scores must be a named numeric vector", call. = TRUE)
  }
  
  ii <- names(scores) %in% genes
  
  if(sum(ii) < min_size){
    return(tibble::tibble(size = sum(ii),
                          median = NA,
                          mean = NA,
                          bg.median = NA,
                          bg.mean = NA,
                          statistic = NA,
                          p.value = NA))
  }
  if(sum(!ii) < 1){
    return(tibble::tibble(size = sum(ii),
                          median = NA,
                          mean = NA,
                          bg.median = NA,
                          bg.mean = NA,
                          statistic = NA,
                          p.value = NA))
  }
  
  if(test == 'wilcoxon'){
    res <- wilcox.test(scores[ii], scores[!ii], alternative = alternative)
  }else if(test == 'ks'){
    res <- ks.test(scores[ii], scores[!ii], alternative = alternative)
  }else{
    stop("ERROR: Invalid test", call. = TRUE)
  }
  
  tibble::tibble(size = sum(ii),
                 median = median(scores[ii]),
                 mean = mean(scores[ii]),
                 bg.median = median(scores),
                 bg.mean = mean(scores),
                 statistic = res$statistic,
                 p.value = res$p.value)
}


Dat
Dat %>% filter(!is.na(terms)) %>% filter(score != 0) %>% head(200)

gsea(dat = Dat %>% filter(!is.na(terms)) %>% filter(score != 0) %>% head(200),
     test = 'wilcoxon', alternative = 'two.sided', min_size = 3) %>%
  filter(!is.na(p.value))
