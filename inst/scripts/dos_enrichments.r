library(tidyverse)
# library(argparser)
library(HMVAR)


setwd("/godot/users/sur/exp/fraserv/2019/today/")
mkres_dir <- "dos/"
annot_dir <- "annots/"
suffix <- ".DoS.table.txt$"
filter_fixed <- FALSE
annot_column <- "GO_terms"
outdir <- "GO"


mkres_files <- list.files(mkres_dir)
mkres_files <- str_subset(mkres_files, pattern = suffix)
mkres_files <- setdiff(mkres_files,"summary.DoS.table.txt")
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
  # i <- 5
  input <- mkres_files[i]
  prefix <- prefixes[i]
  annotations <- str_subset(annot_files, pattern = prefix)
  # input
  # prefix
  # annotations
  
  if(length(annotations) != 1){
    cat("Input:", input, "\n")
    stop("ERROR: missing or extra annotation file")
  }
  cat(input, "\n")
  
  # Read test
  cat("\tReading DoS...\n")
  mkres <- read_tsv(file.path(mkres_dir, input),
                    col_types = cols(gene_id = col_character(),
                                     .default = col_double())) %>%
    filter(!is.na(p.value))
  
  if(filter_fixed){
    mkres <- mkres %>%
      filter(Dn > 0 & Ds > 0)
  }
  # mkres
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
  # annots
  
  # Calculate score
  cat("Calculate signed -log10(p-value)...\n")
  mkres <- mkres %>%
    mutate(score = sign(DoS) * (-log10(p.value)))
  
  # Match genes with annotations
  cat("Match genes with annotations \n")
  gw.test <- annots %>%
    right_join(mkres, by = "gene_id") %>%
    select(gene_id, terms, score) %>%
    filter(!is.na(score)) %>%
    filter(!is.na(terms))
  # gw.test
  cat("\tGenes with annotation: ", sum(!is.na(gw.test$terms)), "\n")
  if(nrow(gw.test %>% filter(!is.na(terms)) %>% filter(score != 0)) > 0){
    cat("Testing\n")
    res <- gsea(gw.test, test = 'wilcoxon', alternative = 'two.sided', min_size = 3) %>%
      filter(!is.na(p.value))
    res <- HMVAR:::annotate_gos(res, colname = 'term')
    
    if(nrow(res) > 0){
      filename <- paste0(str_replace(prefix, '^\\^', ''), ".enrichments.txt")
      filename <- file.path(outdir, filename)
      # filename
      write_tsv(res, path = filename)
    }
  }
  
  Dat <- Dat %>% bind_rows(gw.test)
}

# sign(res$mean - res$bg.mean)
# res
filename <- file.path(outdir, "all_genes_terms_and_scores.txt")
write_tsv(Dat, path = filename)


Dat
Dat %>% filter(!is.na(terms)) %>% filter(score != 0) %>% head(200)

Res <- gsea(dat = Dat,
            test = 'wilcoxon',
            alternative = 'two.sided',
            min_size = 3) %>%
  filter(!is.na(p.value))
Res <- HMVAR:::annotate_gos(Res, colname = 'term')
Res
filename <- file.path(outdir, "summary.enrichments.txt")
write_tsv(Res, path = filename)


Res %>% 
  select(-median, -bg.median, -statistic) %>%
  mutate(sign = sign(mean - bg.mean)) %>%
  select(-mean, -bg.mean) %>%
  # filter(p.adjust(p.value, 'fdr') < 0.5) %>%
  print(n = 100)

