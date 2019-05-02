library(tidyverse)

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

overall_enrich <- function(indir, suffix, count_thres = 3, match_go = TRUE){
  suffix <- paste0(suffix, "$")
    
  files <- list.files(indir)
  bg_files <- str_subset(string = files, pattern = suffix)
  
  # Read
  BG <- indir %>%
    file.path(bg_files) %>%
    map_dfr(read_tsv, col_types = cols(.default = col_character(),
                                       sig_gene = col_logical()))
  
  # Test
  test_res <- test_annotation(BG, count_thres = count_thres, match_go = match_go)
  
  return(test_res)
}


args <- list(indir = "Buccal.mucosa/enrichments/enrichments/",
             outdir = "overall_enrichments/")
args <- list(indir = opts[1],
             outdir = opts[2])

dir.create(args$outdir)

go_res <- overall_enrich(indir = args$indir,
                         suffix = "GO_terms.BG.txt",
                         count_thres = 3,
                         match_go = TRUE)
go_res
filename <- file.path(args$outdir, "go_enrichments.txt")
write_tsv(go_res, filename)

ko_res <- overall_enrich(indir = args$indir,
                         suffix = "KEGG_KOs.BG.txt",
                         count_thres = 3,
                         match_go = FALSE)
ko_res
filename <- file.path(args$outdir, "ko_enrichments.txt")
write_tsv(ko_res, filename)

og_res <- overall_enrich(indir = args$indir,
                         suffix = "OGs.BG.txt",
                         count_thres = 3,
                         match_go = FALSE)
og_res
filename <- file.path(args$outdir, "og_enrichments.txt")
write_tsv(og_res, filename)
