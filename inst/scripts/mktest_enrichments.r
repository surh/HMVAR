library(tidyverse)

process_annotation <- function(mkres, sig_genes,
                               annotation = "GO_terms",
                               count_thres = 3,
                               match_go = TRUE){
  # annotation <- "KEGG_KOs"
  
  BG <- mkres %>%
    select(gene_id, annot = annotation) %>%
    pmap_dfr(expand_annot) %>%
    mutate(sig_gene = gene_id %in% sig_genes) %>%
    filter(!is.na(term))
  
  res <- test_annotation(BG = BG, count_thres = count_thres, match_go = match_go)
  
  return(res)
}

expand_annot <- function(gene_id, annot){
  annot <- str_split(string = annot, pattern = ",") %>% unlist
  tibble(gene_id = gene_id,
         term = annot)
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
                                        columns = c("ONTOLOGY","TERM","DEFINITION"))) %>%
      arrange(p.value)
  }
  
  return(test_res)
}

process_mkres <- function(spec, mktest, annots){
  # i <- 1
  # spec <- files$spec[i]
  # mktest <- files$mktest[i]
  # annots <- files$annots[i]
  
  cat(spec, "\n")
  Mktest <- read_tsv(mktest, 
                     col_types = cols(.default = col_double(),
                                      gene_id = col_character())) %>%
    mutate(kaks = (Dn + Pn) / (Ds + Ps),
           mkratio = ((Dn + 1) / (Ds + 1)) / ((Pn + 1) / (Ps + 1)))
  
  Annot <- read_tsv(annots,
                    comment = "#",
                    col_names = c("query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue",
                                  "seed_ortholog_score",	"predicted_gene_name", "GO_terms",
                                  "KEGG_KOs", "BiGG_reactions", "Annotation_tax_scope", "OGs",
                                  "bestOG|evalue|score", "COG_cat", "eggNOG_annot"),
                    col_types = cols(.default = col_character(),
                                     seed_ortholog_score = col_double(),
                                     seed_ortholog_evalue = col_double())) %>%
    mutate(gene_id = str_replace(string = query_name,
                                 pattern = "\\([+-]\\)_[0-9]",
                                 replacement = "")) %>%
    select(gene_id, everything(), -query_name, -seed_ortholog_evalue, -seed_ortholog_score,
           -Annotation_tax_scope)
  
  # Select columns
  Annot <- Annot %>%
    select(gene_id, predicted_gene_name,
           eggNOG_annot, GO_terms, KEGG_KOs, OGs)
  
  
  # Join
  Mktest <- Mktest %>%
    filter(Ds > 0 & Pn > 0) %>%
    left_join(Annot, by = "gene_id") %>%
    mutate(spec = spec)
  
  return(Mktest)
}


# args <- list(mkdir = "../2019-04-02.hmp_mktest_data/Buccal.mucosa/results/",
#              annotdir = "../2019-04-01.hmp_subsite_annotations/hmp.subsite_annotations/",
#              outdir = "Buccal.mucosa/")
args <- list(mkdir = opts[1],
             annotdir = opts[2],
             outdir = opts[3])

# Get list of files
mktest <- tibble(mktest = file.path(args$mkdir, list.files(args$mkdir)))
mktest <- mktest %>% 
  mutate(spec = str_replace(basename(mktest), pattern = "_mktest.txt$", "")) %>%
  select(spec, mktest)
mktest

annots <-  list.files(args$annotdir) %>%
  str_subset(pattern = ".emapper.annotations$") %>%
  file.path(args$annotdir, .) %>%
  tibble(spec = str_replace(basename(.), pattern = ".emapper.annotations$", ""),
         annots = .)
annots

files <- mktest %>%
  left_join(annots, by = "spec") %>%
  filter(!is.na(annots))
rm(mktest, annots)
files

dir.create(args$outdir)

mkres <- files %>%
  pmap_dfr(process_mkres)
mkres$p.value[ mkres$p.value > 1 ] <- 1
filename <- file.path(args$outdir, "mkres_combined.txt")
write_tsv(mkres, filename)

mkres %>% arrange(p.value)
mkres %>% arrange(desc(Dn))
mkres %>% arrange(desc(OR))
mkres %>% arrange(desc(mkratio))
mkres %>% arrange(desc(kaks))

p1 <- mkres %>%
  ggplot(aes(x = mkratio, y = -log10(p.value))) +
  geom_point(aes(size = Dn), alpha = 0.5) + 
  scale_x_log10() +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "mkres_combined_volcano.png")
ggsave(filename, p1, width = 6, height = 5, dpi = 200)

sig_genes <- mkres %>% filter(OR > 1) %>%
  arrange(p.value)
sig_genes
sig_genes <- sig_genes %>%
  select(gene_id) %>%
  unique() %>%
  unlist

go_res <- process_annotation(mkres = mkres,
                             sig_genes = sig_genes,
                             annotation = "GO_terms",
                             count_thres = 3, match_go = TRUE)
go_res
filename <- file.path(args$outdir, "go_enrichments.txt")
write_tsv(go_res, filename)

ko_res <- process_annotation(mkres = mkres,
                             sig_genes = sig_genes,
                             annotation = "KEGG_KOs",
                             count_thres = 3, match_go = FALSE)
ko_res
filename <- file.path(args$outdir, "ko_enrichments.txt")
write_tsv(ko_res, filename)

og_res <- process_annotation(mkres = mkres,
                             sig_genes = sig_genes,
                             annotation = "OGs",
                             count_thres = 3, match_go = FALSE)
og_res
filename <- file.path(args$outdir, "og_enrichments.txt")
write_tsv(og_res, filename)

