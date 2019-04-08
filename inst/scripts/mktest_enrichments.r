library(tidyverse)

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



args <- list(mkdir = "../2019-04-02.hmp_mktest_data/Buccal.mucosa/results/",
             annotdir = "../2019-04-01.hmp_subsite_annotations/hmp.subsite_annotations/",
             outdir = "Buccal.mucosa/")

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


mkres <- files %>%
  pmap_dfr(process_mkres)

mkres %>% arrange(p.value)
mkres %>% arrange(desc(Dn))
mkres %>% arrange(desc(OR))
mkres %>% arrange(desc(mkratio))
mkres %>% arrange(desc(kaks))

p1 <- mkres %>%
  ggplot(aes(x = mkratio, y = p.value)) +
  geom_point(aes(size = Dn)) +
  scale_y_log10() +
  scale_x_log10() +
  AMOR::theme_blackbox()

