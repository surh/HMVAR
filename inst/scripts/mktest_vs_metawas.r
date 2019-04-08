library(tidyverse)


args <- list(mktest = "../2019-04-02.hmp_mktest_data/Buccal.mucosa/results/Porphyromonas_sp_57899_mktest.txt",
             genes = "../today3/Buccal.mucosa/enrichments/enrichments/Porphyromonas_sp_57899.gene_metawas_counts.txt")


# metawas <- read_tsv(args$lmmres,
#                     col_types = cols(.default = col_double(),
#                                      chr = col_character(),
#                                      rs = col_character(),
#                                      allele1 = col_character(),
#                                      allele0 = col_character(),
#                                      type = col_character())) %>%
#   select(-starts_with("logl_H1"), -starts_with("l_mle"))
# metawas

genes <- read_tsv(args$genes,
                  col_types = cols(.default = col_character(),
                                   start = col_number(),
                                   end = col_number(),
                                   n.variants = col_number(),
                                   n.sig = col_number(),
                                   n.outside = col_number()))
genes

mktest <- read_tsv(args$mktest, 
                   col_types = cols(.default = col_number(),
                                    gene_id = col_character()))
mktest <- mktest %>%
  mutate(kaks = (Dn + Pn) / (Ds + Ps),
         mkratio = ((Dn + 1) / (Ds + 1)) / ((Pn + 1) / (Ps + 1)))
mktest

genes <- genes %>%
  left_join(mktest, by = "gene_id") %>%
  filter(!is.na(Dn)) %>%
  filter(!is.na(kaks)) %>%
  filter(kaks != Inf)


test <- wilcox.test(genes$kaks[genes$n.sig == 0],
                    genes$kaks[genes$n.sig > 0],
                    alternative = "less")


p1 <- genes %>%
  ggplot(aes(x = n.sig > 0, y = kaks)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_point(position = position_jitter(width = 0.2), size = 0.1) +
  scale_y_log10() +
  ggtitle(label = "", subtitle = test$p.value) +
  AMOR::theme_blackbox()
p1

test <- wilcox.test(genes$mkratio[genes$n.sig == 0],
                    genes$mkratio[genes$n.sig > 0],
                    alternative = "less")
p1 <- genes %>%
  ggplot(aes(x = n.sig > 0, y = mkratio)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_point(position = position_jitter(width = 0.2), size = 0.1) +
  scale_y_log10() +
  ggtitle(label = "", subtitle = test$p.value) +
  AMOR::theme_blackbox()
p1

genes %>%
  filter(n.sig > 0) %>%
  ggplot(aes(x = n.sig / n.variants, y = kaks)) +
  geom_point() +
  geom_smooth(method = "lm") +
  AMOR::theme_blackbox()

genes %>%
  filter(n.sig > 0) %>%
  ggplot(aes(x = n.sig / n.variants, y = mkratio)) +
  geom_point() +
  geom_smooth(method = "lm") +
  AMOR::theme_blackbox()











