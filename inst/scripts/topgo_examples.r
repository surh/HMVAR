library(HMVAR)

dat <- read_eggnog("~/micropopgen/exp/2019/2019-04-01.hmp_subsite_annotations/annotations/Catonella_morbi_61904.emapper.annotations")
set.seed(123)
annots <- dat %>%
  dplyr::select(gene_id = query_name, terms = GO_terms) %>%
  dplyr::bind_cols(score = rnorm(nrow(dat)))
annots

# res <- test_go(genes = sig_genes, annots = annots)
# 
# topGO::GenTable(res$topgo_data, fisher = res$topgo_res)
# 

gsea(dat = annots, test = 'wilcoxon', alternative = 'greater', min_size = 10)


set.seed(123)
gene_scores <- rnorm(length(annots$gene_id))
names(gene_scores) <- annots$gene_id
gene_scores[1:200] <- gene_scores[1:200] + rnorm(200, mean = 1)
hist(gene_scores)
gene_scores[1:10]

alternative <- 'greater'
test <- "wilcoxon"
min_size <- 10






d <- annots %>%
  split(.$term) %>%
  purrr::map(~ .x$gene_id) %>%
  purrr::map_int(length) %>% 

sum(d >= 10)

annots %>%
  split(.$term) %>%
  purrr::map(~ .x$gene_id) %>%
  length

D <- annots %>%
  split(.$term) 

d <- D[[8]]
d

ii <- names(gene_scores) %in% d$gene_id

res <- ks.test(1:10, 21:30, alternative = 'greater')

res <- wilcox.test(1:10, 21:30, alternative = 'greater')

