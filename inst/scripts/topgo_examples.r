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



