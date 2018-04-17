library(HMAR)
library(GO.db)
library(qvalue)

setwd("~/micropopgen/exp/2018/today/go_enrichment/")

Tab <- read.table("GO_enrichment_by_comparsions.txt", sep = "\t", header = TRUE)
head(Tab)

Tab$qval <- NULL
for(c in levels(Tab$comparison)){
  cat(c, "\n")
  Tab$qval[ Tab$comparison == c ] <- p.adjust(Tab$pval[ Tab$comparison == c ], 'fdr')
  summary(qvalue::qvalue(Tab$pval[ Tab$comparison == c ]))
}
p1 <- ggplot(Tab, aes(x = pval)) +
  facet_wrap(~ comparison, ncol = 1) +
  geom_histogram(bins = 20) +
  AMOR::theme_blackbox
p1
ggsave("pvalues_go_enrich_by_comparioson.svg", p1, width = 4, height = 10)

head(Tab)

Tab$GO_term <- NA
Tab$GO_desc <- NA
Tab$GO_ontology <- NA
for(i in 1:nrow(Tab)){
  # i <- 1
  
  go <- GO.db::GOTERM[[as.character(Tab$annotation[i])]]
  if(!is.null(go)){
    Tab$GO_term[i] <- go@Term
    Tab$GO_desc[i] <- go@Definition
    Tab$GO_ontology[i] <- go@Ontology
  }
  
  
}
write.table("GO_enrichment_by_comparsions_processed.txt")



###########