# (C) Copyright 2019 Sur Herrera Paredes
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

library(HMAR)
library(GO.db)
library(qvalue)

Sys.time()

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
write.table(Tab, "GO_enrichment_by_comparsions_processed.txt", col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")



########### Repeat for species and comparisons ###########

Tab <- read.table("GO_enrichment_by_comparison_and_species.txt", sep = "\t", header = TRUE)
head(Tab)

Tab$qval <- NULL
for(c in levels(Tab$comparison)){
  cat(c, "\n")
  Tab$qval[ Tab$comparison == c ] <- p.adjust(Tab$pval[ Tab$comparison == c ], 'fdr')
  # summary(qvalue::qvalue(Tab$pval[ Tab$comparison == c ], lambda = seq(0.01, 0.95, 0.05)))
}
p1 <- ggplot(Tab, aes(x = pval)) +
  facet_wrap(~ comparison, ncol = 2) +
  geom_histogram(bins = 20) +
  AMOR::theme_blackbox
p1
ggsave("pvalues_go_enrich_by_comparison_and_species.svg", p1, width = 8, height = 20)

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
write.table(Tab, "GO_enrichment_by_comparisons_and_species_processed.txt", col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")





