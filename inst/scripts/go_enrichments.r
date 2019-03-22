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

library(GO.db)
# library(topGO)
library(HMVAR)
library(tidyverse)

######################
go_enrichment <- function(dat, m, n, threshold = 0.01){
  q <- sum(dat$q.value < threshold)
  k <- nrow(dat)
  pval <- phyper(q = q - 1, 
                 m = m,
                 n = n,
                 k = k,
                 lower.tail = FALSE )
  
  tibble(Annotation = unique(dat$Annotation),
         n.genes = k,
         n.sig = q,
         OR = (q / k) / (m / (m+n)),
         p.value = min(pval, 1 - pval) * 2)
}

get_go_enrichment <- function(Genes, outfile,
                              GOTERM = GO.db::GOTERM,
                              threshold = 0.01){
  
  # Test for enrichment
  go_enrichments <- Genes %>% split(.$Annotation) %>%
    map_dfr(go_enrichment,
            m = length(unique(Genes$gene_id[ Genes$q.value < threshold])),
            n = length(unique(Genes$gene_id[ Genes$q.value >= threshold]))) %>%
    mutate(q.value = p.adjust(p.value, 'fdr'))
  # go_enrichments
  
  
  
  # get GO annotations
  go_annots <- go_enrichments$Annotation %>% map_dfr(function(x, GOTERM){
    term <- GOTERM[[x]]
    
    if(is.null(term)){
      name <- NA
      Ontology <- NA
      Definition <- NA
    }else{
      name <- term@Term
      Ontology <- term@Ontology
      Definition <- term@Definition
    }
    tibble(Annotation = x,
           Term = name,
           Ontology = Ontology,
           Definition = Definition)
  }, GOTERM = GOTERM)
  # go_annots
  
  # Match go annotations
  go_enrichments <- go_enrichments %>%
    left_join(go_annots, by = "Annotation")
  # go_enrichments %>% arrange(desc(OR > 1), q.value)
  # go_enrichments %>% arrange(desc(n.sig))
  
  # Write output
  write_tsv(go_enrichments, outfile)
  
  # Return Genes
  return(Genes)
}

##########################


args <- list(dir.genes = "hmp.vmwa.pvals.genes/",
             dir.annots = "/godot/users/sur/data/genomes/midas_db_v1.2/GO.annotations/",
             threshold = 0.01,
             outdir = "hmp.vmwa.genes.GO/",
             outfile = "hmp.vmwa.genes.GO.all.txt")

args <- list(dir.genes = "../2018-10-19.vmwa_enriched_genes/hmp.vmwa.pvals.genes/",
             dir.annots = "./",
             threshold = 0.01,
             outdir = "hmp.vmwa.genes.GO/",
             outfile = "hmp.vmwa.genes.GO.all.txt")

if(!dir.exists(args$outdir))
  dir.create(args$outdir)

files <- list.files(args$dir.genes)
files

Res <- NULL
for(f in files){
  # f <- files[24]
  name <- str_replace(string = f, pattern = "_vmwa.genes.txt$", replacement = "")
  cat("\t", name, "\n")
  
  # genes_file <- "Porphyromonas_sp_57899_vmwa.genes.txt"
  # go_file <- "Porphyromonas_sp_57899.GO.txt"
  # threshold <- 0.01
  # outfile <- "Porphyromonas_sp_57899_enriched.GO.txt"
  
  # Prepare vars
  genes_file <- paste0(args$dir.genes, "/", f)
  go_file <- paste0(args$dir.annots, "/", name, ".GO.txt")
  outfile <- paste0(args$outdir, "/", name, "_enriched.GO.txt")
  
  # Read and match data
  Genes <- read_tsv(genes_file)
  go <- read_tsv(go_file)
  go <- go %>% select(gene_id = Gene, everything())
  Genes <- Genes %>% left_join(go, by = "gene_id")
  
  # Enrichments
  res <- get_go_enrichment(Genes = Genes, outfile = outfile,
                           GOTERM = GO.db::GOTERM,
                           threshold = args$threshold)
  
  # Store full results
  Res <- Res %>% bind_rows(res)
}

all <- get_go_enrichment(Genes = Res, outfile = args$outfile,
                         GOTERM = GO.db::GOTERM,
                         threshold = args$threshold)

##################

args <- list(dir.genes = "qin2012.vmwa.pvals.genes/",
             dir.annots = "/godot/users/sur/data/genomes/midas_db_v1.2/GO.annotations/",
             threshold = 0.01,
             outdir = "qin2012.vmwa.genes.GO/",
             outfile = "qin2012.vmwa.genes.GO.all.txt")


if(!dir.exists(args$outdir))
  dir.create(args$outdir)

files <- list.files(args$dir.genes)
files

Res <- NULL
for(f in files){
  name <- str_replace(string = f, pattern = "_vmwa.genes.txt$", replacement = "")
  cat("\t", name, "\n")
  
  # genes_file <- "Porphyromonas_sp_57899_vmwa.genes.txt"
  # go_file <- "Porphyromonas_sp_57899.GO.txt"
  # threshold <- 0.01
  # outfile <- "Porphyromonas_sp_57899_enriched.GO.txt"
  
  # Prepare vars
  genes_file <- paste0(args$dir.genes, "/", f)
  go_file <- paste0(args$dir.annots, "/", name, ".GO.txt")
  outfile <- paste0(args$outdir, "/", name, "_enriched.GO.txt")
  
  # Read and match data
  Genes <- read_tsv(genes_file)
  go <- read_tsv(go_file)
  go <- go %>% select(gene_id = Gene, everything())
  Genes <- Genes %>% left_join(go, by = "gene_id")
  
  # Enrichments
  res <- get_go_enrichment(Genes = Genes, outfile = outfile,
                           GOTERM = GO.db::GOTERM,
                           threshold = args$threshold)
  
  # Store full results
  Res <- Res %>% bind_rows(res)
}

all <- get_go_enrichment(Genes = Res, outfile = args$outfile,
                         GOTERM = GO.db::GOTERM,
                         threshold = args$threshold)







