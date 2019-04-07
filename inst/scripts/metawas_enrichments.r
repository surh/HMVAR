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


# Copyright (C) 2019 Sur Herrera Paredes
# Performs functional enrichment on metawas hits from eggnog annotations
# Based on script from /home/sur/micropopgen/exp/2019/2019-01-29.metawas_enrichment/

library(tidyverse)


count_vars <- function(d){
  n.variants <- nrow(d)
  n.sig <- sum(d$type %in% c("int", "both"))
  n.outside <- sum(d$dist > 0)
  tibble(chr = unique(d$chr),
         start = unique(d$start),
         end = unique(d$end),
         n.variants = n.variants,
         n.sig = n.sig,
         n.outside = n.outside)
}

metawas_gene_counts <- function(metawas, closest, annot, outdir = "./", prefix = NULL){
  # Match everything per gene
  all <- metawas %>% full_join(closest, by = c("chr", "ps"))
  
  annot <- annot %>%
    select(gene_id, predicted_gene_name, eggNOG_annot)
  
  
  Genes <- all %>%
    split(.$gene_id) %>%
    map_dfr(~count_vars(.), .id = "gene_id") %>%
    arrange(chr, start) %>%
    left_join(annot, by = "gene_id")
  
  filename <- paste0(c(prefix, "gene_metawas_counts.txt"), collapse = ".")
  filename <- file.path(outdir, filename)
  write_tsv(Genes, filename)
  
  return(filename)
}

#' Sums named vectors by name
#'
#' @param a A named vector.
#' @param b A named vector.
#'
#' @return A named vector.
#' @export
#'
#' @examples
sum_vecs <- function(a, b){
  if(is.null(a)){
    return(b)
  }
  if(is.null(b)){
    return(a)
  }
  
  new <- tibble(term = names(a), count = a) %>%
    full_join(tibble(term = names(b), count = b), by = "term")
  new <- new %>%
    bind_cols(count = new %>% select(-term) %>% rowSums(na.rm = TRUE)) %>%
    select(term, count)
  vec <- new$count
  names(vec) <- new$term
  
  return(vec)
}


# spec <- "Actinomyces_odontolyticus_57475"
# genomes_file <- "hmp_genomes.txt"
# lmm_dir <- "2019-02-19.hmp_metawas/"
# snp_dir <- "metawas_closest/"
# annot_dir <- "annotations/"
# dist_thres <- 500
# count_thres <- 3
# snp_groups <- c("env", "both")
# outdir <- "outptut_env_both"

args <- list(lmmres = "../2019-03-29.hmp_metawas_data/Supragingival.plaque/metawas/lmmpcs/Porphyromonas_sp_57899_lmm.results.txt",
             closest = "../2019-03-29.hmp_metawas_data/Supragingival.plaque/closest/Porphyromonas_sp_57899.closest",
             annotations = "../2019-04-01.hmp_subsite_annotations/hmp.subsite_annotations/Porphyromonas_sp_57899.emapper.annotations",
             dist_thres = 500,
             count_thres = 3,
             outdir = "metawas_enrichments",
             prefix = NULL)

dir.create(args$outdir)

# manhat_dir <- paste0(outdir, "/manhattan")
# dir.create(manhat_dir)
# genes_dir <- paste0(outdir, "/genes/")
# dir.create(genes_dir)
# enrich_dir <- paste0(outdir, "/enrich/")
# dir.create(enrich_dir)
# 
# genomes <- read_tsv(genomes_file, col_names = FALSE, col_types = 'c')
# genomes <- genomes$X1
# 
GENES <- NULL
OG_bg <- NULL
GO_bg <- NULL
KO_bg <- NULL

# Read data
metawas <- read_tsv(args$lmmres,
                    col_types = cols(.default = col_double(),
                                     chr = col_character(),
                                     rs = col_character(),
                                     allele1 = col_character(),
                                     allele0 = col_character(),
                                     type = col_character())) %>%
  select(-starts_with("logl_H1"), -starts_with("l_mle"))
# metawas
closest <- read_tsv(args$closest,
                    col_names = c("chr", "ps", "ps2", "chr2", "start", "end", "gene_id", "dist"),
                    col_types = cols(.default = col_number(),
                                     chr = col_character(),
                                     chr2 = col_character(),
                                     gene_id = col_character())) %>%
  select(-ps2, -chr2)
# closest

# Get genes
genes_tested <- closest %>%
  filter(abs(dist) <= args$dist_thres) %>%
  select(gene_id) %>% unique()
sig_genes <- closest %>%
  left_join(metawas %>%
              filter(type %in% c("int", "both")) %>%
              select(chr, rs, ps),
            by = c("chr", "ps")) %>%
  filter(!is.na(rs)) %>%
  filter(abs(dist) <= args$dist_thres) %>%
  select(gene_id) %>%
  unique()

cat("\tReading annot file...\n")
# annot <- read_tsv(annot_file, comment = "#", col_names = FALSE)
annot <- read_tsv(args$annotations,
                  comment = "#",
                  col_names = c("query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue",
                                "seed_ortholog_score",	"predicted_gene_name", "GO_terms",
                                "KEGG_KOs", "BiGG_reactions", "Annotation_tax_scope", "OGs",
                                "bestOG|evalue|score", "COG_cat", "eggNOG_annot"),
                  col_types = cols(.default = col_character(),
                                   seed_ortholog_score = col_double(),
                                   seed_ortholog_evalue = col_double()))
# annot
# Reformat gene name
annot <- annot %>%
  mutate(gene_id = str_replace(string = query_name,
                               pattern = "\\([+-]\\)_[0-9]",
                               replacement = "")) %>%
  select(gene_id, everything(), -query_name, -seed_ortholog_evalue, -seed_ortholog_score,
         -Annotation_tax_scope)
# annot

# Select onl tested genes
# Only these will be considered in the universe background
annot <- annot %>% filter(gene_id %in% genes_tested$gene_id)

# Metawas gene counts
metawas_gene_counts(metawas = metawas,
                    closest = closest,
                    annot = annot,
                    outdir = args$outdir,
                    prefix = args$prefix)


######
# Create gene-annot table
annotation <- "GO_terms"

expand_annot <- function(gene_id, annot){
  annot <- str_split(string = annot, pattern = ",") %>% unlist
  tibble(gene_id = gene_id,
         term = annot)
}

BG <- annot %>%
  select(gene_id, annot = annotation) %>%
  pmap_dfr(expand_annot) %>%
  filter(!is.na(term)) %>%
  mutate(sig_gene = gene_id %in% sig_genes$gene_id)
BG

filename <- paste0(c(args$prefix, annotation, "BG.txt"), collapse = ".")
filename <- file.path(args$outdir, filename)
write_tsv(BG, filename)

# Test
to_test <- BG %>%
  filter(sig_gene) %>%
  count(term) %>%
  filter(n >= args$count_thres) %>%
  select(term, n.sig = n)

bg_counts <- BG %>%
  filter(term %in% to_test$term) %>%
  count(term) %>%
  select(term, n.bg = n)


test_res <- to_test %>%
  left_join(bg_counts, by = "term") %>%
  pmap_dfr(function(term, n.sig, n.bg, selection.size, bg.size){
    mat <- matrix(c(n.sig, n.bg, selection.size, bg.size), ncol = 2)
    res <- fisher.test(mat)
    tibble(term = term,
           n.sig = n.sig,
           n.bg = n.bg,
           OR = res$estimate,
           p.value = res$p.value)
  }, selection.size = nrow(sig_genes), bg.size = nrow(BG)) %>%
  mutate(q.value = p.adjust(p.value, 'fdr')) %>%
  arrange(q.value)
test_res <- test_res %>%
  bind_cols(test_res %>%
              select(term) %>%
              unlist %>%
              AnnotationDbi::select(GO.db::GO.db, .,
                                    columns = c("ONTOLOGY","TERM","DEFINITION")))
test_res





#####
# Get background
og_bg <- annot$OGs %>% map(str_split, pattern = ",") %>% unlist %>% table
go_bg <- annot$GO_terms %>% map(str_split, pattern = ",") %>% unlist %>% table
ko_bg <- annot$KEGG_KOs %>% map(str_split, pattern = ",") %>% unlist %>% table

OG_bg <- sum_vecs(OG_bg, og_bg)
GO_bg <- sum_vecs(GO_bg, go_bg)
KO_bg <- sum_vecs(KO_bg, ko_bg)

    
    
    
    # Match genes and annotations
    genes <- sig_genes %>% left_join(annot)
    genes
    # filename <- paste0(genes_dir, "/", spec, "_genes.annot")
    # write_tsv(genes, filename)
    GENES <- GENES %>% bind_rows(genes)
    
    # Calculate enrichment
    # ogs <- genes$OGs
    # ogs <- ogs[ !is.na(ogs) ]
    # ogs <- ogs %>% map(str_split, pattern = ',') %>% unlist %>% table
    gos <- genes$GO_terms
    gos <- gos[ !is.na(gos) ]
    gos <- gos %>% map(str_split, pattern = ',') %>% unlist %>% table
    gos <- gos[ gos >= args$count_thres ]
    kos <- genes$KEGG_KOs
    kos <- kos[ !is.na(kos) ]
    kos <- kos %>% map(str_split, pattern = ',') %>% unlist %>% table
    kos <- kos[ kos >= args$count_thres ]
    
    if(length(gos) > 0){
      gos <- tibble(term = names(gos), n = gos)
      gos_res <- gos %>% pmap_dfr(function(term,n,bg = go_bg, total = sum(gos$n)){
        term_bg <- bg[term]
        mat <- matrix(c(n, total, term_bg, sum(bg)), ncol = 2)
        res <- fisher.test(mat)
        p.value <- res$p.value
        OR <- res$estimate
        tibble(term = term, n = n, OR = OR, p.value = p.value)
      }) %>% mutate(q.value = p.adjust(p.value, 'fdr')) %>%
        arrange(q.value)
      
      gos_res <- gos_res %>%
        bind_cols(gos_res %>%
                    select(term) %>%
                    unlist %>%
                    AnnotationDbi::select(GO.db::GO.db, .,
                                          columns = c("ONTOLOGY","TERM","DEFINITION")))
      
      filename <- paste0(c(args$prefix, "go_enrichments.txt"), collapse = ".")
      filename <- file.path(args$outdir, filename)
      write_tsv(gos_res, filename)
    }
    
    if(length(kos) > 0){
      kos <- tibble(term = names(kos), n = kos)
      kos_res <- kos %>% pmap_dfr(function(term,n,bg = ko_bg, total = sum(kos$n)){
        term_bg <- bg[term]
        mat <- matrix(c(n, total, term_bg, sum(bg)), ncol = 2)
        res <- fisher.test(mat)
        p.value <- res$p.value
        OR <- res$estimate
        tibble(term = term, n = n, OR = OR, p.value = p.value)
      }) %>% mutate(q.value = p.adjust(p.value, 'fdr')) %>%
        arrange(q.value)
      filename <- paste0(enrich_dir, "/", spec, "_ko.enrich.txt")
      write_tsv(kos_res, filename)
    }
  }
  
  rm(genes,kos, gos, annot, metawas)
}



rm(og_bg, go_bg, ko_bg)



ogs <- GENES$OGs
ogs <- ogs[ !is.na(ogs) ]
ogs <- ogs %>% map(str_split, pattern = ',') %>% unlist %>% table
ogs <- ogs[ ogs >= count_thres ]
gos <- GENES$GO_terms
gos <- gos[ !is.na(gos) ]
gos <- gos %>% map(str_split, pattern = ',') %>% unlist %>% table
gos <- gos[ gos >= count_thres ]
kos <- GENES$KEGG_KOs
kos <- kos[ !is.na(kos) ]
kos <- kos %>% map(str_split, pattern = ',') %>% unlist %>% table
kos <- kos[ kos >= count_thres ]



if(length(gos) > 0){
  gos <- tibble(term = names(gos), n = gos)
  gos_res <- gos %>% pmap_dfr(function(term,n,bg = GO_bg, total = sum(gos$n)){
    term_bg <- bg[term]
    mat <- matrix(c(n, total, term_bg, sum(bg)), ncol = 2)
    res <- fisher.test(mat)
    p.value <- res$p.value
    OR <- res$estimate
    tibble(term = term, n = n, OR = OR, p.value = p.value)
  }) %>% mutate(q.value = p.adjust(p.value, 'fdr')) %>%
    arrange(q.value)
  filename <- paste0(outdir, "/overall_go.enrich.txt")
  write_tsv(gos_res, filename)
}


if(length(kos) > 0){
  kos <- tibble(term = names(kos), n = kos)
  kos_res <- kos %>% pmap_dfr(function(term,n,bg = KO_bg, total = sum(kos$n)){
    term_bg <- bg[term]
    mat <- matrix(c(n, total, term_bg, sum(bg)), ncol = 2)
    res <- fisher.test(mat)
    p.value <- res$p.value
    OR <- res$estimate
    tibble(term = term, n = n, OR = OR, p.value = p.value)
  }) %>% mutate(q.value = p.adjust(p.value, 'fdr')) %>%
    arrange(q.value)
  filename <- paste0(outdir, "/overall_ko.enrich.txt")
  write_tsv(kos_res, filename)
}


if(length(ogs) > 0){
  ogs <- tibble(term = names(ogs), n = ogs)
  ogs_res <- ogs %>% pmap_dfr(function(term,n,bg = OG_bg, total = sum(ogs$n)){
    term_bg <- bg[term]
    mat <- matrix(c(n, total, term_bg, sum(bg)), ncol = 2)
    res <- fisher.test(mat)
    p.value <- res$p.value
    OR <- res$estimate
    tibble(term = term, n = n, OR = OR, p.value = p.value)
  }) %>% mutate(q.value = p.adjust(p.value, 'fdr')) %>%
    arrange(q.value)
  filename <- paste0(outdir, "/overall_og.enrich.txt")
  write_tsv(ogs_res, filename)
}
