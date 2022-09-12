#!/usr/bin/env Rscript

# (C) Copyright 2022 Sur Herrera Paredes
# This file is part of This program.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with This program.  If not, see <https://www.gnu.org/licenses/>.

# setwd("/cashew/users/sur/exp/fraserv/2022/today3/")



get_gene_variants <- function(x, ref,
                              faa,
                              translate = TRUE,
                              outdir = "output"){
  
  # x <- gene.granges
  # ref <- Genome
  # translate <- TRUE
  # outdir <- args$outdir
  
  
  if(length(x) != 1){
    stop("ERROR: just one genomic range can be used")
  }
  
  
  gene_seq <- seqinr::getFrag(object = ref[ GenomicRanges::seqnames(x) %>% as.character() ],
                              begin = GenomicRanges::start(x),
                              end = GenomicRanges::end(x))
  gene_seq <- seqinr::getSequence(gene_seq)[[1]]
  gene_var <- gene_seq
  
  
  # Replace snps
  snp_pos <- sites_info$ref_pos - GenomicRanges::start(x) + 1
  gene_var[snp_pos] <- sites_info$minor_allele
  
  # Get nt sequence
  strand <- GenomicRanges::strand(x) %>% as.character()
  if(strand == "+"){
    gene_var <- seqinr::as.SeqFastadna(gene_var, name = x$ID, Annot = "ref_alleles")
    gene_seq <- seqinr::as.SeqFastadna(gene_seq, name = x$ID, Annot = "minor_alleles")
  }else if(strand == "-"){
    gene_var <- seqinr::comp(rev(gene_var))
    gene_var <- seqinr::as.SeqFastadna(gene_var, name = x$ID, Annot = "minor_alleles")
    
    gene_seq <- seqinr::comp(rev(gene_seq))
    gene_seq <- seqinr::as.SeqFastadna(gene_seq, name = x$ID, Annot = "minor_alleles")
  }else{
    stop("ERROR")
  }
  
  
  
  gene_nt <- list(gene_seq, gene_var)
  filename <- file.path(outdir, paste0(x$ID, ".fna"))
  seqinr::write.fasta(gene_nt,
                      names = paste(x$ID, c("reference", "minor"), sep = "\t"),
                      file.out = filename)
  
  
  
  if(translate){
    # The strand is already take care of
    # strand <- as.character(GenomicRanges::strand(x))
    # strand <- ifelse(strand == "-", "R",
    #                  ifelse(strand == "+", "F", NA))
    
    gene_aa <- seqinr::translate(gene_var)
    gene_aa <- seqinr::as.SeqFastadna(gene_aa, name = x$ID, Annot = "minor_alleles")
    if( gene_aa[length(gene_aa)] == '*' ){
      gene_aa <- gene_aa[ 1:(length(gene_aa) - 1) ]
    }
    
    
    gene_faa <- list(faa[[1]], gene_aa)
    filename <- file.path(outdir, paste0(x$ID, ".faa"))
    seqinr::write.fasta(gene_faa,
                        names = paste(x$ID, c("reference", "minor"), sep = "\t"),
                        file.out = filename)
  }
  
  NULL
}


library(tidyverse)
library(HMVAR)


# args <- list(midas_dir = "snps_merged/MGYG-HGUT-02438/",
#              pdir = "preHCT_posHCT/MGYG-HGUT-02438.tsv.gz",
#              genome = "uhgg_catalogue/MGYG-HGUT-024/MGYG-HGUT-02438/genome/MGYG-HGUT-02438.fna",
#              proteins = "uhgg_catalogue/MGYG-HGUT-024/MGYG-HGUT-02438/genome/MGYG-HGUT-02438.faa",
#              hits = "gene_hits/MGYG-HGUT-02438.tsv",
#              gff = "uhgg_catalogue/MGYG-HGUT-024/MGYG-HGUT-02438/genome/MGYG-HGUT-02438.gff",
#              outdir = "output/")
args <- list(midas_dir = "snps_merged/",
             pdir = "preHCT_posHCT/",
             genome = "uhgg_catalogue/",
             proteins = "uhgg_catalogue/",
             hits = "gene_hits/",
             gff = "uhgg_catalogue/",
             outdir = "output/")


# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}









Hits <- read_tsv(args$hits,
                 col_types = cols(gene_id = col_character(),
                                  spec = col_character()))
Hits


FAA <- seqinr::read.fasta(args$proteins, seqtype = "AA")
# FAA

Genome <- seqinr::read.fasta(args$genome)

info <- HMVAR::read_midas_info(file.path(args$midas_dir, "snps_info.txt"))
info

Pdir <- read_tsv(args$pdir,
                 col_types = cols(site_id = col_character()))
Pdir

GFF <- plyranges::read_gff(args$gff)




for(gene in Hits$gene_id){
  # gene <- Hits$gene_id[1]
  cat("\t", gene, "\n")
  
  # Get infor of p_directional hits in genes with hits
  sites_info <- info %>%
    filter(gene_id == gene)
  sites_info <- sites_info %>%
    select(site_id, ref_id, ref_pos, ref_allele, major_allele, minor_allele,
           gene_id, amino_acids) %>%
    left_join(Pdir %>%
                select(site_id, p_directional),
              by = "site_id") %>%
    filter(p_directional >= 0.8) %>%
    HMVAR::determine_snp_effect()
  
  # Get gene genomic ranges
  gene.granges <- GFF[ GFF$ID == gene ]
  
  # Get AA
  faa <- FAA[gene]
  
  
  gene_seq <- get_gene_variants(x = gene.granges,
                                ref = Genome,
                                faa = faa,
                                translate = TRUE,
                                outdir = args$outdir)
  
}











