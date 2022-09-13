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
                              sites_info,
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
  gene_var[snp_pos] <- sites_info$alt_allele
  
  # Get nt sequence
  strand <- GenomicRanges::strand(x) %>% as.character()
  if(strand == "+"){
    gene_var <- seqinr::as.SeqFastadna(gene_var, name = x$ID, Annot = "ref_alleles")
    gene_seq <- seqinr::as.SeqFastadna(gene_seq, name = x$ID, Annot = "alt_alleles")
  }else if(strand == "-"){
    gene_var <- seqinr::comp(rev(gene_var))
    gene_var <- seqinr::as.SeqFastadna(gene_var, name = x$ID, Annot = "alt_alleles")
    
    gene_seq <- seqinr::comp(rev(gene_seq))
    gene_seq <- seqinr::as.SeqFastadna(gene_seq, name = x$ID, Annot = "ref_alleles")
  }else{
    stop("ERROR")
  }
  
  
  
  gene_nt <- list(gene_seq, gene_var)
  filename <- file.path(outdir, paste0(x$ID, ".fna"))
  cat("\tWriting ", filename, "\n")
  seqinr::write.fasta(gene_nt,
                      names = paste(x$ID, c("reference", "alternative"), sep = "\t"),
                      file.out = filename)
  
  
  
  if(translate){
    # The strand is already take care of
    # strand <- as.character(GenomicRanges::strand(x))
    # strand <- ifelse(strand == "-", "R",
    #                  ifelse(strand == "+", "F", NA))
    
    gene_aa <- seqinr::translate(gene_var)
    gene_aa <- seqinr::as.SeqFastadna(gene_aa, name = x$ID, Annot = "alt_alleles")
    if( gene_aa[length(gene_aa)] == '*' ){
      gene_aa <- gene_aa[ 1:(length(gene_aa) - 1) ]
    }
    
    
    gene_faa <- list(faa[[1]], gene_aa)
    filename <- file.path(outdir, paste0(x$ID, ".faa"))
    cat("\tWriting ", filename, "\n")
    seqinr::write.fasta(gene_faa,
                        names = paste(x$ID, c("reference", "alternative"), sep = "\t"),
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
             uhgg = "uhgg_catalogue/",
             hits = "gene_hits/",
             outdir = "output/")

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}


list.files(args$hits, full.names = TRUE, recursive = FALSE) %>%
  map(function(hits, uhgg_catalogue, midas_dir, pdir, outdir = "output"){
    # hits <- list.files(args$hits, full.names = TRUE, recursive = FALSE)[1]
    # uhgg_catalogue <- args$uhgg
    # midas_dir <- args$midas_dir
    # pdir <- args$pdir
    
    spec <- basename(hits) %>% str_remove("[.]tsv$")
    cat(spec, "\n")
    
    specdir <- file.path(outdir, spec)
    dir.create(specdir)
    
    # Read data
    Hits <- read_tsv(hits,
                     col_types = cols(gene_id = col_character(),
                                      spec = col_character()))
    
    filename <- file.path(uhgg_catalogue, 
                          str_sub(spec, 1, 13),
                          spec,
                          "genome",
                          paste0(spec, ".faa"))
    FAA <- seqinr::read.fasta(filename, seqtype = "AA")
    
    filename <- file.path(uhgg_catalogue, 
                          str_sub(spec, 1, 13),
                          spec,
                          "genome",
                          paste0(spec, ".fna"))
    Genome <- seqinr::read.fasta(filename)
    
    filename <- file.path(uhgg_catalogue, 
                          str_sub(spec, 1, 13),
                          spec,
                          "genome",
                          paste0(spec, ".gff"))
    GFF <- plyranges::read_gff(filename)
    
    info <- HMVAR::read_midas_info(file.path(midas_dir,
                                             spec,
                                             "snps_info.txt"))
    
    
    Pdir <- read_tsv(file.path(pdir, paste0(spec, ".tsv.gz")),
                     col_types = cols(site_id = col_character()))
  
    
    
    for(gene in Hits$gene_id){
      # gene <- Hits$gene_id[1]
      cat("\t", gene, "\n")
      
      # Get info of p_directional hits in genes with hits
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
      # Get alternative allele
      sites_info <- sites_info %>%
        mutate(alt_allele = ifelse(ref_allele == major_allele, minor_allele, major_allele))
      
      # Get gene genomic ranges
      gene.granges <- GFF[ GFF$ID == gene ]
      
      # Get AA
      faa <- FAA[gene]
      
      
      gene_seq <- get_gene_variants(x = gene.granges,
                                    ref = Genome,
                                    faa = faa,
                                    sites_info = sites_info,
                                    translate = TRUE,
                                    outdir = specdir)
      
    }
    
    spec
  }, uhgg_catalogue = args$uhgg, midas_dir = args$midas_dir, pdir = args$pdir,
  outdir = args$outdir)























