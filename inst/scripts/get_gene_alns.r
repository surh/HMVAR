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

library(tidyverse)
library(seqinr)
library(HMVAR)

#' Checks loaded midas data
#' 
#' Internal
#' 
#' @param Dat A MIDAS data
#'
#' @return TRUE if all tests pass
check_midas_data <- function(Dat){
  
  # Check name of elements
  if(!all(c('info', 'freq', 'depth') %in% names(Dat))){
    stop("ERROR: The midas_data object must have elements c('info', 'freq', 'depth')", call. = TRUE)
  }
  
  # Check element types
  if(!tibble::is_tibble(Dat$info))
    stop("ERROR: info must be a tibble.", call. = TRUE)
  if(!tibble::is_tibble(Dat$freq))
    stop("ERROR: freq must be a tibble.", call. = TRUE)
  if(!tibble::is_tibble(Dat$depth))
    stop("ERROR: depth must be a tibble.", call. = TRUE)
  
  # Check dimensions of freq and depth
  if(any(dim(Dat$freq) != dim(Dat$depth)))
    stop("ERROR: freq and depth must have matching dimmensions.", call. = TRUE)
  
  # Check names of samples and site ids
  if(any(colnames(Dat$freq) != colnames(Dat$depth)))
    stop("ERROR: column names of depth and freq don't match")
  if(any(Dat$freq$site_id != Dat$depth$site_id))
    stop("ERROR: site_id differs between freq and depth.", call. = TRUE)
  if(any(Dat$info$site_id != Dat$freq$site_id))
    stop("ERROR: site_id differens between info and freq")
  
  return(TRUE)
}

#' Assign alleles bases on allele frequency and depth
#' 
#' Internal
#'
#' @param dat A data.frame or tibble
#' @param depth_thres depth_threshold
#' @param freq_thres freq_threshold
#' @param sequence Should the actual nucleodites be returned? If FALSE
#' the allele column will contain 'major' or 'minor'. If TRUE, the
#' allele column will contain the actual nucleotide
#' @param na_rm Should sites from samples that could not be assigned be returned?
#'
#' @return Same tibble as dat but with an allele column
#' 
#' @importFrom magrittr %>%
assign_allele <- function(dat, depth_thres = 1, freq_thres = 0.5, sequence = FALSE, na_rm = TRUE){
  if(!is.data.frame(dat))
    stop("ERROR: dat must be a data.frame.", call. = TRUE)
  if(!all(c('depth', 'freq', 'major_allele', 'minor_allele') %in% colnames(dat)))
    stop("ERROR: dat must have columns c('depth', 'freq', 'major_allele', 'minor_allele').", call. = TRUE)
  
  # Assign allele
  dat <- dat %>%
    dplyr::filter(depth >= depth_thres) %>%
    dplyr::filter(freq != 0.5) %>%
    dplyr::mutate(allele = replace(freq, freq <= freq_thres, 'major')) %>%
    dplyr::mutate(allele = replace(allele, freq >= (1 - freq_thres), 'minor')) %>%
    dplyr::mutate(allele = replace(allele, (freq > freq_thres) & (freq < (1 - freq_thres)), NA))
  
  if(na_rm){
    dat <- dat %>%
      dplyr::filter(!is.na(allele))
  }
  
  if(sequence){
    dat$allele[ dat$allele == 'major' ] <- dat$major_allele[ dat$allele == 'major' ]
    dat$allele[ dat$allele == 'minor' ] <- dat$minor_allele[ dat$allele == 'minor' ]
  }
  
  return(dat)
}


#' Get all MIDAS data from a gene
#' 
#' Returns a list similar to the one by \link{read_midas_data} but
#' only with sites from a given gene.
#'
#' @param Dat A list. See output from \link{read_midas_data}
#' @param gene Name of the gene to be extracted, must exist in
#' 'gene_id' column of the 'info' element from Dat.
#' @param depth_thres Minimum depth trheshold.
#' @param freq_thres Freq threshold.
#'
#' @return A tibble that results of merging the info, dat and freq elements
#' of Dat. Plus allele assignments.
#' 
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
gene_midas_data <- function(Dat, gene, depth_thres = 1, freq_thres = 0.5){
  # gene <- genes[1]
  check_midas_data(Dat)
  
  # Get allele_freq and data from genes
  i <- Dat$info %>%
    dplyr::filter(gene_id == gene)
  f <- Dat$freq %>%
    dplyr::filter(site_id %in% i$site_id) %>%
    arrange(match(site_id, i$site_id))
  d <- Dat$depth %>%
    dplyr::filter(site_id %in% i$site_id) %>%
    arrange(match(site_id, i$site_id))
  all_sites <- match_freq_and_depth(freq = f, depth = d, info = i,
                                    map = NULL, depth_thres = depth_thres)
  all_sites <- all_sites %>%
    dplyr::select(-snp_type, -site_type, -amino_acids, -gene_id, - ref_id)
  
  
  all_sites <- assign_allele(all_sites,
                             depth_thres = depth_thres,
                             freq_thres = freq_thres,
                             sequence = TRUE, na_rm = TRUE)
  
  return(all_sites)
}

#' Get sequence fragment
#' 
#' Internal. Simple wrapper for \link[seqinr]{getFrag}
#'
#' @param seq A sequence fasta object from \link{seqinr}
#' @param ref_id The sequence ID from the fragment,
#' @param start 1-indexed start position.
#' @param end 1-indexed end position (included in fragment).
#'
#' @return A character vector
get_sequence_fragment <- function(seq, ref_id, start, end){
  gene_seq <- seqinr::getFrag(seq[[ ref_id ]], begin = start, end = end)
  gene_seq <- toupper(seqinr::getSequence(gene_seq))
  
  return(gene_seq)
}


#' Get gene SNV alignment
#' 
#' Takes SNP informmation from multiple samples for a sequence
#' fragment, and a reference sequence for that fragment. It produces
#' and 'alignment' object that has one 'aligned' sequnece per sample
#' where the SNPs have been replaced.
#'
#' @param snps A tibble or data.frame. It must have columns
#' 'site_id', 'sample', 'ref_pos', 'ref_allele', 'major_allele',
#' and 'allele'. Column 'minor_allele' is not required at the moment
#' but might be in the future.
#' @param seq A character vector with the sequence of the fragment.
#' @param missing_as How to treat sites in samples when there is no
#' information in snps? Either 'N'(set all those sites to N), 'major'
#' (use the major allele for that site), or 'ref' (use the reference allele). 
#' @param strand If strand is '-', then the sequence will be reverse complemented.
#' @param keep_last_codon This is useful to remove stop codons when gene is a
#' CDS. It will trim the last 3 sites from the alignment.
#'
#' @return An algnment object from \link{seqinr}
#' 
#' @export
#' @importFrom magrittr %>%
#' 
#' 
gene_snv_aln <- function(snps, seq, missing_as = 'major', strand = '+', keep_last_codon = TRUE){
  # Check params
  if(missing(snps) || missing(seq))
    stop("ERROR: snps and seq must be provided", call. = TRUE)
  if(!(strand %in% c('+', '-')))
    stop("ERROR: strand must be '+' or '-'.", call. = TRUE)
  if(!is.character(seq))
    stop("ERROR: seq must be a character vector.", call. = TRUE)
  
  # Get vector of 'reference' alleles. ie. values to use when no data.
  if(missing_as == 'major'){
    alleles <- snps %>%
      split(.$ref_pos) %>%
      purrr:::map_chr(~.x$major_allele[1])
  }else if(missing_as == 'ref'){
    alleles <- snps %>%
      split(.$ref_pos) %>%
      purrr:::map_chr(~.x$ref_allele[1])
  }else if(missing_as == 'N'){
    alleles <- snps %>%
      split(.$ref_pos) %>%
      purrr:::map_chr(~ 'N')
  }else{
    stop("ERROR: missing_as must be 'major', 'ref', or 'N'.", call. = TRUE)
  }
  
  # 
  # dat <- snps %>%
  #   split(.$sample)
  # x <- dat[[1]]
  
  aln <- snps %>%
    split(.$sample) %>%
    purrr::map(function(x, seq, alleles){
      # Find sites with missing data
      x <- tibble::tibble(ref_pos = as.numeric(names(alleles))) %>%
        left_join(x, by = 'ref_pos')
      
      # Prepare vector with sites with data
      changes <- setNames(x$allele, x$ref_pos)
      changes[ is.na(changes) ] <- alleles[ names(changes)[ is.na(changes) ] ]
      
      # # Define changes for missing data
      # if(missing_as == 'N'){
      #   changes[ is.na(changes) ] <- 'N'
      # }else if(missing_as == 'ref'){
      #   changes <- changes[ !is.na(changes) ]
      # }else{
      #   stop("ERROR: missing_as must be 'N' or 'ref'.", call. = TRUE)
      # }
      
      # Replace SNVs in sequence
      new_gene <- seq
      new_gene[ as.numeric(names(changes)) ] <- changes
      
      return(new_gene)
    },
    seq = seq, 
    alleles = alleles)
  
  # Reverse_complement
  if(gene_pos$strand == '-'){
    aln <- aln %>%
      purrr::map(~rev(chartr('TGCAtgca', 'ACGTacgt', .x)))
  }
  
  # Remove last
  if(!keep_last_codon){
    aln <- aln %>%
      purrr::map(~.x[1:(length(.x) - 3)])
  }
  
  # Convert to alignment
  aln <- seqinr::as.alignment(nb = length(aln), nam = names(aln), seq = aln)
  
  return(aln)
}


setwd("~/micropopgen/exp/2019/today2/")

depth_thres <- 1
freq_thres <- 0.5
missing_as <- 'major'
keep_last_codon <- TRUE


map <- read_tsv('map.txt') %>%
  select(sample = ID, Group = Group)
map
Dat <- read_midas_data(midas_dir = "midas_output_small/", map = map)


genome_fasta <- seqinr::read.fasta('genome/genome.fna.gz')
genome_feats <- readr::read_tsv('genome/genome.features.gz',
                                col_types = readr::cols(.default = readr::col_character(),
                                                        start = readr::col_integer(),
                                                        end = readr::col_integer()))

genes <- unique(Dat$info$gene_id)
# gene <- genes[1]
# gene

for(gene in genes){
  # Get ref sequence
  gene_pos <- genome_feats %>%
    filter(gene_id == gene) %>%
    select(start, scaffold_id, end, strand)
  if(nrow(gene_pos) != 1){
    stop("ERROR: no data on genes")
  }
  gene_seq <- get_sequence_fragment(genome_fasta,
                                    ref_id = gene_pos$scaffold_id,
                                    start = gene_pos$start,
                                    end = gene_pos$end)
  
  # Get snp data
  all_sites <- gene_midas_data(Dat = Dat, gene = gene,
                               depth_thres = depth_thres,
                               freq_thres = freq_thres)
  # Set relative positions
  all_sites <- all_sites %>%
    mutate(ref_pos = ref_pos - gene_pos$start + 1)
  
  aln <- gene_snv_aln(snps = all_sites, seq = gene_seq,
                      missing_as = missing_as,
                      strand = gene_pos$strand,
                      keep_last_codon = keep_last_codon)
  
  
  filename <- paste0(gene, '.aln.fasta')
  seqinr::write.fasta(aln$seq, aln$nam, file.out = filename)
}