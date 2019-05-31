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
