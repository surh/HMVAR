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

#' Expand annotation
#' 
#' Internal function
#'
#' @param gene_id Character
#' @param terms Character string of comma-separated terms for the
#' `gene_id`.
#'
#' @return A tibble with one row per term
#' 
#' @importFrom magrittr %>%
expand_annot <- function(gene_id, terms){
  terms <- stringr::str_split(string = terms, pattern = ",") %>%
    unlist %>%
    stringr::str_replace("^ +", "") %>%
    stringr::str_replace(" +$", "")
  tibble::tibble(gene_id = gene_id,
                 term = terms)
}

#' Convert annotation table to gene or annotation list for topGO
#' 
#' Takes a data frame or tibble that maps genes to annotations,
#' and creates either a list of genes or a list of annotation terms
#' that can be ussed with \link[topGO]{annFUN.gene2GO} or
#' \link[topGO]{annFUN.GO2genes}.
#'
#' @param annots A tibble or data frame that has columns 'gene_id' and
#' 'annots'. Column 'annots' must be a comma-separated charachter string
#' with all the terms that annotate a given 'gene_id'.
#' @param direction Either "geneID2GO" or "GO2geneID". Direction of the
#' list that is produced.
#'
#' @return A named list. If `direction = "geneID2GO"`, then the list has
#' one element per 'gene_id' (named after that gene), and each element
#' of the list is a character vector with all the terms that annotate that
#' gene. If `direction = "GO2geneID"`, the the list has one element per
#' annotation term in 'annots' (named after that term), and each element
#' of the list is a character vector with all the genes annotated with
#' that term.
#' 
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' d <- tibble::tibble(gene_id = c('gene1', 'gene2', 'gene3'),
#'                     terms = c('term3', 'term1,term2', 'term2, term3'))
#' annots_to_geneGO(d, direction = "geneID2GO")
#' annots_to_geneGO(d, direction = "GO2geneID")
annots_to_geneGO <- function(annots, direction = "geneID2GO"){
  if(!all(c("gene_id", "terms") %in% colnames(annots))){
    stop("ERROR: missing columns in annots")
  }
  
  if(direction == "geneID2GO"){
    annots <- annots %>%
      purrr::pmap_dfr(expand_annot) %>%
      dplyr::filter(!is.na(term)) %>%
      split(.$gene_id) %>%
      purrr::map(~ .x$term)
    
  }else if(direction == "GO2geneID"){
    annots <- annots %>%
      purrr::pmap_dfr(expand_annot) %>%
      dplyr::filter(!is.na(term)) %>%
      split(.$term) %>%
      purrr::map(~ .x$gene_id)
  }else{
    stop("ERROR: direction must be geneID2GO or GO2geneID", call. = TRUE)
  }
  
  return(annots)
}