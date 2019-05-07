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
#'                     terms = c(NA, 'term1,term2', 'term2, term3'))
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



#' gene selection function
#' 
#' Internal
#'
#' @param thres score threshold to select function
#'
#' @return A function
gene_sel_fun <- function(thres){
  function(x) x < thres
}


#' Gene Ontology enrichment via topGO
#' 
#' Performs Gene Ontology (GO) enrichment analysis via
#' topGO
#'
#' @param genes Either a character vector with the gene
#' identifiers that are 'significant' or a named numeric vector
#' where the vector names are the gene_identifiers of the universe
#' of genes and the numeric values the genes' scores.
#' @param annots Either a data.frame or tibble that has 'gene_id' and
#' 'terms' column or the result of running \link{annots_to_geneGO} on
#' such table.
#' @param ontology Which ontology to test. See help at \link[topGO]{topGOdata-class}.
#' @param description Description for the test. See help at \link[topGO]{topGOdata-class}.
#' @param algorithm Algortithm for test. See help at \link[topGO]{runTest}.
#' @param statistic Statistic to test. See help at \link[topGO]{runTest}.
#' @param node_size Minimum number of genes per term to test that
#' term. See help at \link[topGO]{topGOdata-class}
#' @param score_threshold If genes is a numeric vector, then
#' this should be the 'significance' threshold. E.g. if scores are p-values
#' a common threshold would be 0.05.
#'
#' @return A list with elements topgo_data and topgo_res of class
#' topGOdata and topGOresult respecitveley
#' 
#' @export
#' @importClassesFrom topGO topGOdata
test_go <- function(genes, annots,
                    ontology,
                    description, algorithm, statistic,
                    node_size, score_threshold) UseMethod("test_go")


#' @rdname test_go
#' @method test_go character
#' @export
test_go.character <- function(genes, annots,
                              ontology = "BP",
                              description = '',
                              algorithm = 'classic',
                              statistic = 'fisher',
                              node_size = 3){
  
  # Get annotations as gene -> GO list
  if(is.data.frame(annots)){
    annots <- annots_to_geneGO(annots = annots, direction = "geneID2GO")
  }else if(!is.list(annots)){
    stop("ERROR: annots must be either a data.frame or a list.", call. = TRUE)
  }
  
  # Convert list of significant genes into scores
  gene_scores <- -1*(names(annots) %in% genes)
  names(gene_scores) <- names(annots)
  
  res <- test_go(genes = gene_scores, annots = annots,
                 ontology = ontology, description = description,
                 algorithm = algorithm, statistic = statistic, node_size = node_size,
                 score_threshold = 0)
  
  return(res)
  
}

#' @rdname test_go
#' @method test_go numeric
#' @export
test_go.numeric <- function(genes, annots,
                            ontology = "BP",
                            description = '',
                            algorithm = 'classic',
                            statistic = 'fisher',
                            node_size = 3,
                            score_threshold = 0.05){
  
  # Get annotations as gene -> GO list
  if(is.data.frame(annots)){
    annots <- annots_to_geneGO(annots = annots, direction = "geneID2GO")
  }else if(!is.list(annots)){
    stop("ERROR: annots must be either a data.frame or a list.", call. = TRUE)
  }
  
  topGO::groupGOTerms()
  
  # Create topGO data
  go_data <- new("topGOdata",
                 description = description,
                 ontology = ontology,
                 allGenes = genes,
                 geneSelectionFun = gene_sel_fun(score_threshold),
                 nodeSize = node_size,
                 annot = topGO::annFUN.gene2GO,
                 gene2GO = annots)
  
  # perform topGO test
  go_res <- topGO::runTest(go_data,
                           algorithm = algorithm,
                           statistic = statistic)
  
  return(list(topgo_data = go_data, topgo_res = go_res))
}

#' Title
#'
#' @param genes 
#' @param scores 
#' @param test 
#' @param alternative 
#' @param min_size 
#'
#' @return
#' @export
#'
#' @examples
term_gsea <- function(genes, scores, test = "wilcoxon", alternative = "greater", min_size = 3){
  
  if(length(genes) < min_size){
    return(NULL)
  }
  
  ii <- names(scores) %in% genes
  
  if(test == 'wilcoxon'){
    res <- wilcox.test(scores[ii], scores[!ii], alternative = alternative)
  }else if(test == 'ks'){
    res <- ks.test(scores[ii], scores[!ii], alternative = alternative)
  }else{
    stop("ERROR: Invalid test", call. = TRUE)
  }
  
  tibble::tibble(size = length(genes), statistic = res$statistic, p.value = res$p.value)
}

#' Title
#'
#' @param dat 
#' @param test 
#' @param alternative 
#' @param min_size 
#'
#' @return
#' @export
#'
#' @examples
gsea <- function(dat, test = 'wilcoxon', alternative = 'greater', min_size = 3){
  if(!all(c('gene_id', 'terms', 'score') %in% colnames(dat))){
    stop("ERROR: missing columns", call. = TRUE)
  }
  
  scores <- dat$score
  names(scores) <- dat$gene_id
  dat <- dat %>% dplyr::select(-score)
  
  # Expand annotations
  dat <- dat %>%
    purrr::pmap_dfr(expand_annot) %>%
    dplyr::filter(!is.na(term))
  
  # Clean background
  scores <- scores[ names(scores) %in% unique(dat$gene_id) ]
  
  # Test
  dat <- dat %>%
    split(.$term) %>%
    purrr::map(~ .x$gene_id) %>%
    purrr::map_dfr(term_gsea, scores = scores,
                   test = test,
                   alternative = alternative,
                   min_size = min_size,
                   .id = "term") %>%
    dplyr::arrange(p.value)
  
  return(dat)
}

