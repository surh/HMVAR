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
#' @param ... Other arguments to specific methods
#'
#' @return A list with elements topgo_data and topgo_res of class
#' topGOdata and topGOresult respecitveley. If the topGOdata object
#' cannot be created, it will return a list with two NULL entries.
#' If the topGOdata can be created but the test cannot be performed.
#' The topgo_data will contain the topGOdata object and the topgo_res
#' element will be NULL. It will issue warnings anytime that the
#' topgo_res element is NULL
#' 
#' @export
#' @importClassesFrom topGO topGOdata
test_go <- function(genes, annots,
                    ontology,
                    description, algorithm, statistic,
                    node_size, ...) UseMethod("test_go")


#' @rdname test_go
#' @method test_go character
#' @export
test_go.character <- function(genes, annots,
                              ontology = "BP",
                              description = '',
                              algorithm = 'classic',
                              statistic = 'fisher',
                              node_size = 3, ...){
  
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
#' 
#' @param score_threshold If genes is a numeric vector, then
#' this should be the 'significance' threshold. E.g. if scores are p-values
#' a common threshold would be 0.05.
test_go.numeric <- function(genes, annots,
                            ontology = "BP",
                            description = '',
                            algorithm = 'classic',
                            statistic = 'fisher',
                            node_size = 3,
                            score_threshold = 0.05, ...){
  
  # Get annotations as gene -> GO list
  if(is.data.frame(annots)){
    annots <- annots_to_geneGO(annots = annots, direction = "geneID2GO")
  }else if(!is.list(annots)){
    stop("ERROR: annots must be either a data.frame or a list.", call. = TRUE)
  }

  topGO::groupGOTerms()
  
  go_data <- tryCatch(new("topGOdata",
                          description = description,
                          ontology = ontology,
                          allGenes = genes,
                          geneSelectionFun = gene_sel_fun(score_threshold),
                          nodeSize = node_size,
                          annot = topGO::annFUN.gene2GO,
                          gene2GO = annots),
                      error = function(e) NULL)
  
  if(class(go_data) == "topGOdata"){
    # perform topGO test
    go_res <- tryCatch(topGO::runTest(go_data,
                                      algorithm = algorithm,
                                      statistic = statistic),
                       error = function(e){
                         warning("topGO::runTest failed")
                         return(NULL)})
  }else if(class(go_data) == "NULL"){
    warning("topGOdata object couldn't be created")
    go_res <- NULL
  }else{
    stop("ERROR: unexpected error", call. = TRUE)
  }
  
  return(list(topgo_data = go_data, topgo_res = go_res))
}

#' Gene-set Enrichment analysis on one term.
#' 
#' Performs gene-set enrichment analysis on a group of genes.
#' Tests whether the scores among selected genes differ from
#' the overall score dsitribution.
#' 
#' The function currently doesn't check whether some genes in `genes`
#' are missing from `scores`. It will simply ignore those and test among
#' the `genes` found in scores.
#'
#' @param genes Character vector of gene IDs that belong to
#' a group to test.
#' @param scores Named numeric vector of gene scores to be tested.
#' The 'names' attribute must correspond to the values in `genes`.
#' @param test Which test to perform. Either 'wilcoxon' or 'ks'
#' for Wilcoxon Rank Sum and Kolmogorov-Smirnov tests repsectiveley.
#' Test use R's base \link{wilcox.test} and \link{ks.test} respectiveley.
#' @param alternative The alternative hypothesis to test. Either 'greater',
#' 'less' or 'two.sided'. It  corresponds to option 'alternative' in
#' \link{wilcox.test} or \link{ks.test}. Typically, if scores are p-values
#' one wishes to #' test the hypothesis that p-values within 'genes' are 'less'
#' than expected; while if scores are some other type of value (like
#' fold-change abundance) one is trying to test that those values are
#' 'greater'. Keep in mind that the Kolmogorov-Smirnov test is a test of the
#' maximum difference in cumulative distribution values. Therefore, an
#' alternative 'greater' in this case correspons to cases where
#' score is stochastially smaller than the rest. 
#' @param min_size The minimum number of genes in the group for the test
#' to be performed. Basically if the number of genes that appear in
#' 'scores' is less than 'min_size', the test won't be performed.
#'
#' @return If the test is not performed it returns NULL. If the test
#' is performed it returns a tibble with elements: size (the number 
#' of elements in both 'genes' and 'scores'), statistic (the statistic
#' calculated, depends on the test), and p.value (the p-value of the test).
#' 
#' @export
#'
#' @examples
# Create some fake scores
#' set.seed(123)
#' scores <- rnorm(n = 100)
#' names(scores) <- paste('gene', 1:100, sep = "")
#' 
#' # Select some genes and increase their scores
#' genes <- names(scores)[1:10]
#' scores[genes] <- scores[genes] + rnorm(10, mean = 1)
#' 
#' # Test
#' term_gsea(genes, scores)
#' term_gsea(genes, scores, test = 'ks', alternative = 'less')
term_gsea <- function(genes, scores, test = "wilcoxon",
                      alternative = "greater", min_size = 3){
  
  if(!is.character(genes)){
    stop("ERROR: genes must be a character vector", call. = TRUE)
  }
  if(!is.numeric(scores) || is.null(attr(scores, "names"))){
    stop("ERROR: scores must be a named numeric vector", call. = TRUE)
  }
  
  ii <- names(scores) %in% genes
  
  if(sum(ii) < min_size){
    return(tibble::tibble(size = sum(ii), statistic = NA, p.value = NA))
  }
  if(sum(!ii) < 1){
    return(tibble::tibble(size = sum(ii), statistic = NA, p.value = NA))
  }
  
  if(test == 'wilcoxon'){
    res <- wilcox.test(scores[ii], scores[!ii], alternative = alternative)
  }else if(test == 'ks'){
    res <- ks.test(scores[ii], scores[!ii], alternative = alternative)
  }else{
    stop("ERROR: Invalid test", call. = TRUE)
  }
  
  tibble::tibble(size = sum(ii), statistic = res$statistic, p.value = res$p.value)
}

#' Gene-set enrichment analysis
#' 
#' Performs Gene-set enrichment analysis on all annotation terms
#' for a set of genes.
#'
#' @param dat A data.frame or tibble. It must contain one row per gene
#' and columns 'gene_id', 'terms', and 'score'. Column 'terms' must be
#' of type character and each entry must be a comma-separated character
#' string of all the terms that annotate the corresponding gene.
#' @param test Which test to perform. Either 'wilcoxon' or 'ks'
#' for Wilcoxon Rank Sum and Kolmogorov-Smirnov tests repsectiveley.
#' Test use R's base \link{wilcox.test} and \link{ks.test} respectiveley.
#' @param alternative The alternative hypothesis to test. Either 'greater',
#' 'less' or 'two.sided'. It  corresponds to option 'alternative' in
#' \link{wilcox.test} or \link{ks.test}. Typically, if scores are p-values
#' one wishes to #' test the hypothesis that p-values within 'genes' are 'less'
#' than expected; while if scores are some other type of value (like
#' fold-change abundance) one is trying to test that those values are
#' 'greater'. Keep in mind that the Kolmogorov-Smirnov test is a test of the
#' maximum difference in cumulative distribution values. Therefore, an
#' alternative 'greater' in this case correspons to cases where
#' score is stochastially smaller than the rest. 
#' @param min_size The minimum number of genes in the group for the test
#' to be performed. Basically if the number of genes that appear in
#' 'scores' is less than 'min_size', the test won't be performed.
#'
#' @return A tibble with elements: term (the annotation term ID), size (the number 
#' of elements in both 'genes' and 'scores'), statistic (the statistic
#' calculated, depends on the test), and p.value (the p-value of the test). The
#' tibble is sorted by increasing p-value.
#' 
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' # Make some fake data
#' dat <- tibble::tibble(gene_id = paste('gene', 1:10, sep = ''),
#'                       terms = c('term1,term2,term3',
#'                                 NA,
#'                                 'term2,term3,term4',
#'                                 'term3',
#'                                 'term4,term5',
#'                                 'term6',
#'                                 'term6',
#'                                 'term6,term2',
#'                                 'term6,term7',
#'                                 'term6,term2'),
#'                       score = 1:10)
#' dat
#' 
#' # Test
#' gsea(dat, min_size = 2)
#' gsea(dat, min_size = 3, test = 'ks', alternative = 'less')
gsea <- function(dat, test = 'wilcoxon', alternative = 'greater', min_size = 3){
  if(!is.data.frame(dat)){
    stop("ERROR: dat must be a data.frame or tibble", call. = TRUE)
  }
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
    dplyr::filter(!is.na(statistic)) %>%
    dplyr::arrange(p.value)
  
  return(dat)
}

#' Enrichment analysis of annotation terms
#' 
#' Takes a data.frame or tibble with one gene per row, annotations
#' per gene, and a score per gene. It performs enrichment analysis
#' via \link{gsea} or \link{test_go}.
#'
#' @param dat A data.frame or tibble. It must have columns 'gene_id',
#' 'terms' and 'score'. Values in 'gene_id' must be unique.
#' @param method Either 'gsea' or 'test_go'.
#' @param ... Extra parameters for \link{gsea} or \link{test_go}. If
#' `method = 'test_go'` parameters 'genes', 'scores' and 'ontology'
#' are already provided by this function.
#'
#' @return A tibble with columnns 'term', 'size', 'statistic', and
#' 'p.value'. If at least one term has GO-like IDs, then columns 'ontology'
#' and 'annotation' will be otbained from GO.db.
#' 
#' @export
#' @importFrom magrittr %>%
terms_enrichment <- function(dat, method = 'gsea', ...){
  if(!is.data.frame(dat)){
    stop("ERROR: dat must be a data.frame or tibble", call. = TRUE)
  }
  if(any(!(c('gene_id', 'terms', 'score') %in% colnames(dat)))){
    stop("ERROR: dat must have columns 'gene_id', 'terms' and 'score'.", call. = TRUE)
  }
  if(any(duplicated(dat$gene_id))){
    stop("ERROR: values in 'gene_id' must be unique.", call. = TRUE)
  }
  
  if(method == 'gsea'){
    res <- gsea(dat = dat, ...)
    
    # If terms are GO match with annotations
    if(any(stringr::str_detect(res$term, "^GO:[0-9]{7}"))){
      res <- res %>%
        purrr::pmap_dfr(function(term, size, statistic, p.value){
          t <- GO.db::GOTERM[[term]]
          if(!is.null(t)){
            ontology <- t@Ontology
            annotation <- t@Term
          }else{
            ontology <- NA
            annotation <- NA
          }
          tibble::tibble(term = term,
                         size = size,
                         statistic = statistic,
                         p.value = p.value,
                         ontology = ontology,
                         annotation = annotation)})
    }
    
  }else if(method == 'test_go'){
    genes <- dat$score
    names(genes) <- dat$gene_id
    bp.res <- test_go(genes = genes,
                      annots = dat %>% dplyr::select(gene_id, terms),
                      ontology = 'BP',
                      ...)
    cc.res <- test_go(genes = genes,
                      annots = dat %>% dplyr::select(gene_id, terms),
                      ontology = 'CC',
                      ...)
    mf.res <- test_go(genes = genes,
                      annots = dat %>% dplyr::select(gene_id, terms),
                      ontology = 'MF',
                      ...)
    
    res <- tibble::tibble(term = as.character(NA),
                          size = as.numeric(NA),
                          p.value = as.numeric(NA),
                          ontology = as.character(NA),
                          annotation = as.character(NA)) %>%
      dplyr::filter(!is.na(term))
    
    
    if(class(bp.res$topgo_res) == "topGOresult"){
      res <- res %>%
        dplyr::bind_rows(topGO::GenTable(bp.res$topgo_data,
                                         p.value = bp.res$topgo_res,
                                         topNodes = length(bp.res$topgo_res@score)) %>%
                           dplyr::bind_cols(ontology = rep('BP', length(bp.res$topgo_res@score))) %>%
                           tibble::as_tibble() %>%
                           dplyr::select(term = GO.ID, size = Annotated, p.value, ontology, annotation = Term) %>%
                           dplyr::mutate(p.value = as.numeric(p.value)))
    }
    if(class(cc.res$topgo_res) == "topGOresult"){
      res <- res %>%
        dplyr::bind_rows(topGO::GenTable(cc.res$topgo_data,
                                         p.value = cc.res$topgo_res,
                                         topNodes = length(cc.res$topgo_res@score)) %>%
                           dplyr::bind_cols(ontology = rep('CC', length(cc.res$topgo_res@score))) %>%
                           tibble::as_tibble() %>%
                           dplyr::select(term = GO.ID, size = Annotated, p.value, ontology, annotation = Term) %>%
                           dplyr::mutate(p.value = as.numeric(p.value)))
    }
    if(class(mf.res$topgo_res) == "topGOresult"){
      res <- res %>%
        dplyr::bind_rows(topGO::GenTable(mf.res$topgo_data,
                                         p.value = mf.res$topgo_res,
                                         topNodes = length(mf.res$topgo_res@score)) %>%
                           dplyr::bind_cols(ontology = rep('MF', length(mf.res$topgo_res@score))) %>%
                           tibble::as_tibble() %>%
                           dplyr::select(term = GO.ID, size = Annotated, p.value, ontology, annotation = Term) %>%
                           dplyr::mutate(p.value = as.numeric(p.value)))
    }
    
    # Format
    res <- res %>% dplyr::arrange(p.value)
    
  }else{
    stop("ERROR: method must be 'gsea' or 'test_go'", call. = TRUE)
  }
  
  return(res)
}



#' Sign test
#' 
#' Performs sign test of a score for all annotation terms in a set of genes
#' that occurr a minimum number of times.
#'
#' @param dat A data.frame or tibble. It must contain one row per gene
#' and columns 'gene_id', 'terms', and 'score'. Column 'terms' must be
#' of type character and each entry must be a comma-separated character
#' string of all the terms that annotate the corresponding gene.
#' @param alternative The alternative hypothesis to test. Either 'greater',
#' 'less' or 'two.sided'. It  corresponds to option 'alternative' in
#' \link{binom.test} or \link{ks.test} which is used to determine if the number
#' of positive scores for a given term is more than expected given the overall
#' number of positive scores. Genes with score == 0 are removed prior to
#' any analysis. 
#' @param min_size The minimum number of genes in the group for the test
#' to be performed. Basically if the number of genes that appear in
#' 'scores' is less than 'min_size', the test won't be performed.
#'
#' @return A tibble with elements: term (the annotation term ID), n_scucesses
#' (the number of positive scores in that term), expected (the expected number
#' of positive scores given the number of genes in the term and the overall
#' probability of a positive score), n_trials (the number of genes annotated
#' with the corresponding term), p_success (the overall proability of a
#' positive score among all genes in dat) , and p.value (the p-value of the test).
#' The tibble is sorted by increasing p-value.
#' 
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' # Make some fake data
#' dat <- tibble::tibble(gene_id = paste('gene', 1:10, sep = ''),
#'                       terms = c('term1,term2,term3',
#'                                 NA,
#'                                 'term2,term3,term4',
#'                                 'term3',
#'                                 'term4,term5',
#'                                 'term6',
#'                                 'term6',
#'                                 'term6,term2',
#'                                 'term6,term7',
#'                                 'term6,term2'),
#'                       score = -5:4)
#' dat
#' 
#' # Test
#' sign_test(dat, min_size = 2)
sign_test <- function(dat, alternative = 'two.sided', min_size = 3){
  if(!is.data.frame(dat)){
    stop("ERROR: dat must be a data.frame or tibble", call. = TRUE)
  }
  if(!all(c('gene_id', 'terms', 'score') %in% colnames(dat))){
    stop("ERROR: missing columns", call. = TRUE)
  }
  if(any(duplicated(dat$gene_id))){
    stop("ERROR: values in 'gene_id' must be unique.", call. = TRUE)
  }
  
  # Clean all non-informative scores
  dat <- dat %>%
    dplyr::filter(!is.na(score)) %>%
    dplyr::filter(score != 0)
  
  # Get background prob
  p.success <- sum(dat$score > 0) / nrow(dat)
  
  # Separate scores and annotations
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
    purrr::keep(~length(.x) >= min_size) %>%
    purrr::map_dfr(function(x, scores, p.success, alternative = 'two.sided'){
      successes <- sum(scores[x] > 0)
      trials <- length(x)
      test <- binom.test(x = successes, n = trials,
                         p = p.success, alternative = alternative)
      
      tibble::tibble(n_successes = successes,
                     expected = trials * p.success,
                     n_trials = trials,
                     p_success = p.success,
                     p.value = as.numeric(test$p.value))
    }, scores = scores,
    p.success = p.success,
    alternative = alternative,
    .id = 'term') %>%
    dplyr::arrange(p.value)
  
  return(dat)
}
