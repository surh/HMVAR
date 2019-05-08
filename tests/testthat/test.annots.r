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

context("Everything about annotations")
library(HMVAR)

test_that('annots_to_geneGO',{
  d <- tibble::tibble(gene_id = c('gene1', 'gene2', 'gene3'),
                      terms = c(NA, 'term1,term2', 'term2, term3'))
  e <- list(gene2 = c('term1', 'term2'), gene3 = c('term2', 'term3'))
  expect_identical(annots_to_geneGO(d, 'geneID2GO'), e,
                   info = "geneID to GO")
  
  e <- list(term1 = 'gene2', term2 = c('gene2', 'gene3'), term3 = 'gene3')
  expect_identical(annots_to_geneGO(d, 'GO2geneID'), e,
                   info = "GO to geneID")
  
  expect_error(annots_to_geneGO(d, 'bla'),
               info = 'Bad direction')
})

test_that('test_go', {
  # Helper function
  class_and_dimensions <- function(...){
    res <- test_go(...)
    n <- length(res$topgo_res@score)
    res <- topGO::GenTable(res$topgo_data, p.value = res$topgo_res, topNodes = n)
    
    return(list(class = 'data.frame', dim = dim(res)))
  }
  
  # Basic usage
  a <- tibble::tibble(gene_id = paste('gene', 1:5, sep = ''),
                      terms = c('GO:0005575,GO:0044464,GO:0071944',
                                'GO:0044464',
                                'GO:0071944',
                                'GO:0005623,GO:0005886,GO:0071944',
                                'GO:0005886,GO:0016020,GO:0044464'))
  g <- c('gene1', 'gene2', 'gene5')
  e <- list(class = 'data.frame', dim = as.integer(c(6,6)))
  expect_identical(class_and_dimensions(genes = g, annots = a, ontology = 'CC', node_size = 2),
                   expected = e,
                   info = "Basic usage")
  
  # Pass list of annotations
  a.list <- annots_to_geneGO(a, 'geneID2GO')
  expect_equal(test_go(genes = g, annots = a, ontology = 'CC', node_size = 2),
               test_go(genes = g, annots = a.list, ontology = 'CC', node_size = 2),
               info = "Basic usage, pass list of annotations")
  
  # Pass scores
  g.scores <- c(-1, -1, 0, 0, -1)
  names(g.scores) <- a$gene_id
  expect_equal(test_go(genes = g, annots = a, ontology = 'CC', node_size = 2),
               test_go(genes = g.scores, annots = a, ontology = 'CC', node_size = 2, score_threshold = 0),
               info = "Basic usage, pass scores")
  
  # Badd annot object
  expect_error(test_go(genes = g, annots = c('a', 'b', 'c'), ontology = 'CC', node_size = 2),
               info = "Bad annots")
})