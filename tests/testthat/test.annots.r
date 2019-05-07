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