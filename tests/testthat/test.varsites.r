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

context("varsites")
library(HMVAR)

test_that("Substitution type",{
  i <- dplyr::tibble(major_allele = 'A', minor_allele = 'C')
  e <- i %>% dplyr::bind_cols(substitution = "transversion")
  expect_identical(determine_substitution_type(i), e,
                   info = "Test transversion")
  
  i <- dplyr::tibble(major_allele = 'A', minor_allele = 'G')
  e <- i %>% dplyr::bind_cols(substitution = "transition")
  expect_identical(determine_substitution_type(i), e,
                   info = "Test transition")
  
  i <- dplyr::tibble(major_allele = 'A', minor_allele = 'g')
  e <- i %>% dplyr::bind_cols(substitution = as.character(NA)) %>%
    dplyr::filter(!is.na(substitution))
  expect_identical(determine_substitution_type(i), e,
                   info = "Test NA cleaned")
  
  i <- dplyr::tibble(major_allele = 'A', minor_allele = 'g')
  e <- i %>% dplyr::bind_cols(substitution = as.character(NA))
  expect_identical(determine_substitution_type(i, clean = FALSE), e,
                   info = "Test NA not cleaned")
  
  i <- dplyr::tibble(major_allele = 'A', minor_allele = NA)
  e <- i %>% dplyr::bind_cols(substitution = as.character(NA)) %>%
    dplyr::filter(!is.na(substitution))
  expect_identical(determine_substitution_type(i), e,
                   info = "Ambiguous when NA allele")
})