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

context("DoS")
library(HMVAR)

test_that("calculate_dos", {
  # Test DoS
  d <- tibble::tibble(Dn = c(1, 0),
                      Ds = c(1, 1),
                      Pn = c(1, 1),
                      Ps = c(1, 1))
  e <- d
  dos <- c(0, -0.5)
  e <- e %>% dplyr::bind_cols(DoS = dos,
                              p.value = 2*(1 - pnorm(q = abs(dos) / sd(dos, na.rm = TRUE))))
  expect_identical(calculate_dos(d, test = TRUE, clean = FALSE), e,
                   info = "Test DoS")
  
  
  # Test DoS with NA
  d <- tibble::tibble(Dn = c(1, 0, 0),
                      Ds = c(1, 1, 0),
                      Pn = c(1, 1, 1),
                      Ps = c(1, 1, 1))
  e <- d
  dos <- c(0, -0.5, NaN)
  e <- e %>% dplyr::bind_cols(DoS = dos,
                              p.value = 2*(1 - pnorm(q = abs(dos) / sd(dos, na.rm = TRUE))))
  expect_identical(calculate_dos(d, test = TRUE, clean = FALSE), e,
                   info = "Test DoS with NA")
  
  # Test DoS with NA, clean
  d <- tibble::tibble(Dn = c(1, 0, 0),
                      Ds = c(1, 1, 0),
                      Pn = c(1, 1, 1),
                      Ps = c(1, 1, 1))
  e <- d[1:2, ]
  dos <- c(0, -0.5)
  e <- e %>% dplyr::bind_cols(DoS = dos,
                              p.value = 2*(1 - pnorm(q = abs(dos) / sd(dos, na.rm = TRUE))))
  expect_identical(calculate_dos(d, test = TRUE, clean = TRUE), e,
                   info = "Test DoS with NA, clean")
  
  
  # Test DoS, p.value
  d <- tibble::tibble(Dn = c(1, 0),
                      Ds = c(1, 1),
                      Pn = c(1, 1),
                      Ps = c(1, 1),
                      p.value = c(2, 2))
  e <- d
  dos <- c(0, -0.5)
  e <- e %>% dplyr::bind_cols(DoS = dos,
                              DoS.p.value = 2*(1 - pnorm(q = abs(dos) / sd(dos, na.rm = TRUE))))
  expect_identical(calculate_dos(d, test = TRUE, clean = FALSE), e,
                   info = "Test DoS, p.value")
})

test_that("dos", {
  i <- tibble::tibble(site_id = c('snv1', 'snv2', 'snv3', 'snv4'),
                      minor_allele = rep('A', 4),
                      major_allele = rep('C', 4),
                      gene_id = rep('gene1', 4),
                      amino_acids = c('A,A,A,A', 
                                      'A,A,A,A',
                                      'A,D,E,F',
                                      'A,D.E.F'))
  f <- tibble::tibble(site_id = c('snv1', 'snv2', 'snv3', 'snv4'),
                      s1 = c(1, 1, 1, 1),
                      s2 = c(1, 0, 1, 0),
                      s3 = c(0, 1, 0, 1),
                      s4 = c(0, 0, 0, 0))
  d <- tibble::tibble(site_id = c('snv1', 'snv2', 'snv3', 'snv4'),
                      s1 = c(1, 1, 1, 1),
                      s2 = c(1, 1, 1, 1),
                      s3 = c(1, 1, 1, 1),
                      s4 = c(1, 1, 1, 1))
  m <- tibble::tibble(sample = c('s1', 's2', 's3', 's4'),
                      Group = c('A', 'A', 'B', 'B'))
  e <- tibble::tibble(gene_id = 'gene1',
                      Dn = as.integer(1),
                      Ds = as.integer(1),
                      Pn = as.integer(1),
                      Ps = as.integer(1),
                      DoS = 0)
  expect_identical(dos(info = i,
                       freq = f,
                       depth = d,
                       map = m,
                       depth_thres = 1,
                       freq_thres = 0.5,
                       test = FALSE,
                       clean = FALSE),
                   e, info = "DoS")
  
  # Missing column
  i <- tibble::tibble(site_id = c('snv1', 'snv2', 'snv3', 'snv4'),
                      major_allele = rep('C', 4),
                      gene_id = rep('gene1', 4),
                      amino_acids = c('A,A,A,A', 
                                      'A,A,A,A',
                                      'A,D,E,F',
                                      'A,D.E.F'))
  expect_error(dos(i), info = "Missing column in info")
})