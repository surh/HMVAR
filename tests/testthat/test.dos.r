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