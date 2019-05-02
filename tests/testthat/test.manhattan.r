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

context("manhattan")
library(HMVAR)

test_that("Manhattan plot",{
  # Create data
  dat <- tibble::tibble(ref_id = c('A', 'B'),
                        ref_pos = c(1, 2),
                        p.value = c(0.1, 9))
  # Create some fake ref sizes
  refs <- tibble::tibble(ref_id = c('B', 'A'),
                         start = c(0, 0),
                         end = c(4, 5))
  
  # no refs
  expect_identical(plotgg_manhattan(dat)$data, dat,
                   info = "Manhattan plot, no refs")
  
  # refs
  e <- dat
  e <- e %>% dplyr::bind_cols(start = c(0, 0), end = c(5, 4))
  expect_identical(plotgg_manhattan(dat, refs)$data, e,
                   info = "Manhattan plot, refs")
  
  # extra ref
  refs <- tibble::tibble(ref_id = c('B', 'A', 'C'),
                         start = c(0, 0, 0),
                         end = c(4, 5, 2))
  e <- dat
  e <- e %>%
    dplyr::bind_cols(start = c(0, 0), end = c(5, 4)) %>%
    dplyr::bind_rows(tibble::tibble(ref_id = 'C',
                                    ref_pos = NA,
                                    p.value = NA,
                                    start = 0,
                                    end = 2))
  expect_identical(plotgg_manhattan(dat, refs)$data, e,
                   info = "Manhattan plot, extra ref")
  expect_warning(print(plotgg_manhattan(dat, refs)),
                 info = "Manhattan plot, extra ref, warning")
  
  # missing ref
  refs <- tibble::tibble(ref_id = c('B'),
                         start = c(0),
                         end = c(4))
  e <- dat
  e <- e %>%
    dplyr::bind_cols(start = c(NA, 0), end = c(NA, 4)) 
  expect_identical(plotgg_manhattan(dat, refs)$data, e,
                   info = "Manhattan plot, missing ref")
})


