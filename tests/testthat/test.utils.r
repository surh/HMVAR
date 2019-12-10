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

context("utils")
library(HMVAR)

test_that("Check map", {
  f <- system.file("toy_example/map.txt", package = "HMVAR")
  expect_identical(HMVAR:::check_map(map = f),
                   HMVAR:::check_map(map_file = f),
                   info = "map_file compatibility")
  
  m <- HMVAR:::check_map(f) %>%
    dplyr::select(ID = sample, Group)
  expect_error(HMVAR:::check_map(map_file = m),
               info = "bad column names")
  
  m <- HMVAR:::check_map(f) %>%
    dplyr::mutate(Group = replace(Group,
                                  list = Group == 'AAA',
                                  values = 'non.BBB'))
  expect_identical(HMVAR:::check_map(map = f, focal_group = 'BBB'), m,
                   info = "process column names")
})