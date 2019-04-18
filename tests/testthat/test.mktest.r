context("MK tests")
library(HMVAR)

test_that("SNP distribution (no-missing)",{
  # Set up data
  m <- dplyr::tibble(sample = letters[1:6],
                     Group = rep(LETTERS[1:2], each = 3))
  i <- dplyr::tibble(site_id = 's1')
  d <- dplyr::tibble(site_id = i$site_id,
                     a = 1, b = 1, c = 1 ,
                     d = 1, e = 1, f = 1)
  
  # Fixed
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 0, b = 0, c = 0,
                     d = 1, e = 1, f = 1)
  e <- tibble(site_id = i$site_id,
              distribution = factor("Fixed", levels = c("Fixed", "Invariant", "Polymorphic")))
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0),
                   e,
                   info = "Perfect fixed")
  
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 1, b = 0, c = 0,
                     d = 1, e = 1, f = 0)
  e <- tibble(site_id = i$site_id,
              distribution = factor("Polymorphic", levels = c("Fixed", "Invariant", "Polymorphic")))
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0),
                   e,
                   info = "Perfect polymorphic")
  
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 0, b = 0, c = 0,
                     d = 0, e = 0, f = 0)
  e <- tibble(site_id = i$site_id,
              distribution = factor("Invariant", levels = c("Fixed", "Invariant", "Polymorphic")))
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0),
                   e,
                   info = "Perfect invariant")
  
  
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 1),
                   determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0),
                   info = "Compare freq_thres & (1 - freq_thres)")
})

test_that("SNP distributions (missing data)",{
  # Set up data
  m <- dplyr::tibble(sample = letters[1:6],
                     Group = rep(LETTERS[1:2], each = 3))
  i <- dplyr::tibble(site_id = 's1')
  
  d <- dplyr::tibble(site_id = i$site_id,
                     a = 0, b = 0, c = 1 ,
                     d = 1, e = 1, f = 0)
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 0, b = 0, c = 0,
                     d = 1, e = 1, f = 1)
  e <- tibble(site_id = i$site_id,
              distribution = factor(NA, levels = c("Fixed", "Invariant", "Polymorphic")))
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0),
                   e,
                   info = "Not enough samples in 1 group")
  
  d <- dplyr::tibble(site_id = i$site_id,
                     a = 0, b = 0, c = 1 ,
                     d = 1, e = 0, f = 0)
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 0, b = 0, c = 0,
                     d = 1, e = 1, f = 1)
  e <- tibble(site_id = i$site_id,
              distribution = factor(NA, levels = c("Fixed", "Invariant", "Polymorphic")))
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0),
                   e,
                   info = "Not enough samples in either group")
  
  
  d <- dplyr::tibble(site_id = i$site_id,
                     a = 1, b = 1, c = 1 ,
                     d = 1, e = 1, f = 1)
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 0.5, b = 0.5, c = 0,
                     d = 1, e = 1, f = 1)
  e <- tibble(site_id = i$site_id,
              distribution = factor(NA, levels = c("Fixed", "Invariant", "Polymorphic")))
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0),
                   e,
                   info = "Not enough samples with non 0.5 MAF in one group")
  
  d <- dplyr::tibble(site_id = i$site_id,
                     a = 1, b = 1, c = 1 ,
                     d = 1, e = 1, f = 1)
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 0.2, b = 0.2, c = 0,
                     d = 1, e = 1, f = 1)
  e <- tibble(site_id = i$site_id,
              distribution = factor(NA, levels = c("Fixed", "Invariant", "Polymorphic")))
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0),
                   e,
                   info = "Not enough samples with valid MAF in one group")
  
  d <- dplyr::tibble(site_id = i$site_id,
                     a = 1, b = 0, c = 1 ,
                     d = 1, e = 1, f = 1)
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 0.2, b = 0, c = 0,
                     d = 1, e = 1, f = 1)
  e <- tibble(site_id = i$site_id,
              distribution = factor(NA, levels = c("Fixed", "Invariant", "Polymorphic")))
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0),
                   e,
                   info = "Combination missing valid MAF and depth")
})


# test_that("SNP distributions (fractional)",{
#   f <- dplyr::tibble(site_id = i$site_id,
#                      a = 0.1, b = 0.1, c = 0.1,
#                      d = 0, e = 0, f = 0)
#   e <- tibble(site_id = i$site_id,
#               distribution = factor("Invariant", levels = c("Fixed", "Invariant", "Polymorphic")))
#   expect_identical(determine_snp_dist(info = i,
#                                       freq = f,
#                                       depth = d,
#                                       map = m,
#                                       depth_thres = 1,
#                                       freq_thres = 0.9),
#                    e,
#                    info = "Check that freq_thres uses >= & <=")
# })
 
test_that("SNP distribution errors",{
  # Set up data
  m <- dplyr::tibble(sample = letters[1:6],
                     Group = rep(LETTERS[1:2], each = 3))
  i <- dplyr::tibble(site_id = 's1')
  d <- dplyr::tibble(site_id = i$site_id,
                     a = 1, b = 1, c = 1 ,
                     d = 1, e = 1, f = 0)
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 0, b = 0, c = 0,
                     d = 1, e = 1, f = 1)
  expect_error(determine_snp_dist(info = i,
                                  freq = f,
                                  depth = d,
                                  map = m,
                                  depth_thres = 1,
                                  freq_thres = 2),
                   info = "MAF threshold out of range")
})
 
# test_that("SNP effects",{})
# 
# test_that("MK table")