context("MK tests")
library(HMVAR)

test_that("SNP distribution (no-missing)",{
  # Testing case where depth is always >= than thres
  # Testing also freq_thres
  # Perfect separation
  m <- dplyr::tibble(sample = letters[1:6],
                     Group = rep(LETTERS[1:2], each = 3))
  i <- dplyr::tibble(site_id = 's1')
  d <- dplyr::tibble(site_id = i$site_id,
                     a = 1, b = 1, c = 1 ,
                     d = 1, e = 1, f = 1)
  
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
                                      freq_thres = 0.5),
                   e)
  
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 0, b = 0, c = 0,
                     d = 1, e = 1, f = 0)
  e <- tibble(site_id = i$site_id,
              distribution = factor("Polymorphic", levels = c("Fixed", "Invariant", "Polymorphic")))
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0.5),
                   e)
  
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 0.1, b = 0.1, c = 0.1,
                     d = 0, e = 0, f = 0)
  e <- tibble(site_id = i$site_id,
              distribution = factor("Invariant", levels = c("Fixed", "Invariant", "Polymorphic")))
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0.5),
                   e)
  
  f <- dplyr::tibble(site_id = i$site_id,
                     a = 0.1, b = 0.1, c = 0.1,
                     d = 0, e = 0, f = 0)
  e <- tibble(site_id = i$site_id,
              distribution = factor("Invariant", levels = c("Fixed", "Invariant", "Polymorphic")))
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0.9),
                   e,
                   info = "Check that freq_thres uses >= & <=")
  
  expect_identical(determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0.9),
                   determine_snp_dist(info = i,
                                      freq = f,
                                      depth = d,
                                      map = m,
                                      depth_thres = 1,
                                      freq_thres = 0.1),
                   info = "Compare freq_thres & (1 - freq_thres)")
})
