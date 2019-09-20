library(HMVAR)
library(tidyverse)
setwd("~/micropopgen/exp/2019/today")

site_fst <- function(freq, depth, info, map, depth_thres, method = "Weir-Cockerham"){
  dat <- match_freq_and_depth(freq = freq,
                              depth = depth,
                              info = info,
                              map = map,
                              depth_thres = depth_thres)
  
  # Using Weir-Cockerham 1984 method
  # Get basic quantities
  n_i <- table(dat$Group)
  p_i <- dat %>% split(.$Group) %>% purrr::map_dbl(~mean(.$freq))
  p_i <- p_i[names(n_i)]
  r <- length(n_i)
  n_mean <- nrow(dat) / r
  n_c <- ((r * n_mean) - (sum(n_i ^ 2) / (r * n_mean))) / (r - 1)
  # p_mean <- sum(n_i * p_i) / (r * n_mean)
  p_mean <- mean(dat$freq)
  S_sqrd <- sum(n_i * ((p_i - p_mean) ^ 2)) / ((r-1) * n_mean)
  h_mean <- 0
  
  # Calculate parts
  a <- (n_mean / n_c) * (S_sqrd - (1 / (n_mean - 1)) * ((p_mean * (1 - p_mean)) - ((r - 1) * S_sqrd / r) - (h_mean / 4)))
  b <- (n_mean / (n_mean - 1)) * (p_mean * (1 - p_mean) - ((r - 1) * S_sqrd / r) - ((2 * n_mean - 1) * h_mean / (4 * n_mean)) )
  c <- h_mean / 2
  
  tibble::tibble(r = r, n_mean = n_mean, n_c = n_c,
                 p_mean = p_mean, S_sqrd = S_sqrd,
                 h_mean = h_mean,
                 a = a, b = b,  c = c)
  
}
  

calculate_fst <- function(sites, Dat, map, depth_thres = 1, method = "Weir-Cockerham"){
  res <- sites %>%
    purrr::map_dfr(function(site, Dat, map, depth_thres = 1, method = "Weir-Cockerham"){
      dat <- Dat %>%
        purrr::map(function(d, site){d %>% filter(site_id == site)}, site = site)
      site_fst(freq = dat$freq, depth = dat$depth,
               info = dat$info, map = map, depth_thres = depth_thres)
    }, Dat = Dat, map = map, depth_thres = 1, method = method, .id = "site_id")
  
  return(res)
}

map <- read_tsv("midas/map.txt")
map <- map %>% select(sample=ID, Group)

Dat <- read_midas_data("midas/merged.snps/Veillonella_parvula_57794/", map = map, cds_only = FALSE)
gc()

map <- map %>% filter(sample %in% colnames(Dat$freq))
map$Group %>% table


fst <- calculate_fst(sites = Dat$info$site_id[1:10], Dat = Dat, map = map, depth_thres = 1)
fst

w_size <- 1000
s_size <- 100

refs <- Dat$info %>%
  split(.$ref_id)
ref <- refs[[1]]

ref
Sites <- seq(from = 1, to = max(ref$ref_pos) - w_size + 1, by = s_size) %>%
  purrr::map(function(start, w_size, info, Dat, map, depth_thres = 1){
    sites <- info %>% filter(ref_pos >= start & ref_pos < start + w_size) %>%
      dplyr::select(site_id) %>%
      unlist
    
    res <- calculate_fst(sites = sites, Dat = Dat, map = map, depth_thres = depth_thres)
    
    tibble::tibble(start = start, end = start + w_size,
                   n_loci <- length(sites),
                   Fst = sum(res$a) / (sum(res$a + res$b + res$c)))
    
  }, w_size = w_size, info = ref, Dat = Dat, map = map, depth_thres = 1)

# hist(Sites %>% map_int(length))


Sites %>%
  purrr::map_dfr(function)

