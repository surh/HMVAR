library(HMVAR)
library(tidyverse)
setwd("~/micropopgen/exp/2019/today")

ref_window_fst <- function(ref, w_size, s_size){
  Res <- NULL
  left_ii <- 1
  right_ii <- 2
  for(start in seq(from = 1, to = max(ref$ref_pos) - s_size + 1, by = s_size)){
    end <- start + w_size
    
    for(curr_left in left_ii:nrow(ref)){
      if(ref$ref_pos[curr_left] >= start)
        break
    }
    
    for(curr_right in right_ii:nrow(ref)){
      if(ref$ref_pos[curr_right] > end - 1)
        break
    }
    
    window <- ref[(curr_left):(curr_right - 1), ] %>%
      filter(!is.na(Fst))
    # window %>% print(n = w_size)
    
    res <- tibble(start = start, end = end,
                  ref_id = unique(window$ref_id),
                  n_sites = nrow(window),
                  Fst = sum(window$a / sum(window$a + window$b + window$c)))
    Res <- Res %>%
      bind_rows(res)
    left_ii <- curr_left
    right_ii <- curr_right - 1
  }
  Res <- Res %>%
    mutate(Fst = replace(Fst, Fst < 0, 0))
  # Res
  # ref %>% filter(ref_pos >= 901 & ref_pos < 1901) %>% filter(!is.na(Fst)) %>% nrow
  
  return(Res)
}

site_fst <- function(freq, depth, info, map, depth_thres, method = "Weir-Cockerham"){
  dat <- match_freq_and_depth(freq = freq,
                              depth = depth,
                              info = info,
                              map = map,
                              depth_thres = depth_thres,
                              verbose = FALSE)
  
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

fst <- NULL
start_time <- date()
start_time
for(i in 1:nrow(Dat$info)){
  f <- site_fst(freq = Dat$freq[i,],
                depth = Dat$depth[i,],
                info = Dat$info[i,],
                map = map,
                depth_thres = 1)
  f$site_id <- Dat$info$site_id[i]
  f$ref_id <- Dat$info$ref_id[i]
  f$ref_pos <- Dat$info$ref_pos[i]
  fst <- fst %>%
    bind_rows(f)
  
  if((i %% 1000) == 0)
    cat(i, "\n")
}
end_time <- date()
end_time
fst
fst <- fst %>%
  mutate(Fst = a / (a + b + c)) %>%
  mutate(Fst = replace(Fst, Fst < 0, 0))
fst
summary(fst$Fst)
hist(fst$Fst)

snp_effects <- determine_snp_effect(info = Dat$info %>%
                                      select(site_id, major_allele, minor_allele, amino_acids) %>%
                                      filter(!is.na(amino_acids)))

p1 <- plotgg_manhattan(dat = fst,
                       pval_thres = 0.25,
                       pval_column = 'Fst',
                       log_transform = FALSE)
p1


fst <- fst %>% left_join(snp_effects %>% select(site_id, snp_effect), by = "site_id")
fst <- fst %>%
  mutate(snp_effect = as.character(snp_effect)) %>%
  mutate(snp_effect = replace(snp_effect, is.na(snp_effect), "non-coding"))
fst

fst %>%
  split(.$snp_effect) %>%
  map(~mean(.$Fst, na.rm = TRUE))
fst

ref <- fst %>%
  filter(ref_id == "NZ_ADFU01000001")
ref

ref <- ref %>%
  arrange(ref_pos)
ref

w_fst <- fst %>%
  split(.$ref_id) %>%
  purrr::map_dfr(ref_window_fst, w_size = 1000, s_size = 300)



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

