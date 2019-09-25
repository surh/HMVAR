library(HMVAR)
library(tidyverse)
setwd("~/micropopgen/exp/2019/today")
# devtools::document(pkg = "~/micropopgen/src/HMVAR/")

map <- read_tsv("midas/map.txt")
map <- map %>% select(sample=ID, Group)

Dat <- read_midas_data("midas/merged.snps/Veillonella_parvula_57794/", map = map, cds_only = FALSE)
gc()

map <- map %>% filter(sample %in% colnames(Dat$freq))
map$Group %>% table


# i <- 1
# freq <- Dat$freq[i,]
# support <- Dat$depth[i,]
# info <- Dat$info[i, ]
# map <- map
# support_thres <- 1
# method = "Fst_pool"
# 
# dat <- match_freq_and_depth(freq = freq,
#                             depth = support,
#                             info = info,
#                             map = map,
#                             depth_thres = support_thres,
#                             verbose = FALSE)

calculate_fst <- function(Dat, map,
                          method = "Weir-Cockerham",
                          support_thres = 1,
                          w_size = NULL,
                          s_size = NULL,
                          sorted = FALSE,
                          verbose = TRUE){
  fst <- NULL
  # start_time <- date()
  # start_time
  for(i in 1:nrow(Dat$info)){
    f <- site_fst(freq = Dat$freq[i,],
                  support = Dat$depth[i,],
                  info = Dat$info[i,],
                  map = map,
                  support_thres = support_thres,
                  method = method)
    
    f$site_id <- Dat$info$site_id[i]
    f$ref_id <- Dat$info$ref_id[i]
    f$ref_pos <- Dat$info$ref_pos[i]
    fst <- fst %>%
      dplyr::bind_rows(f)
    
    if(verbose && ((i %% 1000) == 0))
      cat(i, "\n")
  }
  # end_time <- date()
  # end_time
  
  if(method == "Weir-Cockerham"){
    fst <- fst %>%
      dplyr::mutate(Fst = a / (a + b + c)) %>%
      dplyr::mutate(Fst = replace(Fst, Fst < 0, 0))
  }else if(method == "Fstpool"){
    fst <- fst %>%
      dplyr::mutate(Fst_pool = (MSP - MSI) / (MSP + ((n_c - 1) * MSI)) ) %>%
      dplyr::mutate(Fst_pool = replace(Fst_pool, Fst_pool < 0, 0))
  }
    
  if(is.numeric(w_size) && is.numeric(s_size)){
    w_fst <- fst %>%
      split(.$ref_id) %>%
      purrr::map_dfr(ref_window_fst, w_size = 1000, s_size = 300, sorted = sorted) 
  }else{
    return(list(fst = fst, w_fst = NULL))
  }
  
  return(list(fst = fst, w_fst = w_fst))
}

fst <- calculate_fst(Dat = list(info = Dat$info[1:10000, ],
                                depth = Dat$depth[1:10000,],
                                freq = Dat$freq[1:10000,]),
                     method = "Weir-Cockerham",
                     support_thres = 1,
                     map = map,
                     sorted = TRUE,
                     verbose = TRUE)

fst_pool <- calculate_fst(Dat = list(info = Dat$info[1:10000, ],
                                     depth = Dat$depth[1:10000,],
                                     freq = Dat$freq[1:10000,]),
                          method = "Fstpool",
                          support_thres = 1,
                          map = map,
                          sorted = TRUE,
                          verbose = TRUE)




# summary(fst$Fst)
# hist(fst$Fst)
# 
# snp_effects <- determine_snp_effect(info = Dat$info %>%
#                                       select(site_id, major_allele, minor_allele, amino_acids) %>%
#                                       filter(!is.na(amino_acids)))
w_fst <- window_fst(dat = fst$fst, w_size = 1000, s_size = 300, sorted = TRUE)
w_fst_pool <- window_fst(dat = fst_pool$fst, w_size = 1000, s_size = 300, sorted = TRUE)

w_fst
w_fst_pool 


w_fst_pool %>%
  arrange(desc(Fst_pool))



p1 <- ggplot(w_fst, aes(x = (start + end) / 2, y = Fst)) +
  facet_grid(~ ref_id, scales = "free_x", space = "free_x") +
  # geom_point(aes(size = n_sites)) +
  geom_point() +
  geom_hline(yintercept = 0.25, color = "red", size = 2) +
  ggplot2::theme(panel.background = ggplot2::element_blank(), 
                 panel.grid = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(color = "black"), 
                 axis.title = ggplot2::element_text(color = "black", face = "bold"), 
                 axis.line.x.bottom = ggplot2::element_line(),
                 axis.line.y.left = ggplot2::element_line())
p1

p1 <- ggplot(w_fst_pool, aes(x = (start + end) / 2, y = Fst_pool)) +
  facet_grid(~ ref_id, scales = "free_x", space = "free_x") +
  # geom_point(aes(size = n_sites)) +
  geom_point() +
  geom_hline(yintercept = 0.25, color = "red", size = 2) +
  ggplot2::theme(panel.background = ggplot2::element_blank(), 
                 panel.grid = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(color = "black"), 
                 axis.title = ggplot2::element_text(color = "black", face = "bold"), 
                 axis.line.x.bottom = ggplot2::element_line(),
                 axis.line.y.left = ggplot2::element_line())
p1

p1 <- ggplot(w_fst, aes(x = Fst, y = n_sites)) +
  geom_point() +
  geom_smooth(method = "loess") +
  AMOR::theme_blackbox() 
p1

