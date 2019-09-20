library(HMVAR)
library(tidyverse)
setwd("~/micropopgen/exp/2019/today")
devtools::document(pkg = "~/micropopgen/src/HMVAR/")

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

p1 <- ggplot(w_fst, aes(x = Fst, y = n_sites)) +
  geom_point() +
  geom_smooth(method = "loess") +
  AMOR::theme_blackbox() 
p1

