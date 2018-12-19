library(HMVAR)
library(tidyverse)
setwd("/godot/users/sur/exp/fraserv/2018/today3/")

# Parameters
midas_dir <- "../2018-12-14.hmp_mktest/genomes/Granulicatella_adiacens_61980/"
depth_thres <- 1
map_file <- "../2018-12-14.hmp_mktest/hmp_SPvsTD_map.txt"
# genes <- "638301.3.peg.283"


########## This is part of a function in HMVAR #############
########## Should make the function more modular ###########
# Read data
map <- read_tsv(map_file)
info <- read_tsv(paste0(midas_dir, "/snps_info.txt"),col_types = 'ccncccnnnnnccccc',
                 na = 'NA')
depth <- read_midas_abun(paste0(midas_dir, "/snps_depth.txt"))
freq <- read_midas_abun(paste0(midas_dir, "/snps_freq.txt"))

# Process data
# Rename map columns
map <- map %>% select(sample = ID, Group) 
# Clean info
info <- info %>% select(-locus_type, -starts_with("count_"))
# Clean depth and freq
depth <- select_samples_from_abun(depth, map)
freq <- select_samples_from_abun(freq, map)
# Clean map
map <- map %>% filter(sample %in% colnames(depth))

# Select gene data
# info <- info %>% filter(gene_id %in% genes)
info <- info %>% filter(!is.na(gene_id))
freq <- freq %>% filter(site_id %in% info$site_id)
depth <- depth %>% filter(site_id %in% info$site_id)

# Calculate MK parameters
# Calcualate snp effect
info <- determine_snp_effect(info)
# Calculate snp dist
info <- calculate_snp_dist(info = info,
                           freq = freq,
                           depth = depth,
                           map = map,
                           depth_thres = depth_thres)
###########################################################


# Match freqs and depth
depth <- depth %>% gather(key = "sample", value = 'depth', -site_id)
freq <- freq %>% gather(key = "sample", value = 'freq', -site_id)
meta <- info %>% select(site_id, ref_id, ref_pos, snp_effect, distribution)

dat <- depth %>%
  inner_join(freq, by = c("site_id", "sample")) %>%
  left_join(map, by = "sample") %>%
  filter(depth >= depth_thres) %>%
  left_join(meta, by = "site_id") %>%
  mutate(allele = replace(freq, freq < 0.5, 'major')) %>%
  mutate(allele = replace(allele, freq >= 0.5, 'minor')) %>%
  filter(distribution != "Invariant")
dat


# Count number of sites and number of variable per sample
varsites <- dat %>%
  filter(depth >= 1) %>%
  split(.$sample) %>%
  map_dfr(function(d){
    nsites <- nrow(d)
    sites <- d$freq
    nfixed <- sum(sites >= 1 | sites <= 0)
    return(tibble(nsites = nsites,
                  nfixed = nfixed,
                  nvariable = nsites - nfixed))
  }, .id = "sample") %>%
  inner_join(map, by = "sample")
varsites

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*nvariable/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
ggsave("Granulicatella_adiacens_61980_varsites_depth1.png", p1, width = 5, height = 4, dpi = 150)

p1 <- varsites %>%
  ggplot(aes(x = Group, y = nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
ggsave("Granulicatella_adiacens_61980_totsites.png", p1, width = 5, height = 4, dpi = 150)

# varsites %>%
#   ggplot(aes(x = Group, y = nvariable, col = Group)) +
#   geom_boxplot() +
#   geom_point(position = position_jitter(width = 0.2)) +
#   theme(panel.background = element_blank(),
#         panel.grid = element_blank())


varsites <- dat %>%
  filter(depth >= 15) %>%
  split(.$sample) %>%
  map_dfr(function(d){
    nsites <- nrow(d)
    sites <- d$freq
    nfixed <- sum(sites >= 1 | sites <= 0)
    return(tibble(nsites = nsites,
                  nvariable = nsites - nfixed))
  }, .id = "sample") %>%
  inner_join(map, by = "sample")
varsites

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*nvariable/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
p1
ggsave("Granulicatella_adiacens_61980_varsites_depth15.png",
       p1, width = 5, height = 4, dpi = 150)

p1 <- varsites %>%
  ggplot(aes(x = Group, y = nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
p1
ggsave("Granulicatella_adiacens_61980_totsites_depth15.png",
       p1, width = 5, height = 4, dpi = 150)
