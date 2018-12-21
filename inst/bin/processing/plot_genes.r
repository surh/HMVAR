library(HMVAR)
library(tidyverse)
library(argparser)

# setwd("/godot/users/sur/exp/fraserv/2018/today3")

# mktest_file <- "mktest_selected/results/Leptotrichia_shahii_5891_mktest.txt"
# Genes <- read_tsv(mktest_file)
# gene <- Genes$gene[2]
# genes <- gene

# Parameters
args <- list(midas_dir = "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/genomes/Granulicatella_adiacens_61980/",
             depth_thres = 1,
             map_file = "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/hmp_SPvsTD_map.txt",
             gene = "638301.3.peg.283",
             genes = NULL)

# !!!!
genes <- args$gene

# Read map
cat("Reading map...\n")
map <- read_tsv(args$map_file)
map <- map %>%
  select(sample = ID,
         everything())

# Read data
cat("Read MIDAS data...\n")
Dat <- read_midas_data(args$midas_dir, map = map, genes = genes)

# Annotate variants
cat("Annotate variants...")
# Calcualate snp effect
Dat$info <- determine_snp_effect(Dat$info)
# Calculate snp dist
Dat$info <- determine_snp_dist(info = Dat$info,
                               freq = Dat$freq,
                               depth = Dat$depth,
                               map = map,
                               depth_thres = args$depth_thres)






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

dat <- dat %>%
  mutate(site_id = factor(site_id,
                          levels = as.character(unique(sort(as.numeric(dat$site_id))))))


# ggplot(subset(dat, site_id %in% c("68173","68332")),
#        aes(x = freq, group = site_id)) +
#   facet_grid(~ Group) +
#   geom_density(aes(y = ..scaled..))
# 
# ggplot(dat,
#        aes(x = freq, group = site_id)) +
#   facet_grid(~ Group) +
#   geom_density(aes(y = ..scaled..))


# p1 <- ggplot(dat, aes(x = ref_pos, y = sample)) +
#   facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
#   geom_point(aes(size = depth, col = snp_effect)) +
#   theme(axis.text.y = element_blank(),
#         panel.grid = element_blank(),
#         panel.background = element_blank())
# p1
# 
# 
# p1 <- ggplot(dat, aes(x = ref_pos, y = sample)) +
#   facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
#   geom_point(aes(size = freq, col = snp_effect)) +
#   theme(axis.text.y = element_blank(),
#         panel.grid = element_blank(),
#         panel.background = element_blank())
# p1

p1 <- ggplot(dat, aes(x = ref_pos, y = sample)) +
  facet_grid(Group ~ distribution + snp_effect,
             space = "free_y", scales = "free_y") +
  geom_point(aes(col = allele), size = 0.5) +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank())
# p1
ggsave("638301.3.peg.283_allsites.svg", p1, width = 6, height = 4)
ggsave("638301.3.peg.283_allsites.png", p1, width = 12, height = 8, dpi = 150)

# p1 <- ggplot(dat, aes(x = ref_pos, y = sample)) +
#   facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
#   geom_point(aes(size = log2(freq*depth + 1), col = snp_effect)) +
#   theme(axis.text.y = element_blank(),
#         panel.grid = element_blank(),
#         panel.background = element_blank())
# p1

###
# p1 <- ggplot(dat, aes(x = ref_pos, y = sample)) +
#   facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
#   geom_tile(aes(fill = freq, col = snp_effect)) +
#   theme(axis.text.y = element_blank(),
#         panel.grid = element_blank(),
#         panel.background = element_blank())
# p1
# 
# p1 <- ggplot(dat, aes(x = ref_pos, y = sample)) +
#   facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
#   geom_tile(aes(fill = allele)) +
#   theme(axis.text.y = element_blank(),
#         panel.grid = element_blank(),
#         panel.background = element_blank())
# p1


###
p1 <- ggplot(dat, aes(x = ref_pos, y = depth)) +
  geom_line(aes(color = Group, group = sample)) +
  scale_y_log10() +
  theme(panel.background = element_rect(color = "black", fill = NA),
        panel.grid = element_blank())
p1

# p1 <- ggplot(dat, aes(x = ref_pos, y = freq)) +
#   geom_line(aes(color = Group, group = sample)) +
#   theme(panel.background = element_rect(color = "black", fill = NA),
#         panel.grid = element_blank())
# p1

# p1 <- ggplot(dat, aes(x = ref_pos, y = freq*depth)) +
#   geom_line(aes(color = Group, group = sample)) +
#   theme(panel.background = element_rect(color = "black", fill = NA),
#         panel.grid = element_blank())
# p1

p1 <- ggplot(dat, aes(x = site_id, y = allele, fill = allele)) +
  facet_grid(Group ~ .) +
  geom_bar(stat = "identity") +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank())
p1
ggsave("638301.3.peg.283_allsites_bar.png", p1, width = 12, height = 6, dpi = 150)

###
# p1 <- ggplot(dat, aes(x = freq)) +
#   facet_grid(Group ~ snp_effect) +
#   geom_density(aes(fill = Group, y = ..scaled..))
# p1
# 
# p1 <- ggplot(dat, aes(x = freq)) +
#   facet_grid(Group ~ snp_effect) +
#   geom_histogram(aes(fill = Group), bins = 10)
# p1










