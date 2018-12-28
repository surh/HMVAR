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
             freq_thres = 0.5,
             map_file = "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/hmp_SPvsTD_map.txt",
             gene_map = '/home/sur/micropopgen/data/midas_db_v1.2/rep_genomes/Granulicatella_adiacens_61980/genome.features.gz',
             gene = "638301.3.peg.283",
             oudtir = "results/",
             genes = NA,
             heatmap = TRUE,
             lollipop = TRUE,
             depth = TRUE)

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
cat("Annotating variants...")
# Calcualate snp effect
Dat$info <- determine_snp_effect(Dat$info)
# Calculate snp dist
Dat$info <- determine_snp_dist(info = Dat$info,
                               freq = Dat$freq,
                               depth = Dat$depth,
                               map = map,
                               depth_thres = args$depth_thres)

# Match freqs and depth
cat("Matching all data...\n")
depth <- Dat$depth %>% gather(key = "sample", value = 'depth', -site_id)
freq <- Dat$freq %>% gather(key = "sample", value = 'freq', -site_id)

dat <- depth %>%
  inner_join(freq, by = c("site_id", "sample")) %>%
  left_join(map, by = "sample") %>%
  filter(depth >= args$depth_thres) %>%
  left_join(Dat$info, by = "site_id") %>%
  mutate(allele = replace(freq, freq < args$freq_thres, 'major')) %>%
  mutate(allele = replace(allele, freq >= (1 - args$freq_thres), 'minor')) %>%
  mutate(allele = replace(allele, (freq >= args$freq_thres) & (freq < (1 - args$freq_thres)), NA)) %>%
  filter(!is.na(allele)) %>%
  filter(distribution != "Invariant")
dat

#### Plotting
# Prepare outdir
if(!dir.exists(args$outdir)){
  cat("Preparing output directory...\n")
  dir.create(args$outdir)
}

# Heatmap
if(args$heatmap){
  p1 <- ggplot(dat, aes(x = ref_pos, y = sample)) +
    facet_grid(Group ~ distribution + snp_effect,
               space = "free_y", scales = "free_y") +
    geom_point(aes(col = allele), size = 0.5) +
    theme(axis.text.y = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank())
  filename <- paste0(args$outdir, "/", args$gene, "_mkheatmap.png")
  ggsave(filename, p1, width = 12, height = 8, dpi = 200)
}


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
# ggsave("638301.3.peg.283_allsites_bar.png", p1, width = 12, height = 6, dpi = 150)

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



# dat %>%
#   split(.$site_id) %>%
#   map_dfr(function(d){
#     d %>%
#     
#   })


grps <- unique(dat$Group)
# d <- dat %>% filter(site_id == "68173")

res <- dat %>% split(.$site_id) %>%
  map_dfr(function(d, grps){
    g1.major <- d %>% filter(Group == grps[1] & allele == "major") %>% nrow
    g2.major <- d %>% filter(Group == grps[2] & allele == "major") %>% nrow
    g1.minor <- d %>% filter(Group == grps[1] & allele == "minor") %>% nrow
    g2.minor <- d %>% filter(Group == grps[2] & allele == "minor") %>% nrow
    
    tibble(ref_id = unique(d$ref_id),
           ref_pos = unique(d$ref_pos),
           snp_effect = unique(d$snp_effect),
           distribution = unique(d$distribution),
           # allele = c("minor", "major"),
           g1.count.minor = c(g1.minor),
           g2.count.minor = c(g2.minor),
           g1.count.major = c(g1.major),
           g2.count.major = c(g2.major))
  }, grps = grps, .id = "site_id")

res

p1 <- ggplot(res, aes(x = ref_pos)) +
  geom_hline(yintercept = 0) +
  
  # geom_segment(aes(y = g1.count.major,
  #                  yend = g1.count.minor,
  #                  xend = ref_pos,
  #                  col = snp_effect),
  #              size = 0.2) +
  # geom_segment(aes(y = -g2.count.major,
  #                  yend = -g2.count.minor,
  #                  xend = ref_pos,
  #                  col = snp_effect),
  #              size = 0.4) +
  # scale_color_manual(values = c(NA,"black")) +


  # geom_point(aes(y = g1.count.major), color = "red", alpha = 1) +
  # geom_point(aes(y = g1.count.minor), color = "blue", alpha = 1) +
  # 
  # geom_point(aes(y = -g2.count.major), color = "red", alpha = 1) +
  # geom_point(aes(y = -g2.count.minor), color = "blue", alpha = 1) +
  
  
  
  geom_segment(aes(y = 100 * g1.count.minor / (g1.count.major + g1.count.minor),
                   yend = 0,
                   xend = ref_pos,
                   col = snp_effect),
               size = 0.2) +
  geom_segment(aes(y = -100 * g2.count.minor / (g2.count.major + g2.count.minor),
                   yend = 0,
                   xend = ref_pos,
                   col = snp_effect),
               size = 0.2) +
  geom_point(aes(y = 100 * g1.count.minor / (g1.count.major + g1.count.minor)),
             color = "black", alpha = 1,
             fill = "red",
             shape = 21) +
  geom_point(aes(y = -100 * g2.count.minor / (g2.count.major + g2.count.minor)),
             color = "black", alpha = 1,
             fill = "blue",
             shape = 21) +
  
  



  # geom_segment(aes(y = g1.count.minor + g2.count.minor ,
  #                  yend = 0,
  #                  xend = ref_pos),
  #              col = "black",
  #              size = 0.2) +
  # geom_point(aes(y = g1.count.minor + g2.count.minor,
  #                fill = snp_effect),
  #            color = "black", shape = 21) +

  
  
  # scale_x_continuous(limits = c(68173, 68300)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
p1



############################
p1 <- ggplot(res, aes(x = abs(g1.count.major - g1.count.minor))) +
  geom_density(aes(color = snp_effect, fill = snp_effect), alpha = 0.4) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
p1


p1 <- ggplot(res, aes(x = abs(g2.count.major - g2.count.minor))) +
  geom_density(aes(color = snp_effect, fill = snp_effect), alpha = 0.4) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
p1

# res <- aggregate(freq ~ site_id + Group + allele,
#                  data = dat, FUN = length, drop = FALSE) %>%
#   as.tibble %>%
#   left_join(Dat$info, by = "site_id") %>%
#   mutate(freq = replace(freq, is.na(freq), 0))
# res


# site_id ref_pos ref_id group1_count group2_count snp_effect distribution allele

res %>%
  select(site_id, Group, allele, freq, ref_id, ref_pos, gene_id, snp_effect, distribution)

p1 <- ggplot(res, aes(x = ref_pos, group = interaction(site_id,Group))) +
  geom_point(aes(y = freq, color = Group)) +
  # geom_segment(aes(y = freq))
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line())
p1


subset(res, Group == "Tongue.dorsum")
p1 <- ggplot(subset(res, Group == "Tongue.dorsum"),
             aes(x = ref_pos, group = site_id)) +
  geom_point(aes(y = freq, color = allele)) +
  geom_segment(aes(y = freq)) + 
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line())
p1






