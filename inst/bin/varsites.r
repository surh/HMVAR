#!/usr/bin/env Rscript

# (C) Copyright 2018 Sur Herrera Paredes
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

library(HMVAR)
library(tidyverse)
library(argparser)

################ FUNCTIONS ################
process_arguments <- function(){
  p <- arg_parser(paste0("Calculates and plots the number of variable sites ",
                         "per sample within a genome (or a subset of genes in that ",
                         "genome."))
  
  p <- add_argument(p, "midas_dir",
                    help = paste0("Directory with output from merging MIDAS SNPs. ",
                                  "Must contain 'snps_info.txt', 'snps_depth.txt' ",
                                  "and 'snps_freq.txt'."),
                    type = "character")
  p <- add_argument(p, "map_file",
                    help = paste0("Mapping file associating samples with groups. ",
                                  "Must have 'Group' and 'ID' columns."),
                    type = "character")
  p <- add_argument(p, "--genes",
                    help = paste0("File with subset of genes to analyze. ",
                                  "It must be one gene ID per line without headers."),
                    type = "character",
                    default = NULL)
  p <- add_argument(p, "--cds_only",
                    help = paste0("Flag indicating whether only sites at CDSs should ",
                                  "be considered."),
                    flag = TRUE)
  p <- add_argument(p, "--depth_thres",
                    help = paste0("Minumum number of reads at a given site ",
                                  "at a given sample, for that site in that ",
                                  "sample to be included in calculations"),
                    default = 1,
                    type = "integer")
  p <- add_argument(p, "--freq_thres",
                    help = paste0("This value is the threshold that will be ",
                                  "used to determine which allele is present in ",
                                  "a sample. The value from ",
                                  "min(freq_thres, 1 - freq_thres) represents the ",
                                  "distance from 0 or 1, for a site to be assigned ",
                                  "to the major or minor allele respectively. ",
                                  "It must be a value in [0,1]."),
                    default = 0.5,
                    type = "double")
  p <- add_argument(p, "--outdir",
                    help = paste0("Directory to write (and create if needed) to write ",
                                  "the output."),
                    default = "results/",
                    type = "character")
  p <- add_argument(p, "--plot_position",
                    help = paste0("Flag indicating whether to plot the proportion ",
                                  "of variable sites by genomic position. Might be ",
                                  "time consuming."),
                    flag = TRUE)
  p <- add_argument(p, "--prefix",
                    help = paste0("Prefix for output filenames."),
                    default = "out",
                    type = "character")

  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(args$freq_thres < 0 || args$freq_thres > 1)
    stop("ERROR: freq_thres must be a value in [0,1].")
  if(args$depth_thres < 0)
    stop("ERROR: depth_thres must be a non-negative integer.")
  if(!dir.exists(args$midas_dir))
    stop("ERROR: midas_dir does not exist.")
  if(!file.exists(args$map_file))
    stop("ERROR: map_file does not exist.")
  
  return(args)
}

# Parameters
# args <- list(midas_dir = "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/genomes/Granulicatella_adiacens_61980/",
#              map_file = "/godot/users/sur/exp/fraserv/2018/2018-12-14.hmp_mktest/hmp_SPvsTD_map.txt",
#              genes = NA,
#              depth_thres = 1,
#              freq_thres = 0.5,
#              cds_only = FALSE,
#              plot_position = TRUE,
#              prefix = "Granulicatella_adiacens_61980",
#              outdir = "results/")
# args <- process_arguments()

args <- list(midas_dir = system.file("toy_example/merged.snps/", package = "HMVAR"),
             map_file = system.file("toy_example/map.txt", package = "HMVAR"),
             genes = NA,
             depth_thres = 1,
             freq_thres = 0.5,
             cds_only = FALSE,
             plot_position = TRUE,
             prefix = "toyprefix",
             outdir = "results/")

# Read genes
if(is.na(args$genes)){
  args$genes <- NULL
}else{
  cat("Reading list of genes...\n")
  args$genes <- read_tsv(args$genes, col_names = FALSE, col_types = 'c')
  args$genes <- genes$X1
}

# Read map
cat("Reading map...\n")
map <- read_tsv(args$map_file,
                col_types = cols(ID = col_character(),
                                 Group = col_character())) %>%
  select(sample = ID,
         everything())

# # Read MIDAS data
# cat("Reading MIDAS output...\n")
# info <- read_tsv(paste0(args$midas_dir, "/snps_info.txt"),
#                  col_types = 'ccncccnnnnnccccc',
#                  na = 'NA')
# depth <- read_midas_abun(paste0(args$midas_dir, "/snps_depth.txt"))
# freq <- read_midas_abun(paste0(args$midas_dir, "/snps_freq.txt"))
# 
# # Process data
# cat("Processing MIDAS data...\n")
# # Clean info
# info <- info %>% 
#   select(-locus_type, -starts_with("count_"), -ref_allele)
# # Clean depth and freq
# depth <- depth %>%
#   select(site_id, intersect(map$sample, colnames(depth)) )
# freq <- freq %>%
#   select(site_id, intersect(map$sample, colnames(freq)) )
# # Clean map
# map <- map %>% 
#   filter(sample %in% colnames(depth))
# 
# # Select gene/cds data
# if(args$cds_only){
#   cat("Selecting CDS only...\n")
#   info <- info %>% 
#     filter(!is.na(gene_id))
# }
# if(!is.null(genes)){
#   cat("Selecting chosen genes...\n")
#   info <- info %>%
#     filter(gene_id %in% genes)
# }
# freq <- freq %>% 
#   filter(site_id %in% info$site_id)
# depth <- depth %>% 
#   filter(site_id %in% info$site_id)

# Read and process midas sites
midas_data <- read_midas_data(midas_dir = args$midas_dir,
                              map = map, genes = args$genes,
                              cds_only = FALSE)
cat("Determining snp effect...\n")
midas_data$info <- determine_snp_effect(midas_data$info)
cat("Determining snp distribution...\n")
midas_data$info <- determine_snp_dist(info = midas_data$info,
                                      freq = midas_data$freq,
                                      depth = midas_data$depth,
                                      map = map,
                                      depth_thres = args$depth_thres,
                                      freq_thres = args$freq_thres,
                                      clean = FALSE)

# Subsitution type
cat("Determining substitution type...\n")
info <- info %>%
  add_column(substitution = info %>%
               pmap_chr(function(major_allele, minor_allele, ...){
                 base_type <- c(A = "purine", C = "pyrimidine", G = "purine", T = "pyrimidine")
                 if(base_type[major_allele] == base_type[minor_allele]){
                   substitution <- "transition" 
                 }else{
                   substitution <- "transversion"
                 }
                 return(substitution)
               }))


# Match freqs and depth
cat("Matching all data...\n")
depth <- depth %>% gather(key = "sample", value = 'depth', -site_id)
freq <- freq %>% gather(key = "sample", value = 'freq', -site_id)

dat <- depth %>%
  inner_join(freq, by = c("site_id", "sample")) %>%
  left_join(map, by = "sample") %>%
  filter(depth >= args$depth_thres) %>%
  left_join(info, by = "site_id") %>%
  filter(distribution != "Invariant")
# dat

# Prepare output dir
if(!dir.exists(args$outdir)){
  cat("Preparing output directory...\n")
  dir.create(args$outdir)
}

# Count number of sites and number of variable per sample
cat("Calcuating variable sites...\n")
varsites <- dat %>%
  filter(depth >= 1) %>%
  split(.$sample) %>%
  map_dfr(function(d){
    nsites <- nrow(d)
    sites <- d$freq
    nfixed <- sum(sites >= 1 | sites <= 0)
    ntransitions <- sum(d$substitution == "transition")
    nsynonymous <- sum(d$snp_effect == "synonymous", na.rm = TRUE)
    nnonsynonymous <- sum(d$snp_effect == "non-synonymous", na.rm = TRUE)
    return(tibble(nsites = nsites,
                  fixed = nfixed,
                  variable = nsites - nfixed,
                  transitions = ntransitions,
                  transversions = nsites - ntransitions,
                  synonymous = nsynonymous,
                  non.synonymous = nnonsynonymous,
                  other = nsites - nsynonymous - nnonsynonymous))
  }, .id = "sample") %>%
  inner_join(map, by = "sample")
# varsites
cat("\tWriting number of variable sites to file...\n")
filename <- paste0(args$outdir, "/", args$prefix, ".varsites.txt")
write_tsv(varsites, path = filename)


# Plot total. variable and fixed sites
cat("\tPlotting...\n")
p1 <- varsites %>%
  ggplot(aes(x = Group, y = nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line())
# p1
filename <- paste0(args$outdir, "/", args$prefix, ".varsites.nsites.svg")
ggsave(filename, p1, width = 5, height = 4)

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*variable/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_continuous(limits = c(0,100)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line())
# p1
filename <- paste0(args$outdir, "/", args$prefix, ".varsites.percvarsites.svg")
ggsave(filename, p1, width = 5, height = 4)

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*fixed/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_continuous(limits = c(0,100)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line())
# p1
filename <- paste0(args$outdir, "/", args$prefix, ".varsites.percfixsites.svg")
ggsave(filename, p1, width = 5, height = 4)

# Plot transitions and transversions
# p1 <- varsites %>%
#   ggplot(aes(x = Group, y = transitions, col = Group)) +
#   geom_boxplot() +
#   geom_point(position = position_jitter(width = 0.2)) +
#   theme(panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA))
# p1
# 
# p1 <- varsites %>%
#   ggplot(aes(x = Group, y = transversions, col = Group)) +
#   geom_boxplot() +
#   geom_point(position = position_jitter(width = 0.2)) +
#   theme(panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA))
# p1

p1 <- varsites %>%
  ggplot(aes(x = Group, y = transitions/transversions, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_log10(breaks = function(x){seq(from = x[1], to = x[2], length.out = 5)}) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line())
# p1
filename <- paste0(args$outdir, "/", args$prefix, ".varsites.subs_type.svg")
ggsave(filename, p1, width = 5, height = 4)

# Plot synonymous and non-synonymous
p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*synonymous/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_continuous(limits = c(0,100)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line())
# p1
filename <- paste0(args$outdir, "/", args$prefix, ".varsites.percsynom.svg")
ggsave(filename, p1, width = 5, height = 4)

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*non.synonymous/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_continuous(limits = c(0,100)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line())
# p1
filename <- paste0(args$outdir, "/", args$prefix, ".varsites.percnonsynom.svg")
ggsave(filename, p1, width = 5, height = 4)

p1 <- varsites %>%
  ggplot(aes(x = Group, y = 100*other/nsites, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_continuous(limits = c(0,100)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line())
# p1
filename <- paste0(args$outdir, "/", args$prefix, ".varsites.percother.svg")
ggsave(filename, p1, width = 5, height = 4)

p1 <- varsites %>%
  ggplot(aes(x = Group, y = non.synonymous/synonymous, col = Group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_log10(breaks = function(x){seq(from = x[1], to = x[2], length.out = 5)}) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line())
p1
filename <- paste0(args$outdir, "/", args$prefix, ".varsites.dnds.svg")
ggsave(filename, p1, width = 5, height = 4)

#### Plot position

if(args$plot_position){
  cat("Calculating variable samples per position...\n")
  pos.varsites <- dat %>%
    split(.$site_id) %>%
    map_dfr(function(d){
      res <- d %>%
        split(.$Group) %>% map_dfr(function(d){
          nsamples <- nrow(d)
          nfixed <- sum(d$freq == 0 | d$freq == 1)
          return(tibble(site_id = unique(d$site_id),
                        nsamples = nsamples,
                        fixed = nfixed))
        }, .id = "Group")
    }) %>% left_join(info, by = "site_id")
  # pos.varsites
  cat("\tWriting variable samples per position to file...\n")
  filename <- paste0(args$outdir, "/", args$prefix, ".posvarsites.txt")
  write_tsv(pos.varsites, path = filename)
  
  cat("\tPlotting...\n")
  p1 <- ggplot(pos.varsites, aes(x = ref_pos, y = fixed / nsamples)) +
    facet_grid(~ ref_id, space = "free_x", scales = "free_x") +
    # geom_line(aes(color = Group)) +
    geom_point(aes(color = Group), size = 0.05, alpha = 0.05) +
    geom_smooth(aes(color = Group), se = FALSE,
                method = "gam", formula = y ~ s(x, bs = "cs")) +
    scale_y_continuous(limits = c(0,1)) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90),
          axis.line.x.bottom = element_line(),
          axis.line.y.left = element_line())
  # p1
  filename <- paste0(args$outdir, "/", args$prefix, ".posvarsites.propfixed.png")
  ggsave(filename, p1, width = 12, height = 4, dpi = 200)
}