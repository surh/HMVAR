#!/usr/bin/env Rscript

# (C) Copyright 2018-2019 Sur Herrera Paredes
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
midas_data$info <- determine_substitution_type(midas_data$info, clean = FALSE)
midas_data$info

# Reformat data
dat <- match_freq_and_depth(freq = midas_data$freq,
                            depth = midas_data$depth,
                            info = midas_data$info,
                            map = map,
                            depth_thres = args$depth_thres)

# For each site determine if it is homogeneous or heterogenous
dat <- determine_sample_dist(dat)
dat

# Count number of sites and number of variable per sample
cat("Calcuating variable sites...\n")
varsites <- dat %>%
  filter(depth >= 1) %>%
  split(.$sample) %>%
  map_dfr(function(d){
    nsites <- nrow(d)
    # sites <- d$freq
    nhomogeneous <- sum(d$sample_dist == "homogeneous", na.rm = TRUE)
    nheterogeneous <- sum(d$sample_dist == "heterogeneous", na.rm = TRUE)
    ntransitions <- sum(d$substitution == "transition", na.rm = TRUE)
    ntransvertions <- sum(d$substitution == "transversion", na.rm = TRUE)
    nsynonymous <- sum(d$snp_effect == "synonymous", na.rm = TRUE)
    nnonsynonymous <- sum(d$snp_effect == "non-synonymous", na.rm = TRUE)
    return(tibble(nsites = nsites,
                  n_homogeneous = nhomogeneous,
                  n_heterogeneous = nheterogeneous,
                  n_transitions = ntransitions,
                  n_transversions = ntransvertions,
                  n_synonymous = nsynonymous,
                  n_non.synonymous = nnonsynonymous))
  }, .id = "sample") %>%
  inner_join(map, by = "sample")
varsites

# Plot synonymous & non-synonymous
p1 <- plotgg_stacked_columns(dat = varsites,
                       x = "sample",
                       columns = c("n_synonymous", "n_non.synonymous"),
                       facet_formula = ~ Group,
                       gather_key = "type",
                       gather_value = "nloci")
p1

# Plot homogeneous & heterogeneous
p1 <- plotgg_stacked_columns(dat = varsites,
                             x = "sample",
                             columns = c("n_homogeneous", "n_heterogeneous"),
                             facet_formula = ~ Group,
                             gather_key = "type",
                             gather_value = "nloci")
p1


#### Plot position
#' Size of groups
#' 
#' Internal function
#'
#' @param d A data frame or tibble
#' @param columns Names of columns to group
#'
#' @return A tibble
#' @importFrom magrittr %>%
group_size <- function(d, columns){
  d %>% split(.[,columns]) %>%
    purrr:::map_int(~nrow(.)) %>%
    t %>%
    tibble::as_tibble()
}

group <- c('Group')
group
cat("Calculating variable samples per position...\n")
dat
D <- dat %>%
  split(.$site_id)
d <- D[[1]]
d



dat %>%
  split(.$site_id) %>%
  map_dfr(group_size, columns = 'sample_dist', .id = 'site_id')




variable_dist_per_site <- function(dat, variable, group = NULL){
  group <- "Group"
  variable <- 'sample_dist'
  
  if(!all(c("site_id", "ref_id", "ref_pos", group, variable) %in% colnames(dat))){
    stop("ERROR: missing columns in dat")
  }
  
  if(is.null(group)){
    res <- dat %>%
      split(.$site_id) %>%
      map_dfr(group_size, columns = 'sample_dist', .id = 'site_id')
  }else{
    res <- dat %>%
      split(.$site_id) %>%
      map_dfr(function(d, column){
        d %>% split(.[,column]) %>%
          map_dfr(group_size, columns = 'sample_dist', .id = column)},
        column = 'Group',
        .id = 'site_id')
  }
  dat
  res
  res <- dat %>% select(site_id, ref_id, ref_pos) %>% filter(!duplicated(.)) %>% full_join(res, by = "site_id")
}





group_size(d = d, columns = 'sample_dist')


pos.varsites <- dat %>%
  split(.$site_id) %>%
  map_dfr(function(d){
    res <- d %>%
      split(.$Group) %>% map_dfr(function(d){
        nsamples <- nrow(d)
        n_homogeneous <- sum(d$sample_dist == "homogeneous")
        return(tibble(site_id = unique(d$site_id),
                      nsamples = nsamples,
                      n_homogeneous = n_homogeneous))
      }, .id = "Group")
  }) %>% left_join(midas_data$info, by = "site_id")
# pos.varsites
cat("\tWriting variable samples per position to file...\n")
filename <- paste0(args$outdir, "/", args$prefix, ".posvarsites.txt")
write_tsv(pos.varsites, path = filename)

cat("\tPlotting...\n")
p1 <- ggplot(pos.varsites, aes(x = ref_pos, y = n_homogeneous / nsamples)) +
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
p1
filename <- paste0(args$outdir, "/", args$prefix, ".posvarsites.propfixed.png")
ggsave(filename, p1, width = 12, height = 4, dpi = 200)


# # Prepare output dir
# if(!dir.exists(args$outdir)){
#   cat("Preparing output directory...\n")
#   dir.create(args$outdir)
# }