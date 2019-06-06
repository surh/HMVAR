# (C) Copyright 2019 Sur Herrera Paredes
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

library(tidyverse)
# library(seqinr)
library(HMVAR)

# setwd("~/micropopgen/exp/2019/today4/")
# setwd("/godot/users/sur/exp/fraserv/2019/today5/")

depth_thres <- 1
freq_thres <- 0.5

keep_last_codon <- TRUE
outdir <- 'results/'

# mkdir <- "../2019-04-02.hmp_mktest_data/Buccal.mucosa/results/"
# genome_dir <- "~/micropopgen/data/genomes/midas_db_v1.2/hmp.subsite/"
# closest_dir <- "../today3/hmp.subsite.closest/"

# mkdir <- "../2019-04-02.hmp_subsite_mktest/Buccal.mucosa/results/"
# genome_dir <- "~/micropopgen/data/midas_db_v1.2/rep_genomes/"
# closest_dir <- "hmp.subsite.closest/"
# map_file <- "../2019-04-02.hmp_subsite_mktest/Buccal.mucosa/map.txt"
# midas_dir <- "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d.5/"
# group <- "Buccal.mucosa"
# missing_as <- 'major'
# filter_sample <- TRUE

mkdir <- "mkres/"
genome_dir <- "genomes/"
closest_dir <- "closest/"
map_file <- "map.txt"
midas_dir <- "midas/"
outdir <- 'results/'
group <- opts[1]
missing_as <- opts[2]
filter_samples <- opts[3]

if(!dir.exists(outdir)){
  dir.create(outdir)
}
map <- read_tsv(map_file,
                col_types = cols(.default = col_character()))
if(as.logical(filter_samples)){
  map <- map %>%
    filter(Group == group)
}
map <- map %>%
  select(sample = ID, Group)


mkres_files <- list.files(mkdir)
for(mkres_file in mkres_files){
  # mkres_file <- mkres_files[6]
  spec <- str_replace(mkres_file, "_mktest.txt$", '')
  cat(spec, "\n")
  mkres_file <- file.path(mkdir, mkres_file)
  mkres <- read_tsv(mkres_file,
                    col_types = cols(.default = col_double(),
                                     gene_id = col_character()))
  mkres <- mkres %>% 
    filter(p.value < 0.1)
  # mkres
  if(nrow(mkres) > 0){
    cat("\tYES\n")
    
    # Read genome data
    genome_fasta <- seqinr::read.fasta(file.path(genome_dir, spec, 'genome.fna.gz'))
    genome_feats <- readr::read_tsv(file.path(genome_dir, spec, 'genome.features.gz'),
                                    col_types = readr::cols(.default = readr::col_character(),
                                                            start = readr::col_integer(),
                                                            end = readr::col_integer()))
    closest <- read_tsv(file.path(closest_dir, paste0(spec, ".closest.txt")), col_names = FALSE,
                        col_types = cols(.default = col_character(),
                                         X2 = col_number(),
                                         X3 = col_number(),
                                         X8 = col_number(),
                                         X9 = col_number(),
                                         X13 = col_number()))
    
    # Select closest
    closest <- closest %>%
      filter(X4 %in% mkres$gene_id)
    
    
    # Select genes
    genes <- union(closest$X10, mkres$gene_id)
    genome_feats <- genome_feats %>%
      filter(gene_id %in% genes)
    
    # Read midas snv data
    Dat <- read_midas_data(midas_dir = file.path(midas_dir,spec),
                           map = map,
                           genes = genes)
    
    dir.create(file.path(outdir, spec))
    # Get every gene
    for(gene in genome_feats$gene_id){
      # gene <- genome_feats$gene_id[1]
      # Get ref sequence
      gene_pos <- genome_feats %>%
        filter(gene_id == gene) %>%
        select(start, scaffold_id, end, strand)
      if(nrow(gene_pos) != 1){
        stop("ERROR: no data on genes")
      }
      gene_seq <- HMVAR:::get_sequence_fragment(genome_fasta,
                                                ref_id = gene_pos$scaffold_id,
                                                start = gene_pos$start,
                                                end = gene_pos$end)
      
      # Get snp data
      all_sites <- gene_midas_data(Dat = Dat, gene = gene,
                                   depth_thres = depth_thres,
                                   freq_thres = freq_thres)
      # Set relative positions
      all_sites <- all_sites %>%
        mutate(ref_pos = ref_pos - gene_pos$start + 1)
      
      aln <- gene_snv_aln(snps = all_sites, seq = gene_seq,
                          missing_as = missing_as,
                          strand = gene_pos$strand,
                          keep_last_codon = keep_last_codon)
      
      
      filename <- file.path(outdir,spec, paste0(gene, '.aln.fasta'))
      if(file.exists(filename))
        stop(paste("ERROR: file", filename, "exists"), call. = TRUE)
      seqinr::write.fasta(aln$seq, aln$nam, file.out = filename)
    }
    
  }
  
}
# 
# map <- read_tsv('map.txt') %>%
#   select(sample = ID, Group = Group)
# map
# Dat <- read_midas_data(midas_dir = "midas_output_small/", map = map)
# 
# 
# genome_fasta <- seqinr::read.fasta('genome/genome.fna.gz')
# genome_feats <- readr::read_tsv('genome/genome.features.gz',
#                                 col_types = readr::cols(.default = readr::col_character(),
#                                                         start = readr::col_integer(),
#                                                         end = readr::col_integer()))
# 
# genes <- unique(Dat$info$gene_id)
# # gene <- genes[1]
# # gene
# 
# for(gene in genes){
#   # Get ref sequence
#   gene_pos <- genome_feats %>%
#     filter(gene_id == gene) %>%
#     select(start, scaffold_id, end, strand)
#   if(nrow(gene_pos) != 1){
#     stop("ERROR: no data on genes")
#   }
#   gene_seq <- get_sequence_fragment(genome_fasta,
#                                     ref_id = gene_pos$scaffold_id,
#                                     start = gene_pos$start,
#                                     end = gene_pos$end)
#   
#   # Get snp data
#   all_sites <- gene_midas_data(Dat = Dat, gene = gene,
#                                depth_thres = depth_thres,
#                                freq_thres = freq_thres)
#   # Set relative positions
#   all_sites <- all_sites %>%
#     mutate(ref_pos = ref_pos - gene_pos$start + 1)
#   
#   aln <- gene_snv_aln(snps = all_sites, seq = gene_seq,
#                       missing_as = missing_as,
#                       strand = gene_pos$strand,
#                       keep_last_codon = keep_last_codon)
#   
#   
#   filename <- file.path(outdir, paste0(gene, '_aln.fasta'))
#   seqinr::write.fasta(aln$seq, aln$nam, file.out = filename)
# }