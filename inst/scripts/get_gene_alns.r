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
library(seqinr)
library(HMVAR)

setwd("~/micropopgen/exp/2019/today2/")

depth_thres <- 1
freq_thres <- 0.5
missing_as <- 'major'
keep_last_codon <- TRUE


map <- read_tsv('map.txt') %>%
  select(sample = ID, Group = Group)
map
Dat <- read_midas_data(midas_dir = "midas_output_small/", map = map)


genome_fasta <- seqinr::read.fasta('genome/genome.fna.gz')
genome_feats <- readr::read_tsv('genome/genome.features.gz',
                                col_types = readr::cols(.default = readr::col_character(),
                                                        start = readr::col_integer(),
                                                        end = readr::col_integer()))

genes <- unique(Dat$info$gene_id)
# gene <- genes[1]
# gene

for(gene in genes){
  # Get ref sequence
  gene_pos <- genome_feats %>%
    filter(gene_id == gene) %>%
    select(start, scaffold_id, end, strand)
  if(nrow(gene_pos) != 1){
    stop("ERROR: no data on genes")
  }
  gene_seq <- get_sequence_fragment(genome_fasta,
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
  
  
  filename <- paste0(gene, '.aln.fasta')
  seqinr::write.fasta(aln$seq, aln$nam, file.out = filename)
}