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


library(HMVAR)
library(tidyverse)
######################
gene_enrichment <- function(dat, m, n){
  q <- sum(dat$P < 0.01)
  k <- nrow(dat)
  pval <- phyper(q = q - 1, 
                 m = m,
                 n = n,
                 k = k,
                 lower.tail = FALSE )
  
  tibble(gene_id = unique(dat$gene_id),
         ref_id = unique(dat$ref_id),
         mid = round(mean(dat$ref_pos)),
         n.snps = k,
         n.sig = q,
         OR = (q / k) / (m / (m+n)),
         p.value = min(pval, 1 - pval) * 2)
}


match_vmwa_info <- function(vmwa_file, info_file, name, outfile, outdir.plots, threshold = 0.01){
  # Read and process data
  vmwares <- read_tsv(vmwa_file)
  vmwares <- vmwares %>% filter(!is.na(beta)) %>%
    filter(!is.na(SNP)) %>%
    arrange(SNP)
  # vmwares
  
  info <- read_tsv(info_file, col_types = 'icncccnnnnnccccc', na = 'NA')
  info <- info %>%
    select(SNP = site_id, everything(), -starts_with("count_"), -snp_type, -ends_with("_allele"), -amino_acids) %>%
    filter(SNP %in% vmwares$SNP)
  # info
  
  # Match info and vmwa
  vmwares <- vmwares %>% inner_join(info, by = "SNP")
  # vmwares
  
  # Calculate gene stats
  gene_aggregate <- vmwares %>% split(.$gene_id) %>%
    map_dfr(gene_enrichment, m = sum(vmwares$P < threshold), n = sum(vmwares$P >= threshold)) %>%
    mutate(q.value = p.adjust(p.value, 'fdr'))
  # gene_aggregate %>% arrange(q.value) %>% filter(OR > 1) %>% head(20)
  write_tsv(gene_aggregate, outfile)
  
  ### Plot
  # p1 <- ggplot(vmwares, aes(x = ref_pos)) +
  #   facet_grid(~ref_id, scales = "free_x", space = "free_x") +
  #   geom_point(aes(y = -log10(P), col = P < threshold)) +
  #   theme_blackbox() +
  #   theme(strip.text = element_text(angle = 90),
  #         axis.text.x = element_blank(),
  #         panel.background = element_rect(color = NA))
  # p1
  
  p1 <- ggplot(gene_aggregate, aes(x = p.value)) +
    geom_histogram(bins = 20) +
    theme_blackbox()
  filename <- paste0(outdir.plots, "/", name, "_pval.hist.png")
  ggsave(filename, p1, width = 6, height = 4, dpi = 300)
  
  filename <- paste0(outdir.plots, "/", name, "_pval.qqplot.png")
  png(filename, width = 6, height = 6, units = 'in', res = 300)
  pval_qqplot(gene_aggregate$p.value)
  dev.off()
  
  p1 <- ggplot(gene_aggregate, aes(x = mid, y = -log10(p.value) * sign(log(OR)))) +
    facet_grid(~ ref_id, scales = "free_x", space = "free_x") +
    geom_point(aes(color = q.value < threshold)) +
    theme(strip.text = element_text(angle = 90),
          axis.text.x = element_blank(),
          panel.background = element_rect(color = NA))
  filename <- paste0(outdir.plots, "/", name, "_gen.manhat.png")
  ggsave(filename, p1, width = 25, height = 5, dpi = 300)
}
#########################

# args <- list(dir.vmwa = "/godot/users/sur/exp/fraserv/2018/today4/qin2012.vmwa/",
#              dir.info = "/godot/shared_data/metagenomes/qin2012/midas/merge/2018-09-03.merged.snps/",
#              threshold = 0.01,
#              outdir.plots = "qin2012.vmwa.pvals.plots",
#              outdir.genes = "qin2012.vmwa.pvals.genes")

args <- list(dir.vmwa = "/godot/users/sur/exp/fraserv/2018/today4/hmp.vmwa/",
             dir.info = "/godot/shared_data/metagenomes/hmp/midas/merge/2018-02-07.merge.snps.d1/",
             threshold = 0.01,
             outdir.plots = "hmp.vmwa.pvals.plots",
             outdir.genes = "hmp.vmwa.pvals.genes")


if(!dir.exists(args$outdir.plots))
  dir.create(args$outdir.plots)

if(!dir.exists(args$outdir.genes))
  dir.create(args$outdir.genes)


files <- list.files(args$dir.vmwa)
files

for (f in files){
  name <- str_replace(string = f, pattern = "_associations.txt$", replacement = "")
  cat("\t", name, "\n")
  
  # vmwa_file <- "Porphyromonas_sp_57899_associations.txt"
  # info_file <- "snps_info.txt"
  # gff_file <- "genome.features.gz"
  # outfile <- "Porphyromonas_sp_57899_vmwa.genes.txt"
  # name <- "Porphyromonas_sp_57899"
  # outdir.plots <- "gene_aggregate_plots"
  
  vmwa_file <- paste0(args$dir.vmwa, "/", f)
  info_file <- paste0(args$dir.info, "/", name,"/snps_info.txt")
  outfile <- paste0(args$outdir.genes, "/" , name, "_vmwa.genes.txt")
  
  match_vmwa_info(vmwa_file = vmwa_file,
                  info_file = info_file,
                  name = name,
                  outfile = outfile,
                  outdir.plots = args$outdir.plots,
                  threshold = args$threshold)
}
