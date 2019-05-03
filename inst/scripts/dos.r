#!/usr/bin/env Rscript

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
library(argparser)
library(HMVAR)

process_arguments <- function(){
  p <- arg_parser(paste("Calculate DoS statistics on a single genome or a
                        set of pre-processed genomes."))
  
  # Positional arguments
  p <- add_argument(p, "input",
                    help = paste("Input can be either a single file or a directory.",
                                 "If a single file is passed, it hould be a tab-delimited",
                                 "file where each row corresponds to a gene, and it must",
                                 "include columns Dn, Ds, Pn, and Ps corresponding to the",
                                 "values in the McDonald-Kreitman contingency table.\n",
                                 "If a directory is passed, it should either contain a",
                                 "set of tab-delimited files as above, or a set of",
                                 "sub-directories, each one of which is the ouptut from",
                                 "<midas_merge.py snps>. The option --dir_type allows",
                                 "you to specify which type of directory it is, and",
                                 "the option --suffix lets you specify the suffix of",
                                 "the tab-delimited files."),
                    type = "character")
  p <- add_argument(p, "map_file",
                    help = paste0("Mapping file that associates samples ",
                                  "with groups. It must be a tab-delimited ",
                                  "file with headers 'ID' and 'Group'"),
                    type = "character")
  p <- add_argument(p, "outdir",
                    help = paste("This should be the outpu directory for the output"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--dir_type",
                     help = paste("Either 'tabs' or 'midas' to indicate whether",
                                  "tab-delimited files or midas_merge.py directories",
                                  "within the input directory",
                                  "are to be processed. If <input> is a file, then",
                                  "this option will be ignored."),
                     type = "character",
                     default = "files")
  p <- add_argument(p, "--genes",
                    help = paste0("File with list of gene IDs to test. ",
                                  "It must contain one gene id per line, ",
                                  "without header, and the IDs must ",
                                  "correspond to the gene_id column from ",
                                  "the snps_info.txt file if --dir_type=midas", 
                                  "or to the same column in the tab-delimited files",
                                  "if --dir_type=tabs or if <input> was a file.",
                                  "If NULL, all ",
                                  "genes will be tested."),
                    default = NULL,
                    type = "character")
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
                                  " distance from 0 or 1, for a site to be assigned ",
                                  "to the major or minor allele respectively. ",
                                  "It must be a value in [0,1]."),
                    default = 0.5,
                    type = "double")
  p <- add_argument(p, "--focal_group",
                    help = paste("Group to use as a focal group to compare against all other.",
                                 "If passed, the script will convert all values in the Group",
                                 "column of map into 'non.<focal_group>', effectiveley",
                                 "ensuring a dichotomous grouping. If not passed, the script",
                                 "will expect only two levels in the Group column of map and it",
                                 "will fail if more than one group is found."),
                    default = NULL,
                    type = "character")
  p <- add_argument(p, "--suffix",
                    help = paste("The suffix of the tab-delimited files. It will be used",
                                 "to identify which files to process, and also, to",
                                 "determine the output files. Basically if an input",
                                 "file is of the form /path/<prefix><suffix>",
                                 "the ouptut will be of the form <outdir>/<prefix>*").
                    type = "character",
                    default = ".txt")
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  
  return(args)
}

# args <- list(mktest = "../2019-04-02.hmp_mktest_data/Buccal.mucosa/results/",
#              outdir = "results/")
args <- list(mktest = opts[1],
             outdir = opts[2])

# files <- list.files(args$mktest)
# f <- files[2]
# f
# outdir <- args$outdir


dir.create(args$outdir)
dos <- list.files(args$mktest) %>%
  map_dfr(function(f, outdir){
    spec <- str_replace(f, "_mktest.txt$", "")
    cat(spec, "\n")
    f <- file.path(args$mktest, f)
    d <- read_tsv(f,
                  col_types = cols(.default = col_double(),
                                   gene_id = col_character()))
    d <- d %>%
      mutate(DoS = (Dn / (Dn + Ds)) - (Pn / (Pn + Ps)))  %>%
      filter(!is.na(DoS)) %>%
      # mutate(DoS.zscore = DoS / sd(DoS)) %>%
      mutate(DoS.pvalue = 2*(1 - pnorm(q = abs(DoS) / sd(DoS)))) %>%
      mutate(spec = spec)
    # d %>% arrange(DoS.pvalue) %>%
    #   mutate(q.value = p.adjust(DoS.pvalue, 'fdr'))
    filename <- paste0(c(spec, "DoS.table.txt"), collapse = ".")
    filename <- file.path(outdir, filename)
    write_tsv(d, filename)
 
    if(nrow(d) > 1){
      p1 <- ggplot(d, aes(x = DoS)) +
        geom_histogram(bins = 15) +
        ggtitle(label = spec) +
        AMOR::theme_blackbox()
      filename <- paste0(c(spec, "DoS.histogram.png"), collapse = ".")
      filename <- file.path(outdir, filename)
      ggsave(filename, p1, width = 6, height = 5, dpi = 200)
      
      p1 <- ggplot(d, aes(x = DoS.pvalue)) +
        geom_histogram(bins = 20) +
        ggtitle(label = spec) +
        AMOR::theme_blackbox()
      filename <- paste0(c(spec, "DoS.pval.histogram.png"), collapse = ".")
      filename <- file.path(outdir, filename)
      ggsave(filename, p1, width = 6, height = 5, dpi = 200)
      
      filename <- paste0(c(spec, "DoS.pval.qqplot.png"), collapse = ".")
      filename <- file.path(outdir, filename)
      png(filename = filename, width = 6, height = 5, units = "in", res = 200)
      HMVAR::pval_qqplot(d$DoS.pvalue)
      dev.off()
    }
    
    d
  }, outdir = args$outdir, .id = "file")

dos
dos %>% arrange(desc(DoS))
dos %>% arrange(DoS.pvalue)
write_tsv(dos, "summary.dos.txt")

p1 <- ggplot(dos, aes(x = DoS)) +
  geom_histogram(bins = 15) +
  AMOR::theme_blackbox()
p1
ggsave("summary.dos.hist.png", width = 6, height = 5, dpi = 200)

p1 <- ggplot(dos, aes(x = DoS)) +
  geom_density(aes(fill = spec), alpha = 0.3) +
  scale_fill_discrete(guide =  FALSE) +
  AMOR::theme_blackbox()
p1
ggsave("summary.dos.specdensity.png", width = 6, height = 5, dpi = 200)

p1 <- ggplot(dos, aes(x = DoS.pvalue)) +
  geom_histogram(bins = 20) +
  AMOR::theme_blackbox()
p1
ggsave("summary.dospval.hist.png", width = 6, height = 5, dpi = 200)

png("summary.dos.pvalqqplot.png", width = 6, height = 6, units = "in", res = 200)
HMVAR::pval_qqplot(dos$DoS.pvalue)
dev.off()

