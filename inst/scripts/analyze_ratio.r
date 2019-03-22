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

library(HMAR)
library(ggplot2)
library(qvalue)

Sys.time()

setwd("~/micropopgen/exp/2018/today/")
# devtools::document("~/micropopgen/src/HMAR/")

# dir <- opts[1]
# outdir <- opts[2]
# # which <- opts[3]
# plot <- as.logical(opts[4])
# pattern <- opts[5]

dir <- "/home/sur/micropopgen/exp/2018/today/mktest_fraserv/"
outdir <- "pvals_res/"
which <- "ratio"
plot <- 1
pattern <- "^mk_results"

if(is.na(outdir))
  stop("ERROR: must provide outfile name")

if(is.na(pattern))
  pattern <- "^mk_results"

cat("dir:", dir, "\n")
cat("outdir:", outdir, "\n")
cat("which:", which, "\n")
cat("plot:", plot, "\n")
cat("pattern:", pattern, "\n")

dir.create(outdir)

species_dirs <- list.dirs(dir, recursive = FALSE)

min.nperm <- 75
qval_thres <- 0.1

RES <- NULL
for(d in species_dirs){
  # d <- species_dirs[1]
  # d <- species_dirs[2]
  species <- basename(d)
  cat(species, "\n")
  
  files <- get_mk_results_files(d, pattern)
  print(files)
  
  for(f in files){
    # f <- "mktest/Streptococcus_mitis_60474/mk_results.Buccal.mucosa_Supragingival.plaque.txt"
    # f <- files[1]
    cat("\t", f, "\n")
    
    comparison <- strsplit(x = sub(pattern = ".txt$", replacement = "",
                                   x = sub(pattern = "^mk_results.", replacement = "",
                                           x = basename(f))),
                           split = "_")[[1]]
    
    Tab <- read.table(f, header = TRUE, sep = "\t")
    Tab <- Tab[,c("gene","contig","start","end","Dn","Ds","Pn","Ps","ratio","ratio.pval","ratio.perm")]
    head(Tab)
    
    cat("\tsubsetting\n")
    # Remove not defined
    Tab <- subset(Tab, Ds > 0 & Pn > 0)
    
    # Remove few permutations
    Tab <- subset(Tab, ratio.perm >= min.nperm)
    
    # Remove low ratio, becaus ethe way permutation was done
    # Tab <- subset(Tab, ratio < 1)
    
    if (nrow(Tab) == 0){
      next
    }
    
    cat("\thistogram\n")
    outfile <- paste(outdir,"/", species, ".", comparison[1], "_", comparison[2], "_ratio.pval.hist.svg", sep = "")
    # summary(qvalue::qvalue(Tab$ratio.pval))
    p1 <- ggplot(Tab, aes(x = ratio.pval)) +
      geom_histogram(bins = 20) +
      AMOR::theme_blackbox
    # p1
    ggsave(outfile, p1, width = 6, height = 4) 
    
    # head(Tab)
    # p1 <- ggplot(Tab, aes(x = log2((Dn * Ps) / (Ds * Pn)), y = -log10(ratio.pval) )) +
    #   geom_point()
    # p1
    
    cat("\tratiodos\n")
    Tab$DoS <- (Tab$Dn / (Tab$Dn + Tab$Ds)) - (Tab$Pn / (Tab$Pn + Tab$Ps))
    # p1 <- ggplot(Tab, aes(x = DoS)) +
    #   geom_histogram(bins = 20) +
    #   AMOR::theme_blackbox
    # p1
    p1 <- ggplot(Tab, aes(x = ratio, y = DoS, col = -log10(ratio.pval))) +
      geom_point() +
      geom_smooth() +
      AMOR::theme_blackbox
    # p1
    outfile <- paste(outdir,"/", species, ".", comparison[1], "_", comparison[2], "_ratio.dos.svg", sep = "")
    ggsave(outfile, p1, width = 6, height = 4) 
    # ks.test(Tab$DoS, y = 'norm')
    

    # FDR
    Tab$ratio.qval <- p.adjust(Tab$ratio.pval, 'fdr')
    # head(Tab[ order(Tab$ratio.qval, decreasing = FALSE), ])
    
    # head(Tab[ order(Tab$ratio.pval, decreasing = FALSE), ], 50)
    # head(Tab[ order(Tab$ratio.pval, decreasing = TRUE), ], 50)
    
    Tab$Species <- species
    Tab$A <- comparison[1]
    Tab$B <- comparison[2]
    
    cat("\tqvals\n")
    # Select passing threshold
    Tab <- subset(Tab, ratio.qval <= qval_thres)
    
    # head(Tab[ order(Tab$ratio, decreasing = TRUE), ], 20)
    RES <- rbind(RES, Tab)
  }
}
write.table(RES, file = "signifincant_mktest", sep = "\t", quote = FALSE)

# p1 <- ggplot(RES, aes(x = ratio, y = DoS, col = -log10(ratio.pval))) +
#   geom_point() +
#   geom_smooth() +
#   AMOR::theme_blackbox
# p1
