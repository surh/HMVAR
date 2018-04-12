library(HMAR)
library(ggplot2)
library(qvalue)

# This must be run with stitch_file.r
# setwd("~/micropopgen/exp/2018/today3/")
# devtools::document("~/micropopgen/src/HMAR/")

dir <- opts[1]
outfile <- opts[2]
pattern <- opts[3]

# dir <- "~/micropopgen/exp/2018/today3/mktest/"
# outfile <- "results.txt"
# pattern <- NA

if(is.na(outfile))
  stop("ERROR: must provide outfile name")

if(is.na(pattern))
  pattern <- "^mk_results"

species_dirs <- list.dirs(dir, recursive = FALSE)

Res <- NULL
for(d in species_dirs){
  species <- basename(d)
  cat(species, "\n")
  
  files <- get_mk_results_files(d, pattern)
  
  if(length(files) > 0){
    stats <- sapply(files, calculate_genome_wide_ni)
    
    res <- data.frame(Species = species, NI = stats[1,], alpha = stats[2,], row.names = NULL)
    Res <- rbind(Res, res)
  }else{
    cat("\t", species, " had no result files.\n")
  }
}

write.table(Res, outfile, sep = "\t", col.names = TRUE, row.names = FALSE)