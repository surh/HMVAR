library(HMAR)
library(ggplot2)
library(qvalue)

# setwd("~/micropopgen/exp/2018/today3/")
# devtools::document("~/micropopgen/src/HMAR/")


dir <- opts[1]
outdir <- opts[2]
which <- opts[3]
plot <- as.logical(opts[4])
pattern <- opts[5]

# dir <- "/home/sur/micropopgen/exp/2018/today3/mktest/"
# outdir <- "pvals_res/"
# which <- "all"
# plot <- 1
# pattern <- "^mk_results"

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

Pi0 <- NULL
for(d in species_dirs){
  # d <- species_dirs[7]
  # d <- species_dirs[2]
  species <- basename(d)
  cat(species, "\n")
  
  files <- get_mk_results_files(d, pattern)
  # print(files)
  
  for(f in files){
    # f <- "mktest/Streptococcus_mitis_60474/mk_results.Buccal.mucosa_Supragingival.plaque.txt"
    # f <- files[1]
    cat("\t", f, "\n")
    
    # Res <- lapply(files, check_pvals_in_file, which = which, plot = plot)
    res <- check_pvals_in_file(f, which = which, plot = plot)
    if(length(res) > 0){
      pi0 <- sapply(res, function(x) x$pi0)
      Pi0 <- rbind(Pi0, data.frame(Species = species, t(pi0), file = f))
      
      if(plot){
        prefix <- basename(f)
        # prefix <- sub(pattern = pattern, replacement = "", x = prefix)
        prefix <- paste(species, ".", prefix, sep = "")
        prefix <- sub(pattern = ".txt$", replacement = "", x = prefix)
        prefix <- paste(outdir, "/", prefix, sep = "")
        prefix <- sub(pattern = ".$", replacement = "", x = prefix)
        
        sapply(names(res), function(name, prefix, res){
          p1 <- res[[name]]$p1
          p1 <- p1 + ggtitle(label = name)
          outfile <- paste(prefix, ".", name,".pvals.svg", sep = "")
          # print(outfile)
          ggplot2::ggsave(outfile, p1, width = 6, height = 4)
          return(NULL)
        }, prefix = prefix, res = res)
        
      }
      
      
    }
    rm(res)
  }
}
outfile <- paste(outdir, "/Pi0_summmary.txt", sep = "")
write.table(Pi0, outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
