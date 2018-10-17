library(HMAR)
# library(ggplot2)
# library(qvalue)
# setwd("~/micropopgen/exp/2018/today3/")

Sys.time()

# setwd("~/micropopgen/exp/2018/today/")

dir <- "/home/sur/micropopgen/exp/2018/today/mktest_fraserv/"
outdir <- "all_mktests/"
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
    
    # Get CDS file
    cat("\tGetting CDS")
    cds_file <- paste("~/micropopgen/data/genomes/midas_db_v1.2/tested_species/FAA/", species,
                      ".CDS.faa", sep = "")
    comp_cds <- paste(outdir,"/",paste(c(comparison,"CDS.faa"), collapse = "_"), sep = "")
    cmd <- paste("cat ", cds_file, " >> ", comp_cds, sep = "")
    system(cmd)
    
    cat("\tratiodos\n")
    Tab$DoS <- (Tab$Dn / (Tab$Dn + Tab$Ds)) - (Tab$Pn / (Tab$Pn + Tab$Ps))
    Tab$ratio.qval <- p.adjust(Tab$ratio.pval, 'fdr')

    Tab$Species <- species
    Tab$A <- comparison[1]
    Tab$B <- comparison[2]
    
    RES <- rbind(RES, Tab)
  }
}
write.table(RES, file = "all_mktest", sep = "\t", quote = FALSE)