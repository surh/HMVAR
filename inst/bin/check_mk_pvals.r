library(HMAR)
library(ggplot2)
library(qvalue)

# setwd("~/micropopgen/exp/2018/today3/")
# devtools::document("~/micropopgen/src/HMAR/")


which <- 'ratio'
plot <- TRUE


dir <- opts[1]
outdir <- opts[2]
pattern <- opts[3]

if(is.na(outdir))
  stop("ERROR: must provide outfile name")

if(is.na(pattern))
  pattern <- "^mk_results"

cat("dir:", dir, "\n")
cat("outfile:", outfile, "\n")
cat("pattern:", pattern, "\n")

dir.create(outdir)

species_dirs <- list.dirs(dir, recursive = FALSE)

Res <- NULL
for(d in species_dirs){
  species <- basename(d)
  cat(species, "\n")
  
  files <- get_mk_results_files(d, pattern)
  
  file <- "mktest/Streptococcus_mitis_60474/mk_results.Buccal.mucosa_Supragingival.plaque.txt"
}

  

check_pvals_in_file <- function(file, which, plot=TRUE){
  Tab <- read.table(file, header = TRUE, sep = "\t")
  default <- c("gene", "contig", "start", "end", "Dn", "Ds", "Pn", "Ps")
  
  # Get pval columns
  if(which == "all"){
    tests <- setdiff(colnames(Tab), default)
    tests <- tests[ grep(pattern = ".pval", x = tests, invert = TRUE) ]
    tests <- tests[ grep(pattern = ".perm", x = tests, invert = TRUE) ]
  }else{
    tests <- which
  }
  
  Res <- list()
  for(t in tests){
    # t <- tests[8]
    res <- check_pvalues(Tab[,t], Tab[,paste(t,".pval",sep="")], plot = plot)
    Res[[t]] <- res
  }
  
  return(Res)
}


check_pvalues <- function(estimates, pvals, plot = TRUE){
  pvals <- pvals[ !is.na(estimates) ]
  qvals <- qvalue::qvalue(pvals)
  # qvals.sum <- summary(qvals)
  pi0 <- qvals$pi0
  
  p1 <- NULL
  if(plot){
    p1 <- ggplot(data.frame(pvals),aes(x = pvals)) +
      geom_histogram(bins = 20) +
      AMOR::theme_blackbox
  }
  
  res <- list(pi0 = pi0, p1 = p1)
  
  return(res)
}


