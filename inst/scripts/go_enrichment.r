Sys.time()

file <- "~/micropopgen/exp/2018/today/signifincant_mktest"
Tab <- read.table(file, header = TRUE, sep = "\t")
head(Tab)

go.dir <- "/home/sur/micropopgen/data/genomes/midas_db_v1.2/GO.annotations/"
go.dir <- "/home/sur/micropopgen/exp/2018/today/GO.annotations/"


comparisons <- levels(interaction(Tab$A, Tab$B, sep = "_", drop = TRUE))

RES <- NULL
for(c in comparisons){
  c <- comparisons[1]
  cat(c, "\n")
  
  # Get data from comparison
  dat <- droplevels(subset(Tab, A == strsplit(x = c, split = "_")[[1]][1] & B == strsplit(x = c, split = "_")[[1]][2]))
  specs <- levels(dat$Species)
  
  ANNOTS <- NULL
  for(s in specs){
    s <- specs[1]
    cat("\t", s, "\n")
    
    go_file <- paste(go.dir, "/", s, ".GO.txt", sep = "")
    cat("\t", go_file, "\n")
    
    go <- read.table(go_file, header = TRUE)
    
    ANNOTS <- rbind(ANNOTS, go)
    rm(go)
  }
  
  bg_counts <- table(ANNOTS$Annotation)
  
  sig_annots <- droplevels(ANNOTS[ ANNOTS$Gene %in% levels(dat$gene), ])
  sig_counts <- table(sig_annots$Annotation)
  
  for(i in 1:length(sig_counts)){
    # i <- 1
    a <- names(sig_counts)[i]
    cat("\t", a, "\n")
  
    q <- sig_counts[i]
    m <- bg_counts[a]
    n <- nrow(ANNOTS) - m
    k <- nrow(sig_annots)
    
    pval <- 1 - phyper(q = q - 1, m = m, n = n, k = k)
    
    res <- data.frame(comparison = c, annotation = a, nsig = q, nbg = m, pval = pval,
                      row.names = NULL)
    RES <- rbind(RES, res)
  }
}