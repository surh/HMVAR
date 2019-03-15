library(tidyverse)

indir <- "../today5/output/"
outdir <- "manhattans/"
dir.create(outdir)

files <- list.files(indir)
for(lmm_file in files){
  # lmm_file <- "Actinomyces_odontolyticus_57475_lmm.results.txt"
  genome <- basename(lmm_file) %>% str_replace("_lmm.results.txt$", "")
  cat(genome, "\n")
  lmm <- read_tsv(file.path(indir,lmm_file),
                  col_types = 'ccnnnnnnnc')
  # lmm
  
  p1 <- ggplot(lmm,
               aes(x = ps, y = -log10(p_lrt.lmm))) +
    facet_grid(~chr, scales = "free_x", space = "free_x") +
    geom_point(aes(color = type)) +
    geom_hline(yintercept = -log10(1e-6), color = "red", size = 2) +
    theme_classic()
  # print(p1)
  filename <- paste0(outdir, "/", genome, ".manhattan.lmm.png")
  ggsave(filename, p1, width = 16, height = 4)
  
  
  p1 <- ggplot(lmm,
               aes(x = ps, y = -log10(p_lrt.lmmpcs))) +
    facet_grid(~chr, scales = "free_x", space = "free_x") +
    geom_point(aes(color = type)) +
    geom_hline(yintercept = -log10(1e-6), color = "red", size = 2) +
    theme_classic()
  # print(p1)
  filename <- paste0(outdir, "/", genome, ".manhattan.lmmpcs.png")
  ggsave(filename, p1, width = 16, height = 4)
}
