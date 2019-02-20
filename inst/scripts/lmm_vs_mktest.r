library(tidyverse)

fisher <- function(gene, contig, start, end, Dn, Ds, Pn, Ps){
  mat <- matrix(c(Dn, Ds, Pn, Ps), ncol = 2)
  res <- fisher.test(mat)
  tibble(gene = gene, contig = contig, start = start,
         end = end, Dn = Dn, Ds = Ds, Pn = Pn, Ps = Ps,
         OR = res$estimate,
         p.value = res$p.value,
         mkratio = ((Dn+1)/(Ds+1)) / ((Pn+1)/(Ps+1)))
}


# lmmres <- read_tsv("../2019-01-28.hmp_metawas/output/Actinomyces_viscosus_57672_lmm.assoc.txt",
#                    col_types = 'ccnnnnn')

mkdir <- "/godot/users/sur/exp/sherlock2/2018-10-11.hmp_mktest/results/"
lmmclosedir <- "/godot/users/sur/exp/fraserv/2019/2019-01-28.hmp_metawas_closest/output/"
genomes_file <- "genomes.txt"
outdir <- "output/"

genomes <- read_tsv(genomes_file, col_names = FALSE, col_types = 'c')
genomes <- genomes$X1

dir.create(outdir)


Res <- NULL
for(genome in genomes){
  # genome <- "Actinomyces_viscosus_57672"
  cat(genome,"\n")
  mkfile <- paste0(mkdir, "/", genome, "_mktest.txt")
  lmmclosefile <- paste0(lmmclosedir, "/", genome, ".closest")
  
  if(!file.exists(mkfile) || !file.exists(lmmclosefile)){
    cat("\tskipping\n")
    next
  }
  
  mkres <- read_tsv(mkfile,
                    col_types = 'ccnnnnnnnnn',
                    na = 'nan')
  lmmclose <- read_tsv(lmmclosefile,
                       col_types = 'cnncnncn', col_names = FALSE)
  if(nrow(lmmclose) == 0){
    cat("\tskipping")
    next
  }
  
  mkres <- mkres %>%
    select(gene, contig, start, end, Dn, Ds, Pn, Ps) %>%
    pmap_dfr(fisher)
  mkres %>% arrange(desc(mkratio))
  
  # Plot
  mkres <- mkres %>% mutate(metawas = gene %in% lmmclose$X7)
  mkres <- mkres %>% filter(Ds > 0 & Pn > 0)
  
  if(sum(mkres$metawas == TRUE) == 0){
    cat("\tSkipping\n")
    next
  }
  
  p1 <- ggplot(mkres %>% filter(Ds > 0 & Pn > 0),
               aes(x = log2(mkratio), y = -log10(p.value))) +
    geom_point(aes(col = metawas)) +
    AMOR::theme_blackbox()
  p1
  filename <- paste0(outdir, "/", genome, ".volcano.png")
  ggsave(filename, p1, width = 6, height = 4)
  
  p1 <- ggplot(mkres %>% filter(Ds > 0 & Pn > 0),
               aes(x = metawas, y = mkratio)) +
    # geom_boxplot(outlier.color = NA) +
    geom_violin() +
    geom_point(aes(col = metawas), position = position_jitter(width = 0.3)) +
    scale_y_log10() +
    AMOR::theme_blackbox()
  p1
  filename <- paste0(outdir, "/", genome, ".violin.png")
  ggsave(filename, p1, width = 6, height = 4)
  
  m1 <- lm(log10(mkratio) ~ metawas, data = mkres %>% filter(Ds > 0 & Pn > 0))
  res <- summary(m1)$coefficients[2,]
  
  Res <- bind_rows(Res,
                   tibble(genome = genome,
                          estimate = res[1],
                          SE = res[2],
                          t.value = res[3],
                          p.value = res[4]))
}
Res
write_tsv(Res, "hmp_mkres_vs_lmmclosest.txt")
