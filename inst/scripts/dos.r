library(tidyverse)

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

