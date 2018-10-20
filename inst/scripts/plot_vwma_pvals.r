library(HMVAR)
library(tidyverse)

# args <- list(dir = "~/micropopgen/exp/2018/2018-10-17.vmwa_mktest_data_qin2012_hmp/qin2012.vmwa.DvsND/",
#              threshold = 0.01,
#              outdir = "qin2012.vmwa.pvals.plots",
#              outfile = "qin2012.vmwa.pvals.summary.txt")

args <- list(dir = "~/micropopgen/exp/2018/2018-10-17.vmwa_mktest_data_qin2012_hmp/hmp.vmwa.SPvsTD/",
             threshold = 0.01,
             outdir = "hmp.vmwa.pvals.plots",
             outfile = "hmp.vmwa.pvals.summary.txt")


if(!dir.exists(args$outdir))
  dir.create(args$outdir)

files <- list.files(args$dir)
files

Res <- NULL
for (f in files){
  # f <- files[47]
  name <- str_replace(string = f, pattern = "_associations.txt$", replacement = "")
  cat("\t", name, "\n")
  
  vmwares <- read_tsv(paste0(args$dir, "/", f))
  vmwares <- vmwares %>% filter(!is.na(beta)) %>%
    filter(!is.na(SNP))
  vmwares
  
  q.vals <- qvalue::qvalue(vmwares$p.value)
  print(summary(q.vals))
  
  res <- tibble(name = name, ntests = nrow(vmwares),
                t.Pi0 = q.vals$pi0, t.sig = sum(vmwares$p.value < args$threshold))
  
  Q.vals <- qvalue::qvalue(vmwares$P)
  print(summary(Q.vals))
  
  res <- res %>% mutate(P.Pi0 = Q.vals$pi0, P.sig = sum(vmwares$P < args$threshold))
  Res <- Res %>% bind_rows(res)
  
  # plot(qvalue::qvalue(vmwares$p.value))
  
  # p1 <- ggplot(vmwares, aes(x = p.value)) +
  #   geom_histogram(bins = 20) +
  #   theme_blackbox()
  # p1
  # 
  # pval_qqplot(vmwares$p.value)
  # 
  # p1 <- ggplot(vmwares, aes(x = beta, y = -log10(p.value))) +
  #   geom_point(aes(col = p.value < args$threshold)) +
  #   scale_x_continuous(limits = c(-1,1))
  # p1
  # 
  # p1 <- ggplot(vmwares, aes(x = SNP)) +
  #   geom_point(aes(y = -log10(p.value), col = p.value < args$threshold)) +
  #   theme_blackbox() +
  #   theme(axis.text.x = element_blank(),
  #         panel.background = element_rect(color = NA))
  # p1
  
  ### Compare
  p1 <- ggplot(vmwares, aes(x = -log10(p.value), y = -log10(P))) +
    geom_point() +
    theme_blackbox()
  filename <- paste0(args$outdir, "/", name, "_compare.pvals.png")
  ggsave(filename, p1, width = 6, height = 6, dpi = 300)
  
  ###
  p1 <- ggplot(vmwares, aes(x = P)) +
    geom_histogram(bins = 20) +
    theme_blackbox()
  filename <- paste0(args$outdir, "/", name, "_pval.hist.png")
  ggsave(filename, p1, width = 6, height = 4, dpi = 300)
  
  filename <- paste0(args$outdir, "/", name, "_pval.qqplot.png")
  png(filename = filename, units = 'in', width = 6, height = 6, res = 300)
  pval_qqplot(vmwares$P)
  dev.off()
  
  p1 <- ggplot(vmwares, aes(x = beta, y = -log10(P))) +
    geom_point(aes(col = p.value < args$threshold)) +
    scale_x_continuous(limits = c(-1, 1)) +
    theme_blackbox()
  filename <- paste0(args$outdir, "/", name, "_volcano.png")
  ggsave(filename, p1, width = 6, height = 4, dpi = 300)
  
  # p1 <- ggplot(vmwares, aes(x = SNP)) +
  #   geom_point(aes(y = -log10(p.value), col = P < args$threshold)) +
  #   theme_blackbox() +
  #   theme(axis.text.x = element_blank(),
  #         panel.background = element_rect(color = NA))
  # p1
}

Res <- Res %>% mutate(exp = ntests * args$threshold,
                      t.extra = t.sig - exp,
                      P.extra = P.sig - exp)
Res
write_tsv(Res, args$outfile)
