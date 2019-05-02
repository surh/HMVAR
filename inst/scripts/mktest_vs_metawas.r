library(tidyverse)

match_metawas_mktest <- function(spec, mktest, genes, outdir){
  # i <- 8
  # spec <- files$spec[i]
  # mktest <- files$mktest[i]
  # genes <- files$genes[i]
  # outdir <- "results/"

  cat(spec, "\n")
  cat("\tRead gene counts...\n")
  Genes <- read_tsv(genes,
                    col_types = cols(.default = col_character(),
                                     start = col_number(),
                                     end = col_number(),
                                     n.variants = col_number(),
                                     n.sig = col_number(),
                                     n.outside = col_number()))
  # Genes
  
  cat("\tRead mktest...\n")
  Mktest <- read_tsv(mktest, 
                     col_types = cols(.default = col_double(),
                                      gene_id = col_character())) %>%
    mutate(kaks = (Dn + Pn) / (Ds + Ps),
           mkratio = ((Dn + 1) / (Ds + 1)) / ((Pn + 1) / (Ps + 1)))
  # Mktest
  
  cat("\tJoin...\n")
  Genes <- Genes %>%
    left_join(Mktest, by = "gene_id") %>%
    filter(!is.na(Dn)) %>%
    filter(!is.na(kaks)) %>%
    filter(kaks != Inf)
  
  if(nrow(Genes) == 0){
    return(Genes)
  }
  
  cat("\tPlotting...\n")
  if(sum(Genes$n.sig > 0) > 0){
    test <- wilcox.test(Genes$kaks[Genes$n.sig == 0],
                        Genes$kaks[Genes$n.sig > 0],
                        alternative = "less")
  }else{
    test <- list(p.value = 1)
  }
  p1 <- Genes %>%
    ggplot(aes(x = n.sig > 0, y = kaks)) +
    geom_violin(draw_quantiles = 0.5) +
    geom_point(position = position_jitter(width = 0.2), size = 0.1) +
    scale_y_log10() +
    ggtitle(label = "", subtitle = test$p.value) +
    AMOR::theme_blackbox()
  # p1
  filename <- paste0(c(spec,"nsig_vs_kaks.png"), collapse = ".")
  filename <- file.path(outdir, filename)
  ggsave(filename, p1, width = 6, height = 5, dpi = 200)
  
  if(sum(Genes$n.sig > 0) > 0){
    test <- wilcox.test(Genes$mkratio[Genes$n.sig == 0],
                        Genes$mkratio[Genes$n.sig > 0],
                        alternative = "less")
  }else{
    test <- list(p.value = 1)
  }
  p1 <- Genes %>%
    ggplot(aes(x = n.sig > 0, y = mkratio)) +
    geom_violin(draw_quantiles = 0.5) +
    geom_point(position = position_jitter(width = 0.2), size = 0.1) +
    scale_y_log10() +
    ggtitle(label = "", subtitle = test$p.value) +
    AMOR::theme_blackbox()
  # p1
  filename <- paste0(c(spec,"nsig_vs_mkratio.png"), collapse = ".")
  filename <- file.path(outdir, filename)
  ggsave(filename, p1, width = 6, height = 5, dpi = 200)
  
  p1 <- Genes %>%
    filter(n.sig > 0) %>%
    ggplot(aes(x = n.sig / n.variants, y = kaks)) +
    geom_point() +
    geom_smooth(method = "lm") +
    AMOR::theme_blackbox()
  filename <- paste0(c(spec,"sigprop_vs_kaks.png"), collapse = ".")
  filename <- file.path(outdir, filename)
  ggsave(filename, p1, width = 6, height = 5, dpi = 200)
  
  p1 <- Genes %>%
    filter(n.sig > 0) %>%
    ggplot(aes(x = n.sig / n.variants, y = mkratio)) +
    geom_point() +
    geom_smooth(method = "lm") +
    AMOR::theme_blackbox()
  filename <- paste0(c(spec,"sigprop_vs_mkratio.png"), collapse = ".")
  filename <- file.path(outdir, filename)
  ggsave(filename, p1, width = 6, height = 5, dpi = 200)
  
  Genes$spec <- spec
  
  return(Genes)
}

# args <- list(mkdir = "../2019-04-02.hmp_mktest_data/Buccal.mucosa/results/",
#              metawasdir = "../today3/Buccal.mucosa/enrichments/enrichments/",
#              outdir = "Buccal.mucosa/results",
#              overall_out = "Buccal.mucosa/overall/")
args <- list(mkdir = opts[1],
             metawasdir = opts[2],
             outdir = opts[3],
             overall_out = opts[4])

mktest <- tibble(mktest = file.path(args$mkdir, list.files(args$mkdir)))
mktest <- mktest %>% 
  mutate(spec = str_replace(basename(mktest), pattern = "_mktest.txt$", "")) %>%
  select(spec, mktest)
mktest

metawas <- list.files(args$metawasdir) %>%
  str_subset(pattern = ".gene_metawas_counts.txt$") %>%
  file.path(args$metawasdir, .) %>%
  tibble(spec = str_replace(basename(.), pattern = ".gene_metawas_counts.txt$", ""),
         genes = .)
metawas

files <- mktest %>%
  left_join(metawas, by = "spec") %>%
  filter(!is.na(genes))
rm(mktest, metawas)
files

# spec <- files$spec[1]
# mktest <- files$mktest[1]
# genes <- files$genes[1]
# outdir <- "results/"

dir.create(args$outdir)

Full <- files %>%
  pmap_dfr(match_metawas_mktest,
           outdir = args$outdir)

# Overall results
dir.create(args$overall_out)

test <- wilcox.test(Full$kaks[Full$n.sig == 0],
                    Full$kaks[Full$n.sig > 0],
                    alternative = "two.sided")
p1 <- Full %>%
  ggplot(aes(x = n.sig > 0, y = kaks)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_point(aes(col = spec), 
             position = position_jitter(width = 0.2),
             size = 0.1) +
  scale_color_discrete(guide = FALSE) +
  scale_y_log10() +
  ggtitle(label = "", subtitle = test$p.value) +
  AMOR::theme_blackbox()
p1
filename <- paste0(c("overall","nsig_vs_kaks.png"), collapse = ".")
filename <- file.path(args$overall_out, filename)
ggsave(filename, p1, width = 8, height = 5, dpi = 200)

test <- wilcox.test(Full$mkratio[Full$n.sig == 0],
                    Full$mkratio[Full$n.sig > 0],
                    alternative = "two.sided")
p1 <- Full %>%
  ggplot(aes(x = n.sig > 0, y = mkratio)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_point(aes(col = spec), 
             position = position_jitter(width = 0.2),
             size = 0.1) +
  scale_color_discrete(guide = FALSE) +
  scale_y_log10() +
  ggtitle(label = "", subtitle = test$p.value) +
  AMOR::theme_blackbox()
p1
filename <- paste0(c("overall","nsig_vs_mkratio.png"), collapse = ".")
filename <- file.path(args$overall_out, filename)
ggsave(filename, p1, width = 6, height = 5, dpi = 200)

p1 <- Full %>%
  filter(n.sig > 0) %>%
  ggplot(aes(x = n.sig / n.variants, y = kaks)) +
  geom_point(aes(col = spec)) +
  scale_color_discrete(guide = FALSE) +
  geom_smooth(method = "lm") +
  AMOR::theme_blackbox()
p1
filename <- paste0(c("overall","sigprop_vs_kaks.png"), collapse = ".")
filename <- file.path(args$overall_out, filename)
ggsave(filename, p1, width = 6, height = 5, dpi = 200)

p1 <- Full %>%
  filter(n.sig > 0) %>%
  ggplot(aes(x = n.sig / n.variants, y = mkratio)) +
  geom_point(aes(col = spec)) +
  scale_color_discrete(guide = FALSE) +
  geom_smooth(method = "lm") +
  AMOR::theme_blackbox()
p1
filename <- paste0(c("overall","sigprop_vs_mkratio.png"), collapse = ".")
filename <- file.path(args$overall_out, filename)
ggsave(filename, p1, width = 6, height = 5, dpi = 200)


filename <- paste0(c("overall","mktest_metawas.txt"), collapse = ".")
filename <- file.path(args$overall_out, filename)
write_tsv(Full, filename)

# Select ones that can be tested
Dat <- Full %>% filter(Ds > 0 & Pn > 0) %>%
  arrange(desc(OR))
Dat

test <- ks.test(Full$kaks[Full$n.sig == 0],
                    Full$kaks[Full$n.sig > 0],
                    alternative = "greater")
p1 <- Dat %>%
  ggplot(aes(x = p.value)) +
  geom_histogram(aes(fill = n.sig > 0), bins = 20) +
  ggtitle(label = "", subtitle = test$p.value) +
  AMOR::theme_blackbox()
p1
filename <- paste0(c("overall","pvals_hist.png"), collapse = ".")
filename <- file.path(args$overall_out, filename)
ggsave(filename, p1, width = 8, height = 5, dpi = 200)

p1 <- Dat %>%
  ggplot(aes(x = mkratio, y = -log10(p.value))) +
  geom_point(aes(size = Dn)) +
  scale_x_log10() +
  AMOR::theme_blackbox()
p1
filename <- paste0(c("overall","volcano.png"), collapse = ".")
filename <- file.path(args$overall_out, filename)
ggsave(filename, p1, width = 6, height = 5, dpi = 200)