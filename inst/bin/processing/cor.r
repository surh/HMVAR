args <- list(vmwa.genes.dir = "hmp.vmwa.pvals.genes/",
             mkres.dir = "../2018-12-14.hmp_mktest/results/",
             q.value.thres = 0.01)

# List of files from mktest
mkres_files <- list.files(args$mkres.dir)







f <- mkres_files[12]
species <- str_replace(string = f, pattern = "_mktest.txt$", replacement = "")
vmwa_file <- paste0(species, "_vmwa.genes.txt")
cat(species, "\n", "\t", f, "\n\t", vmwa_file, "\n")

# Read data
vmwa <- read_tsv(paste0(args$vmwa.genes.dir, "/", vmwa_file))
vmwa <- vmwa %>% select(everything(),vmwa.OR = OR, vmwa.p.value = p.value, vmwa.q.value = q.value)
mkres <- read_tsv(paste0(args$mkres.dir, "/", f))
mkres <- mkres %>% select(gene_id = gene, everything()) %>% select(-ratio.pval)

# MK test
mkpval <- mkres %>%
  pmap_dbl(.f = function(Dn, Ds, Pn, Ps, ...){
    mat <- matrix(c(Dn, Ds, Pn, Ps), ncol = 2)
    res <- fisher.test(mat)
    return(res$p.value)})
mkres <- mkres %>% add_column(mk.p.value = mkpval)

# Join
Full <- vmwa %>% full_join(mkres, by = "gene_id") %>% select(-ref_id, -start, -end)
Full <- Full %>%
  add_column(significant = 1*(Full$vmwa.p.value < 0.01) + 2*(Full$mk.p.value < 0.01)) %>%
  mutate(significant = factor(significant)) %>%
  mutate(significant = recode_factor(significant,
                                     `0` = "none",
                                     `1` = "vmwa",
                                     `2` = "mktest",
                                     `3` = "both"))

p1 <- ggplot(Full, aes(x = vmwa.OR, y = ratio)) +
  geom_point(aes(col = significant)) +
  scale_y_log10() +
  scale_x_log10() +
  AMOR::theme_blackbox()
p1

p1 <- ggplot(Full,
             aes(x = -log10(vmwa.p.value),
                 y = -log10(mk.p.value))) +
  geom_point(aes(col = significant)) +
  # scale_y_log10() +
  # scale_x_log10() +
  AMOR::theme_blackbox()
p1
ggsave("wvma_mk_cor.png", p1, width = 5, height = 5, dpi = 150)
