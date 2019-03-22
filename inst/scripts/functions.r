# (C) Copyright 2019 Sur Herrera Paredes
# 
# This file is part of HMVAR.
# 
# HMVAR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# HMVAR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with HMVAR.  If not, see <http://www.gnu.org/licenses/>.

##########################
myfun <- function(gene, contig, start, end, Dn, Ds, Pn, Ps, ratio, ratio.pval, ratio.perm, pos){
  test <- fisher.test(matrix(c(Dn,Pn,Ds,Ps), ncol = 2))
  pval <- test$p.value
  or <- test$estimate
  tibble(p.value = pval, OR = or)
}

process_mkres <- function(f, dir){
  mkres <- read_tsv(paste0(dir, '/', f), col_types = 'ccnnnnnnnnn', na = 'nan')
  mkres <- mkres %>%
    filter(Ds > 0 & Pn > 0) %>%
    mutate(pos = (start + end) / 2)
  fisher <- mkres %>% pmap_dfr(myfun)
  mkres <- bind_cols(mkres, fisher)
  
  # fisher exact can return numbers just above 1
  mkres$p.value[ mkres$p.value > 1 ] <- 1
  
  return(mkres)
}

plot_mkres <- function(mkres, name, dir, threshold = 0.01){
  # # Compare naive ratio with ML
  # p1 <- ggplot(mkres, aes(x = ratio, y = OR)) +
  #   geom_point()
  # p1
  
  ### P-value histogram
  p1 <- ggplot(mkres, aes(x = p.value)) +
    geom_histogram(bins = 20) +
    theme_blackbox()
  filename <- paste0(dir, "/", name, "_pval.hist.svg")
  ggsave(filename, p1, width = 6, height = 4)
  
  ### p-value qqplot
  svglite::svglite(filename <- paste0(dir, "/", name, "_pval.qqplot.svg"))
  pval_qqplot(mkres$p.value)
  dev.off()
  
  ### p-value as a function of number of diff subs
  p1 <- ggplot(mkres, aes(x = Dn, y = -log(p.value))) +
    geom_point() +
    geom_smooth()
  filename <- paste0(dir, "/", name, "_pval.Dn.svg")
  ggsave(filename, p1, width = 6, height = 4)
  
  p1 <- ggplot(mkres, aes(x = Ds, y = -log(p.value))) +
    geom_point() +
    geom_smooth()
  filename <- paste0(dir, "/", name, "_pval.Ds.svg")
  ggsave(filename, p1, width = 6, height = 4)
  
  p1 <- ggplot(mkres, aes(x = Pn, y = -log(p.value))) +
    geom_point() +
    geom_smooth()
  filename <- paste0(dir, "/", name, "_pval.Pn.svg")
  ggsave(filename, p1, width = 6, height = 4)
  
  p1 <- ggplot(mkres, aes(x = Ps, y = -log(p.value))) +
    geom_point() +
    geom_smooth()
  filename <- paste0(dir, "/", name, "_pval.Ps.svg")
  ggsave(filename, p1, width = 6, height = 4)
  
  ### p-value function of all subs
  p1 <- ggplot(mkres, aes(x = Dn + Ds + Pn + Ps, y = -log(p.value))) +
    geom_point() +
    geom_smooth()
  filename <- paste0(dir, "/", name, "_pval.subs.svg")
  ggsave(filename, p1, width = 6, height = 4)
  
  p1 <- ggplot(mkres, aes(x = Dn + Ds + Pn + Ps, fill = p.value == 1)) +
    geom_density(alpha = 0.3)
  filename <- paste0(dir, "/", name, "_pval.subs.dens.svg")
  ggsave(filename, p1, width = 6, height = 4)
  
  ### plot ratio and pval across genome
  # p1 <- ggplot(mkres, aes(x = pos)) +
  #   facet_grid(~contig, scales = "free_x") +
  #   geom_line(aes(y = ratio)) +
  #   geom_point(aes(y = -log10(p.value), col = p.value < threshold)) +
  #   theme_blackbox() 
  # p1
  
  # p1 <- ggplot(mkres, aes(x = pos)) +
  #   facet_grid(~contig, scales = "free_x") +
  #   geom_line(aes(y = ratio)) +
  #   theme_blackbox() 
  # p1
  
  p1 <- ggplot(mkres, aes(x = pos)) +
    facet_grid(~contig, scales = "free_x", space = "free_x") +
    geom_point(aes(y = -log10(p.value), col = p.value < threshold)) +
    theme_blackbox() +
    theme(strip.text = element_text(angle = 90),
          axis.text.x = element_blank(),
          panel.background = element_rect(color = NA))
  # filename <- paste0(dir, "/", name, "_pval.manhat.svg")
  filename <- paste0(dir, "/", name, "_pval.manhat.png")
  ggsave(filename, p1, width = 25, height = 4, dpi = 300)
  
  p1 <- ggplot(mkres, aes(x = log2(ratio), y = -log10(p.value))) +
    geom_point(aes(color = p.value < threshold)) +
    theme_blackbox()
  filename <- paste0(dir, "/", name, "_volcano.svg")
  ggsave(filename, p1, width = 6, height = 4)
}

#######################
