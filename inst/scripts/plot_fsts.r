
p1 <- ggplot(fst_pool$fst, aes(x= Fst_pool, y = FSTpool)) +
  geom_point() +
  AMOR::theme_blackbox()
p1

p1 <- ggplot(fst_pool$fst, aes(x= Fst, y = Fst_pool)) +
  geom_point() +
  AMOR::theme_blackbox()
p1

pairs(fst_pool$fst[,9:11])

fst_pool$fst %>%
  gather(key = "Estimator", value = "Fst", Fst, Fst_pool, FSTpool) %>%
  ggplot(aes(x = ref_pos, y = Fst)) +
  facet_grid(Estimator ~ .) +
  geom_point()

p1 <- ggplot(fst_pool$fst, aes(x = ref_pos, y = Fst)) +
  geom_point() +
  AMOR::theme_blackbox()
p1

summary(fst_pool$fst$FSTpool)
hist(log10(fst_pool$fst$FSTpool))

# summary(fst$Fst)
# hist(fst$Fst)
# 
# snp_effects <- determine_snp_effect(info = Dat$info %>%
#                                       select(site_id, major_allele, minor_allele, amino_acids) %>%
#                                       filter(!is.na(amino_acids)))
w_fst <- window_fst(dat = fst$fst, w_size = 1000, s_size = 300, sorted = TRUE)
w_fst_pool <- window_fst(dat = fst_pool$fst, w_size = 1000, s_size = 300, sorted = TRUE)

w_fst
w_fst_pool 


w_fst_pool %>%
  arrange(desc(Fst_pool))



p1 <- ggplot(w_fst, aes(x = (start + end) / 2, y = Fst)) +
  facet_grid(~ ref_id, scales = "free_x", space = "free_x") +
  # geom_point(aes(size = n_sites)) +
  geom_point() +
  geom_hline(yintercept = 0.25, color = "red", size = 2) +
  ggplot2::theme(panel.background = ggplot2::element_blank(), 
                 panel.grid = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(color = "black"), 
                 axis.title = ggplot2::element_text(color = "black", face = "bold"), 
                 axis.line.x.bottom = ggplot2::element_line(),
                 axis.line.y.left = ggplot2::element_line())
p1

p1 <- ggplot(w_fst_pool, aes(x = (start + end) / 2, y = Fst_pool)) +
  facet_grid(~ ref_id, scales = "free_x", space = "free_x") +
  # geom_point(aes(size = n_sites)) +
  geom_point() +
  geom_hline(yintercept = 0.25, color = "red", size = 2) +
  ggplot2::theme(panel.background = ggplot2::element_blank(), 
                 panel.grid = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(color = "black"), 
                 axis.title = ggplot2::element_text(color = "black", face = "bold"), 
                 axis.line.x.bottom = ggplot2::element_line(),
                 axis.line.y.left = ggplot2::element_line())
p1

p1 <- ggplot(w_fst, aes(x = Fst, y = n_sites)) +
  geom_point() +
  geom_smooth(method = "loess") +
  AMOR::theme_blackbox() 
p1

