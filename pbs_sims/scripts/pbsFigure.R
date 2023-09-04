# Valeria Anorve-Garibay
# MUC19 #
# August 2023

# Figure S30 and S31

library(ggplot2)
library(ggExtra)
library(ggpubr)

source("sim-pbs-functions.R")

all_snps_st = 0.1408079128199867
den_snps_st = 0.27948744669683

# S30
# vmodel
vneutral = select_sims("pbs/pbs_vneutral_uniform.txt")

# no introgression
vneutral_noIntro = ggplot(vneutral, aes(x = position, y = pbs, color = factor(mt))) +
  geom_point(size = 2, shape = 20) + labs(y = "", x = "", title = "") + theme_classic() +
  scale_color_manual(values = c("gray70", "gray70")) +
  scale_x_continuous(breaks = c(0, 218562, 437124), labels = c("0KB", "218KB", "437KB")) +
  ylim(-1.1,1) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16)) +
  geom_hline(yintercept = all_snps_st, color = "red",  linetype = "dashed")

vneutral_noIntro = ggMarginal(vneutral_noIntro, type="density", margins = "y",
           colour = "black", fill = "black", size = 10)

vneutral_noIntro_hist = ggplot(vneutral, aes(x = pbs)) +
  geom_histogram() + labs(y = "", x = "", title = "") + theme_classic() +
  geom_vline(xintercept = all_snps_st, color = "red",  linetype = "dashed") +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16))

length(vneutral$pbs[which(vneutral$pbs >= all_snps_st)])/length(vneutral$pbs)
length(vneutral$pbs[which(vneutral$pbs <= all_snps_st)])/length(vneutral$pbs)

# introgression
vneutral_DEN = vneutral[which(vneutral$pO == 1),]
vneutral_intro = ggplot(vneutral_DEN, aes(x = position, y = pbs, color = factor(mt))) +
  geom_point(size = 2, shape = 20) + labs(y = "", x = "", title = "") + theme_classic() +
  scale_color_manual(values = c("gray70", "gray70")) +
  scale_x_continuous(breaks = c(0, 218562, 437124), labels = c("0KB", "218KB", "437KB")) +
  ylim(-1.1,1) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16)) +
  geom_hline(yintercept = den_snps_st, color = "red",  linetype = "dashed")
vneutral_intro = ggMarginal(vneutral_intro, type="density", margins = "y",
                              colour = "black", fill = "black", size = 10)

vneutral_intro_hist = ggplot(vneutral_DEN, aes(x = pbs)) +
  geom_histogram() + labs(y = "", x = "", title = "") + theme_classic() +
  geom_vline(xintercept = den_snps_st, color = "red",  linetype = "dashed") +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16))

length(vneutral_DEN$pbs[which(vneutral_DEN$pbs >= den_snps_st)])/length(vneutral_DEN$pbs)
length(vneutral_DEN$pbs[which(vneutral_DEN$pbs <= den_snps_st)])/length(vneutral_DEN$pbs)

vneutral_gg = ggarrange(vneutral_noIntro, vneutral_intro)
vneutral_gg_hist = ggarrange(vneutral_noIntro_hist, vneutral_intro_hist)

vnegative = select_sims("pbs/pbs_vnegative_uniform.txt")

# no introgression
vnegative_noIntro = ggplot(vnegative, aes(x = position, y = pbs, color = factor(mt))) +
  geom_point(size = 2, shape = 20) + labs(y = "", x = "", title = "") + theme_classic() +
  scale_color_manual(values = c("gray70", "firebrick3","gray70")) +
  scale_x_continuous(breaks = c(0, 218562, 437124), labels = c("0KB", "218KB", "437KB")) +
  ylim(-1.1,1) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16)) +
  geom_hline(yintercept = all_snps_st, color = "red",  linetype = "dashed")

vnegative_noIntro = ggMarginal(vnegative_noIntro, type="density", margins = "y",
           colour = "black", fill = "black", size = 10)

vnegative_noIntro_hist = ggplot(vnegative, aes(x = pbs)) +
  geom_histogram() + labs(y = "", x = "", title = "") + theme_classic() +
  geom_vline(xintercept = all_snps_st, color = "red",  linetype = "dashed") +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16))

length(vnegative$pbs[which(vnegative$pbs >= all_snps_st)])/length(vnegative$pbs)
length(vnegative$pbs[which(vnegative$pbs <= all_snps_st)])/length(vnegative$pbs)

# introgression
vnegative_DEN = vnegative[which(vnegative$pO == 1),]
vnegative_intro = ggplot(vnegative_DEN, aes(x = position, y = pbs, color = factor(mt))) +
  geom_hline(yintercept = den_snps_st, color = "red",  linetype = "dashed") +
  geom_point(size = 2, shape = 20) + labs(y = "", x = "", title = "") + theme_classic() +
  scale_color_manual(values = c("gray70", "firebrick3","gray70")) +
  scale_x_continuous(breaks = c(0, 218562, 437124), labels = c("0KB", "218KB", "437KB")) +
  ylim(-1.1,1) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16))
vnegative_intro = ggMarginal(vnegative_intro, type="density", margins = "y",
           colour = "black", fill = "black", size = 10)

vnegative_intro_hist = ggplot(vnegative_DEN, aes(x = pbs)) +
  geom_histogram() + labs(y = "", x = "", title = "") + theme_classic() +
  geom_vline(xintercept = den_snps_st, color = "red",  linetype = "dashed") +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16))

length(vnegative_DEN$pbs[which(vnegative_DEN$pbs >= den_snps_st)])/length(vnegative_DEN$pbs)
length(vnegative_DEN$pbs[which(vnegative_DEN$pbs <= den_snps_st)])/length(vnegative_DEN$pbs)

vnegative_gg = ggarrange(vnegative_noIntro, vnegative_intro)
vnegative_gg_hist = ggarrange(vnegative_noIntro_hist, vnegative_intro_hist)

pdf("figures/vmodel_map.pdf", paper = "USr", width = 14, height = 8)
ggarrange(vneutral_gg, vnegative_gg, ncol = 1, nrow = 2)
dev.off()

pdf("figures/vmodel_uniform_hist.pdf", paper = "USr", width = 14, height = 8)
ggarrange(vneutral_gg_hist, vnegative_gg_hist, ncol = 1, nrow = 2)
dev.off()

# S31
# smodel
sneutral = select_sims("pbs/pbs_sneutral_uniform.txt")

# no introgression
sneutral_noIntro = ggplot(sneutral, aes(x = position, y = pbs, color = factor(mt))) +
  geom_point(size = 2, shape = 20) + labs(y = "", x = "", title = "") + theme_classic() +
  scale_color_manual(values = c("gray70", "gray70")) +
  scale_x_continuous(breaks = c(0, 218562, 437124), labels = c("0KB", "218KB", "437KB")) +
  ylim(-1.1,1) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16)) +
  geom_hline(yintercept = all_snps_st, color = "red",  linetype = "dashed")

sneutral_noIntro = ggMarginal(sneutral_noIntro, type="density", margins = "y",
                              colour = "black", fill = "black", size = 10)

sneutral_noIntro_hist = ggplot(sneutral, aes(x = pbs)) +
  geom_histogram() + labs(y = "", x = "", title = "") + theme_classic() +
  geom_vline(xintercept = all_snps_st, color = "red",  linetype = "dashed") +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16))

length(sneutral$pbs[which(sneutral$pbs >= all_snps_st)])/length(sneutral$pbs)
length(sneutral$pbs[which(sneutral$pbs <= all_snps_st)])/length(sneutral$pbs)

# introgression
sneutral_DEN = sneutral[which(sneutral$pO == 3),]
sneutral_intro = ggplot(sneutral_DEN, aes(x = position, y = pbs, color = factor(mt))) +
  geom_point(size = 2, shape = 20) + labs(y = "", x = "", title = "") + theme_classic() +
  scale_color_manual(values = c("gray70", "gray70")) +
  scale_x_continuous(breaks = c(0, 218562, 437124), labels = c("0KB", "218KB", "437KB")) +
  ylim(-1.1,1) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16)) +
  geom_hline(yintercept = den_snps_st, color = "red",  linetype = "dashed")

sneutral_intro = ggMarginal(sneutral_intro, type="density", margins = "y",
                            colour = "black", fill = "black", size = 10)

sneutral_intro_hist = ggplot(sneutral_DEN, aes(x = pbs)) +
  geom_histogram() + labs(y = "", x = "", title = "") + theme_classic() +
  geom_vline(xintercept = den_snps_st, color = "red",  linetype = "dashed") +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16)) + xlim(min(sneutral_DEN$pbs), 0.31)

length(sneutral_DEN$pbs[which(sneutral_DEN$pbs >= den_snps_st)])/length(sneutral_DEN$pbs)
length(sneutral_DEN$pbs[which(sneutral_DEN$pbs <= den_snps_st)])/length(sneutral_DEN$pbs)

sneutral_gg = ggarrange(sneutral_noIntro, sneutral_intro)
sneutral_gg_hist = ggarrange(sneutral_noIntro_hist, sneutral_intro_hist)

snegative = select_sims("pbs/pbs_snegative_uniform.txt")

# no introgression
snegative_noIntro = ggplot(snegative, aes(x = position, y = pbs, color = factor(mt))) +
  geom_point(size = 2, shape = 20) + labs(y = "", x = "", title = "") + theme_classic() +
  scale_color_manual(values = c("gray70", "firebrick3","gray70")) +
  scale_x_continuous(breaks = c(0, 218562, 437124), labels = c("0KB", "218KB", "437KB")) +
  ylim(-1.1,1) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16)) +
  geom_hline(yintercept = all_snps_st, color = "red",  linetype = "dashed")

snegative_noIntro = ggMarginal(snegative_noIntro, type="density", margins = "y",
                               colour = "black", fill = "black", size = 10)

snegative_noIntro_hist = ggplot(snegative, aes(x = pbs)) +
  geom_histogram() + labs(y = "", x = "", title = "") + theme_classic() +
  geom_vline(xintercept = all_snps_st, color = "red",  linetype = "dashed") +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16))

length(snegative$pbs[which(snegative$pbs >= all_snps_st)])/length(snegative$pbs)
length(snegative$pbs[which(snegative$pbs <= all_snps_st)])/length(snegative$pbs)

# introgression
snegative_DEN = snegative[which(snegative$pO == 3),]
snegative_intro = ggplot(snegative_DEN, aes(x = position, y = pbs, color = factor(mt))) +
  geom_hline(yintercept = den_snps_st, color = "red",  linetype = "dashed") +
  geom_point(size = 2, shape = 20) + labs(y = "", x = "", title = "") + theme_classic() +
  scale_color_manual(values = c("gray70", "firebrick3","gray70")) +
  scale_x_continuous(breaks = c(0, 218562, 437124), labels = c("0KB", "218KB", "437KB")) +
  ylim(-1.1,1) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16))

snegative_intro = ggMarginal(snegative_intro, type="density", margins = "y",
                             colour = "black", fill = "black", size = 10)

snegative_intro_hist = ggplot(snegative_DEN, aes(x = pbs)) +
  geom_histogram() + labs(y = "", x = "", title = "") + theme_classic() +
  geom_vline(xintercept = den_snps_st, color = "red",  linetype = "dashed") +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 16),
        axis.title.x =element_text(hjust = 0.5, size = 16),
        axis.title.y =element_text(hjust = 0.5, size = 16)) + xlim(min(snegative_DEN$pbs),0.31)

length(snegative_DEN$pbs[which(snegative_DEN$pbs >= den_snps_st)])/length(snegative_DEN$pbs)
length(snegative_DEN$pbs[which(snegative_DEN$pbs <= den_snps_st)])/length(snegative_DEN$pbs)

snegative_gg = ggarrange(snegative_noIntro, snegative_intro)
snegative_gg_hist = ggarrange(snegative_noIntro_hist, snegative_intro_hist)

pdf("figures/smodel_uniform.pdf", paper = "USr", width = 14, height = 8)
ggarrange(sneutral_gg, snegative_gg, ncol = 1, nrow = 2)
dev.off()

pdf("figures/smodel_uniform_hist.pdf", paper = "USr", width = 14, height = 8)
ggarrange(sneutral_gg_hist, snegative_gg_hist, ncol = 1, nrow = 2)
dev.off()

