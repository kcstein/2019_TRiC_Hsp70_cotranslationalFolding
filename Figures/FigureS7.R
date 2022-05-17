### Stalling at max positions
plot <- ggplot(data = tric_translatome_dt_max[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))]) +
  #geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(data = tric_translatome_dt_max[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = tric_translatome_dt_max[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, occupancy), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#1F78B4') +
  xlim(-25, 25) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of binding") +
  theme(legend.position = c(.75,1.1), plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(), axis.text.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/A.TRiCpausing_max_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = tric_translatome_dt_max[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, binding = movingAverage(tric_odds1, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, binding), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#A6CEE3') +
  stat_summary(aes(adjusted_start, binding), fun.y = "mean", geom = "line", size = 1.25, color = '#A6CEE3') +
  xlim(-25, 25) + scale_y_continuous(position = "right", breaks = c(0,1,2,3,4,5), labels = c("0.0","1.0","2.0","3.0","4.0","5.0")) + coord_cartesian(ylim = c(0,5))
G <- plot + theme_classic(20) + labs(y = "TRiC enrichment (odds ratio)", x = "Codons from start of binding") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/A.TRiCbinding_max_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(data = ssb_translatome_dt_max[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))]) +
  #geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(data = ssb_translatome_dt_max[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = ssb_translatome_dt_max[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, occupancy), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#33A02C') +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#33A02C') +
  xlim(-25, 25) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of binding") +
  theme(legend.position = c(.75,1.1), plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(), axis.text.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/A.SSBpausing_max_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = ssb_translatome_dt_max[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, binding = movingAverage(ssb_odds1, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, binding), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#B2DF8A') +
  stat_summary(aes(adjusted_start, binding), fun.y = "mean", geom = "line", size = 1.25, color = '#B2DF8A') +
  xlim(-25, 25) + scale_y_continuous(position = "right", breaks = c(0,1,2,3,4,5), labels = c("0.0","1.0","2.0","3.0","4.0","5.0")) + coord_cartesian(ylim = c(0,5))
G <- plot + theme_classic(20) + labs(y = "Ssb enrichment (odds ratio)", x = "Codons from start of binding") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/A.SSBbinding_max_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)


plot <- ggplot(data = bukau_translatome_dt_max[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))]) +
  #geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(data = bukau_translatome_dt_max[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = bukau_translatome_dt_max[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, occupancy), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#6A3D9A') +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#6A3D9A') +
  xlim(-25, 25) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of binding") +
  theme(legend.position = c(.75,1.1), plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(), axis.text.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/A.Bukaupausing_max_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = bukau_translatome_dt_max[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, binding = movingAverage(ssb_odds1, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, binding), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#CAB2D6') +
  stat_summary(aes(adjusted_start, binding), fun.y = "mean", geom = "line", size = 1.25, color = '#CAB2D6') +
  xlim(-25, 25) + scale_y_continuous(position = "right", breaks = c(0,1,2,3,4,5), labels = c("0.0","1.0","2.0","3.0","4.0","5.0")) + coord_cartesian(ylim = c(0,5))
G <- plot + theme_classic(20) + labs(y = "Ssb enrichment (odds ratio)", x = "Codons from start of binding") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/A.Bukaubinding_max_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)


### Bukau stalling in translatome
plot <- ggplot(data = bukau_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))]) +
  #geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(data = bukau_translatome_dt[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = bukau_translatome_dt[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, occupancy), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#6A3D9A') +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#6A3D9A') +
  xlim(-25, 25) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of binding") +
  theme(legend.position = c(.75,1.1), plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(), axis.text.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/B.Bukaupausing_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = bukau_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, binding = movingAverage(ssb_odds1, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, binding), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#CAB2D6') +
  stat_summary(aes(adjusted_start, binding), fun.y = "mean", geom = "line", size = 1.25, color = '#CAB2D6') +
  xlim(-25, 25) + scale_y_continuous(position = "right", breaks = c(0,1,2,3,4), labels = c("0.0","1.0","2.0","3.0","4.0")) + coord_cartesian(ylim = c(0,4))
G <- plot + theme_classic(20) + labs(y = "Ssb enrichment (odds ratio)", x = "Codons from start of binding") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/B.Bukaubinding_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)


### All translatomes, all positions
plot <- ggplot(data = bukau_translatome_dt[C_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(C_WT_adjusted_norm, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#BCBDDC') + 
  stat_summary(data = bukau_translatome_dt[L_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(L_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#9E9AC8') + 
  stat_summary(data = bukau_translatome_dt[G_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(G_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#807DBA') + 
  stat_summary(data = bukau_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#6A51A3') + 
  #stat_summary(data = bukau_translatome_dt[B_WT_rpc >= 20, .(adjusted_start, occupancy = movingAverage(B_WT_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#54278F') + 
  scale_x_continuous(limits = c(-35, 25), breaks = c(-30,-20,-10,0,10,20)) + coord_cartesian(ylim = c(0.8,1.3))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of Ssb binding") +
  theme(axis.text = element_text(size = 20, color = "black"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/C.Bukaubinding_AllTranslatomes.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = ssb_translatome_dt[C_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(C_WT_adjusted_norm, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#A1D99B') + 
  stat_summary(data = ssb_translatome_dt[L_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(L_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#74C476') + 
  stat_summary(data = ssb_translatome_dt[G_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(G_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#41AB5D') + 
  stat_summary(data = ssb_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#238B45') + 
  #stat_summary(data = ssb_translatome_dt[B_WT_rpc >= 20, .(adjusted_start, occupancy = movingAverage(B_WT_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#238B45') + 
  scale_x_continuous(limits = c(-35, 25), breaks = c(-30,-20,-10,0,10,20)) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of Ssb binding") +
  theme(axis.text = element_text(size = 20, color = "black"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/C.SSBbinding_AllTranslatomes.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = tric_translatome_dt[C_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(C_WT_adjusted_norm, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#9ECAE1') + 
  stat_summary(data = tric_translatome_dt[L_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(L_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#6BAED6') + 
  stat_summary(data = tric_translatome_dt[G_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(G_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#4292C6') + 
  stat_summary(data = tric_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#2171B5') + 
  #stat_summary(data = tric_translatome_dt[B_WT_rpc >= 20, .(adjusted_start, occupancy = movingAverage(B_WT_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#238B45') + 
  scale_x_continuous(limits = c(-35, 25), breaks = c(-30,-20,-10,0,10,20)) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of TRiC binding") +
  theme(axis.text = element_text(size = 20, color = "black"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/C.TRiCbinding_AllTranslatomes.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)


### All translatomes, max positions
plot <- ggplot(data = bukau_translatome_dt_max[C_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(C_WT_adjusted_norm, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#BCBDDC') + 
  stat_summary(data = bukau_translatome_dt_max[L_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(L_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#9E9AC8') + 
  stat_summary(data = bukau_translatome_dt_max[G_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(G_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#807DBA') + 
  stat_summary(data = bukau_translatome_dt_max[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#6A51A3') + 
  #stat_summary(data = bukau_translatome_dt_max[B_WT_rpc >= 20, .(adjusted_start, occupancy = movingAverage(B_WT_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#54278F') + 
  scale_x_continuous(limits = c(-35, 25), breaks = c(-30,-20,-10,0,10,20)) + coord_cartesian(ylim = c(0.8,1.3))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of Ssb binding") +
  theme(axis.text = element_text(size = 20, color = "black"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/D.BukaubindingMax_AllTranslatomes.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = ssb_translatome_dt_max[C_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(C_WT_adjusted_norm, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#A1D99B') + 
  stat_summary(data = ssb_translatome_dt_max[L_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(L_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#74C476') + 
  stat_summary(data = ssb_translatome_dt_max[G_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(G_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#41AB5D') + 
  stat_summary(data = ssb_translatome_dt_max[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#238B45') + 
  #stat_summary(data = ssb_translatome_dt_max[B_WT_rpc >= 20, .(adjusted_start, occupancy = movingAverage(B_WT_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#238B45') + 
  scale_x_continuous(limits = c(-35, 25), breaks = c(-30,-20,-10,0,10,20)) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of Ssb binding") +
  theme(axis.text = element_text(size = 20, color = "black"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/D.SSBbindingMax_AllTranslatomes.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = tric_translatome_dt_max[C_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(C_WT_adjusted_norm, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#9ECAE1') + 
  stat_summary(data = tric_translatome_dt_max[L_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(L_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#6BAED6') + 
  stat_summary(data = tric_translatome_dt_max[G_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(G_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#4292C6') + 
  stat_summary(data = tric_translatome_dt_max[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#2171B5') + 
  #stat_summary(data = tric_translatome_dt_max[B_WT_rpc >= 20, .(adjusted_start, occupancy = movingAverage(B_WT_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#238B45') + 
  scale_x_continuous(limits = c(-35, 25), breaks = c(-30,-20,-10,0,10,20)) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of TRiC binding") +
  theme(axis.text = element_text(size = 20, color = "black"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/D.TRiCbindingMax_AllTranslatomes.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

### Response letter
plot <- ggplot(data = doring_translatome_dt_B1[C_WT_adjusted_rpc >= 0.5 & peak_startaa > 75, .(adjusted_start, occupancy = movingAverage(C_WT_adjusted_norm, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = 'gray10') + 
  stat_summary(data = doring_translatome_dt_B1[L_WT_adjusted_rpc >= 0.5 & peak_startaa > 75, .(adjusted_start, occupancy = movingAverage(L_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = 'gray30') + 
  stat_summary(data = doring_translatome_dt_B1[G_WT_adjusted_rpc >= 0.5 & peak_startaa > 75, .(adjusted_start, occupancy = movingAverage(G_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(data = doring_translatome_dt_B1[WT_adjusted_rpc >= 0.5 & peak_startaa > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = 'gray70') + 
  #stat_summary(data = doring_translatome_dt_B1[B_WT_adjusted_rpc >= 0.5 & peak_startaa > 75, .(adjusted_start, occupancy = movingAverage(B_WT_adjusted_norm, n=5, center=T))], aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#238B45') + 
  scale_x_continuous(limits = c(-35, 25), breaks = c(-30,-20,-10,0,10,20)) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of Ssb binding") +
  theme(axis.text = element_text(size = 20, color = "black"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank())
ggsave("/Users/KevinStein/Desktop/DoringSitesAndORFsubset_AllTranslatomes.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)


### SRP
plot <- ggplot(data = srp_translatome_dt[WT_adjusted_rpc >= 0.5, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))]) +
  #geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(data = srp_translatome_dt[WT_adjusted_rpc_random >= 0.5, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = srp_translatome_dt[WT_adjusted_rpc_random >= 0.5, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, occupancy), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#FF7F00') +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#FF7F00') +
  xlim(-25, 25) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of SRP binding") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/E.SRPpausing_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)


### Codon optimality
temp <- as.data.table(unique(bukau_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75,name]))
colnames(temp) <- c("orf")
plot <- ggplot(data = bukau_tAI[bukau_tAI$name %in% temp$orf]) + xlim(-25, 25) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_random, tAI_ma), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(aes(adjusted_random, tAI_ma), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, tAI_ma), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#6A3D9A') +
  stat_summary(aes(adjusted_start, tAI_ma), fun.y = "mean", geom = "line", size = 1.25, color = '#6A3D9A') +
  stat_summary(data = bukau_tAI[bukau_tAI$name %in% temp$orf & adjusted_start <= -1 & adjusted_start >= -3], 
               aes(adjusted_start, tAI_ma), fun.y = "mean", geom = "line", size = 1.25, color = '#FF7F00') +
  stat_summary(data = bukau_tAI[bukau_tAI$name %in% temp$orf & adjusted_start >= 1 & adjusted_start <= 3], 
               aes(adjusted_start, tAI_ma), fun.y = "mean", geom = "line", size = 1.25, color = '#FF7F00') + 
  coord_cartesian(ylim = c(0.4,0.5))
G <- plot + theme_classic(20) + labs(y = "Codon optimality (sTAI)", x = "Codons from start of Ssb binding") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS7/F.Bukaubinding_CodonOptimality.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)
wilcox.test(bukau_tAI[bukau_tAI$name %in% temp$orf & adjusted_start <= -1 & adjusted_start >= -3]$tAI, bukau_tAI[bukau_tAI$name %in% temp$orf & adjusted_random <= -1 & adjusted_random >= -3]$tAI, alternative = 't')
wilcox.test(bukau_tAI[bukau_tAI$name %in% temp$orf & adjusted_start >= 1 & adjusted_start <= 3]$tAI, bukau_tAI[bukau_tAI$name %in% temp$orf & adjusted_random >= 1 & adjusted_random <= 3]$tAI, alternative = 't')

