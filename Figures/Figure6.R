temp <- tric_peaks5[orf == "YGR145W", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 65), color = "gray70", size = 2) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 29), color = "#41B6C4", size = 2) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 29, yend = 56), color = "#CAB2D6", size = 2) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 56, yend = 64), color = "#FDBF6F", size = 2) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), plot.background = element_blank(), panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/A.DomainSchematic.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

temp <- tric_peaks5[orf == "YGR145W", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 51), color = "gray70", size = 2) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 18), color = "#CAB2D6", size = 2) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 18, yend = 26), color = "#FDBF6F", size = 2) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), plot.background = element_blank(), panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/A.DomainSchematic.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Stalling in translatome
plot <- ggplot(data = tric_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))]) +
  #geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(data = tric_translatome_dt[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = tric_translatome_dt[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, occupancy), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#1F78B4') +
  xlim(-25, 25) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of binding") +
  theme(legend.position = c(.75,1.1), plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(), axis.text.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/A.TRiCpausing_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = tric_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, binding = movingAverage(tric_odds1, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, binding), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#A6CEE3') +
  stat_summary(aes(adjusted_start, binding), fun.y = "mean", geom = "line", size = 1.25, color = '#A6CEE3') +
  xlim(-25, 25) + scale_y_continuous(position = "right", breaks = c(0,1,2,3,4), labels = c("0.0","1.0","2.0","3.0","4.0")) + coord_cartesian(ylim = c(0,4))
G <- plot + theme_classic(20) + labs(y = "TRiC enrichment (odds ratio)", x = "Codons from start of binding") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/A.TRiCbinding_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = ssb_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(WT_adjusted_norm, n=5, center=T))]) +
  #geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(data = ssb_translatome_dt[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = ssb_translatome_dt[WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, occupancy), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#33A02C') +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#33A02C') +
  xlim(-25, 25) + coord_cartesian(ylim = c(0.8,1.4))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codons from start of binding") +
  theme(legend.position = c(.75,1.1), plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(), axis.text.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/B.SSBpausing_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = ssb_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, binding = movingAverage(ssb_odds1, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_start, binding), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#B2DF8A') +
  stat_summary(aes(adjusted_start, binding), fun.y = "mean", geom = "line", size = 1.25, color = '#B2DF8A') +
  xlim(-25, 25) + scale_y_continuous(position = "right", breaks = c(0,1,2,3,4), labels = c("0.0","1.0","2.0","3.0","4.0")) + coord_cartesian(ylim = c(0,4))
G <- plot + theme_classic(20) + labs(y = "Ssb enrichment (odds ratio)", x = "Codons from start of binding") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/B.SSBbinding_WTtranslatome.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)


### Codon optimality
temp <- as.data.table(unique(tric_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75,name]))
colnames(temp) <- c("orf")
plot <- ggplot(data = tric_tAI[tric_tAI$name %in% temp$orf]) + xlim(-25, 25) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_random, tAI_ma), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(aes(adjusted_random, tAI_ma), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, tAI_ma), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(adjusted_start, tAI_ma), fun.y = "mean", geom = "line", size = 1.25, color = '#1F78B4') +
  coord_cartesian(ylim = c(0.35,0.45))
G <- plot + theme_classic(20) + labs(y = "Codon optimality (sTAI)", x = "Codons from start of TRiC binding") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/C.TRiCbinding_CodonOptimality.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

temp <- as.data.table(unique(ssb_translatome_dt[WT_adjusted_rpc >= 0.5 & peak > 75,name]))
colnames(temp) <- c("orf")
plot <- ggplot(data = ssb_tAI[ssb_tAI$name %in% temp$orf]) + xlim(-25, 25) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_random, tAI_ma), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(aes(adjusted_random, tAI_ma), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, tAI_ma), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#33A02C') +
  stat_summary(aes(adjusted_start, tAI_ma), fun.y = "mean", geom = "line", size = 1.25, color = '#33A02C') +
  coord_cartesian(ylim = c(0.35,0.45))
G <- plot + theme_classic(20) + labs(y = "Codon optimality (sTAI)", x = "Codons from start of Ssb binding") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/D.SSBbinding_CodonOptimality.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)


### Tubulins
# TUB1
meta_cols <- c("SSB" = "#4DAF4A", "TRiC" = "#377EB8")
temp <- tric_dt[orf == "YML085C"]
temp <- tric_fishers[orf == "YFL037W"]
temp[, tric_sum_ma := movingAverage(tric_sum, n=5, center=T)]
temp[, ribo_sum_ma := movingAverage(ribo_sum, n=5, center=T)]
temp[, ssb_Schx_sum_ma := movingAverage(ssb_Schx_sum, n=5, center=T)]
temp[, ssb_Rchx_sum_ma := movingAverage(ssb_Rchx_sum, n=5, center=T)]
for (i in 1:nrow(temp)) {
  tric <- matrix(c(temp[i]$tric_sum_ma, temp[i]$ribo_sum_ma, (temp[i]$tric_total - temp[i]$tric_sum_ma), (temp[i]$ribo_total - temp[i]$ribo_sum_ma)), nrow = 2)
  tric_test <- fisher.test(tric)
  temp[i, tric_odds := tric_test$estimate]
  temp[i, tric_pvalue := tric_test$p.value]
  #ssb <- matrix(c(temp[i]$ssb_Schx_sum_ma, temp[i]$ssb_Rchx_sum_ma, (temp[i]$ssb_Schx_total - temp[i]$ssb_Schx_sum_ma), (temp[i]$ssb_Rchx_total - temp[i]$ssb_Rchx_sum_ma)), nrow = 2)
  #ssb_test <- fisher.test(ssb)
  #temp[i, ssb_odds := ssb_test$estimate]
  #temp[i, ssb_pvalue := ssb_test$p.value]
}
temp[, tric_padj := p.adjust(tric_pvalue, method = "BH"), by = orf]
temp[, ssb_padj := p.adjust(ssb_pvalue, method = "BH"), by = orf]
temp[, tric_odds_ma := movingAverage(tric_odds, n=5, center=T), by = orf]
temp[, ssb_odds_ma := movingAverage(ssb_odds, n=5, center=T), by = orf]
plot <- ggplot(aes(position, occupancy), data=temp[orf == "YML085C" & tric_odds_ma < Inf, .(position, occupancy = movingAverage(tric_odds, n=5, center=T))]) + 
  geom_line(aes(color = "TRiC"), size = 1.25) +
  geom_line(data = temp[orf == "YML085C" & ssb_odds_ma < Inf, .(position, occupancy2 = movingAverage(ssb_odds, n=5, center=T))], aes(position, occupancy2, color = "SSB"), size = 1.25) + 
  scale_color_manual(labels = c("Ssb", "TRiC"), values = meta_cols, name = "") +
  scale_x_continuous(breaks = c(0,200,400))
G <- plot + theme_classic(20) + labs(y = "Enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = c(.01,.99), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("TUB1")
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/E.TUB1_odds.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

temp <- tric_peaks5[orf == "YML085C", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = tric_peaks5[orf == "YML085C"], aes(x = orf, xend = orf, y = 0, yend = length), color = "gray70", size = 2) +
  geom_segment(data = tric_peaks5[orf == "YML085C"], aes(x = orf, xend = orf, y = 0, yend = 335), color = "gray30", size = 2) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/E.TUB1_domains.pdf", plot, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=translatome_dt[orf == "YML085C", .(position, occupancy = movingAverage(WT_norm, n=5, center=T))]) + xlim(145,195) +
  geom_line(color = "gray30", size = 1.25) + scale_y_continuous(breaks = c(0,1,2), labels = c("0.0","1.0","2.0"))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codon position") +
  theme(legend.position = c(.75,1.1), plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(), axis.text.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/E.TUB1_elongationrate1.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=temp[orf == "YML085C" & tric_odds_ma < Inf, .(position, occupancy = movingAverage(tric_odds, n=5, center=T))]) + 
  geom_vline(xintercept = 170, color = "gray70", linetype = 'dashed', size = 1.25) +
  geom_line(aes(color = "TRiC"), size = 1.25) +
  scale_color_manual(labels = c("TRiC"), values = meta_cols, name = "") +
  xlim(145,195) + scale_y_continuous(position = "right")
G <- plot + theme_classic(20) + labs(y = "TRiC enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/E.TUB1_elongationrate1_binding.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=translatome_dt[orf == "YML085C", .(position, occupancy = movingAverage(WT_norm, n=5, center=T))]) + xlim(332,382) +
  geom_line(color = "#E31A1C", size = 1.25) + scale_y_continuous(breaks = c(0,1,2), labels = c("0.0","1.0","2.0"))
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codon position") +
  theme(legend.position = c(.75,1.1), plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(), axis.text.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/E.TUB1_elongationrate2.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=temp[orf == "YML085C" & tric_odds_ma < Inf, .(position, occupancy = movingAverage(tric_odds, n=5, center=T))]) + 
  geom_vline(xintercept = 357, color = "gray70", linetype = 'dashed', size = 1.25) +
  geom_line(aes(color = "TRiC"), size = 1.25) +
  scale_color_manual(labels = c("TRiC"), values = meta_cols, name = "") +
  xlim(332,382) + scale_y_continuous(position = "right")
G <- plot + theme_classic(20) + labs(y = "TRiC enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/E.TUB1_elongationrate2_binding.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)


# TUB2
plot <- ggplot(aes(position, occupancy), data=temp[orf == "YFL037W" & tric_odds_ma < Inf & R2X > 7, .(position, occupancy = movingAverage(tric_odds, n=5, center=T))]) + 
  geom_line(data = temp[orf == "YFL037W" & ssb_odds_ma < Inf, .(position, occupancy2 = movingAverage(ssb_odds, n=5, center=T))], aes(position, occupancy2, color = "SSB"), size = 1.25) + 
  geom_line(aes(color = "TRiC"), size = 1.25) + 
  scale_color_manual(labels = c("Ssb", "TRiC"), values = meta_cols, name = "") +
  scale_x_continuous(breaks = c(0,200,400))+coord_cartesian(ylim = c(0,3))
G <- plot + theme_classic(20) + labs(y = "Enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = c(.01,.99), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("TUB2")
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/F.TUB2_odds.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

temp <- tric_peaks5[orf == "YFL037W", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = tric_peaks5[orf == "YFL037W"], aes(x = orf, xend = orf, y = 0, yend = length), color = "gray70", size = 2) +
  geom_segment(data = tric_peaks5[orf == "YFL037W"], aes(x = orf, xend = orf, y = 0, yend = 146), color = "gray30", size = 2) +
  geom_segment(data = tric_peaks5[orf == "YFL037W"], aes(x = orf, xend = orf, y = 146, yend = 332), color = "gray50", size = 2) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/F.TUB2_domains.pdf", plot, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=translatome_dt[orf == "YFL037W", .(position, occupancy = movingAverage(WT_norm, n=5, center=T))]) + xlim(145,195) +
  geom_line(color = "#E31A1C", size = 1.25)
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codon position") +
  theme(legend.position = c(.75,1.1), plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(), axis.text.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/F.TUB2_elongationrate1.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=temp[orf == "YFL037W" & tric_odds_ma < Inf, .(position, occupancy = movingAverage(tric_odds, n=5, center=T))]) + 
  geom_vline(xintercept = 170, color = "gray70", linetype = 'dashed', size = 1.25) +
  geom_line(aes(color = "TRiC"), size = 1.25) +
  scale_color_manual(labels = c("TRiC"), values = meta_cols, name = "") +
  xlim(145,195) + scale_y_continuous(position = "right", breaks = c(0,1,2,3), labels = c("0.0","1.0","2.0","3.0"))
G <- plot + theme_classic(20) + labs(y = "TRiC enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/F.TUB2_elongationrate1_binding.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=translatome_dt[orf == "YFL037W", .(position, occupancy = movingAverage(WT_norm, n=5, center=T))]) + xlim(323,373) +
  geom_line(color = "#E31A1C", size = 1.25)
G <- plot + theme_classic(20) + labs(y = "Norm. translatome reads", x = "Codon position") +
  theme(legend.position = c(.75,1.1), plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"), axis.line.x = element_blank(), axis.text.x = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/F.TUB2_elongationrate2.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=temp[orf == "YFL037W" & tric_odds_ma < Inf, .(position, occupancy = movingAverage(tric_odds, n=5, center=T))]) + 
  geom_vline(xintercept = 348, color = "gray70", linetype = 'dashed', size = 1.25) +
  geom_line(aes(color = "TRiC"), size = 1.25) +
  scale_color_manual(labels = c("TRiC"), values = meta_cols, name = "") +
  xlim(323,373) + scale_y_continuous(position = "right", breaks = c(0,1,2,3), labels = c("0.0","1.0","2.0","3.0"))
G <- plot + theme_classic(20) + labs(y = "TRiC enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = "none", axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig6/F.TUB2_elongationrate2_binding.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

