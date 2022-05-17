### ENP2 with adding 1 to each position and using moving average before calculating Fishers
temp <- tric_peaks5[orf == "YGR145W", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  #geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = position, yend = (position+5)), color = "#1F78B4", size = 15) +
  #geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 214, yend = 219), color = "#33A02C", size = 15) +
  #geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 260, yend = 265), color = "#33A02C", size = 15) +
  #geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 301, yend = 306), color = "#33A02C", size = 15) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = length), color = "gray70", size = 2) +
  geom_segment(data = wd_repeats[orf == "YGR145W"], aes(x = orf, xend = orf, y = start, yend = end), color = "#CAB2D6", size = 8) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 98, yend = 136), color = "#CAB2D6", size = 8) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 145, yend = 170), color = "#CAB2D6", size = 8) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig4/A.ENP2_domains.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=temp[orf == "YGR145W" & ssb_odds_ma < Inf, .(position, occupancy = movingAverage(ssb_odds, n=5, center=T))]) + 
  geom_vline(xintercept = 84, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = 126, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 175, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 208, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 256, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 299, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = 342, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_line(aes(color = "1"),size = 1.25) +
  scale_color_manual(values = c("1" = "#33A02C"), labels = c("Ssb"), name = "")
G <- plot + theme_classic(20) + labs(y = "Ssb enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = c(.55,.9), legend.justification = c("left", "top"), axis.text = element_text(size = 16, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("ENP2")
ggsave("/Users/KevinStein/Desktop/Figures/Fig4/A.ENP2_comparisonSsb.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=temp[orf == "YGR145W" & tric_odds_ma < Inf, .(position, occupancy = movingAverage(tric_odds, n=5, center=T))]) + 
  geom_vline(xintercept = 84, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = 126, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 175, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 208, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 256, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 299, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = 342, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_line(aes(color = "1"),size = 1.25) +
  scale_color_manual(values = c("1" = "#1F78B4"), labels = c("TRiC"), name = "")
G <- plot + theme_classic(20) + labs(y = "TRiC enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = c(.55,.9), legend.justification = c("left", "top"), axis.text = element_text(size = 16, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("ENP2")
ggsave("/Users/KevinStein/Desktop/Figures/Fig4/A.ENP2_comparisonTRiC.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


### Metagene of WD
temp1 <- wd_repeats[wd_repeats$orf %in% tric_substrates5$orf]
setkeyv(temp1, c("orf"))
#setkeyv(substrates_fishers, c("orf"))
#setkeyv(temp, c("orf"))
#wd_repeats_dt <- substrates_fishers[temp1]
wd_repeats_dt1 <- substrates_fishers[temp1, allow.cartesian = TRUE]
wd_repeats_dt1[, wdposition_start30 := position - (start+30)]
plot <- ggplot(data = wd_repeats_dt1[(wd_repeats_dt1$orf %in% tric_substrates5$orf) & repeat. == "WD1", .(wdposition_start30, tric = movingAverage(tric_odds, n=15, center=T),
                                                                                                          ssb = movingAverage(ssb_odds, n=15, center=T))]) +
  geom_hline(yintercept = 1, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = temp1[repeat. == "WD1"]$start_avg30, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = temp1[repeat. == "WD2"]$start_avg30, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = temp1[repeat. == "WD3"]$start_avg30, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = temp1[repeat. == "WD4"]$start_avg30, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = temp1[repeat. == "WD5"]$start_avg30, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = temp1[repeat. == "WD6"]$start_avg30, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = temp1[repeat. == "WD7"]$start_avg30, color = 'gray50', linetype = 'dashed', size = 1) +
  stat_summary(aes(wdposition_start30, ssb), fun.data = "median_hilow", geom = "ribbon", alpha = 0.5,
               fun.args=list(conf.int=0.5), fill = '#33A02C') +
  stat_summary(aes(wdposition_start30, tric), fun.data = "median_hilow", geom = "ribbon", alpha = 0.5,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(wdposition_start30, ssb, color = '#33A02C'), fun.y = "median", geom = "line", size = 1.25) +
  stat_summary(aes(wdposition_start30, tric, color = '#1F78B4'), fun.y = "median", geom = "line", size = 1.25) + 
  scale_color_manual(labels = c("TRiC", "Ssb"), values = c("#1F78B4", "#33A02C"), name = "") +
  coord_cartesian(ylim = c(0,2.4), xlim = c(0, 375))
G <- plot + theme_classic(20) + labs(y = "Enrichment (odds ratio)", x = "Codons after 1st WD repeat emerges") +
  theme(legend.position = c(.85,1.1), legend.justification = c("left", "top"), legend.text = element_text(size=20), legend.background = element_blank(),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig4/C.WD_firstrepeat_median3.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

temp <- tric_peaks5[orf == "YGR145W", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 375), color = "gray70", size = 2) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 37), color = "#CAB2D6", size = 8) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 61, yend = 100), color = "#CAB2D6", size = 8) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 123, yend = 161), color = "#CAB2D6", size = 8) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 180, yend = 218), color = "#CAB2D6", size = 8) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 234, yend = 272), color = "#CAB2D6", size = 8) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 281, yend = 319), color = "#CAB2D6", size = 8) +
  geom_segment(data = tric_peaks5[orf == "YGR145W"], aes(x = orf, xend = orf, y = 335, yend = 371), color = "#CAB2D6", size = 8) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
        axis.text = element_blank(), plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig4/C.WD_domains.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
