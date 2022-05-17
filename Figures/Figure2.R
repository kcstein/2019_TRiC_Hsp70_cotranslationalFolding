### Example substrates
plot <- ggplot(aes(position, occupancy), data=temp[orf == "YLL018C" & tric_odds < Inf, .(position, occupancy = movingAverage(tric_odds, n=5, center=T))]) + 
  geom_line(data = temp[orf == "YLL018C" & ssb_odds < Inf, .(position, occupancy2 = movingAverage(ssb_odds, n=5, center=T))], aes(position, occupancy2, color = "SSB"), size = 1.25) + 
  geom_line(aes(color = "TRiC"), size = 1.25) + 
  scale_color_manual(labels = c("SSB", "TRiC"), values = meta_cols, name = "") +
  scale_x_continuous(breaks = c(0,200,400))
G <- plot + theme_classic(20) + labs(y = "Enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = c(.01,.99), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=16),
        axis.text = element_text(size = 16, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("DPS1")
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/A.DPS1_YLL018C_odds.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

temp <- tric_peaks_all[orf == "YLL018C", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = tric_peaks_all[orf == "YLL018C"], aes(x = orf, xend = orf, y = 0, yend = length), color = "gray70", size = 2) +
  geom_segment(data = domains_cath[orf == "YLL018C" & Domain == "G3DSA:2.40.50.140"], aes(x = orf, xend = orf, y = Start, yend = End), color = "#CAB2D6", size = 8) +
  geom_segment(data = domains_cath[orf == "YLL018C" & Domain == "G3DSA:3.30.930.10"], aes(x = orf, xend = orf, y = Start, yend = End), color = "#CAB2D6", size = 8) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/A.DPS1_domains.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(aes(position, occupancy), data=temp[orf == "YHL015W" & tric_odds < Inf, .(position, occupancy = movingAverage(tric_odds, n=3, center=T))]) + 
  geom_line(data = temp[orf == "YHL015W" & ssb_odds < Inf, .(position, occupancy2 = movingAverage(ssb_odds, n=3, center=T))], aes(position, occupancy2, color = "SSB"), size = 1.25) + 
  geom_line(aes(color = "TRiC"), size = 1.25) + 
  scale_color_manual(labels = c("SSB", "TRiC"), values = meta_cols, name = "") + xlim(0,136)
G <- plot + theme_classic(20) + labs(y = "Enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = c(.01,.99), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=16),
        axis.text = element_text(size = 16, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("RPS20")
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/A.RPS20_YHL015W_odds.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

temp <- ssb_peaks5[orf == "YHL015W", .SD[which.min(peak)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = ssb_peaks5[orf == "YHL015W"], aes(x = orf, xend = orf, y = 0, yend = length), color = "gray70", size = 2) +
  geom_segment(data = domains_cath[orf == "YHL015W" & Domain == "G3DSA:3.30.70.600"], aes(x = orf, xend = orf, y = Start, yend = End), color = "#CAB2D6", size = 8) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/A.RPS20_domains.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Number of recruitment sites per orf
localization_tric[, peaks_group := ifelse(PeaksPerOrf == 1, 1, 0)]
localization_tric[, peaks_group := ifelse(PeaksPerOrf == 2, 2, peaks_group)]
localization_tric[, peaks_group := ifelse(PeaksPerOrf == 3, 3, peaks_group)]
localization_tric[, peaks_group := ifelse(PeaksPerOrf == 4, 4, peaks_group)]
localization_tric[, peaks_group := ifelse(PeaksPerOrf >= 5, 5, peaks_group)]
localization_ssb[, peaks_group := ifelse(PeaksPerOrf == 1, 1, 0)]
localization_ssb[, peaks_group := ifelse(PeaksPerOrf == 2, 2, peaks_group)]
localization_ssb[, peaks_group := ifelse(PeaksPerOrf == 3, 3, peaks_group)]
localization_ssb[, peaks_group := ifelse(PeaksPerOrf == 4, 4, peaks_group)]
localization_ssb[, peaks_group := ifelse(PeaksPerOrf >= 5, 5, peaks_group)]

plot <- ggplot(localization_tric[localization == "Cytonuclear"], aes(x="1.TRiC")) + geom_bar(aes(fill = factor(peaks_group)), color = "black", size = 0.2, position = position_fill(reverse = T), show.legend = T) +
  geom_bar(data = localization_ssb[localization == "Cytonuclear"], aes(x="2.SSB", fill = factor(peaks_group)), color = "black", size = 0.2, position = position_fill(reverse = T), show.legend = T) +
  scale_fill_brewer(type = "seq", palette = "Greys", direction = -1, labels = c("1","2","3","4",">5"), name = "", guide = guide_legend(keywidth = 2, keyheight = 2)) +
  scale_x_discrete(labels = c("TRiC", "SSB")) + guides(fill = guide_legend(reverse = TRUE))
plot <- plot + theme_classic(20) + labs(y = "Fraction of substrates\nby # of recruitment sites", x = "") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/B.PeakNumber_Cytonuclear.pdf", plot, width = 5, height = 4, dpi = 300, useDingbats = F)


### Plot traces from start codon
plot <- ggplot(data = tric_substrates5_dt[length >= 300, .(position, tric = movingAverage(log2(tric_odds_ma), n=15, center=T))]) + xlim(0, 300) +
  geom_hline(yintercept = 0, color = 'gray50', linetype = 'dashed', size = 1.25) + geom_vline(xintercept = 30, color = 'gray50', linetype = 'dashed', size = 1.25) +
  stat_summary(aes(position, tric), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = '#1F78B4',
               fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(position, tric, color = "#1F78B4"), fun.y = "median", geom = "line", size = 1.25) +
  scale_color_manual(labels = c("TRiC substrates"), values = c("#1F78B4"), name = "") +
  coord_cartesian(ylim=c(-0.6,0.5)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5))
G <- plot + theme_classic(20) + labs(y = "TRiC enrichment\n(log2 odds ratio)", x = "Codons from start") +
  theme(legend.position = c(.35,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/C.TRiCsubstrates_start.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = ssb_substrates5_dt[length >= 300, .(position, ssb = movingAverage(log2(ssb_odds_ma), n=15, center=T))]) + xlim(0, 300) +
  geom_hline(yintercept = 0, color = 'gray50', linetype = 'dashed', size = 1.25) + geom_vline(xintercept = 30, color = 'gray50', linetype = 'dashed', size = 1.25) +
  stat_summary(aes(position, ssb), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = "#33A02C",
               fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(position, ssb, color = "#33A02C"), fun.y = "median", geom = "line", size = 1.25) +
  scale_color_manual(labels = c("Ssb substrates"), values = c("#33A02C"), name = "") +
  coord_cartesian(ylim=c(-2,0.25))
G <- plot + theme_classic(20) + labs(y = "Ssb enrichment\n(log2 odds ratio)", x = "Codons from start") +
  theme(legend.position = c(.45,.35), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/D.SSBsubstrates_start.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


### Plot traces from stop codon
plot <- ggplot(data = tric_substrates5_dt[length >= 300, .(stopdist, tric = movingAverage(log2(tric_odds_ma), n=15, center=T))]) + xlim(-300, 0) +
  geom_hline(yintercept = 0, color = 'gray50', linetype = 'dashed', size = 1.25) +
  stat_summary(aes(stopdist, tric), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = '#1F78B4',
               fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(stopdist, tric), color = "#1F78B4", fun.y = "median", geom = "line", size = 1.25) +
  coord_cartesian(ylim=c(-0.6,0.5)) + scale_y_continuous(position = "right", breaks = c(-0.5,0,0.5))
G <- plot + theme_classic(20) + labs(y = "", x = "Codons from stop") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/C.TRiCSubstrates_stop.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = ssb_substrates5_dt[length >= 300, .(stopdist, ssb = movingAverage(log2(ssb_odds_ma), n=15, center=T))]) + xlim(-300, 0) +
  geom_hline(yintercept = 0, color = 'gray50', linetype = 'dashed', size = 1.25) + 
  stat_summary(aes(stopdist, ssb), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = "#33A02C",
               fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(stopdist, ssb), color = "#33A02C", fun.y = "median", geom = "line", size = 1.25) +
  coord_cartesian(ylim=c(-2,0.25)) + scale_y_continuous(position = "right")
G <- plot + theme_classic(20) + labs(y = "", x = "Codons from stop") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/D.SSBsubstrates_stop.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


### Significant residues by normalized position
i <- cbind(match(tric_peaks5all$position1, tric_substrates5_dt$position1))
tric_peaks5all <- cbind(tric_peaks5all, position_norm = tric_substrates5_dt[i]$position_norm)
tric_peaks5all <- cbind(tric_peaks5all, tric_odds_ma = tric_substrates5_dt[i]$tric_odds_ma)
plot <- ggplot(tric_substrates5_dt[(!tric_substrates5_dt$position1 %in% tric_peaks5all$position1) & position_norm >= 0 & position_norm <= 1 & tric_odds_ma < Inf], aes(position_norm, tric_odds_ma)) + geom_point(color = 'gray75', size = 2) +
  geom_vline(xintercept = 0.5, color = 'gray50', linetype = 'dashed', size = 1.25) +
  geom_point(data = tric_peaks5all[position_norm >= 0 & position_norm <= 1 & tric_odds_ma < Inf], aes(position_norm, tric_odds_ma, color = tric_padj), size = 2) + 
  scale_color_distiller(type = "seq", palette = "Blues", limits = c(0,0.05), name = "Padj", guide = guide_colorbar(barwidth = 1, barheight = 5)) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,1), labels = c("0" = "Start", "1" = "Stop")) +
  ylim(0,50)
plot <- plot + theme_classic(20) + labs(y = "TRiC enrichment (odds ratio)", x = "Norm. codon position") +
  theme(legend.text = element_text(size=20), legend.title = element_text(size=20),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/E.TRiC_allsubstrates.tiff", plot, width = 6, height = 4, dpi = 300)

i <- cbind(match(ssb_peaks5all$position1, ssb_substrates5_dt$position1))
ssb_peaks5all <- cbind(ssb_peaks5all, position_norm = ssb_substrates5_dt[i]$position_norm)
ssb_peaks5all <- cbind(ssb_peaks5all, ssb_odds_ma = ssb_substrates5_dt[i]$ssb_odds_ma)
plot <- ggplot(ssb_substrates5_dt[(!ssb_substrates5_dt$position1 %in% ssb_peaks5all$position1) & position_norm >= 0 & position_norm <= 1 & ssb_odds_ma < Inf], aes(position_norm, ssb_odds_ma)) + geom_point(color = 'gray75', size = 2) +
  geom_vline(xintercept = 0.5, color = 'gray50', linetype = 'dashed', size = 1.25) +
  geom_point(data = ssb_peaks5all[position_norm >= 0 & position_norm <= 1 & ssb_odds_ma < Inf], aes(position_norm, ssb_odds_ma, color = ssb_padj), size = 2) +
  scale_color_distiller(type = "seq", palette = "Greens", limits = c(0,0.05), name = "Adjusted\np-value", guide = guide_colorbar(barwidth = 1, barheight = 5)) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,1), labels = c("0" = "Start", "1" = "Stop")) +
  ylim(0,20)
plot <- plot + theme_classic(20) + labs(y = "Ssb enrichment (odds ratio)", x = "Norm. codon position") +
  theme(legend.text = element_text(size=20), legend.title = element_text(size=20),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/F.SSB_allsubstrates.tiff", plot, width = 6, height = 4, dpi = 300)


### Positional enrichment by localization
localization_tric1 <- localization[tric_peaks5]
localization_ssb1 <- localization[ssb_peaks5]

plot <- ggplot(localization_tric1[localization != "Exception" & localization != "SSTMD" & localization != "TA"], aes(localization, position_norm)) + 
  geom_hline(yintercept = 0.5, color = 'gray50', linetype = 'dashed', size = 1.25) +
  geom_violin(aes(fill = factor(localization1))) +
  scale_fill_manual(values = c("Cytonuclear"="#377EB8", "ER-targeted"="#377EB8", "Mitochondria"="#377EB8"), name = "") +
  coord_flip() + scale_x_discrete(limits = rev(unique(sort(localization_tric1[localization != "Exception" & localization != "SSTMD" & localization != "TA"]$localization)))) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(20) + labs(y = "Norm. codon position", x = "") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text = element_text(color = "black", size = 20),
        legend.position = "none")
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/G.TRiCAllpeaks_localization1_pvalueSScyto0.0113.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(localization_tric1[localization == "Mitochondria"]$position_norm, localization_tric1[localization == "TMD"]$position_norm, alternative = 't')

plot <- ggplot(localization_ssb1[localization != "Exception" & localization != "SSTMD" & localization != "TA"], aes(localization, position_norm)) + 
  geom_hline(yintercept = 0.5, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_violin(aes(fill = factor(localization1))) +
  scale_fill_manual(values = c("Cytonuclear"="#4DAF4A", "ER-targeted"="#4DAF4A", "Mitochondria"="#4DAF4A"), name = "") +
  coord_flip() + scale_x_discrete(limits = rev(unique(sort(localization_ssb1[localization != "Exception" & localization != "SSTMD" & localization != "TA"]$localization)))) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(20) + labs(y = "Norm. codon position", x = "") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text = element_text(color = "black", size = 20),
        legend.position = "none")
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/H.SsbAllpeaks_localization1.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(localization_ssb1[localization == "TMD"]$position_norm, localization_ssb1[localization == "Mitochondria"]$position_norm, alternative = 't')


### Heat map of binding sites in shared substrates
plot <- ggplot(ssb_peaks5_max[ssb_peaks5_max$orf %in% tric_peaks5_max$orf], aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray90", color = "gray90") +
  geom_hline(yintercept = 30, color = 'gray50', linetype = 'dashed', size = 1.5) +
  geom_point(data = ssb_peaks5[ssb_peaks5$orf %in% tric_peaks5_max$orf], aes(x = reorder(orf, -length), y = peak, color = "2SSB1"), alpha = 0.6, size = 1.75) +
  geom_point(data = tric_peaks5[tric_peaks5$orf %in% ssb_peaks5_max$orf], aes(x = reorder(orf, -length), y = peak, color = "1TRiC1"), size = 1.75) + 
  coord_flip(ylim = c(0,1000)) +
  scale_color_manual(values = allcolors, labels = c("1TRiC1"="TRiC association","2SSB1"="Ssb association"), name = "") +
  scale_y_continuous(breaks = c(0, 500, 1000), labels = c(0, 500, 1000))
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.35,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 1, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/I.SharedSubstrates_AllBindingSites.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(ssb_peaks5[ssb_peaks5$orf %in% tric_peaks5_max$orf], aes(position_norm, fill = "2SSB1", color = "2SSB1")) + 
  geom_vline(xintercept = 0.5, color = 'gray50', linetype = 'dashed', size = 1.25) +
  geom_density(size = 1.25, alpha = 0.3) +
  geom_density(data = tric_peaks5[tric_peaks5$orf %in% ssb_peaks5_max$orf], aes(position_norm, fill = '1TRiC1', color = "1TRiC1"), size = 1.25, alpha = 0.3) + 
  scale_fill_manual(values = allcolors, labels = c("TRiC association","SSB association"), name = "", guide = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
  scale_color_manual(values = allcolors, labels = c("TRiC association","SSB association"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(30) + labs(x = "Norm. codon position", y = "Density") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text.x = element_text(color = "black", size = 30),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig2/J.SharedSubstrates_AllBindingSites_density.pdf", plot, width = 5, height = 4, dpi = 300, useDingbats = F)
wilcox.test(tric_peaks5[tric_peaks5$orf %in% ssb_peaks5_max$orf]$position_norm, ssb_peaks5[ssb_peaks5$orf %in% tric_peaks5_max$orf]$position_norm)
wilcox.test(tric_peaks5_max[tric_peaks5_max$orf %in% ssb_peaks5_max$orf]$position_norm, ssb_peaks5_max[ssb_peaks5_max$orf %in% tric_peaks5_max$orf]$position_norm)

