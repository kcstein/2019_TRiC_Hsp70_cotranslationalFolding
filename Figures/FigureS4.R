### Peak number by subset
plot <- ggplot(localization_tric[localization != "Exception" & localization != "TA" & localization != "SSTMD"], aes(x=localization)) + geom_bar(aes(fill = factor(peaks_group)), color = "black", size = 0.2, position = position_fill(reverse = T), show.legend = T) +
  scale_fill_brewer(type = "seq", palette = "Blues", direction = -1, labels = c("1","2","3","4",">5"), name = "") +
  scale_x_discrete(labels = c("Cytonuc", "Mito", "SS", "TMD")) + guides(fill = guide_legend(reverse = TRUE))
plot <- plot + theme_classic(20) + labs(y = "Fraction of substrates\nby # of recruitment sites", x = "") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/A.PeakNumber_TRiC.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(localization_ssb[localization != "Exception" & localization != "TA" & localization != "SSTMD"], aes(x=localization)) + geom_bar(aes(fill = factor(peaks_group)), color = "black", size = 0.2, position = position_fill(reverse = T), show.legend = T) +
  scale_fill_brewer(type = "seq", palette = "Greens", direction = -1, labels = c("1","2","3","4",">5"), name = "") +
  scale_x_discrete(labels = c("Cytonuc", "Mito", "SS", "TMD")) + guides(fill = guide_legend(reverse = TRUE))
plot <- plot + theme_classic(20) + labs(y = "Fraction of substrates\nby # of recruitment sites", x = "") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/A.PeakNumber_SSB.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(localization_tric[localization == "SS"]$PeaksPerOrf, localization_tric[localization == "Cytonuclear"]$PeaksPerOrf, alternative = 't')
wilcox.test(localization_ssb[localization == "Mitochondria"]$PeaksPerOrf, localization_ssb[localization == "Cytonuclear"]$PeaksPerOrf, alternative = 't')


### Binding by substrate localization
localization_ssb_cyto <- localization_ssb[localization == "Cytonuclear"]
localization_ssb_mito <- localization_ssb[localization == "Mitochondria"]
localization_ssb_tmd <- localization_ssb[localization == "TMD"]
localization_ssb_ss <- localization_ssb[localization == "SS"]
localization_ssb1_cyto <- localization_ssb1[localization == "Cytonuclear"]
localization_ssb1_mito <- localization_ssb1[localization == "Mitochondria"]
localization_ssb1_tmd <- localization_ssb1[localization == "TMD"]
localization_ssb1_ss <- localization_ssb1[localization == "SS"]
localization_tric_cyto <- localization_tric[localization == "Cytonuclear"]
localization_tric_mito <- localization_tric[localization == "Mitochondria"]
localization_tric_tmd <- localization_tric[localization == "TMD"]
localization_tric_ss <- localization_tric[localization == "SS"]
localization_tric1_cyto <- localization_tric1[localization == "Cytonuclear"]
localization_tric1_mito <- localization_tric1[localization == "Mitochondria"]
localization_tric1_tmd <- localization_tric1[localization == "TMD"]
localization_tric1_ss <- localization_tric1[localization == "SS"]

### Heat maps and density for substrates by localization
# TRiC
plot <- ggplot(localization_tric_cyto, aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray90", color = "gray90") +
  geom_point(data = localization_tric1_cyto, aes(x = reorder(orf, -length), y = peak, color = "1TRiC1"), size = 1.75) + 
  coord_flip(ylim = c(0,1000)) +
  scale_color_manual(values = allcolors, labels = c("1TRiC1"="TRiC association"), name = "") +
  scale_y_continuous(breaks = c(0, 500, 1000), labels = c(0, 500, 1000))
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.35,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 1, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/B.TRiCcyto_AllBindingSites.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(localization_tric1_cyto, aes(position_norm, fill = "1TRiC1", color = "1TRiC1")) + geom_density(size = 1, alpha = 0.3) +
  scale_fill_manual(values = allcolors, labels = c("1TRiC1"="TRiC association"), name = "") +
  scale_color_manual(values = allcolors, labels = c("1TRiC1"="TRiC association"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(24) + labs(x = "Norm. codon position", y = "Density") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text.x = element_text(color = "black", size = 24),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/B.TRiCcyto_AllBindingSites_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(localization_tric_mito, aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray90", color = "gray90") +
  geom_point(data = localization_tric1_mito, aes(x = reorder(orf, -length), y = peak, color = "1TRiC1"), size = 2) + 
  coord_flip(ylim = c(0,1000)) +
  scale_color_manual(values = allcolors, labels = c("1TRiC1"="TRiC association"), name = "") +
  scale_y_continuous(breaks = c(0, 500, 1000), labels = c(0, 500, 1000))
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.35,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 6, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/C.TRiCmito_AllBindingSites.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(localization_tric1_mito, aes(position_norm, fill = "1TRiC1", color = "1TRiC1")) + geom_density(size = 1, alpha = 0.3) +
  scale_fill_manual(values = allcolors, labels = c("1TRiC1"="TRiC association"), name = "") +
  scale_color_manual(values = allcolors, labels = c("1TRiC1"="TRiC association"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(24) + labs(x = "Norm. codon position", y = "Density") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text.x = element_text(color = "black", size = 24),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/C.TRiCmito_AllBindingSites_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

test1 <- secretory.mito_protein_topology[secretory.mito_protein_topology$orf %in% localization_tric_ss$orf]
test1 <- test1[!is.na(start)]
plot <- ggplot(localization_tric_ss[localization_tric_ss$orf %in% test1$orf], aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray90", color = "gray90") +
  geom_segment(data = test1[(topology == "SS")], aes(x = reorder(orf, -length), xend = reorder(orf, -length), y = start, yend = end, color = "lO"), alpha = 0.8, size = 2.8, show.legend = F) +
  geom_point(data = localization_tric1_ss[localization_tric1_ss$orf %in% test1$orf], aes(x = reorder(orf, -length), y = peak, color = "1TRiC1"), size = 2) + 
  coord_flip() +
  scale_color_manual(values = allcolors, labels = c("1TRiC1"="TRiC association", "lO"="Signal sequence"), name = "") +
  scale_y_continuous(breaks = c(0, 500, 1000), labels = c(0, 500, 1000))
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.35,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 6, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/D.TRiCss_AllBindingSites.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(localization_tric1_ss, aes(position_norm, fill = "1TRiC1", color = "1TRiC1")) + geom_density(size = 1, alpha = 0.3) +
  scale_fill_manual(values = allcolors, labels = c("1TRiC1"="TRiC association"), name = "") +
  scale_color_manual(values = allcolors, labels = c("1TRiC1"="TRiC association"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(24) + labs(x = "Norm. codon position", y = "Density") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text.x = element_text(color = "black", size = 24),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/D.TRiCss_AllBindingSites_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

test1 <- secretory.mito_protein_topology[secretory.mito_protein_topology$orf %in% localization_tric_tmd$orf]
test1 <- test1[!is.na(start)]
plot <- ggplot(localization_tric_tmd[localization_tric_tmd$orf %in% test1$orf], aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray90", color = "gray90") +
  geom_segment(data = test1[(topology == "Transmembrane")], aes(x = reorder(orf, -length), xend = reorder(orf, -length), y = start, yend = end, color = "lO"), alpha = 0.6, size = 2.8, show.legend = F) +
  geom_segment(data = test1[(topology == "Cytoplasmic")], aes(x = reorder(orf, -length), xend = reorder(orf, -length), y = start, yend = end, color = "lP"), alpha = 0.8, size = 2.8, show.legend = F) +
  geom_point(data = localization_tric1_tmd[localization_tric1_tmd$orf %in% test1$orf], aes(x = reorder(orf, -length), y = peak, color = "1TRiC1"), size = 2) + 
  coord_flip() +
  scale_color_manual(values = allcolors, labels = c("1TRiC1"="TRiC association", "lP"="Cytoplasmic domain", "lO"="Transmembrane domain", name = "")) +
  scale_y_continuous(breaks = c(0, 500, 1000), labels = c(0, 500, 1000))
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.35,1.04), legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 6, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 20), legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/E.TRiCtmd_AllBindingSites.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(localization_tric1_tmd, aes(position_norm, fill = "1TRiC1", color = "1TRiC1")) + geom_density(size = 1, alpha = 0.3) +
  scale_fill_manual(values = allcolors, labels = c("1TRiC1"="TRiC association"), name = "") +
  scale_color_manual(values = allcolors, labels = c("1TRiC1"="TRiC association"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(24) + labs(x = "Norm. codon position", y = "Density") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text.x = element_text(color = "black", size = 24),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/E.TRiCtmd_AllBindingSites_density.pdf", plot, width = 6, height = 4, dpi = 300)

# SSB
plot <- ggplot(localization_ssb_cyto, aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray90", color = "gray90") +
  geom_point(data = localization_ssb1_cyto, aes(x = reorder(orf, -length), y = peak, color = "2SSB1"), size = 1.5) + 
  coord_flip(ylim = c(0,1000)) +
  scale_color_manual(values = allcolors, labels = c("2SSB1"="Ssb association"), name = "") +
  scale_y_continuous(breaks = c(0, 500, 1000), labels = c(0, 500, 1000))
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.35,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 1, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/B.SSBcyto_AllBindingSites.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(localization_ssb1_cyto, aes(position_norm, fill = "2SSB1", color = "2SSB1")) + geom_density(size = 1, alpha = 0.3) +
  scale_fill_manual(values = allcolors, labels = c("2SSB1"="Ssb association"), name = "") +
  scale_color_manual(values = allcolors, labels = c("2SSB1"="Ssb association"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(24) + labs(x = "Norm. codon position", y = "Density") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text.x = element_text(color = "black", size = 24),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/B.SSBcyto_AllBindingSites_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(localization_ssb_mito, aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray90", color = "gray90") +
  geom_point(data = localization_ssb1_mito, aes(x = reorder(orf, -length), y = peak, color = "2SSB1"), size = 2) + 
  coord_flip(ylim = c(0,1000)) +
  scale_color_manual(values = allcolors, labels = c("2SSB1"="Ssb association"), name = "") +
  scale_y_continuous(breaks = c(0, 500, 1000), labels = c(0, 500, 1000))
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.35,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 2, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/C.SSBmito_AllBindingSites.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(localization_ssb1_mito, aes(position_norm, fill = "2SSB1", color = "2SSB1")) + geom_density(size = 1, alpha = 0.3) +
  scale_fill_manual(values = allcolors, labels = c("2SSB1"="Ssb association"), name = "") +
  scale_color_manual(values = allcolors, labels = c("2SSB1"="Ssb association"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(24) + labs(x = "Norm. codon position", y = "Density") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text.x = element_text(color = "black", size = 24),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/C.SSBmito_AllBindingSites_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

test1 <- secretory.mito_protein_topology[secretory.mito_protein_topology$orf %in% localization_ssb_ss$orf]
test1 <- test1[!is.na(start)]
plot <- ggplot(localization_ssb_ss[localization_ssb_ss$orf %in% test1$orf], aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray90", color = "gray90") +
  geom_segment(data = test1[(topology == "SS")], aes(x = reorder(orf, -length), xend = reorder(orf, -length), y = start, yend = end, color = "lO"), alpha = 0.8, size = 2, show.legend = F) +
  geom_point(data = localization_ssb1_ss[localization_ssb1_ss$orf %in% test1$orf], aes(x = reorder(orf, -length), y = peak, color = "2SSB1"), size = 2) + 
  coord_flip() +
  scale_color_manual(values = allcolors, labels = c("2SSB1"="Ssb association", "lO"="Signal sequence"), name = "") +
  scale_y_continuous(breaks = c(0, 500, 1000), labels = c(0, 500, 1000))
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.35,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 4, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/D.SSBss_AllBindingSites.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(localization_ssb1_ss, aes(position_norm, fill = "2SSB1", color = "2SSB1")) + geom_density(size = 1, alpha = 0.3) +
  scale_fill_manual(values = allcolors, labels = c("2SSB1"="Ssb association"), name = "") +
  scale_color_manual(values = allcolors, labels = c("2SSB1"="Ssb association"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(24) + labs(x = "Norm. codon position", y = "Density") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text.x = element_text(color = "black", size = 24),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/D.SSBss_AllBindingSites_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

test1 <- secretory.mito_protein_topology[secretory.mito_protein_topology$orf %in% localization_ssb_tmd$orf]
test1 <- test1[!is.na(start)]
plot <- ggplot(localization_ssb_tmd[localization_ssb_tmd$orf %in% test1$orf], aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray90", color = "gray90") +
  geom_segment(data = test1[(topology == "Transmembrane")], aes(x = reorder(orf, -length), xend = reorder(orf, -length), y = start, yend = end, color = "lO"), alpha = 0.6, size = 2, show.legend = F) +
  geom_segment(data = test1[(topology == "Cytoplasmic")], aes(x = reorder(orf, -length), xend = reorder(orf, -length), y = start, yend = end, color = "lP"), alpha = 0.8, size = 2, show.legend = F) +
  geom_point(data = localization_ssb1_tmd[localization_ssb1_tmd$orf %in% test1$orf], aes(x = reorder(orf, -length), y = peak, color = "2SSB1"), size = 2) + 
  coord_flip() +
  scale_color_manual(values = allcolors, labels = c("2SSB1"="Ssb association", "lP"="Cytoplasmic domain", "lO"="Transmembrane domain", name = "")) +
  scale_y_continuous(breaks = c(0, 500, 1000), labels = c(0, 500, 1000))
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.35,1.04), legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 4, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 16), legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/E.SSBtmd_AllBindingSites.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(localization_ssb1_tmd, aes(position_norm, fill = "2SSB1", color = "2SSB1")) + geom_density(size = 1, alpha = 0.3) +
  scale_fill_manual(values = allcolors, labels = c("2SSB1"="Ssb association"), name = "") +
  scale_color_manual(values = allcolors, labels = c("2SSB1"="Ssb association"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(24) + labs(x = "Norm. codon position", y = "Density") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text.x = element_text(color = "black", size = 24),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/E.SSBtmd_AllBindingSites_density.pdf", plot, width = 6, height = 4, dpi = 300)


### Metagene from start
plot <- ggplot(data = substrates_fishers[(substrates_fishers$orf %in% localization_ssb_cyto$orf) & length >= 300,
                                         .(position, ssb = movingAverage(log2(ssb_odds_ma), n=15, center=T))]) + xlim(0, 300) +
  geom_hline(yintercept = 0, color = 'gray50', linetype = 'dashed', size = 1) + geom_vline(xintercept = 30, color = 'gray50', linetype = 'dashed', size = 1) +
  stat_summary(aes(position, ssb), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = '#33A02C',
               fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(position, ssb, color = "SSB"), fun.y = "median", geom = "line", size = 1.25) +
  stat_summary(data = substrates_fishers[(substrates_fishers$orf %in% localization_tric_cyto$orf) & length >= 300,
                                         .(position, tric = movingAverage(log2(tric_odds_ma), n=15, center=T))],
               aes(position, tric), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = '#1F78B4',
               fun.args=list(conf.int=0.5)) + 
  stat_summary(data = substrates_fishers[(substrates_fishers$orf %in% localization_tric_cyto$orf) & length >= 300,
                                         .(position, tric = movingAverage(log2(tric_odds_ma), n=15, center=T))],
               aes(position, tric, color = "TRiC"), fun.y = "median", geom = "line", size = 1.25) +
  scale_color_manual(labels = c("Ssb substrates", "TRiC substrates"), values = cols, name = "") +
  guides(color = guide_legend(reverse = TRUE))
G <- plot + theme_classic(20) + labs(y = "Enrichment (odds ratio)", x = "Codons from start") +
  theme(legend.position = c(.55,.35), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/B.AllSubstrates_cyto_start.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = substrates_fishers[(substrates_fishers$orf %in% localization_ssb_mito$orf) & length >= 300,
                                         .(position, ssb = movingAverage(log2(ssb_odds_ma), n=15, center=T))]) + xlim(0, 300) +
  geom_hline(yintercept = 0, color = 'gray50', linetype = 'dashed', size = 1) + geom_vline(xintercept = 30, color = 'gray50', linetype = 'dashed', size = 1) +
  stat_summary(aes(position, ssb), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = '#33A02C',
               fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(position, ssb, color = "SSB"), fun.y = "median", geom = "line", size = 1.25) +
  stat_summary(data = substrates_fishers[(substrates_fishers$orf %in% localization_tric_mito$orf) & length >= 300,
                                         .(position, tric = movingAverage(log2(tric_odds_ma), n=15, center=T))],
               aes(position, tric), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = '#1F78B4',
               fun.args=list(conf.int=0.5)) + 
  stat_summary(data = substrates_fishers[(substrates_fishers$orf %in% localization_tric_mito$orf) & length >= 300,
                                         .(position, tric = movingAverage(log2(tric_odds_ma), n=15, center=T))],
               aes(position, tric, color = "TRiC"), fun.y = "median", geom = "line", size = 1.25) +
  scale_color_manual(labels = c("Ssb substrates", "TRiC substrates"), values = cols, name = "") +
  guides(color = guide_legend(reverse = TRUE))
G <- plot + theme_classic(20) + labs(y = "Enrichment (odds ratio)", x = "Codons from start") +
  theme(legend.position = c(.55,.35), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/C.AllSubstrates_mito_start.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = substrates_fishers[(substrates_fishers$orf %in% localization_ssb_ss$orf) & length >= 300,
                                         .(position, ssb = movingAverage(log2(ssb_odds_ma), n=15, center=T))]) + xlim(0, 300) +
  geom_hline(yintercept = 0, color = 'gray50', linetype = 'dashed', size = 1) + geom_vline(xintercept = 30, color = 'gray50', linetype = 'dashed', size = 1) +
  stat_summary(aes(position, ssb), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = '#33A02C',
               fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(position, ssb, color = "SSB"), fun.y = "median", geom = "line", size = 1.25) +
  stat_summary(data = substrates_fishers[(substrates_fishers$orf %in% localization_tric_ss$orf) & length >= 300,
                                         .(position, tric = movingAverage(log2(tric_odds_ma), n=15, center=T))],
               aes(position, tric), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = '#1F78B4',
               fun.args=list(conf.int=0.5)) + 
  stat_summary(data = substrates_fishers[(substrates_fishers$orf %in% localization_tric_ss$orf) & length >= 300,
                                         .(position, tric = movingAverage(log2(tric_odds_ma), n=15, center=T))],
               aes(position, tric, color = "TRiC"), fun.y = "median", geom = "line", size = 1.25) +
  scale_color_manual(labels = c("Ssb substrates", "TRiC substrates"), values = cols, name = "") +
  guides(color = guide_legend(reverse = TRUE))
G <- plot + theme_classic(20) + labs(y = "Enrichment (odds ratio)", x = "Codons from start") +
  theme(legend.position = c(.55,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/D.AllSubstrates_ss_start.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = substrates_fishers[(substrates_fishers$orf %in% localization_ssb_tmd$orf) & length >= 300,
                                         .(position, ssb = movingAverage(log2(ssb_odds_ma), n=15, center=T))]) + xlim(0, 300) +
  geom_hline(yintercept = 0, color = 'gray50', linetype = 'dashed', size = 1) + geom_vline(xintercept = 30, color = 'gray50', linetype = 'dashed', size = 1) +
  stat_summary(aes(position, ssb), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = '#33A02C',
               fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(position, ssb, color = "SSB"), fun.y = "median", geom = "line", size = 1.25) +
  stat_summary(data = substrates_fishers[(substrates_fishers$orf %in% localization_tric_tmd$orf) & length >= 300,
                                         .(position, tric = movingAverage(log2(tric_odds_ma), n=15, center=T))],
               aes(position, tric), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3, fill = '#1F78B4',
               fun.args=list(conf.int=0.5)) + 
  stat_summary(data = substrates_fishers[(substrates_fishers$orf %in% localization_tric_tmd$orf) & length >= 300,
                                         .(position, tric = movingAverage(log2(tric_odds_ma), n=15, center=T))],
               aes(position, tric, color = "TRiC"), fun.y = "median", geom = "line", size = 1.25) +
  scale_color_manual(labels = c("Ssb substrates", "TRiC substrates"), values = cols, name = "") +
  guides(color = guide_legend(reverse = TRUE))
G <- plot + theme_classic(20) + labs(y = "Enrichment (odds ratio)", x = "Codons from start") +
  theme(legend.position = c(.55,.35), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS4/E.AllSubstrates_tmd_start.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

