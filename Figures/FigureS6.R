### Heat map of substrates having domains of interest
# Rossmann (3.40.50 or SSF52540); WD (2.130.10 or SSF50978)
temp <- as.data.table(unique(domains_cath_tric5_max[Gene3D_topo == "2.130.10"]$orf))
temp4 <- tric_peaks5_max[tric_peaks5_max$orf %in% temp$V1]
temp5 <- ssb_peaks5_max[ssb_peaks5_max$orf %in% temp$V1]
plot <- ggplot(temp4, aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray80") +
  geom_segment(data = domains_cath_tric5_max[(domains_cath_tric5_max$orf %in% temp$V1) & Gene3D_topo == "2.130.10"], aes(x = reorder(orf, -length), xend = reorder(orf, -length), y = Start, yend = End, color = "lO"), alpha = 0.6, show.legend = F, size = 3) +
  geom_point(data = ssb_peaks5[ssb_peaks5$orf %in% temp5$orf], aes(x = reorder(orf, -length), y = peak, color = "2SSB1"), alpha = 0.6, size = 2) + 
  geom_point(data = tric_peaks5[tric_peaks5$orf %in% temp4$orf], aes(x = reorder(orf, -length), y = peak, color = "1TRiC1"), size = 2) + 
  coord_flip(ylim = c(0,1000)) + 
  scale_color_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC association", "2SSB1"="Ssb association","lO"="WD40 repeats"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.45,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 6, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/FigS6/B.WD_AllBindingSites_CATH.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

temp <- as.data.table(unique(domains_cath_tric5_max[Gene3D_topo == "3.40.50"]$orf))
temp4 <- tric_peaks5_max[tric_peaks5_max$orf %in% temp$V1]
temp5 <- ssb_peaks5_max[ssb_peaks5_max$orf %in% temp$V1]
plot <- ggplot(temp4, aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray90", color = "gray90") +
  geom_segment(data = domains_cath_tric5_max[(domains_cath_tric5_max$orf %in% temp$V1) & Gene3D_topo == "3.40.50"], aes(x = reorder(orf, -length), xend = reorder(orf, -length), y = Start, yend = End, color = "lP"), show.legend = F, alpha = 0.8, size = 1) +
  geom_point(data = ssb_peaks5[ssb_peaks5$orf %in% temp5$orf], aes(x = reorder(orf, -length), y = peak, color = "2SSB1"), alpha = 0.6, size = 2) + 
  geom_point(data = tric_peaks5[tric_peaks5$orf %in% temp4$orf], aes(x = reorder(orf, -length), y = peak, color = "1TRiC1"), size = 2) + 
  coord_flip(ylim = c(0,1000)) + 
  scale_color_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC association", "2SSB1"="Ssb association", "lP"="Rossmann fold"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.45,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 2, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/FigS6/C.Rossmann_AllBindingSites_CATH.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


# Significance of WD vs Rossmann
temp <- as.data.table(unique(domains_cath_tric5_max[Gene3D_topo == "3.40.50"]$orf))
temp1 <- as.data.table(unique(domains_cath_tric5_max[Gene3D_topo == "2.130.10"]$orf))
temp2 <- as.data.table(unique(domains_cath_ssb5_max[Gene3D_topo == "3.40.50"]$orf))
temp3 <- as.data.table(unique(domains_cath_ssb5_max[Gene3D_topo == "2.130.10"]$orf))
wilcox.test(tric_peaks5[(tric_peaks5$orf %in% temp$V1)]$position_norm, tric_peaks5[(tric_peaks5$orf %in% temp1$V1)]$position_norm)
wilcox.test(tric_peaks5_max[(tric_peaks5_max$orf %in% temp$V1)]$position_norm, tric_peaks5_max[(tric_peaks5_max$orf %in% temp1$V1)]$position_norm)
wilcox.test(ssb_peaks5[(ssb_peaks5$orf %in% temp2$V1)]$position_norm, ssb_peaks5[(ssb_peaks5$orf %in% temp3$V1)]$position_norm)
wilcox.test(ssb_peaks5_max[(ssb_peaks5_max$orf %in% temp2$V1)]$position_norm, ssb_peaks5_max[(ssb_peaks5_max$orf %in% temp3$V1)]$position_norm)
wilcox.test(tric_peaks5[(tric_peaks5$orf %in% temp$V1)]$position_norm, ssb_peaks5[(ssb_peaks5$orf %in% temp2$V1)]$position_norm)
wilcox.test(tric_peaks5_max[(tric_peaks5_max$orf %in% temp$V1)]$position_norm, ssb_peaks5_max[(ssb_peaks5_max$orf %in% temp2$V1)]$position_norm)
wilcox.test(tric_peaks5[(tric_peaks5$orf %in% temp1$V1)]$position_norm, ssb_peaks5[(ssb_peaks5$orf %in% temp3$V1)]$position_norm)
wilcox.test(tric_peaks5_max[(tric_peaks5_max$orf %in% temp1$V1)]$position_norm, ssb_peaks5_max[(ssb_peaks5_max$orf %in% temp3$V1)]$position_norm)
plot <- ggplot(ssb_peaks5[(ssb_peaks5$orf %in% temp3$V1)], aes(position_norm, fill = "lO", color = "lO")) + geom_density(size = 1, alpha = 0.7) +
  geom_density(data = ssb_peaks5[(ssb_peaks5$orf %in% temp2$V1)], aes(position_norm, fill = "lP", color = "lP"), size = 1, alpha = 0.7) + 
  geom_density(data = tric_peaks5[(tric_peaks5$orf %in% temp1$V1)], aes(position_norm, fill = "dO", color = "dO"), size = 1, alpha = 0.3) +
  geom_density(data = tric_peaks5[(tric_peaks5$orf %in% temp$V1)], aes(position_norm, fill = 'dP', color = "dP"), size = 1, alpha = 0.3) +
  scale_fill_manual(values = allcolors, labels = c("dO"="TRiC WD40","dP"="TRiC Rossmann","lO"="Ssb WD40","lP"="Ssb Rossmann"), name = "", guide = guide_legend(keywidth = 1, keyheight = 1)) +
  scale_color_manual(values = allcolors, labels = c("dO"="TRiC WD40","dP"="TRiC Rossmann","lO"="Ssb WD40","lP"="Ssb Rossmann"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(20) + labs(x = "Normalized codon position", y = "Density") +
  theme(legend.position = c(.01,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(color = "black", size = 20), legend.text = element_text(size = 20),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS6/D.TRiCandSSB_AllBindingSites_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Significance of early versus late TRiC binding for Rossmann and WD40
temp <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/domains_freq.csv", stringsAsFactors = T, header = T)
temp <- as.data.table(temp)
temp$Description <- factor(temp$Description,levels = c("Rossmann fold","WD40"))
plot <- ggplot(temp[(Chaperone == "Early" | Chaperone == "Late") & !is.na(Description)], aes(x = Description,y=log2, fill = factor(Chaperone))) + geom_col(position = "dodge", color = "black", size = 0.5) + 
  scale_fill_manual(limits = c("Early","Late"), labels = c("Early", "Late"), values = c("#1F78B4", "#A6CEE3"), name = "", guide = guide_legend(keywidth = 1, keyheight = 1)) +
  scale_x_discrete(limits = c("Rossmann fold", "WD40"), labels = c("Rossmann", "WD40"))
plot <- plot + theme_classic(20) + labs(y = "log2 domain enrichment", x = "") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS6/E.EarlyLateDomains.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Rossmann domains heat map and density
temp <- domains_cath_tric5_start[Gene3D_topo == "3.40.50"]
setkeyv(temp, c("domain_alias"))
temp1 <- temp[, .SD[which.min(Start)], by = domain_alias]
temp1[, new.length := length - Start + 1]
plot <- ggplot(temp1, aes(x = reorder(domain_alias, -DomainLength), y = new.length)) + geom_col(fill = "gray90") +
  geom_col(data = temp1, aes(x = reorder(domain_alias, -DomainLength), y = DomainLength), fill = "#CAB2D6", alpha = 0.6, show.legend = F) +
  geom_point(data = domains_cath_tric5_start[domains_cath_tric5_start$domain_alias %in% temp1$domain_alias], aes(x = reorder(domain_alias, -DomainLength), y = peak.start, color = "1TRiC1"), size = 2) + 
  coord_flip(ylim = c(0,600)) + 
  scale_color_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC association"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Distance from domain start (codons)", x = "Domain") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS6/F.TRiC_RossmannDomains_AllBindingSites_CATH.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(domains_cath_ssb5_start[Gene3D_topo == "3.40.50"], aes(length.ratio, fill = '2SSB1', color = "2SSB1")) + 
  geom_vline(xintercept = 1, color = 'gray70', linetype = 'dashed', size = 1.25) +
  geom_density(size = 1, alpha = 0.3) +
  geom_density(data = domains_cath_tric5_start[Gene3D_topo == "3.40.50"], aes(length.ratio, fill = "1TRiC1", color = "1TRiC1"), size = 1, alpha = 0.3) + 
  scale_fill_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC","2SSB1"="Ssb"), name = "", guide = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
  scale_color_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC","2SSB1"="Ssb"), name = "") +
  scale_x_continuous(limits = c(0,2), breaks = c(0, 1), labels = c("0" = "Domain start", "1" = "Domain end"))
plot <- plot + theme_classic(24) + labs(x = "Normalized domain position", y = "Density") +
  theme(legend.position = c(.65,.95), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text = element_text(color = "black", size = 20), legend.text = element_text(size = 20), legend.background = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS6/F.RossmannDomainsCATH_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(domains_cath_tric5_start[Gene3D_topo == "3.40.50"]$length.ratio, domains_cath_ssb5_start[Gene3D_topo == "3.40.50"]$length.ratio)

temp <- tric_peaks_all[orf == "YGR145W", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 600), color = "gray70", size = 2) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 175), color = "#CAB2D6", alpha = 0.7, size = 8) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 30), color = "#6A3D9A", alpha = 0.7, size = 8) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
        axis.text = element_blank(), plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS6/F.DomainSchematic.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### WD domains heat map and density
temp <- domains_cath_tric5_start[Gene3D_topo == "2.130.10"]
setkeyv(temp, c("domain_alias"))
temp1 <- temp[, .SD[which.min(Start)], by = domain_alias]
temp1[, new.length := length - Start + 1]
plot <- ggplot(temp1, aes(x = reorder(domain_alias, -DomainLength), y = new.length)) + geom_col(fill = "gray90") +
  geom_col(data = temp1, aes(x = reorder(domain_alias, -DomainLength), y = DomainLength), fill = "#FDBF6F", alpha = 0.6, show.legend = F) +
  geom_point(data = domains_cath_tric5_start[domains_cath_tric5_start$domain_alias %in% temp1$domain_alias], aes(x = reorder(domain_alias, -DomainLength), y = peak.start, color = "1TRiC1"), size = 2) + 
  coord_flip(ylim = c(0,600)) + 
  scale_color_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC association"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Distance from domain start (codons)", x = "Domain") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS6/G.TRiC_WD_AllBindingSites_CATH.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(domains_cath_ssb5_start[Gene3D_topo == "2.130.10"], aes(length.ratio, fill = '2SSB1', color = "2SSB1")) + 
  geom_vline(xintercept = 1, color = 'gray70', linetype = 'dashed', size = 1.25) +
  geom_density(size = 1, alpha = 0.3) +
  geom_density(data = domains_cath_tric5_start[Gene3D_topo == "2.130.10"], aes(length.ratio, fill = "1TRiC1", color = "1TRiC1"), size = 1, alpha = 0.3) + 
  scale_fill_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC","2SSB1"="Ssb"), name = "", guide = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
  scale_color_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC","2SSB1"="Ssb"), name = "") +
  scale_x_continuous(limits = c(0,2), breaks = c(0, 1), labels = c("0" = "Domain start", "1" = "Domain end"))
plot <- plot + theme_classic(24) + labs(x = "Normalized domain position", y = "Density") +
  theme(legend.position = c(.65,.95), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text = element_text(color = "black", size = 20), legend.text = element_text(size = 20), legend.background = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS6/G.WD_CATH_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(domains_cath_tric5_start[Gene3D_topo == "2.130.10"]$length.ratio, domains_cath_ssb5_start[Gene3D_topo == "2.130.10"]$length.ratio)

temp <- tric_peaks_all[orf == "YGR145W", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 600), color = "gray70", size = 2) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 175), color = "#FDBF6F", alpha = 0.7, size = 8) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 30), color = "#FF7F00", alpha = 0.7, size = 8) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigSS6/G.DomainSchematic.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### ENP2 in Bukau dataset
plot <- ggplot(aes(position, occupancy), data=temp[orf == "YGR145W" & ssb_odds_ma < Inf, .(position, occupancy = movingAverage(ssb_odds, n=5, center=T))]) + 
  geom_vline(xintercept = 84, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = 126, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 175, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 208, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 256, color = 'gray50', linetype = 'dashed', size = 1) + 
  geom_vline(xintercept = 299, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_vline(xintercept = 342, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_line(aes(color = "1"),size = 1.25) +
  scale_color_manual(values = c("1" = "#6A3D9A"), labels = c("Ssb"), name = "")
G <- plot + theme_classic(20) + labs(y = "Ssb enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = c(.55,.9), legend.justification = c("left", "top"), axis.text = element_text(size = 16, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("ENP2")
ggsave("/Users/KevinStein/Desktop/Figures/FigS6/H.ENP2_Bukau.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)
