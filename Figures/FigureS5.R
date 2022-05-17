### Enrichment around domains
plot <- ggplot(domains_cath_tric5_align[tric_odds_ma < Inf], aes(adjusted_end, tric_odds)) + 
  geom_vline(xintercept = 0, color = 'gray50', linetype = "dashed", size = 1.25) +
  geom_point(color = 'gray85') +
  geom_point(data = domains_cath_tric5_align[(domains_cath_tric5_align$position1 %in% tric_peaks5all$position1)], aes(adjusted_end, tric_odds, color = tric_padj), size = 2) +
  scale_color_distiller(type = "seq", palette = "Blues", limits = c(0,0.05), name = "Padj", guide = guide_colorbar(barwidth = 1, barheight = 5)) +
  xlim(-200,200) + ylim(0,75)
plot <- plot + theme_classic(20) + labs(y = "TRiC enrichment (odds ratio)", x = "Distance from\ndomain end (codons)") +
  theme(legend.text = element_text(size=20), legend.title = element_text(size=20),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS5/A.TRiC_CATHdomaincentered_All.tiff", plot, width = 6, height = 4, dpi = 300)
plot <- ggplot(domains_cath_ssb5_align[ssb_odds_ma < Inf], aes(adjusted_end, ssb_odds)) + 
  geom_vline(xintercept = 0, color = 'gray50', linetype = "dashed", size = 1.25) +
  geom_point(color = 'gray85') +
  geom_point(data = domains_cath_ssb5_align[(domains_cath_ssb5_align$position1 %in% ssb_peaks5all$position1)], aes(adjusted_end, ssb_odds, color = ssb_padj), size = 2) +
  scale_color_distiller(type = "seq", palette = "Greens", limits = c(0,0.05), name = "Padj", guide = guide_colorbar(barwidth = 1, barheight = 5)) +
  xlim(-200,200) + ylim(0,75)
plot <- plot + theme_classic(20) + labs(y = "Ssb enrichment (odds ratio)", x = "Distance from\ndomain end (codons)") +
  theme(legend.text = element_text(size=20), legend.title = element_text(size=20),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS5/B.SSB_CATHdomaincentered_All.tiff", plot, width = 6, height = 4, dpi = 300)


### Domain significance
plot <- ggplot(domains_cath_tric5_start, aes(1, length.ratio, fill = "1TRiC1"))  + geom_boxplot(notch = T) +
  geom_boxplot(data = domains_cath_tric5_start, aes(2, random.length.ratio, fill = "1TRiC2"), notch = T) +
  geom_boxplot(data = domains_cath_ssb5_start, aes(3, length.ratio, fill = "2SSB1"), notch = T) +
  geom_boxplot(data = domains_cath_ssb5_start, aes(4, random.length.ratio, fill = "2SSB2"), notch = T) +
  geom_boxplot(data = domains_cath_tric5_max_start, aes(5, length.ratio, fill = "1TRiC1"), notch = T) +
  geom_boxplot(data = domains_cath_tric5_max_start, aes(6, random.length.ratio, fill = "1TRiC2"), notch = T) +
  geom_boxplot(data = domains_cath_ssb5_max_start, aes(7, length.ratio, fill = "2SSB1"), notch = T) +
  geom_boxplot(data = domains_cath_ssb5_max_start, aes(8, random.length.ratio, fill = "2SSB2"), notch = T) +
  geom_boxplot(data = domains_scop_tric5_start, aes(9, length.ratio, fill = "1TRiC1"), notch = T) +
  geom_boxplot(data = domains_scop_tric5_start, aes(10, random.length.ratio, fill = "1TRiC2"), notch = T) +
  geom_boxplot(data = domains_scop_ssb5_start, aes(11, length.ratio, fill = "2SSB1"), notch = T) +
  geom_boxplot(data = domains_scop_ssb5_start, aes(12, random.length.ratio, fill = "2SSB2"), notch = T) +
  geom_boxplot(data = domains_scop_tric5_max_start, aes(13, length.ratio, fill = "1TRiC1"), notch = T) +
  geom_boxplot(data = domains_scop_tric5_max_start, aes(14, random.length.ratio, fill = "1TRiC2"), notch = T) + coord_cartesian(ylim = c(0,2.5)) +
  geom_boxplot(data = domains_scop_ssb5_max_start, aes(15, length.ratio, fill = "2SSB1"), notch = T) +
  geom_boxplot(data = domains_scop_ssb5_max_start, aes(16, random.length.ratio, fill = "2SSB2"), notch = T) + coord_cartesian(ylim = c(0,2.5)) +
  scale_fill_manual(values = allcolors, labels = c("1TRiC1"="TRiC association","1TRiC2"="TRiC random sampling",
                                                   "2SSB1"="SSB association","2SSB2"="SSB random sampling"), name = "") +
  scale_x_continuous(breaks = c(2.5,6.5,10.5,14.5), labels = c("CATH\nall sites","CATH\nZmax","SCOP\nall sites","SCOP\nZmax"))
plot <- plot + theme_classic(20) + labs(y = "Norm. domain position", x = "") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text = element_text(color = "black", size = 20), legend.text = element_text(size = 10))
ggsave("/Users/KevinStein/Desktop/Figures/FigS5/C.Domains_boxplot.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Significance
wilcox.test(domains_cath_tric5_start$length.ratio, domains_cath_ssb5_start$length.ratio, alternative = 't')
wilcox.test(domains_cath_tric5_max_start$length.ratio, domains_cath_ssb5_max_start$length.ratio, alternative = 't')
wilcox.test(domains_scop_tric5_start$length.ratio, domains_scop_ssb5_start$length.ratio, alternative = 't')
wilcox.test(domains_scop_tric5_max_start$length.ratio, domains_scop_ssb5_max_start$length.ratio, alternative = 't')
wilcox.test(domains_cath_tric5_start$length.ratio, domains_cath_tric5_start$random.length.ratio, alternative = 't')
wilcox.test(domains_cath_tric5_max_start$length.ratio, domains_cath_tric5_max_start$random.length.ratio, alternative = 't')
wilcox.test(domains_scop_tric5_start$length.ratio, domains_scop_tric5_start$random.length.ratio)
wilcox.test(domains_scop_tric5_max_start$length.ratio, domains_scop_tric5_max_start$random.length.ratio)
wilcox.test(domains_cath_ssb5_start$length.ratio, domains_cath_ssb5_start$random.length.ratio)
wilcox.test(domains_cath_ssb5_max_start$length.ratio, domains_cath_ssb5_max_start$random.length.ratio)
wilcox.test(domains_scop_ssb5_start$length.ratio, domains_scop_ssb5_start$random.length.ratio)
wilcox.test(domains_scop_ssb5_max_start$length.ratio, domains_scop_ssb5_max_start$random.length.ratio)

temp <- wilcox.test(domains_cath_tric5_start$length.ratio, domains_cath_bukau5_start$length.ratio, alternative = 't')
temp$p.value
wilcox.test(domains_cath_tric5_max_start$length.ratio, domains_cath_bukau5_max_start$length.ratio, alternative = 't')
wilcox.test(domains_cath_ssb5_start$length.ratio, domains_cath_bukau5_start$length.ratio)
wilcox.test(domains_cath_ssb5_max_start$length.ratio, domains_cath_bukau5_max_start$length.ratio)


### Comparison to Bukau
plot <- ggplot(domains_cath_tric5_start, aes(1, length.ratio, fill = "1TRiC1"))  + geom_boxplot(notch = T) +
  geom_boxplot(data = domains_cath_bukau5_start, aes(2, length.ratio, fill = "dP"), notch = T) +
  geom_boxplot(data = domains_cath_tric5_max_start, aes(3, length.ratio, fill = "1TRiC1"), notch = T) +
  geom_boxplot(data = domains_cath_bukau5_max_start, aes(4, length.ratio, fill = "dP"), notch = T) +
  coord_cartesian(ylim = c(0,2.5)) +
  scale_fill_manual(values = allcolors, labels = c("1TRiC1"="TRiC association","2SSB1"="SSB association","dP"="Doring"), name = "") +
  scale_x_continuous(breaks = c(2,4), labels = c("CATH\nall sites","CATH\nZmax"))
plot <- plot + theme_classic(20) + labs(y = "Norm. domain position", x = "") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text = element_text(color = "black", size = 20), legend.text = element_text(size = 10))
ggsave("/Users/KevinStein/Desktop/Figures/FigS5/D.Domains_boxplot_Doring.pdf", plot, width = 4, height = 4, dpi = 300, useDingbats = F)


### Binding in "early" domains
temp <- domains_cath_tric5_max_start[domainend_norm < 0.75]
setkeyv(temp, c("domain_alias"))
temp1 <- temp[, .SD[which.min(Start)], by = domain_alias]
temp1[, new.length := length - Start + 1]
plot <- ggplot(temp1, aes(x = reorder(domain_alias, -DomainLength), y = new.length)) + geom_col(fill = "gray90") +
  geom_col(data = temp1, aes(x = reorder(domain_alias, -DomainLength), y = DomainLength), fill = "#E31A1C", alpha = 0.2, show.legend = F) +
  geom_point(data = domains_cath_tric5_start[domains_cath_tric5_start$domain_alias %in% temp1$domain_alias], aes(x = reorder(domain_alias, -DomainLength), y = peak.start, color = "1TRiC1"), size = 2) + 
  coord_flip(ylim = c(0,600)) + 
  scale_color_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC association"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Distance from domain start (codons)", x = "Domain") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS5/E.TRiC_EarlyDomains_AllBindingSites_CATH.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(domains_cath_ssb5_start[domainend_norm < 0.75], aes(length.ratio, fill = '2SSB1', color = "2SSB1")) + 
  geom_vline(xintercept = 1, color = 'gray70', linetype = 'dashed', size = 1.25) +
  geom_density(size = 1, alpha = 0.3) +
  geom_density(data = domains_cath_tric5_start[domainend_norm < 0.75], aes(length.ratio, fill = "1TRiC1", color = "1TRiC1"), size = 1, alpha = 0.3) + 
  scale_fill_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC","2SSB1"="Ssb","2SSB2"="Ssb Zmax"), name = "", guide = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
  scale_color_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC","2SSB1"="Ssb","2SSB2"="Ssb Zmax"), name = "") +
  scale_x_continuous(limits = c(0,2), breaks = c(0, 1), labels = c("0" = "Domain start", "1" = "Domain end"))
plot <- plot + theme_classic(24) + labs(x = "Normalized domain position", y = "Density") +
  theme(legend.position = c(.65,.95), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text = element_text(color = "black", size = 20), legend.text = element_text(size = 20), legend.background = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS5/F.EarlyDomainsCATH_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(domains_cath_tric5_start[domainend_norm < 0.75]$length.ratio, domains_cath_ssb5_start[domainend_norm < 0.75]$length.ratio)

temp <- tric_peaks_all[orf == "YGR145W", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 600), color = "gray70", size = 2) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 175), color = "#FB9A99", alpha = 0.7, size = 8) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 30), color = "#E31A1C", alpha = 0.7, size = 8) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
        axis.text = element_blank(), plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS5/E.DomainSchematic.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Secondary structure
plot <- ggplot(data = psipred_ssb_dt[, .(adjusted_start, ssb = movingAverage(sheet, n=25, center=T))]) + 
  stat_summary(data = psipred_tric_dt[, .(adjusted_start, tric = movingAverage(sheet, n=25, center=T))], aes(adjusted_start, tric), 
               fun.data = "mean_cl_boot", geom = "ribbon", size = 1.25, alpha = 0.3, fill = '#1F78B4', fun.args=list(conf.int=0.5)) + 
  stat_summary(data = psipred_bukau_dt[, .(adjusted_start, bukau = movingAverage(sheet, n=25, center=T))], aes(adjusted_start, bukau), 
               fun.data = "mean_cl_boot", geom = "ribbon", size = 1.25, alpha = 0.3, fill = '#6A3D9A', fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(adjusted_start, ssb), fun.data = "mean_cl_boot", geom = "ribbon", size = 1.25, alpha = 0.3, fill = '#33A02C', fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(adjusted_start, ssb, color = '2SSB1'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(data = psipred_tric_dt[, .(adjusted_start, tric = movingAverage(sheet, n=25, center=T))], aes(adjusted_start, tric, color = '1TRiC1'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(data = psipred_bukau_dt[, .(adjusted_start, bukau = movingAverage(sheet, n=25, center=T))], aes(adjusted_start, bukau, color = 'dP'), fun.y = "mean", geom = "line", size = 1.25) +
  scale_color_manual(limits = c("1TRiC1", "2SSB1", "dP"), labels = c("TRiC", "Ssb", "Ssb Doring"), values = allcolors, name = "") +
  coord_cartesian(ylim = c(.19,.24), xlim = c(-200, -20)) +
  scale_y_continuous(breaks = c(0.20,0.22,0.24))
G <- plot + theme_classic(20) + labs(y = "Beta sheet propensity", x = "Codons from\nchaperone association") +
  theme(legend.position = c(.01,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS5/G.SecondaryStructure_fromstart.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


### Composition
composition_log2$identity <- factor(composition_log2$region,levels = c("TRiC_nascent", "Ssb_nascent", "Bukau_nascent"))
plot <- ggplot(composition_log2[region == "TRiC_nascent" | region == "Ssb_nascent" | region == "Bukau_nascent"], aes(x="fA", y=A, fill = identity)) + 
  geom_hline(yintercept = 0, color = "black") +
  geom_col(position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="rC", y=C, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="aD", y=D, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="bE", y=E, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="kF", y=F, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="sG", y=G, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="eH", y=H, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="gI", y=I, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="cK", y=K, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="hL", y=L, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="iM", y=M, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="nN", y=N, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="tP", y=P, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="oQ", y=Q, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="dR", y=R, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="pS", y=S, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="qT", y=T, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="jV", y=V, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="lW", y=W, fill = identity), position = "dodge", size = 0.3, color = "black") +
  geom_col(aes(x="mY", y=Y, fill = identity), position = "dodge", size = 0.3, color = "black") +
  scale_fill_manual(values = c("#1F78B4", "#33A02C","#6A3D9A"),limits = c("TRiC_nascent", "Ssb_nascent", "Bukau_nascent"), labels = c("TRiC","Ssb","Ssb, Doring"), name = "",
                    guide = guide_legend(keywidth = 1, keyheight = 1))
plot <- plot + theme_classic(20) + labs(y = "log2 fold change", x = "") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(size = 14),
        legend.position = c(.5,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), axis.text = element_text(color = "black", size = 14))
ggsave("/Users/KevinStein/Desktop/Figures/FigS5/H.Nascent.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
