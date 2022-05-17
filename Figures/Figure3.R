### Binding in reference to domain
# Heat maps and density plots
temp <- tric_peaks_all[orf == "YGR145W", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 600), color = "gray70", size = 2) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 0, yend = 190), color = "#CAB2D6", alpha = 0.7, size = 8) +
  geom_segment(data = tric_peaks_all[orf == "YGR145W"], aes(x = orf, xend = orf, y = 30, yend = 35), color = "gray50", size = 10) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
        axis.text = element_blank(), plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig3/B.DomainSchematic.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(domains_cath_ssb5_start, aes(length.ratio, fill = '2SSB1', color = "2SSB1")) + 
  geom_vline(xintercept = 1, color = 'gray70', linetype = 'dashed', size = 1.25) +
  geom_density(size = 1.25, alpha = 0.3) +
  geom_density(data = domains_cath_tric5_start, aes(length.ratio, fill = "1TRiC1", color = "1TRiC1"), size = 1.25, alpha = 0.3) + 
  scale_fill_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC","2SSB1"="Ssb"), name = "", guide = guide_legend(keywidth = 1, keyheight = 1)) +
  scale_color_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC","2SSB1"="Ssb"), name = "") +
  scale_x_continuous(limits = c(0,2), breaks = c(0, 1), labels = c("0" = "Domain start", "1" = "Domain end"))
plot <- plot + theme_classic(24) + labs(x = "Normalized domain position", y = "Density") +
  theme(legend.position = c(.75,1.1), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text = element_text(color = "black", size = 20), legend.text = element_text(size = 20), legend.background = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig3/B.DomainsCATH_AllBindingSites_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(domains_cath_tric5_start$length.ratio, domains_cath_ssb5_start$length.ratio, alternative = 't')

temp <- domains_cath_ssb5_start
setkeyv(temp, c("domain_alias"))
temp1 <- temp[, .SD[which.min(Start)], by = domain_alias]
temp1[, new.length := length - Start + 1]
plot <- ggplot(temp1, aes(x = reorder(domain_alias, -DomainLength), y = new.length)) + geom_col(fill = "gray90") +
  geom_col(data = temp1, aes(x = reorder(domain_alias, -DomainLength), y = DomainLength), fill = "#CAB2D6", alpha = 0.6, show.legend = F) +
  geom_hline(yintercept = 30, color = 'gray50', linetype = 'dashed', size = 1.5) +
  geom_point(data = domains_cath_ssb5_start, aes(x = reorder(domain_alias, -DomainLength), y = peak.start, color = "2SSB1"), alpha = 0.6, size = 1.5) + 
  coord_flip(ylim = c(0,600)) + 
  scale_color_manual(values = allcolors, labels = c("2SSB1"="Ssb association"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Distance from domain start (codons)", x = "Domain") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/Fig3/C.SSB_AllDomains_AllBindingSites_CATH.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

temp <- domains_cath_tric5_start
setkeyv(temp, c("domain_alias"))
temp1 <- temp[, .SD[which.min(Start)], by = domain_alias]
temp1[, new.length := length - Start + 1]
plot <- ggplot(temp1, aes(x = reorder(domain_alias, -DomainLength), y = new.length)) + geom_col(fill = "gray90") +
  geom_col(data = temp1, aes(x = reorder(domain_alias, -DomainLength), y = DomainLength), fill = "#CAB2D6", alpha = 0.7, show.legend = F) +
  geom_hline(yintercept = 30, color = 'gray50', linetype = 'dashed', size = 1.5) +
  geom_point(data = domains_cath_tric5_start, aes(x = reorder(domain_alias, -DomainLength), y = peak.start, color = "1TRiC1"), size = 1.5) + 
  coord_flip(ylim = c(0,600)) + 
  scale_color_manual(values = allcolors, labels = c("1TRiC1"="TRiC association"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Distance from domain start (codons)", x = "Domain") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/Fig3/D.TRiC_AllDomains_AllBindingSites_CATH.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Secondary structure
secondarystructure <- c("Sheet"="red", "Helix"="cyan", "Coil"="magenta")
plot <- ggplot(data = psipred_ssb_dt[, .(adjusted_start, ssb = movingAverage(sheet, n=25, center=T))]) + 
  stat_summary(data = psipred_ssb_dt[, .(adjusted_random, ssbR = movingAverage(sheet, n=25, center=T))], aes(adjusted_random, ssbR), fun.data = "mean_cl_boot", geom = "ribbon", size = 1.25, alpha = 0.3, fill = '#B2DF8A', fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(adjusted_start, ssb), fun.data = "mean_cl_boot", geom = "ribbon", size = 1.25, alpha = 0.3, fill = '#33A02C', fun.args=list(conf.int=0.5)) + 
  stat_summary(data = psipred_ssb_dt[, .(adjusted_random, ssbR = movingAverage(sheet, n=25, center=T))], aes(adjusted_random, ssbR, color = '2SSB2'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(aes(adjusted_start, ssb, color = '2SSB1'), fun.y = "mean", geom = "line", size = 1.25) +
  scale_color_manual(limits = c("2SSB1","2SSB2"), labels = c("Ssb", "Random sampling"), values = allcolors, name = "") +
  coord_cartesian(ylim = c(.19,.24), xlim = c(-200, -20)) +
  scale_y_continuous(breaks = c(0.20,0.22,0.24))
G <- plot + theme_classic(20) + labs(y = "Beta sheet propensity", x = "Codons from Ssb association") +
  theme(legend.position = c(.01,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig3/E.SecondaryStructure_Ssb_juststart.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = psipred_tric_dt[, .(adjusted_start, tric = movingAverage(sheet, n=25, center=T))]) + 
  stat_summary(data = psipred_tric_dt[, .(adjusted_random, tricR = movingAverage(sheet, n=25, center=T))], aes(adjusted_random, tricR), 
               fun.data = "mean_cl_boot", geom = "ribbon", size = 1.25, alpha = 0.3, fill = '#A6CEE3', fun.args=list(conf.int=0.5)) + 
  stat_summary(aes(adjusted_start, tric), fun.data = "mean_cl_boot", geom = "ribbon", size = 1.25, alpha = 0.3, fill = '#1F78B4', fun.args=list(conf.int=0.5)) + 
  stat_summary(data = psipred_tric_dt[, .(adjusted_random, tricR = movingAverage(sheet, n=25, center=T))], aes(adjusted_random, tricR, color = '1TRiC2'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(aes(adjusted_start, tric, color = '1TRiC1'), fun.y = "mean", geom = "line", size = 1.25) +
  scale_color_manual(limits = c("1TRiC1", "1TRiC2"), labels = c("TRiC", "Random sampling"), values = allcolors, name = "") +
  coord_cartesian(ylim = c(.19,.24), xlim = c(-200, -20)) +
  scale_y_continuous(breaks = c(0.20,0.22,0.24))
G <- plot + theme_classic(20) + labs(y = "Beta sheet propensity", x = "Codons from TRiC association") +
  theme(legend.position = c(.01,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig3/F.SecondaryStructure_TRiC_juststart.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


### Charge
plot <- ggplot(data = hydro_tric) + xlim(1,40) +
  stat_summary(aes(position, charge), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(position, charge), fun.y = "mean", geom = "line", size = 1, color = '#1F78B4') +
  stat_summary(data = hydro_ssb, aes(position, charge), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#33A02C') +
  stat_summary(data = hydro_ssb, aes(position, charge), fun.y = "mean", geom = "line", size = 1, color = '#33A02C') +
  #stat_summary(data = hydro_bukau, aes(position, charge), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
  #            fun.args=list(conf.int=0.5), fill = '#6A3D9A') +
  #stat_summary(data = hydro_bukau, aes(position, charge), fun.y = "mean", geom = "line", size = 1, color = '#6A3D9A') +
  stat_summary(data = hydro_random, aes(position, charge), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = hydro_random, aes(position, charge), fun.y = "mean", geom = "line", size = 1, color = 'gray50')
G <- plot + theme_classic(20) + labs(y = "Charge", x = "Codons from chaperone association") +
  theme(legend.position = c(.01,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/Fig3/G.Charge_fromstart.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)


### Significance of particular CATH topology
temp <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/domains_freq.csv", stringsAsFactors = T, header = T)
temp <- as.data.table(temp)
temp$Description <- factor(temp$Description,levels = c("Nucleotide-binding alpha-beta plait", "Rossmann fold", "TIM barrel", "WD40 repeat beta-propeller"))
plot <- ggplot(temp[Chaperone == "1TRiC" | Chaperone == "2SSB"], aes(x = Description,y=-log10(pvalue), fill = factor(Chaperone))) + geom_col(position = "dodge", color = "black", size = 0.2) + scale_y_continuous(limits = c(0,4)) +
  geom_hline(yintercept = -log10(0.05), color = '#CAB2D6', linetype = 'dashed', size = 1) +
  scale_fill_manual(limits = c("1TRiC","2SSB"), labels = c("TRiC", "Ssb"), values = c("#1F78B4", "#33A02C"), name = "", guide = guide_legend(keywidth = 1, keyheight = 1)) +
  scale_x_discrete(labels = c("Nucleotide-binding\nalpha-beta plait", "Rossmann fold", "TIM barrel", "WD40 repeat\nbeta-propeller"))
plot <- plot + theme_classic(18) + labs(y = "-log10 (p-value)", x = "") +
  theme(legend.position = c(.05,1.2), legend.background = element_blank(), legend.justification = c("left", "top"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 15), axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 15, color = "black"),
        plot.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig3/H.Domains.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

