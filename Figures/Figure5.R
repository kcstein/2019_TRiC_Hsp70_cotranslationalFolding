### TIM barrels
plot <- ggplot(domains_cath_tric5_start[Gene3D_topo == "3.40.50"], aes(length.ratio, fill = "dP", color = "dP")) + 
  geom_vline(xintercept = 1, color = 'gray70', linetype = 'dashed', size = 1.25) +
  geom_density(size = 1, alpha = 0.3) +
  geom_density(data = domains_cath_tric5_start[Gene3D_topo == "2.130.10"], aes(length.ratio, fill = "dO", color = "dO"), size = 1, alpha = 0.3) +
  geom_density(data = domains_cath_tric5_start[Gene3D_topo == "3.20.20"], aes(length.ratio, fill = "dG", color = "dG"), size = 1, alpha = 0.3) +
  scale_fill_manual(values = c(allcolors), labels = c("dO"="WD","dP"="Ross","dG"="TIM"), name = "", guide = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
  scale_color_manual(values = c(allcolors), labels = c("dO"="WD","dP"="Ross","dG"="TIM"), name = "") +
  scale_x_continuous(limits = c(0,2), breaks = c(0, 1), labels = c("0" = "Domain start", "1" = "Domain end"))
plot <- plot + theme_classic(24) + labs(x = "Normalized domain position", y = "Density") +
  theme(legend.position = c(.65,.95), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text = element_text(color = "black", size = 20), legend.text = element_text(size = 20), legend.background = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig5/A.WD_Ross_TIM_CATH_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(domains_cath_tric5_start[Gene3D_topo == "3.20.20"]$length.ratio, domains_cath_tric5_start[Gene3D_topo == "3.40.50"]$length.ratio)
wilcox.test(domains_cath_tric5_start[Gene3D_topo == "3.20.20"]$length.ratio, domains_cath_tric5_start[Gene3D_topo == "2.130.10"]$length.ratio)

temp <- domains_cath_tric5_start[Gene3D_topo == "3.20.20"]
setkeyv(temp, c("domain_alias"))
temp1 <- temp[, .SD[which.min(Start)], by = domain_alias]
temp1[, new.length := length - Start + 1]
plot <- ggplot(temp1, aes(x = reorder(domain_alias, -DomainLength), y = new.length)) + geom_col(fill = "gray90") +
  geom_col(data = temp1, aes(x = reorder(domain_alias, -DomainLength), y = DomainLength), fill = "gray70", alpha = 0.6, show.legend = F) +
  geom_point(data = domains_cath_tric5_start[domains_cath_tric5_start$domain_alias %in% temp1$domain_alias], aes(x = reorder(domain_alias, -DomainLength), y = peak.start, color = "1TRiC1"), size = 2) + 
  coord_flip(ylim = c(0,600)) + 
  scale_color_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC association"), name = "")
plot <- plot + theme_classic(20) + labs(y = "Distance from domain start (codons)", x = "Domain") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/Fig5/B.TRiC_TIM_AllBindingSites_CATH.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(domains_cath_ssb5_start[Gene3D_topo == "3.20.20"], aes(length.ratio, fill = '2SSB1', color = "2SSB1")) + 
  geom_vline(xintercept = 1, color = 'gray70', linetype = 'dashed', size = 1.25) +
  geom_density(size = 1, alpha = 0.3) +
  geom_density(data = domains_cath_tric5_start[Gene3D_topo == "3.20.20"], aes(length.ratio, fill = "1TRiC1", color = "1TRiC1"), size = 1, alpha = 0.3) + 
  scale_fill_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC","2SSB1"="Ssb"), name = "", guide = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
  scale_color_manual(values = c(allcolors), labels = c("1TRiC1"="TRiC","2SSB1"="Ssb"), name = "") +
  scale_x_continuous(limits = c(0,2), breaks = c(0, 1), labels = c("0" = "Domain start", "1" = "Domain end"))
plot <- plot + theme_classic(24) + labs(x = "Normalized domain position", y = "Density") +
  theme(legend.position = c(.65,.95), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text = element_text(color = "black", size = 20), legend.text = element_text(size = 20), legend.background = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig5/B.TIM_CATH_density.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(domains_cath_tric5_start[Gene3D_topo == "3.20.20"]$length.ratio, domains_cath_ssb5_start[Gene3D_topo == "3.20.20"]$length.ratio)


### ATP2
# qPCR
temp <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/qPCR.csv", stringsAsFactors = T, header = T)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
temp1 <- summarySE(temp, measurevar="log2", groupvars=c("sample"))
plot <- ggplot(temp1, aes(x=sample, y=log2, fill=sample)) + geom_bar(stat="identity", size = 0.5, color = "black") +
  geom_point(data = temp, aes(sample,log2), color = "black", size = 1.5) +
  geom_errorbar(aes(ymin=log2-se, ymax=log2+se), width=.2) + ylim(0,3) +
  scale_fill_manual(limits = c("ATP2","atp2-PA"), labels = c("ATP2","atp2-PA"), values = c("gray50", "#CB181D"), name = "") +
  scale_x_discrete(labels = c("", ""))
plot <- plot + theme_classic(20) + labs(y = "", x = "") + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/Fig5/G.Atp2_qPCR.pdf", plot, width = 2, height = 4, dpi = 300, useDingbats = F)

# ribosome profiling
plot <- ggplot(aes(position, occupancy), data=tric_substrates5_dt[orf == "YJR121W" & R2X > 7 & tric_odds_ma < Inf, .(position, occupancy = movingAverage(tric_odds, n=5, center=T))]) + 
  geom_line(data = atp2_fishers[orf == "YJR121W" & Atp2_R1 > 7 & atp2_odds_ma < Inf, .(position, occupancy2 = movingAverage(atp2_odds, n=5, center=T))], aes(position, occupancy2, color = "red"), size = 1.25) + 
  geom_line(aes(color = "gray40"), size = 1.25) + 
  scale_color_manual(labels = c("ATP2", "atp2-P353,355A"), values = c("black", "red"), name = "")
G <- plot + theme_classic(20) + labs(y = "TRiC enrichment (odds ratio)", x = "Codon position") +
  theme(legend.position = c(.01,.99), legend.background = element_blank(), legend.justification = c("left", "top"), legend.text = element_text(size=20, face = "italic"),
        axis.text = element_text(size = 20, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("ATP2")
ggsave("/Users/KevinStein/Desktop/Figures/Fig5/H.ATP2.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

# ATP2 domain schematic
temp <- tric_dt[orf == "YJR121W", .SD[which.min(position)], by = orf]
plot <- ggplot(temp, aes(x = orf, y = length)) +
  geom_segment(data = temp[orf == "YJR121W"], aes(x = orf, xend = orf, y = 0, yend = length), color = "gray70", size = 2) +
  geom_segment(data = domains_cath[orf == "YJR121W" & Domain == "G3DSA:2.40.10.170"], aes(x = orf, xend = orf, y = Start, yend = End), color = "#41B6C4", size = 2) +
  geom_segment(data = domains_cath[orf == "YJR121W" & Domain == "G3DSA:3.40.50.300"], aes(x = orf, xend = orf, y = Start, yend = End), color = "#CAB2D6", size = 2) +
  geom_segment(data = domains_cath[orf == "YJR121W" & Domain == "G3DSA:1.10.1140.10"], aes(x = orf, xend = orf, y = Start, yend = End), color = "#FDBF6F", size = 2) +
  coord_flip() 
plot <- plot + theme_classic(20) + labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), plot.background = element_blank(), panel.background = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/Fig5/H.ATP2_domains.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
