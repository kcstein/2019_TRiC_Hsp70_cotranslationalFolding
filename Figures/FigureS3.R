### Significance of positional recruitment
# All substrates
plot <- ggplot(tric_peaks5, aes("4.TRiC", position_norm)) + 
  geom_hline(yintercept = 0.5, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_boxplot(fill="#1F78B4", notch = T) +
  geom_boxplot(data = ssb_peaks5, aes("3.SSB", position_norm), fill = "#33A02C", notch = T) +
  geom_boxplot(data = tric_peaks5_max, aes("2.TRiC", position_norm), fill = "#A6CEE3", notch = T) +
  geom_boxplot(data = ssb_peaks5_max, aes("1.SSB", position_norm), fill = "#B2DF8A", notch = T) +
  scale_x_discrete(labels = c("4.TRiC" = "All TRiC","3.SSB" = "All Ssb","2.TRiC" = "Max TRiC","1.SSB" = "Max Ssb")) + coord_flip() +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(27) + labs(y = "Norm. codon position", x = "") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text = element_text(color = "black", size = 27))
ggsave("/Users/KevinStein/Desktop/Figures/FigS3/A.PositionalRecruitment_Allsubstrates_p4.918926e-24,9.190131e-19.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
temp <- wilcox.test(tric_peaks5$position_norm, ssb_peaks5$position_norm, alternative = 't')
temp$p.value
temp <- wilcox.test(tric_peaks5_max$position_norm, ssb_peaks5_max$position_norm, alternative = 't')
temp$p.value


# Significance of positional enrichment
plot <- ggplot(tric_peaks5[position_norm <= (1/3)], aes("1.tric", tric_odds_ma)) + geom_boxplot(fill = "#C6DBEF", notch = T) + 
  geom_boxplot(data = tric_peaks5[position_norm > (1/3) & position_norm <= (2/3)], aes("2.tric", tric_odds_ma), fill = "#6BAED6", notch = T) +
  geom_boxplot(data = tric_peaks5[position_norm > (2/3)], aes("3.tric", tric_odds_ma), fill = "#2171B5", notch = T) +
  scale_x_discrete(labels = c("N-terminus", "Middle", "C-terminus")) +
  coord_cartesian(ylim = c(0,10))
plot <- plot + theme_classic(20) + labs(y = "TRiC enrichment (odds ratio)", x = "") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS3/B.PositionalEnrichment_TRiC_p2.668e-06,0.007724,1.14e-10.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(tric_peaks5[position_norm > (1/3) & position_norm <= (2/3)]$tric_odds_ma,
            tric_peaks5[position_norm > (2/3)]$tric_odds_ma)
wilcox.test(tric_peaks5[position_norm <= (1/3)]$tric_odds_ma,
            tric_peaks5[position_norm > (2/3)]$tric_odds_ma)

plot <- ggplot(ssb_peaks5[position_norm <= (1/3)], aes("1.ssb", ssb_odds_ma)) + geom_boxplot(fill = "#C7E9C0", notch = T) + 
  geom_boxplot(data = ssb_peaks5[position_norm > (1/3) & position_norm <= (2/3)], aes("2.ssb", ssb_odds_ma), fill = "#74C476", notch = T) +
  geom_boxplot(data = ssb_peaks5[position_norm > (2/3)], aes("3.ssb", ssb_odds_ma), fill = "#238B45", notch = T) +
  scale_x_discrete(labels = c("N-terminus", "Middle", "C-terminus")) +
  coord_cartesian(ylim = c(0,10))
plot <- plot + theme_classic(20) + labs(y = "Ssb enrichment (odds ratio)", x = "") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS3/C.PositionalEnrichment_SSB_p0.0001101,NS,0.01169.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(ssb_peaks5[position_norm > (1/3) & position_norm <= (2/3)]$ssb_odds_ma,
            ssb_peaks5[position_norm > (2/3)]$ssb_odds_ma)
wilcox.test(ssb_peaks5[position_norm <= (1/3)]$ssb_odds_ma,
            ssb_peaks5[position_norm > (2/3)]$ssb_odds_ma)


# Heat map of shared substrates max
plot <- ggplot(ssb_peaks5_max[ssb_peaks5_max$orf %in% tric_peaks5_max$orf], aes(x = reorder(orf, -length), y = length)) + geom_col(fill = "gray90", color = "gray90") +
  geom_point(data = ssb_peaks5_max[ssb_peaks5_max$orf %in% tric_peaks5_max$orf], aes(x = reorder(orf, -length), y = peak, color = '2SSB1'), alpha = 0.6, size = 2) +
  geom_point(data = tric_peaks5_max[tric_peaks5_max$orf %in% ssb_peaks5_max$orf], aes(x = reorder(orf, -length), y = peak, color = '1TRiC1'), size = 2) + 
  coord_flip(ylim = c(0,1000)) +
  scale_color_manual(values = allcolors, labels = c("1TRiC1"="TRiC association","2SSB1"="Ssb association"), name = "") +
  scale_y_continuous(breaks = c(0, 500, 1000), labels = c(0, 500, 1000))
plot <- plot + theme_classic(20) + labs(y = "Codon position", x = "Gene") +
  theme(legend.position = c(.35,1.1), legend.background = element_blank(), legend.justification = c("left", "top"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 1, color = "black"), axis.text.x = element_text(color = "black"), legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/FigS3/D.SharedSubstrates_maxBindingSites.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
plot <- ggplot(ssb_peaks5_max[ssb_peaks5_max$orf %in% tric_peaks5_max$orf], aes(position_norm, fill = "2SSB1", color = "2SSB1")) + geom_density(size = 1, alpha = 0.3) +
  geom_density(data = tric_peaks5_max[tric_peaks5_max$orf %in% ssb_peaks5_max$orf], aes(position_norm, fill = '1TRiC1', color = "1TRiC1"), size = 1, alpha = 0.3) + 
  scale_fill_manual(values = allcolors, labels = c("TRiC association","SSB association"), name = "", guide = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
  scale_color_manual(values = allcolors, labels = c("TRiC association","SSB association"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(30) + labs(x = "Norm. codon position", y = "Density") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), panel.background = element_blank(), axis.text.x = element_text(color = "black", size = 30),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS3/D.SharedSubstrates_maxBindingSites_density.pdf", plot, width = 5, height = 4, dpi = 300, useDingbats = F)
wilcox.test(tric_peaks5_max[tric_peaks5_max$orf %in% ssb_peaks5_max$orf]$position_norm, ssb_peaks5_max[ssb_peaks5_max$orf %in% tric_peaks5_max$orf]$position_norm)

# Response letter
temp <- data.table(tric = tric_peaks5_max[tric_peaks5_max$orf %in% ssb_peaks5_max$orf]$peak_start,
                   tric_norm = tric_peaks5_max[tric_peaks5_max$orf %in% ssb_peaks5_max$orf]$position_norm,
                   ssb = ssb_peaks5_max[ssb_peaks5_max$orf %in% tric_peaks5_max$orf]$peak_start,
                   ssb_norm = ssb_peaks5_max[ssb_peaks5_max$orf %in% tric_peaks5_max$orf]$position_norm)
temp[, diff := tric - ssb]
temp[, diff_norm := tric_norm - ssb_norm]
length(temp$diff)
length(temp[diff > -1]$diff)
ggplot(temp, aes(diff))+geom_histogram()


### Comparison to Bukau
temp1 <- tric_substrates5[tric_substrates5$orf %in% ssb_substrates5$orf]
temp2 <- temp1[temp1$orf %in% bukau_substrates5$orf]
plot <- ggplot(tric_peaks5[(tric_peaks5$orf %in% temp2$orf)], aes("3.TRiC", position_norm)) + 
  geom_hline(yintercept = 0.5, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_boxplot(fill="#1F78B4", notch = T) +
  geom_boxplot(data = ssb_peaks5[(ssb_peaks5$orf %in% temp2$orf)], aes("2.SSB", position_norm), fill = "#33A02C", notch = T) +
  geom_boxplot(data = bukau_peaks5[(bukau_peaks5$orf %in% temp2$orf)], aes("1.Doring", position_norm), fill = "#6A3D9A", notch = T) +
  scale_x_discrete(labels = c("Doring", "Ssb", "TRiC")) + coord_flip() +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(27) + labs(y = "Norm. codon position", x = "") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text = element_text(color = "black", size = 27))
ggsave("/Users/KevinStein/Desktop/Figures/FigS3/E.PositionalRecruitment_Sharedsubstrates_p2.28e-07,5.997e-06,NS.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
wilcox.test(tric_peaks5[(tric_peaks5$orf %in% temp2$orf)]$position_norm, ssb_peaks5[(ssb_peaks5$orf %in% temp2$orf)]$position_norm)
wilcox.test(tric_peaks5[(tric_peaks5$orf %in% temp2$orf)]$position_norm, bukau_peaks5[(bukau_peaks5$orf %in% temp2$orf)]$position_norm)
wilcox.test(ssb_peaks5[(ssb_peaks5$orf %in% temp2$orf)]$position_norm, bukau_peaks5[(bukau_peaks5$orf %in% temp2$orf)]$position_norm)


### Pulse-chase
temp <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/pulsechase.csv", stringsAsFactors = T, header = T)
temp <- as.data.table(temp)
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
temp1 <- summarySE(temp, measurevar="ratio", groupvars=c("sample"))
plot <- ggplot(temp1, aes(x=sample, y=ratio, fill=sample)) + geom_bar(stat="identity", size = 0.5, color = "black") +
  geom_point(data = temp, aes(sample,ratio), color = "black", size = 1.5, position = "jitter") +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), width=.2) + coord_cartesian(ylim = c(0,0.35)) +
  scale_fill_manual(limits = c("1Control","2WT","3deltaSSB"), labels = c("Control","WT","deltaSSB"), values = c("gray80", "gray50", "#CB181D"), name = "") +
  scale_x_discrete(labels = c("", "",""))
plot <- plot + theme_classic(20) + labs(y = "", x = "") + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS3/F.pulsechase.pdf", plot, width = 2, height = 4, dpi = 300, useDingbats = F)
t.test(temp[sample == "2WT"]$ratio, temp[sample == "3deltaSSB"]$ratio, alternative = 'l')


### Different protein length of TRiC substrates
plot <- ggplot(tric_peaks5[length <= 600], aes(position_norm, fill = '1TRiC2', color = "1TRiC2")) + 
  geom_vline(xintercept = 0.5, color = 'gray50', linetype = 'dashed', size = 1.25) +
  geom_density(size = 1.25, alpha = 0.3) +
  geom_density(data = tric_peaks5[length > 600], aes(position_norm, fill = '1TRiC1', color = "1TRiC1"), size = 1.25, alpha = 0.3) +
  scale_fill_manual(values = allcolors, labels = c("> 600 aa","< 600 aa"), name = "", guide = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
  scale_color_manual(values = allcolors, labels = c("> 600 aa","< 600 aa"), name = "") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 1), labels = c("0" = "Start", "1" = "Stop"))
plot <- plot + theme_classic(20) + labs(x = "Norm. codon position", y = "Density") +
  theme(legend.position = c(.05,1.1), legend.background = element_blank(), legend.justification = c("left", "top"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        plot.background = element_blank(), axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS3/G.LongVSshort_TRiCsubstrates_AllBindingSites_density.pdf", plot, width = 5, height = 4, dpi = 300, useDingbats = F)
wilcox.test(tric_peaks5[length <= 600]$position_norm, tric_peaks5[length > 600]$position_norm)
