### Replicates
plot <- ggplot(expression_tric[R1X_tpm != 0 & R2X_tpm != 0], aes(R1X_tpm, R2X_tpm)) + geom_point(color = 'black', fill = "black", alpha = 0.5, size = 2) + 
  scale_x_log10(limits = c(1e-2,1.5e5), breaks = c(1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-2", "-1", "0","1", "2", "3","4","5")) + 
  scale_y_log10(limits = c(1e-2,1.5e5), breaks = c(1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-2", "-1", "0","1", "2", "3","4","5")) + annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
plot <- plot + theme_classic(20) + labs(y = "Translatome replicate 2\n(log10 TPM)", x = "Translatome replicate 1\n(log10 TPM)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
cor(expression_tric[R1X_tpm != 0 & R2X_tpm != 0]$R1X_tpm, expression_tric[R1X_tpm != 0 & R2X_tpm != 0]$R2X_tpm, method = "pearson")
ggsave("/Users/KevinStein/Desktop/Figures/FigS1/A.TRiC_RiboReplicates_Pearson0.99.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(expression_tric[T1X_tpm != 0 & T2X_tpm != 0], aes(T1X_tpm, T2X_tpm)) + geom_point(color = 'black', fill = "black", alpha = 0.5, size = 2) + 
  scale_x_log10(limits = c(1e-2,1.5e5), breaks = c(1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-2", "-1", "0","1", "2", "3","4","5")) + 
  scale_y_log10(limits = c(1e-2,1.5e5), breaks = c(1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-2", "-1", "0","1", "2", "3","4","5")) + annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
plot <- plot + theme_classic(20) + labs(y = "TRiC-bound replicate 2\n(log10 TPM)", x = "TRiC-bound replicate 1\n(log10 TPM)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
cor(expression_tric[T1X_tpm != 0 & T2X_tpm != 0]$T1X_tpm, expression_tric[T1X_tpm != 0 & T2X_tpm != 0]$T2X_tpm, method = "pearson")
ggsave("/Users/KevinStein/Desktop/Figures/FigS1/A.TRiC_TRiCReplicates_Pearson0.99.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(expression_tric[ssb_Rchx1_tpm != 0 & ssb_Rchx2_tpm != 0], aes(ssb_Rchx1_tpm, ssb_Rchx2_tpm)) + geom_point(color = 'black', fill = "black", alpha = 0.5, size = 2) + 
  scale_x_log10(limits = c(1e-2,1.5e5), breaks = c(1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-2", "-1", "0","1", "2", "3","4","5")) + 
  scale_y_log10(limits = c(1e-2,1.5e5), breaks = c(1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-2", "-1", "0","1", "2", "3","4","5")) + annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
plot <- plot + theme_classic(20) + labs(y = "Translatome replicate 2\n(log10 TPM)", x = "Translatome replicate 1\n(log10 TPM)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
cor(expression_tric[ssb_Rchx1_tpm != 0 & ssb_Rchx2_tpm != 0]$ssb_Rchx1_tpm, expression_tric[ssb_Rchx1_tpm != 0 & ssb_Rchx2_tpm != 0]$ssb_Rchx2_tpm, method = "pearson")
ggsave("/Users/KevinStein/Desktop/Figures/FigS1/B.SSB_RiboReplicates_Pearson0.99.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(expression_tric[ssb_Schx1_tpm != 0 & ssb_Schx2_tpm != 0], aes(ssb_Schx1_tpm, ssb_Schx2_tpm)) + geom_point(color = 'black', fill = "black", alpha = 0.5, size = 2) + 
  scale_x_log10(limits = c(1e-2,1.5e5), breaks = c(1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-2", "-1", "0","1", "2", "3","4","5")) + 
  scale_y_log10(limits = c(1e-2,1.5e5), breaks = c(1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-2", "-1", "0","1", "2", "3","4","5")) + annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
plot <- plot + theme_classic(20) + labs(y = "Ssb-bound replicate 2\n(log10 TPM)", x = "Ssb-bound replicate 1\n(log10 TPM)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
cor(expression_tric[ssb_Schx1_tpm != 0 & ssb_Schx2_tpm != 0]$ssb_Schx1_tpm, expression_tric[ssb_Schx1_tpm != 0 & ssb_Schx2_tpm != 0]$ssb_Schx2_tpm, method = "pearson")
ggsave("/Users/KevinStein/Desktop/Figures/FigS1/B.SSB_SSBReplicates_Pearson0.94.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### TRiC subunit enrichment
plot <- ggplot(expression_tric[orf == "YDR188W" | orf  == "YDR212W" | orf == "YIL142W" | orf == "YDL143W" | orf == "YJL111W" | orf == "YDL229W" | orf == "YNL209W"], 
               aes(name1, ratio1)) + geom_point(size = 3, color = "black", alpha = 0.6, fill = "gray40", shape = 21) +
  geom_point(data = expression_tric[orf == "YDR188W" | orf  == "YDR212W" | orf == "YIL142W" | orf == "YDL143W" | orf == "YJL111W" | orf == "YDL229W" | orf == "YNL209W"],
             aes(name1, ratio2), size = 3, color = "black", alpha = 0.6, fill = "gray40", shape = 21) +
  geom_point(data = expression_tric[orf == "YJL014W" | orf == "YJR064W" | orf == "YJL008C" | orf == "YML124C" | orf == "YML085C"],
             aes(name1, ratio1), size = 3, color = "#1F78B4", alpha = 0.9, fill = "#1F78B4") +
  geom_point(data = expression_tric[orf == "YJL014W" | orf == "YJR064W" | orf == "YJL008C" | orf == "YML124C" | orf == "YML085C"],
             aes(name1, ratio2), size = 3, color = "#1F78B4", alpha = 0.9, fill = "#1F78B4") +
  scale_y_log10(limits = c(0.8,25), breaks = c(1,5,10,20)) + annotation_logticks(sides = "l", short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
plot <- plot + theme_classic(20) + labs(y = "Enrichment\n(TRiC tpm / Total tpm)", x = "") +
  theme(axis.text.y = element_text(size = 20, color = "black"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(),
        axis.text.x = element_text(angle=45, vjust=0.5, hjust=1, size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS1/C.TRiCenrichment.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Puromycin
plot <- ggplot(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & group == 1]) +
  geom_point(aes(Rpur_tpm, Tpur_tpm), color = "#E31A1C", fill = "#E31A1C", alpha = 0.8, size = 3) + 
  geom_point(aes(ribo_tpm, tric_tpm), color = "#1F78B4", fill = "#1F78B4", alpha = 0.8, size = 3) + 
  geom_point(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & tricsubunit == 1], aes(Rpur_tpm, Tpur_tpm), color="#E31A1C", fill = "#E31A1C", alpha = 0.8, size = 3) +
  geom_point(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & tricsubunit == 1], aes(ribo_tpm, tric_tpm), color="#1F78B4", fill = "#1F78B4", alpha = 0.8, size = 3) +
  geom_text(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & tricsubunit == 1], aes(Rpur_tpm, Tpur_tpm,label=tricname), hjust=1.1,vjust=-.1, size=3) +
  geom_text(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & tricsubunit == 1], aes(ribo_tpm, tric_tpm,label=tricname), hjust=1.1,vjust=-.1, size=3) +
  scale_x_log10(limits = c(1e0,1e3), breaks = c(1e0,1e1,1e2,1e3), labels = c("0","1", "2", "3")) + 
  scale_y_log10(limits = c(1e0,4e3), breaks = c(1e0,1e1,1e2,1e3), labels = c("0","1", "2", "3")) + annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
plot <- plot + theme_classic(20) + labs(y = "TRiC-bound (log10 TPM)", x = "Translatome (log10 TPM)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), 
        axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS1/D.TRiCpuromycin.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Bukau comparison
# Example
plot <- ggplot(aes(position, occupancy), data=bukau_dt["YDL195W", .(position, occupancy = movingAverage(B_ssb1_I_rpm, n=3, center=T))]) + 
  geom_vline(xintercept = bukau_dataset_dt[orf == "YDL195W"]$peak, color = "gray50", linetype = "dashed") +
  geom_line(aes(color = "1"),size = 1.25) + 
  geom_line(aes(position, occupancy2, color = "2"), data=bukau_dt["YDL195W", .(position, occupancy2 = movingAverage(B_ssb1_T_rpm, n=3, center=T))], size = 1.25) +
  #geom_vline(xintercept = bukau_dataset_dt[orf == "YDL195W"]$peak, color = "green", linetype = "dashed") +
  scale_color_manual(values = c("1" = "red", "2" = "black"), labels = c("Ssb IP", "Total"), name = "")
G <- plot + theme_classic(20) + labs(y = "Ribosome occupancy (RPM)", x = "Codon position") +
  theme(axis.text = element_text(size = 16, color = "black"), plot.title = element_text(size = 20, color = "black", face = "italic", hjust = 0.5)) +
  ggtitle("SEC31")
ggsave("/Users/KevinStein/Desktop/Figures/FigS1/F.SEC31.pdf", G, width = 6, height = 4, dpi = 300, useDingbats = F)

# Doring sites pvalues
bukau_dataset_dt[, peak_startaa := ceiling((peak_start / 3))]
bukau_dataset_dt[, peak_endaa := floor(((peak_start + width) / 3))]
doring_fishers[, position1 := as.character(base::paste(orf, position, sep = "_"))]
temp <- doring_fishers[doring_fishers$position1 %in% doring_peaks_all$position1]
plot <- ggplot(temp, aes(ssb1_padj)) + geom_histogram(color = 'black')
G <- plot + theme_classic(20) + labs(y = "Count", x = "Adjusted p-value") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/Figures/FigS1/G.DoringPeaks_pvalue.pdf", G, width = 5, height = 4, dpi = 300, useDingbats = F)

# Venn Diagram
length(bukau_dataset_dt[bukau_dataset_dt$orf %in% doring_fishers$orf]$orf)
length(bukau_peaks5$orf)
length(bukau_peaks5[bukau_peaks5$position1 %in% doring_peaks_all$position1]$orf)
#length(bukau_peaks_all[bukau_peaks_all$position1 %in% doring_peaks_all$position1]$orf)
library(VennDiagram)
venn.diagram(
  x = list(A = 1:1821, B= 377:11333), # 1821 peaks from my analysis; 10957 Doring peaks since removed 4 overlapping orfs; 1445 overlap
  category.names = c("Ours", "Doring"), 
  filename = 'Doring_comparison.png',
  lwd = 2,
  lty = 'solid',
  col = c('#33A02C', '#6A3D9A'),
  fill = c('#33A02C', '#6A3D9A'),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.03, 0.03),
  cat.fontfamily = "sans", inverted = T)

# no numbers
venn.diagram(
  x = list(A = 1:1821, B= 377:11333),
  category.names = c("Ours", "Doring"),
  filename = 'Doring_comparison1.png',
  lwd = 2,
  lty = 'solid',
  col = c('#33A02C', '#6A3D9A'),
  fill = c('#33A02C', '#6A3D9A'),
  ext.line.lwd = 0,
  cex = 0,
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(45, -45),
  cat.dist = c(0.06, 0.06),
  cat.fontfamily = "sans", inverted = T)


### Crosslinker
expression1_noX[, tricsubunit := ifelse(orf == "YJL014W" | orf == "YJR064W" | orf == "YJL008C", 1, 0)]
expression1_noX[, tricname := ifelse(orf == "YJL014W", "CCT3", 0)] 
expression1_noX[, tricname := ifelse(orf == "YJR064W", "CCT5", tricname)]
expression1_noX[, tricname := ifelse(orf == "YJL008C", "CCT8", tricname)]
expression_tric[, ratio1 := T1X_tpm / R1X_tpm]
expression_tric[, ratio2 := T2X_tpm / R2X_tpm]
expression1_noX[, ratio1 := T1O_tpm / R1O_tpm]
expression1_noX[, ratio2 := T2O_tpm / R2O_tpm]
expression_tric[, pur_ratio1 := T1pur_tpm / R1pur_tpm]
expression_tric[, pur_ratio2 := T2pur_tpm / R2pur_tpm]

expression_tric[, name1 := ifelse(orf == "YJL014W", "a.CCT3", 0)] 
expression_tric[, name1 := ifelse(orf == "YJR064W", "b.CCT5", name1)]
expression_tric[, name1 := ifelse(orf == "YJL008C", "c.CCT8", name1)]
expression_tric[, name1 := ifelse(orf == "YDR212W", "d.CCT1", name1)]
expression_tric[, name1 := ifelse(orf == "YIL142W", "e.CCT2", name1)]
expression_tric[, name1 := ifelse(orf == "YDL143W", "f.CCT4", name1)]
expression_tric[, name1 := ifelse(orf == "YDR188W", "g.CCT6", name1)]
expression_tric[, name1 := ifelse(orf == "YJL111W", "h.CCT7", name1)]
expression_tric[, name1 := ifelse(orf == "YDL229W", "i.SSB1", name1)]
expression_tric[, name1 := ifelse(orf == "YNL209W", "j.SSB2", name1)]
expression_tric[, name1 := ifelse(orf == "YML085C", "k.TUB1", name1)]
expression_tric[, name1 := ifelse(orf == "YML124C", "l.TUB3", name1)]

plot <- ggplot(expression_tric[tricsubunit == 1], aes("+\nRep1", ratio1)) + geom_point(size = 2.5) +
  geom_text(data = expression_tric[tricsubunit == 1], aes("+\nRep1", ratio1,label=tricname), hjust=1.1,vjust=-.1, size=4) +
  geom_point(data = expression_tric[tricsubunit == 1], aes("+\nRep2", ratio2), size = 2.5) +
  geom_text(data = expression_tric[tricsubunit == 1], aes("+\nRep2", ratio2,label=tricname), hjust=1.1,vjust=-.1, size=4) +
  geom_point(data = expression1_noX[tricsubunit == 1], aes("-\nRep1", ratio1), size = 2.5) +
  geom_text(data = expression1_noX[tricsubunit == 1], aes("-\nRep1", ratio1,label=tricname), hjust=1.1,vjust=-.1, size=4) +
  geom_point(data = expression1_noX[tricsubunit == 1], aes("-\nRep2", ratio2), size = 2.5) +
  geom_text(data = expression1_noX[tricsubunit == 1], aes("-\nRep2", ratio2,label=tricname), hjust=1.1,vjust=-.1, size=4) +
  scale_y_log10(limits = c(0.8,25), breaks = c(1,5,10,20)) + annotation_logticks(sides = "l", short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"))
plot <- plot + theme_classic(20) + labs(y = "Enrichment\n(TRiC tpm / Total tpm)", x = "") +
  theme(axis.text = element_text(size = 20, color = "black"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank())
ggsave("/Users/KevinStein/Desktop/Figures/FigS1/J.CrosslinkerComparison.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


temp <- as.data.table(cbind(as.character(tric_peaks5$orf), tric_peaks5$peak, tric_peaks5$position1))
colnames(temp) <- c("orf", "tric.position", "tric.name")
setkeyv(temp, c("orf"))
temp2 <- as.data.table(cbind(as.character(tric_noX_peaks5$orf), tric_noX_peaks5$peak, tric_noX_peaks5$position1))
colnames(temp2) <- c("orf", "noX.position", "noX.name")
setkeyv(temp2, c("orf"))
temp <- temp[temp$orf %in% temp2$orf]
temp2 <- temp2[temp2$orf %in% temp$orf]
temp3 <- temp2[temp, allow.cartesian = T]
temp3$tric.position <- as.integer(temp3$tric.position)
temp3$noX.position <- as.integer(temp3$noX.position)
temp3[, diff := tric.position - noX.position]
setkeyv(temp3, c("tric.name"))
temp4 <- temp3[, .SD[which.min(abs(diff))], by = tric.name]
setkeyv(temp3, c("noX.name"))
temp5 <- temp3[, .SD[which.min(abs(diff))], by = noX.name]

plot <- ggplot(temp5, aes(diff, color = "#A6CEE3")) + 
  geom_freqpoly(data = temp4, aes(diff, color = "#1F78B4"), size = 1, binwidth = 15) +
  geom_freqpoly(size = 1, binwidth = 10) +
  scale_fill_manual(labels = c("Crosslinked sites mapped to\nNot crosslinked", "Not crosslinked sites mapped to\ncrosslinked"), values = c("#1F78B4", "#A6CEE3"), name = "") +
  scale_color_manual(labels = c("Crosslinked sites mapped to\nNot crosslinked", "Not crosslinked sites mapped to\ncrosslinked"), values = c("#1F78B4", "#A6CEE3"), name = "") +
  coord_cartesian(xlim = c(-200,200)) 
plot <- plot + theme_classic(20) + labs(y = "Count", x = "Distance between binding sites (codons)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"),
        legend.position = c(.6,.99), legend.justification = c("left", "top"), legend.background = element_blank(), 
        legend.text = element_text(size = 10))
ggsave("/Users/KevinStein/Desktop/Figures/FigS1/K.ComparingSitesWithAndWithoutCrosslinker.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Response letter
temp8 <- tric_peaks5[tric_peaks5$orf %in% tric_noX_peaks5$orf, .SD[which.max(odds)], by = orf]
temp9 <- tric_noX_peaks5[tric_noX_peaks5$orf %in% tric_peaks5$orf, .SD[which.max(odds)], by = orf]
plot <- ggplot(temp8, aes(peak, temp9$peak)) + geom_point() + xlim(0,1000) + ylim(0,1000)
plot <- plot + theme_classic(20) + labs(y = "Not crosslinked sample (codon)", x = "Crosslinked sample (codons)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"))
ggsave("/Users/KevinStein/Desktop/ComparingSitesWithAndWithoutCrosslinker_scatter_R0.82.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
cor(temp8$peak, temp9$peak)


### Comparison of our sites to Bukau sites identified using DESeq, all sites
temp <- as.data.table(cbind(ssb_peaks5$orf, ssb_peaks5$peak, ssb_peaks5$position1))
colnames(temp) <- c("orf", "ssb.position", "ssb.name")
setkeyv(temp, c("orf"))
temp2 <- as.data.table(cbind(bukau_peaks5$orf, bukau_peaks5$peak, bukau_peaks5$position1))
colnames(temp2) <- c("orf", "bukau.position", "bukau.name")
setkeyv(temp2, c("orf"))
temp <- temp[temp$orf %in% temp2$orf]
temp2 <- temp2[temp2$orf %in% temp$orf]
temp3 <- temp2[temp, allow.cartesian = T]
#temp4 <- temp[temp2, allow.cartesian = T]
temp3$ssb.position <- as.integer(temp3$ssb.position)
temp3$bukau.position <- as.integer(temp3$bukau.position)
temp3[, diff := ssb.position - bukau.position]
setkeyv(temp3, c("ssb.name"))
temp4 <- temp3[, .SD[which.min(abs(diff))], by = ssb.name]
setkeyv(temp3, c("bukau.name"))
temp5 <- temp3[, .SD[which.min(abs(diff))], by = bukau.name]

plot <- ggplot(temp5, aes(diff, color = "#6A3D9A")) + geom_freqpoly(size = 1, binwidth = 15) + 
  geom_freqpoly(data = temp4, aes(diff, color = "#33A02C"), size = 1, binwidth = 15) +
  scale_fill_manual(labels = c("This study mapped to\nDoring, et al", "Doring, et al mapped\nto this study"), values = c("#33A02C", "#6A3D9A"), name = "") +
  scale_color_manual(labels = c("This study mapped to\nDoring, et al", "Doring, et al mapped\nto this study"), values = c("#33A02C", "#6A3D9A"), name = "") +
  coord_cartesian(xlim = c(-250,250)) 
plot <- plot + theme_classic(20) + labs(y = "Count", x = "Distance between binding sites (codons)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), axis.text = element_text(size = 20, color = "black"),
        legend.position = c(.6,.99), legend.justification = c("left", "top"), legend.background = element_blank(), 
        legend.text = element_text(size = 10))
ggsave("/Users/KevinStein/Desktop/Figures/FigS1/L.OurSitesMappedToBukau.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)
