### Gene-level expression ###
# TRiC
expression_tric <- expression1[tric_deseq_gene]
expression_tric[, group := ifelse(padj < 0.05 & log2FoldChange >= 1 & ribo_tpm > 1, 1, 0)]
expression_tric[, tricsubunit := ifelse(orf == "YJL014W" | orf == "YJR064W" | orf == "YJL008C", 1, 0)]
expression_tric[, tricname := ifelse(orf == "YJL014W", "CCT3", 0)] 
expression_tric[, tricname := ifelse(orf == "YJR064W", "CCT5", tricname)]
expression_tric[, tricname := ifelse(orf == "YJL008C", "CCT8", tricname)]
i <- cbind(match(expression_tric$orf, domains1$orf))
expression_tric <- cbind(expression_tric, orfname = domains1[i]$name)

plot <- ggplot(expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & group == 0]) +
  geom_point(aes(ribo_tpm, tric_tpm), fill = "gray75", color = "gray75", alpha = 0.5, size = 2) + 
  geom_point(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & group == 1], aes(ribo_tpm, tric_tpm, color = factor(group)), fill = "#1F78B4", alpha = 0.8, size = 2) + 
  geom_point(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & tricsubunit == 1], aes(ribo_tpm, tric_tpm), color="#E31A1C", size = 2) +
  geom_text(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & group == 1], aes(ribo_tpm, tric_tpm,label=orfname), hjust=1.1,vjust=-.1, size=3) +
  geom_text(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & tricsubunit == 1], aes(ribo_tpm, tric_tpm,label=tricname), hjust=1.1,vjust=-.1, size=3) +
  scale_x_log10(limits = c(1e-1,1.5e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
  scale_y_log10(limits = c(1e-1,1.5e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm")) +
  scale_color_manual(limits = c("1"), labels = c("Substrates"), values = c("#1F78B4"), name = "Pearson's r = 0.99")
plot <- plot + theme_classic(20) + labs(y = "TRiC-bound (log10 TPM)", x = "Translatome (log10 TPM)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), legend.title=element_text(size=20), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))
cor(expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & orf != "YJL014W" & orf != "YJR064W" & orf != "YDR188W" & orf != "YJL008C"]$ribo_tpm,
    expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & orf != "YJL014W" & orf != "YJR064W" & orf != "YDR188W" & orf != "YJL008C"]$tric_tpm, method = "pearson")
ggsave("/Users/KevinStein/Desktop/Figures/Fig1/C.TRiC_deseq.gene.expression.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Ssb
expression_ssb <- expression1[ssb_deseq_gene]
expression_ssb[, group := ifelse(padj < 0.05 & log2FoldChange >= 1 & ssb_Rchx_tpm > 1, 1, 0)]
expression_ssb[, ssbsubunit := ifelse(orf == "YDL229W" | orf == "YNL209W", 1, 0)]
expression_ssb[, ssbname := ifelse(orf == "YDL229W", "SSB1", 0)] 
expression_ssb[, ssbname := ifelse(orf == "YNL209W", "SSB2", ssbname)]
i <- cbind(match(expression_ssb$orf, domains1$orf))
expression_ssb <- cbind(expression_ssb, orfname = domains1[i]$name)
plot <- ggplot(expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Schx1_total+ssb_Schx2_total) > 100 & group == 0]) +
  geom_point(aes(ssb_Rchx_tpm, ssb_Schx_tpm), color = "gray75", fill = "gray75", alpha = 0.5, size = 2) + 
  geom_point(data = expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Schx1_total+ssb_Schx2_total) > 100 & group == 1 & codon_group == 1], aes(ssb_Rchx_tpm, ssb_Schx_tpm, color = factor(group)), fill = "#33A02C", alpha = 0.8, size = 2) + 
  geom_point(data = expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Schx1_total+ssb_Schx2_total) > 100 & ssbsubunit == 1], aes(ssb_Rchx_tpm, ssb_Schx_tpm), color="#6A3D9A", size = 2) +
  geom_text(data = expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Schx1_total+ssb_Schx2_total) > 100 & ssbname == "SSB1"], aes(ssb_Rchx_tpm, ssb_Schx_tpm,label=ssbname), hjust=-.1,vjust=.5, size=3) +
  geom_text(data = expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Schx1_total+ssb_Schx2_total) > 100 & ssbname == "SSB2"], aes(ssb_Rchx_tpm, ssb_Schx_tpm,label=ssbname), hjust=1.1,vjust=.5, size=3) +
  scale_x_log10(limits = c(1e-1,1.5e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
  scale_y_log10(limits = c(1e-1,1.5e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm")) +
  scale_color_manual(limits = c("1"), labels = c("Substrates"), values = c("#33A02C"), name = "Pearson's r = 0.95")
plot <- plot + theme_classic(20) + labs(y = "Ssb-bound (log10 TPM)", x = "Translatome (log10 TPM)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), legend.title=element_text(size=20), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))
cor(expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Rchx1_total+ssb_Rchx2_total) > 100 & orf != "YDL229W" & orf != "YNL209W"]$ssb_Rchx_tpm,
    expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Rchx1_total+ssb_Rchx2_total) > 100 & orf != "YDL229W" & orf != "YNL209W"]$ssb_Schx_tpm, method = "pearson")
ggsave("/Users/KevinStein/Desktop/Figures/Fig1/C.SSB_deseq.gene.expression.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Codon-level expression
# TRiC
expression_tric[, codon_group := ifelse(expression_tric$orf %in% tric_substrates5$orf, 1, 0)]
plot <- ggplot(expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100]) +
  geom_point(aes(ribo_tpm, tric_tpm), color = "gray75", size = 2) + 
  geom_point(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & group == 1 & codon_group == 1], aes(ribo_tpm, tric_tpm, color = factor(codon_group)), size = 2) + 
  geom_point(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & tricsubunit == 1], aes(ribo_tpm, tric_tpm), color="#E31A1C", size = 2) +
  geom_text(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & group == 1 & log2FoldChange > 1.2], aes(ribo_tpm, tric_tpm,label=orfname), hjust=1.1,vjust=-.1, size=3) +
  geom_text(data = expression_tric[orf == "YJR121W" | orf == "YFL039C"], aes(ribo_tpm, tric_tpm,label=orfname), hjust=1.1,vjust=-.1, size=3) +
  geom_text(data = expression_tric[(R1X_total+R2X_total) > 100  & (T1X_total+T2X_total) > 100 & tricsubunit == 1], aes(ribo_tpm, tric_tpm,label=tricname), hjust=1.1,vjust=-.1, size=3) +
  scale_x_log10(limits = c(1e-1,1.5e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
  scale_y_log10(limits = c(1e-1,1.5e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm")) +
  scale_color_manual(limits = c("1"), labels = c("Substrates"), values = c("#1F78B4"), name = "Pearson's r = 0.99")
plot <- plot + theme_classic(20) + labs(y = "TRiC-bound (log10 TPM)", x = "Translatome (log10 TPM)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), legend.title=element_text(size=20), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/Fig1/C.TRiC_deseq.codon.expression.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

# Ssb
expression_ssb[, codon_group := ifelse(expression_ssb$orf %in% ssb_substrates5$orf, 1, 0)]
plot <- ggplot(expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Schx1_total+ssb_Schx2_total) > 100]) +
  geom_point(aes(ssb_Rchx_tpm, ssb_Schx_tpm), color = "gray75", size = 2) + 
  geom_point(data = expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Schx1_total+ssb_Schx2_total) > 100 & group == 1 & codon_group == 1], aes(ssb_Rchx_tpm, ssb_Schx_tpm, color = factor(codon_group)), size = 2) + 
  geom_point(data = expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Schx1_total+ssb_Schx2_total) > 100 & ssbsubunit == 1], aes(ssb_Rchx_tpm, ssb_Schx_tpm), color="#6A3D9A", size = 2) +
  geom_text(data = expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Schx1_total+ssb_Schx2_total) > 100 & group == 1 & log2FoldChange > 1.3], aes(ssb_Rchx_tpm, ssb_Schx_tpm,label=orfname), hjust=1.1,vjust=-.1, size=3) +
  geom_text(data = expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Schx1_total+ssb_Schx2_total) > 100 & ssbname == "SSB1"], aes(ssb_Rchx_tpm, ssb_Schx_tpm,label=ssbname), hjust=-.1,vjust=.5, size=3) +
  geom_text(data = expression_ssb[(ssb_Rchx1_total+ssb_Rchx2_total) > 100  & (ssb_Schx1_total+ssb_Schx2_total) > 100 & ssbname == "SSB2"], aes(ssb_Rchx_tpm, ssb_Schx_tpm,label=ssbname), hjust=1.1,vjust=.5, size=3) +
  scale_x_log10(limits = c(1e-1,1.5e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + 
  scale_y_log10(limits = c(1e-1,1.5e5), breaks = c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5), labels = c("-1", "0","1", "2", "3","4","5")) + annotation_logticks(short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm")) +
  scale_color_manual(limits = c("1"), labels = c("Substrates"), values = c("#33A02C"), name = "Pearson's r = 0.95")
plot <- plot + theme_classic(20) + labs(y = "Ssb-bound (log10 TPM)", x = "Translatome (log10 TPM)") +
  theme(legend.position = c(.01,.99), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), legend.title=element_text(size=20), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/Fig1/C.SSB_deseq.codon.expression.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)


### Volcano plots
plot <- ggplot(tric_fishers[orf != "YJL014W" & orf != "YJR064W" & orf != "YDR188W" & orf != "YJL008C" & tric_odds < Inf & tric_odds > 0]) + 
  #geom_rect(aes(xmin=0, xmax=Inf, ymin=-log(0.05,10), ymax=Inf), fill = "#DEEBF7", alpha = 0.3) +
  geom_point(aes(log2(tric_odds), -log10(tric_padj)), color = "gray75", fill = "gray75", alpha = 0.5, size = 2) +
  geom_point(data = tric_peaks5all[orf != "YJL014W" & orf != "YJR064W" & orf != "YDR188W" & orf != "YJL008C"], aes(log2(tric_odds), -log10(tric_padj), color = "#1F78B4"), fill = "#1F78B4", alpha = 0.8, size = 2) +
  geom_point(data = tric_peaks5all[orf == "YFL039C"], aes(log2(tric_odds), -log10(tric_padj)), color = "#E31A1C") +
  #geom_text(data = tric_peaks5all[position1 == "YFL039C_321" | position1 == "YFL039C_339" | position1 == "YFL039C_352"], aes(log2(tric_odds), -log10(tric_padj), label=position1), hjust=1.1,vjust=-.1, size=3) +
  scale_x_continuous() + scale_y_continuous() + 
  scale_color_manual(labels = c("Substrate residues"), values = c("#1F78B4"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "TRiC enrichment (log2 odds ratio)") +
  theme(legend.position = "none", legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), legend.title=element_text(size=20), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/Fig1/D.TRiC_deseq.codon.residues_volcano.tiff", plot, width = 6, height = 4, dpi = 300)

plot <- ggplot(ssb_fishers[orf != "YDL229W" & orf != "YNL209W" & ssb_odds < Inf & ssb_odds > 0]) + 
  geom_point(aes(log2(ssb_odds), -log10(ssb_padj)), color = "gray75", fill = "gray75", alpha = 0.5, size = 2) +
  geom_point(data = ssb_peaks5all[orf != "YDL229W" & orf != "YNL209W"], aes(log2(ssb_odds), -log10(ssb_padj), color = "#33A02C"), fill = "#33A02C", alpha = 0.8, size = 2) +
  scale_x_continuous() + scale_y_continuous() + 
  scale_color_manual(labels = c("Substrate residues"), values = c("#33A02C"), name = "")
plot <- plot + theme_classic(20) + labs(y = "-log10 (adjusted p-value)", x = "Ssb enrichment (log2 odds ratio)") +
  theme(legend.position = "none", legend.background = element_blank(), legend.justification = c("left", "top"), panel.border = element_rect(color = "black", fill = NA, size = 1), axis.line = element_blank(), legend.title=element_text(size=20), legend.text = element_text(size=20),
        axis.text = element_text(size = 20, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))
ggsave("/Users/KevinStein/Desktop/Figures/Fig1/D.SSB_deseq.codon.residues_volcano.tiff", plot, width = 6, height = 4, dpi = 300)


### Localization
localization <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Annotation/localization_chartron.csv", stringsAsFactors = T, header = T)
localization <- as.data.table(localization)
setkeyv(localization, c("orf"))
localization[, localization := ifelse(localization == "mito", "Mitochondria", as.character(localization))]
localization[, localization1 := ifelse(localization == "SS" | localization == "SSTMD" |
                                         localization == "TA" | localization == "TMD" |
                                         localization == "Exception", "ER-targeted", as.character(localization))]
localization$localization <- as.factor(localization$localization)
localization$localization1 <- as.factor(localization$localization1)
localization_tric <- localization[tric_peaks5_max]
localization_ssb <- localization[ssb_peaks5_max]
sum(localization_tric$localization == "Cytonuclear") / 565
sum(localization_tric$localization == "Mitochondria") / 565
sum(localization_tric[, localization == "Exception" | localization == "TA"]) / 565
sum(localization_tric[, localization == "SS" | localization == "SSTMD" | localization == "TMD"]) / 565
sum(localization_ssb$localization == "Cytonuclear") / 1343
sum(localization_ssb$localization == "Mitochondria") / 1343
sum(localization_ssb[, localization == "Exception" | localization == "TA"]) / 1343
sum(localization_ssb[, localization == "SS" | localization == "SSTMD" | localization == "TMD"]) / 1343
localization[, localization2 := ifelse(localization == "SS" | localization == "SSTMD" | localization == "TMD", "ER-targeted", as.character(localization))]
localization[, localization2 := ifelse(localization == "TA" | localization == "Exception", "Other", as.character(localization2))]
localization_tric[, localization2 := ifelse(localization == "SS" | localization == "SSTMD" | localization == "TMD", "ER-targeted", as.character(localization))]
localization_tric[, localization2 := ifelse(localization == "TA" | localization == "Exception", "Other", as.character(localization2))]
localization_ssb[, localization2 := ifelse(localization == "SS" | localization == "SSTMD" | localization == "TMD", "ER-targeted", as.character(localization))]
localization_ssb[, localization2 := ifelse(localization == "TA" | localization == "Exception", "Other", as.character(localization2))]

plot <- ggplot(localization_tric, aes(x="1.TRiC")) + geom_bar(aes(fill = factor(localization2)), color = "black", size = 0.2, position = position_fill(reverse = T), show.legend = T) +
  geom_bar(data = localization_ssb, aes(x="2.SSB", fill = factor(localization2)), color = "black", size = 0.2, position = position_fill(reverse = T), show.legend = F) +
  geom_bar(data = localization, aes(x="3.Translatome", fill = factor(localization2)), color = "black", size = 0.2, position = position_fill(reverse = T), show.legend = F) +
  scale_fill_brewer(type = "seq", palette = "Greys", direction = -1, name = "", guide = guide_legend(keywidth = 2, keyheight = 2)) +
  scale_x_discrete(labels = c("TRiC", "Ssb", "Translatome")) + guides(fill = guide_legend(reverse = TRUE))
plot <- plot + theme_classic(20) + labs(y = "Fraction of dataset", x = "") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/Fig1/E.SubstrateProportion_tricP0.002_ssbP1.05e-10.pdf", plot, width = 5, height = 4, dpi = 300, useDingbats = F)


### Venn Diagram
library(VennDiagram)
# substrates
length(tric_substrates5[tric_substrates5$orf %in% ssb_substrates5$orf]$orf)
length(tric_substrates5$orf)
length(ssb_substrates5$orf)
venn.diagram(
  x = list(A = 1:565, B= 93:1435), # 1343 ssb substrates, 565 tric substrates (92 tric only), 473 overlap, 1435 total
  category.names = c("TRiC", "SSB"),
  filename = 'substrates_venn.png',
  lwd = 2,
  lty = 'solid',
  col = c('#1F78B4', '#33A02C'),
  fill = c('#1F78B4', '#33A02C'),
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
  x = list(A = 1:565, B= 93:1435),
  category.names = c("TRiC", "SSB"),
  filename = 'substrates_venn1.png',
  lwd = 2,
  lty = 'solid',
  col = c('#1F78B4', '#33A02C'),
  fill = c('#1F78B4', '#33A02C'),
  cex = 0,
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(45, -45),
  cat.dist = c(0.03, 0.03),
  cat.fontfamily = "sans", inverted = T)

# residues in shared substrates
temp <- tric_peaks5[tric_peaks5$orf %in% ssb_substrates5$orf]
temp1 <- ssb_peaks5[ssb_peaks5$orf %in% tric_substrates5$orf]
length(temp$orf)
length(temp1$orf)
length(temp[temp$position1 %in% temp1$position1]$orf)
venn.diagram(
  x = list(A = 1:788, B= 676:1813), # shared substrates: 1138 ssb residues, 788 tric peaks (675 tric only), 113 overlap, 1813 total
  category.names = c("TRiC", "SSB"), 
  filename = 'residues_venn.png',
  lwd = 2,
  lty = 'solid',
  col = c('#1F78B4', '#33A02C'),
  fill = c('#1F78B4', '#33A02C'),
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
  x = list(A = 1:788, B= 676:1813),
  category.names = c("TRiC", "SSB"),
  filename = 'residues_venn1.png',
  lwd = 2,
  lty = 'solid',
  col = c('#1F78B4', '#33A02C'),
  fill = c('#1F78B4', '#33A02C'),
  ext.line.lwd = 0,
  cex = 0,
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(45, -45),
  cat.dist = c(0.03, 0.03),
  cat.fontfamily = "sans", inverted = T)

