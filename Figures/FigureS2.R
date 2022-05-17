# Properties
plot <- ggplot(data = properties_tricsubstrates, aes("1.TRiC", sheet)) + geom_boxplot(fill = "#1F78B4", notch = T) + 
  geom_boxplot(data = properties_ssbsubstrates, aes("2.SSB", sheet), fill = "#33A02C", notch = T) +
  geom_boxplot(data = properties_notsubstrates, aes("3.Not", sheet), fill = "gray50", notch = T) +
  geom_boxplot(data = properties_background, aes("4.Translatome", sheet), fill = "gray80", notch = T) +
  scale_x_discrete(labels = c("TRiC","Ssb","Not\nbound","Transla-\ntome")) +
  scale_y_continuous(limits = c(0,0.5), breaks = c(0,0.1,0.2,0.3,0.4,0.5))
plot <- plot + theme_classic(20) + labs(y = "Beta-sheet\npropensity (a.u.)", x = "") + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS2/A.Sheet.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = properties_tricsubstrates, aes("1.TRiC", helix)) + geom_boxplot(fill = "#1F78B4", notch = T) + 
  geom_boxplot(data = properties_ssbsubstrates, aes("2.SSB", helix), fill = "#33A02C", notch = T) +
  geom_boxplot(data = properties_notsubstrates, aes("3.Not", helix), fill = "gray50", notch = T) +
  geom_boxplot(data = properties_background, aes("4.Translatome", helix), fill = "gray80", notch = T) +
  scale_x_discrete(labels = c("TRiC","Ssb","Not\nbound","Transla-\ntome")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1))
plot <- plot + theme_classic(20) + labs(y = "Alpha-helix\npropensity (a.u.)", x = "") + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS2/A.Helix.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = properties_tricsubstrates, aes("1.TRiC", PPI)) + geom_boxplot(fill = "#1F78B4", notch = T) + 
  geom_boxplot(data = properties_ssbsubstrates, aes("2.SSB", PPI), fill = "#33A02C", notch = T) +
  geom_boxplot(data = properties_notsubstrates, aes("3.Not", PPI), fill = "gray50", notch = T) +
  geom_boxplot(data = properties_background, aes("4.Translatome", PPI), fill = "gray80", notch = T) +
  scale_x_discrete(labels = c("TRiC","Ssb","Not\nbound","Transla-\ntome")) +
  scale_y_continuous(limits = c(0,250), breaks = c(0,50,100,150,200,250))
plot <- plot + theme_classic(20) + labs(y = "# of protein-protein\ninteractions", x = "") + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS2/B.PPI.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = properties_tricsubstrates, aes("1.TRiC", ID)) + geom_boxplot(fill = "#1F78B4", notch = T) + 
  geom_boxplot(data = properties_ssbsubstrates, aes("2.SSB", ID), fill = "#33A02C", notch = T) +
  geom_boxplot(data = properties_notsubstrates, aes("3.Not", ID), fill = "gray50", notch = T) +
  geom_boxplot(data = properties_background, aes("4.Translatome", ID), fill = "gray80", notch = T) +
  scale_x_discrete(labels = c("TRiC","Ssb","Not\nbound","Transla-\ntome")) +
  scale_y_continuous(limits = c(0,0.8), breaks = c(0,0.25,0.5,0.75))
plot <- plot + theme_classic(20) + labs(y = "% Intrinsic disorder", x = "") + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS2/C.IntrinsicDisorder.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = properties_tricsubstrates, aes("1.TRiC", Protein.Length)) + geom_boxplot(fill = "#1F78B4", notch = T) + 
  geom_boxplot(data = properties_ssbsubstrates, aes("2.SSB", Protein.Length), fill = "#33A02C", notch = T) +
  geom_boxplot(data = properties_notsubstrates, aes("3.Not", Protein.Length), fill = "gray50", notch = T) +
  geom_boxplot(data = properties_background, aes("4.Translatome", Protein.Length), fill = "gray80", notch = T) +
  scale_x_discrete(labels = c("TRiC","Ssb","Not\nbound","Transla-\ntome")) +
  scale_y_continuous(limits = c(0,1500), breaks = c(0,500,1000,1500))
plot <- plot + theme_classic(20) + labs(y = "Protein length (aa)", x = "") + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS2/D.ProteinLength.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

plot <- ggplot(data = properties_tricsubstrates, aes("1.TRiC", DomainLength)) + geom_boxplot(fill = "#1F78B4", notch = T) + 
  geom_boxplot(data = properties_ssbsubstrates, aes("2.SSB", DomainLength), fill = "#33A02C", notch = T) +
  geom_boxplot(data = properties_notsubstrates, aes("3.Not", DomainLength), fill = "gray50", notch = T) +
  geom_boxplot(data = properties_background, aes("4.Translatome", DomainLength), fill = "gray80", notch = T) +
  scale_x_discrete(labels = c("TRiC","Ssb","Not\nbound","Transla-\ntome")) +
  scale_y_continuous(limits = c(0,800), breaks = c(0,200,400,600,800))
plot <- plot + theme_classic(20) + labs(y = "Length of longest\ndomain (aa)", x = "") + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black", size = 20))
ggsave("/Users/KevinStein/Desktop/Figures/FigS2/E.DomainLength.pdf", plot, width = 6, height = 4, dpi = 300, useDingbats = F)

