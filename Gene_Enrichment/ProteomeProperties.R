### Load table with properties
properties <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Yeast/YeastProteomeCategorization/protein_properties.txt", header=T)
properties <- as.data.table(properties)

### Add longest domain
# Yes: SUPERFAMILY, Gene3D; Pfam if gene doesn't have domain in one of the first databases
# No: PANTHER, PIRSF, ProSitePatterns, ProSiteProfiles, SignalP_GRAM_POSITIVE, SignalP_GRAM_NEGATIVE, TMHMM,
#  SignalP_EUK, Coils, TIGRFAM, SMART, Hamap, PRINTS, CDD, ProDom, SFLD, Pfam

setkeyv(domains, c("SystematicGene"))
temp <- domains[Method == "Gene3D" | Method == "SUPERFAMILY", .SD[which.max(DomainLength)], by = SystematicGene]
temp1 <- domains[Method == "Pfam", .SD[which.max(DomainLength)], by = SystematicGene]
temp2 <- temp1[!temp1$SystematicGene %in% temp$SystematicGene]
temp3 <- rbind(temp,temp2)
i <- cbind(match(properties$orf, temp3$SystematicGene))
properties <- cbind(properties, DomainLength = temp3[i]$DomainLength)

setkeyv(domains1, c("SystematicGene"))
temp <- domains1[Method == "GENE3D" | Method == "SUPERFAMILY", .SD[which.max(DomainLength)], by = SystematicGene]
temp1 <- domains1[Method == "Pfam", .SD[which.max(DomainLength)], by = SystematicGene]
temp2 <- temp1[!temp1$SystematicGene %in% temp$SystematicGene]
temp3 <- rbind(temp,temp2)
i <- cbind(match(properties$orf, temp3$SystematicGene))
properties <- cbind(properties, DomainLength1 = temp3[i]$DomainLength)

### Subset by TRiC or SSB substrates
properties <- as.data.table(properties)
setkeyv(properties, c("orf"))
setkeyv(tric_substrates5, c("orf"))
setkeyv(ssb_substrates5, c("orf"))
properties_background <- properties[as.character(unique(tric_dt$orf))]
setkeyv(properties_background, c("orf"))
properties_tricsubstrates <- properties_background[tric_substrates5$orf]
properties_nottric <- properties_background[!tric_substrates5$orf]
properties_ssbsubstrates <- properties_background[ssb_substrates5$orf]
properties_notssb <- properties_background[!ssb_substrates5$orf]
properties_notsubstrates <- properties_background[!substrates$orf]

### Plot properties of substrates and non-substrates with violin plot
property1 <- NULL
tric_sig1 <- NULL
ssb_sig1 <- NULL
tric_ssb1 <- NULL
tric_med1 <- NULL
ssb_med1 <- NULL
nottric_med1 <- NULL
notssb_med1 <- NULL
for (i in colnames(properties)) {
  if (i != "orf") {
    filename <- paste("/Users/KevinStein/Desktop/ProteinProperties/", i, ".pdf", sep = "")  
    plot <- ggplot(data = properties_background, aes("1.Translatome", properties_background[[i]])) + geom_boxplot(fill = "gray75") + 
      geom_boxplot(data = properties_tricsubstrates, aes("3.TRiC", properties_tricsubstrates[[i]]), fill = "#1F78B4") +
      geom_boxplot(data = properties_nottric, aes("2.Not TRiC", properties_nottric[[i]]), fill = "#A6CEE3") +
      geom_boxplot(data = properties_ssbsubstrates, aes("5.SSB", properties_ssbsubstrates[[i]]), fill = "#33A02C") +
      geom_boxplot(data = properties_notssb, aes("4.Not SSB", properties_notssb[[i]]), fill = "#B2DF8A") +
      scale_x_discrete(labels = c("Translatome","Not TRiC","TRiC","Not SSB","SSB"))
    G <- plot + theme_classic(14) + labs(y = i, x = "") +
      theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
    ggsave(filename, G, width = 6, height = 4, dpi = 300)
    property <- i
    property1 <- c(property1, property)
    tric_sig <- wilcox.test(na.omit(properties_tricsubstrates[[i]]), na.omit(properties_nottric[[i]]))$p.value
    tric_sig1 <- c(tric_sig1, tric_sig)
    ssb_sig <- wilcox.test(na.omit(properties_ssbsubstrates[[i]]), na.omit(properties_notssb[[i]]))$p.value
    ssb_sig1 <- c(ssb_sig1, ssb_sig)
    tric_ssb <- wilcox.test(na.omit(properties_tricsubstrates[[i]]), na.omit(properties_ssbsubstrates[[i]]))$p.value
    tric_ssb1 <- c(tric_ssb1, tric_ssb)
    tric_med <- median(na.omit(properties_tricsubstrates[[i]]))
    tric_med1 <- c(tric_med1, tric_med)
    nottric_med <- median(na.omit(properties_nottric[[i]]))
    nottric_med1 <- c(nottric_med1, nottric_med)
    ssb_med <- median(na.omit(properties_ssbsubstrates[[i]]))
    ssb_med1 <- c(ssb_med1, ssb_med)
    notssb_med <- median(na.omit(properties_notssb[[i]]))
    notssb_med1 <- c(notssb_med1, notssb_med)
  }
}
properties_sig <- data.table(property = property1, tric_pvalue = tric_sig1,
                                ssb_pvalue = ssb_sig1, tric_ssb_pvalue = tric_ssb1,
                             tric_med = tric_med1, nottric_med = nottric_med1,
                             ssb_med = ssb_med1, notssb_med = notssb_med1)
properties_sig[, diff_med := tric_med - ssb_med]
View(properties_sig)

