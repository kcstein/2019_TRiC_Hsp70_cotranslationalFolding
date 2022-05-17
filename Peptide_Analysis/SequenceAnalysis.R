### Composition
composition_ssb <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/SequenceAnalysis/Composition/IndividualFiles/composition_ssb.csv", header = T)
composition_bukau <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/SequenceAnalysis/Composition/IndividualFiles/composition_bukau.csv", header = T)
composition_random <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/SequenceAnalysis/Composition/IndividualFiles/composition_random.csv", header = T)
composition_difference <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/SequenceAnalysis/Composition/composition_difference.csv", header = T)
composition_tric <- as.data.table(composition_tric)
composition_ssb <- as.data.table(composition_ssb)
composition_bukau <- as.data.table(composition_bukau)
composition_random <- as.data.table(composition_random)
composition_difference <- as.data.table(composition_difference)

composition_difference$identity <- factor(composition_difference$region,levels = c("TRiC_nascent", "SSB_nascent", "Bukau_nascent"))
plot <- ggplot(composition_difference[region == "TRiC_nascent" | region == "SSB_nascent" | region == "Bukau_nascent"], aes(x="A", y=A, fill = identity)) + geom_col(position = "dodge") + 
  geom_col(aes(x="C", y=C, fill = identity), position = "dodge") +
  geom_col(aes(x="D", y=D, fill = identity), position = "dodge") +
  geom_col(aes(x="E", y=E, fill = identity), position = "dodge") +
  geom_col(aes(x="F", y=F, fill = identity), position = "dodge") +
  geom_col(aes(x="G", y=G, fill = identity), position = "dodge") +
  geom_col(aes(x="H", y=H, fill = identity), position = "dodge") +
  geom_col(aes(x="I", y=I, fill = identity), position = "dodge") +
  geom_col(aes(x="K", y=K, fill = identity), position = "dodge") +
  geom_col(aes(x="L", y=L, fill = identity), position = "dodge") +
  geom_col(aes(x="M", y=M, fill = identity), position = "dodge") +
  geom_col(aes(x="N", y=N, fill = identity), position = "dodge") +
  geom_col(aes(x="P", y=P, fill = identity), position = "dodge") +
  geom_col(aes(x="Q", y=Q, fill = identity), position = "dodge") +
  geom_col(aes(x="R", y=R, fill = identity), position = "dodge") +
  geom_col(aes(x="S", y=S, fill = identity), position = "dodge") +
  geom_col(aes(x="T", y=T, fill = identity), position = "dodge") +
  geom_col(aes(x="V", y=V, fill = identity), position = "dodge") +
  geom_col(aes(x="W", y=W, fill = identity), position = "dodge") +
  geom_col(aes(x="Y", y=Y, fill = identity), position = "dodge") +
  scale_fill_manual(values = c("#1F78B4", "#33A02C","#6A3D9A"),limits = c("TRiC_nascent", "SSB_nascent", "Bukau_nascent"), labels = c("TRiC","SSB","Doring, et al"), name = "",
                    guide = guide_legend(keywidth = 0.75, keyheight = 0.75))
plot <- plot + theme_classic(12) + labs(y = "log2(fold change)", x = "") + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
  legend.position = c(.80,.35), legend.justification = c("left", "top"), axis.text = element_text(color = "black", size = 10))
ggsave("/Users/KevinStein/Desktop/Nascent.pdf", plot, width = 6, height = 4, dpi = 300)

composition_difference$identity <- factor(composition_difference$region,levels = c("TRiC_tunnel", "SSB_tunnel", "Bukau_tunnel"))
plot <- ggplot(composition_difference[region == "TRiC_tunnel" | region == "SSB_tunnel" | region == "Bukau_tunnel"], aes(x="A", y=A, fill = identity)) + geom_col(position = "dodge") + 
  geom_col(aes(x="C", y=C, fill = identity), position = "dodge") +
  geom_col(aes(x="D", y=D, fill = identity), position = "dodge") +
  geom_col(aes(x="E", y=E, fill = identity), position = "dodge") +
  geom_col(aes(x="F", y=F, fill = identity), position = "dodge") +
  geom_col(aes(x="G", y=G, fill = identity), position = "dodge") +
  geom_col(aes(x="H", y=H, fill = identity), position = "dodge") +
  geom_col(aes(x="I", y=I, fill = identity), position = "dodge") +
  geom_col(aes(x="K", y=K, fill = identity), position = "dodge") +
  geom_col(aes(x="L", y=L, fill = identity), position = "dodge") +
  geom_col(aes(x="M", y=M, fill = identity), position = "dodge") +
  geom_col(aes(x="N", y=N, fill = identity), position = "dodge") +
  geom_col(aes(x="P", y=P, fill = identity), position = "dodge") +
  geom_col(aes(x="Q", y=Q, fill = identity), position = "dodge") +
  geom_col(aes(x="R", y=R, fill = identity), position = "dodge") +
  geom_col(aes(x="S", y=S, fill = identity), position = "dodge") +
  geom_col(aes(x="T", y=T, fill = identity), position = "dodge") +
  geom_col(aes(x="V", y=V, fill = identity), position = "dodge") +
  geom_col(aes(x="W", y=W, fill = identity), position = "dodge") +
  geom_col(aes(x="Y", y=Y, fill = identity), position = "dodge") +
  scale_fill_manual(values = c("#1F78B4", "#33A02C","#6A3D9A"),limits = c("TRiC_tunnel", "SSB_tunnel", "Bukau_tunnel"), labels = c("TRiC","SSB","Doring, et al"), name = "",
                    guide = guide_legend(keywidth = 0.75, keyheight = 0.75))
plot <- plot + theme_classic(12) + labs(y = "log2(fold change)", x = "") + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
                                                                                 legend.position = c(.80,.95), legend.justification = c("left", "top"), axis.text = element_text(color = "black", size = 10))
ggsave("/Users/KevinStein/Desktop/Tunnel.pdf", plot, width = 6, height = 4, dpi = 300)

composition_difference$identity <- factor(composition_difference$region,levels = c("TRiC_active", "SSB_active", "Bukau_active"))
plot <- ggplot(composition_difference[region == "TRiC_active" | region == "SSB_active" | region == "Bukau_active"], aes(x="A", y=A, fill = identity)) + geom_col(position = "dodge") + 
  geom_col(aes(x="C", y=C, fill = identity), position = "dodge") +
  geom_col(aes(x="D", y=D, fill = identity), position = "dodge") +
  geom_col(aes(x="E", y=E, fill = identity), position = "dodge") +
  geom_col(aes(x="F", y=F, fill = identity), position = "dodge") +
  geom_col(aes(x="G", y=G, fill = identity), position = "dodge") +
  geom_col(aes(x="H", y=H, fill = identity), position = "dodge") +
  geom_col(aes(x="I", y=I, fill = identity), position = "dodge") +
  geom_col(aes(x="K", y=K, fill = identity), position = "dodge") +
  geom_col(aes(x="L", y=L, fill = identity), position = "dodge") +
  geom_col(aes(x="M", y=M, fill = identity), position = "dodge") +
  geom_col(aes(x="N", y=N, fill = identity), position = "dodge") +
  geom_col(aes(x="P", y=P, fill = identity), position = "dodge") +
  geom_col(aes(x="Q", y=Q, fill = identity), position = "dodge") +
  geom_col(aes(x="R", y=R, fill = identity), position = "dodge") +
  geom_col(aes(x="S", y=S, fill = identity), position = "dodge") +
  geom_col(aes(x="T", y=T, fill = identity), position = "dodge") +
  geom_col(aes(x="V", y=V, fill = identity), position = "dodge") +
  geom_col(aes(x="W", y=W, fill = identity), position = "dodge") +
  geom_col(aes(x="Y", y=Y, fill = identity), position = "dodge") +
  scale_fill_manual(values = c("#1F78B4", "#33A02C","#6A3D9A"),limits = c("TRiC_active", "SSB_active", "Bukau_active"), labels = c("TRiC","SSB","Doring, et al"), name = "",
                    guide = guide_legend(keywidth = 0.75, keyheight = 0.75))
plot <- plot + theme_classic(12) + labs(y = "log2(fold change)", x = "") + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
                                                                                 legend.position = c(.80,.95), legend.justification = c("left", "top"), axis.text = element_text(color = "black", size = 10))
ggsave("/Users/KevinStein/Desktop/active.pdf", plot, width = 6, height = 4, dpi = 300)

composition_difference$identity <- factor(composition_difference$region,levels = c("TRiC_downstream", "SSB_downstream", "Bukau_downstream"))
plot <- ggplot(composition_difference[region == "TRiC_downstream" | region == "SSB_downstream" | region == "Bukau_downstream"], aes(x="A", y=A, fill = identity)) + geom_col(position = "dodge") + 
  geom_col(aes(x="C", y=C, fill = identity), position = "dodge") +
  geom_col(aes(x="D", y=D, fill = identity), position = "dodge") +
  geom_col(aes(x="E", y=E, fill = identity), position = "dodge") +
  geom_col(aes(x="F", y=F, fill = identity), position = "dodge") +
  geom_col(aes(x="G", y=G, fill = identity), position = "dodge") +
  geom_col(aes(x="H", y=H, fill = identity), position = "dodge") +
  geom_col(aes(x="I", y=I, fill = identity), position = "dodge") +
  geom_col(aes(x="K", y=K, fill = identity), position = "dodge") +
  geom_col(aes(x="L", y=L, fill = identity), position = "dodge") +
  geom_col(aes(x="M", y=M, fill = identity), position = "dodge") +
  geom_col(aes(x="N", y=N, fill = identity), position = "dodge") +
  geom_col(aes(x="P", y=P, fill = identity), position = "dodge") +
  geom_col(aes(x="Q", y=Q, fill = identity), position = "dodge") +
  geom_col(aes(x="R", y=R, fill = identity), position = "dodge") +
  geom_col(aes(x="S", y=S, fill = identity), position = "dodge") +
  geom_col(aes(x="T", y=T, fill = identity), position = "dodge") +
  geom_col(aes(x="V", y=V, fill = identity), position = "dodge") +
  geom_col(aes(x="W", y=W, fill = identity), position = "dodge") +
  geom_col(aes(x="Y", y=Y, fill = identity), position = "dodge") +
  scale_fill_manual(values = c("#1F78B4", "#33A02C","#6A3D9A"),limits = c("TRiC_downstream", "SSB_downstream", "Bukau_downstream"), labels = c("TRiC","SSB","Doring, et al"), name = "",
                    guide = guide_legend(keywidth = 0.75, keyheight = 0.75))
plot <- plot + theme_classic(12) + labs(y = "log2(fold change)", x = "") + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
                                                                                 legend.position = c(.80,.95), legend.justification = c("left", "top"), axis.text = element_text(color = "black", size = 10))
ggsave("/Users/KevinStein/Desktop/Downstream.pdf", plot, width = 6, height = 4, dpi = 300)


ssbnascent <- NULL
ssbtunnel <- NULL
ssbactive <- NULL
ssbdownstream <- NULL
tricnascent <- NULL
trictunnel <- NULL
tricactive <- NULL
tricdownstream <- NULL
bukaunascent <- NULL
bukautunnel <- NULL
bukauactive <- NULL
bukaudownstream <- NULL
for (i in 2:21) {
  tricnascent1 <- wilcox.test(composition_tric[region == "nascent"][[i]], composition_random[region == "30mer"][[i]])$p.value
  trictunnel1 <- wilcox.test(composition_tric[region == "tunnel"][[i]], composition_random[region == "30mer"][[i]])$p.value
  tricactive1 <- wilcox.test(composition_tric[region == "active"][[i]], composition_random[region == "5mer"][[i]])$p.value
  tricdownstream1 <- wilcox.test(composition_tric[region == "downstream"][[i]], composition_random[region == "4mer"][[i]])$p.value
  ssbnascent1 <- wilcox.test(composition_ssb[region == "nascent"][[i]], composition_random[region == "30mer"][[i]])$p.value
  ssbtunnel1 <- wilcox.test(composition_ssb[region == "tunnel"][[i]], composition_random[region == "30mer"][[i]])$p.value
  ssbactive1 <- wilcox.test(composition_ssb[region == "active"][[i]], composition_random[region == "5mer"][[i]])$p.value
  ssbdownstream1 <- wilcox.test(composition_ssb[region == "downstream"][[i]], composition_random[region == "4mer"][[i]])$p.value
  bukaunascent1 <- wilcox.test(composition_bukau[region == "nascent"][[i]], composition_random[region == "30mer"][[i]])$p.value
  bukautunnel1 <- wilcox.test(composition_bukau[region == "tunnel"][[i]], composition_random[region == "30mer"][[i]])$p.value
  bukauactive1 <- wilcox.test(composition_bukau[region == "active"][[i]], composition_random[region == "5mer"][[i]])$p.value
  bukaudownstream1 <- wilcox.test(composition_bukau[region == "downstream"][[i]], composition_random[region == "4mer"][[i]])$p.value
  tricnascent <- c(tricnascent, tricnascent1)
  trictunnel <- c(trictunnel, trictunnel1)
  tricactive <- c(tricactive, tricactive1)
  tricdownstream <- c(tricdownstream, tricdownstream1)
  ssbnascent <- c(ssbnascent, ssbnascent1)
  ssbtunnel <- c(ssbtunnel, ssbtunnel1)
  ssbactive <- c(ssbactive, ssbactive1)
  ssbdownstream <- c(ssbdownstream, ssbdownstream1)
  bukaunascent <- c(bukaunascent, bukaunascent1)
  bukautunnel <- c(bukautunnel, bukautunnel1)
  bukauactive <- c(bukauactive, bukauactive1)
  bukaudownstream <- c(bukaudownstream, bukaudownstream1)
}
aa <- c(colnames(composition_tric[, c(2:21)]))

composition1 <- data.table(aa = aa, tricnascent = tricnascent,
                            trictunnel = trictunnel,
                            tricactive = tricactive,
                            tricdownstream = tricdownstream,
                            ssbnascent = ssbnascent,
                            ssbtunnel = ssbtunnel,
                            ssbactive = ssbactive,
                            ssbdownstream = ssbdownstream,
                            bukautunnel = bukautunnel,
                            bukaunascent = bukaunascent,
                            bukauactive = bukauactive,
                            bukaudownstream = bukaudownstream)


### Codon optimality and codon usage
nTE_dt <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/CodonUsage/normalizedTE_scale.csv", header = T)
cTE_dt <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/CodonUsage/classicalTE_scale.csv", header = T)
tAI_dt <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/CodonUsage/tAI_stAIcalc.csv", header = T)
codonusage_dt <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/CodonUsage/CodonUsage_yeast.csv", header = T)
nTE_dt <- as.data.table(nTE_dt)
cTE_dt <- as.data.table(cTE_dt)
tAI_dt <- as.data.table(tAI_dt)
codonusage_dt <- as.data.table(codonusage_dt)

temp <- data.table(orf = tric_peaks5_max[,1],
                   name = tric_peaks5_max[,2],
                   length = tric_peaks5_max[,12],
                   stop_peak = tric_peaks5_max[,13],
                   peak = tric_peaks5_max[,9],
                   peak_start = tric_peaks5_max[,7])
names(temp) <- c("orf", "name", "length", "stop_peak", "peak", "peak_start")

random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(30:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

setkeyv(temp, c("orf"))
setkeyv(tric_fishers, c("orf"))
tric_tAI <- tric_fishers[temp, allow.cartesian = TRUE]
i <- cbind(match(tric_tAI$position1, tric_dt$position1))
tric_tAI <- cbind(tric_tAI, codon = tric_dt[i]$codon)
i <- cbind(match(tric_tAI$codon, codonusage_dt$codon))
tric_tAI <- cbind(tric_tAI, per1000 = codonusage_dt[i]$per1000)
tric_tAI <- cbind(tric_tAI, fraction = codonusage_dt[i]$fraction)
i <- cbind(match(tric_tAI$codon, cTE_dt$codons))
tric_tAI <- cbind(tric_tAI, cTE = cTE_dt[i]$Scer)
i <- cbind(match(tric_tAI$codon, nTE_dt$codons))
tric_tAI <- cbind(tric_tAI, nTE = nTE_dt[i]$Scer)
i <- cbind(match(tric_tAI$codon, tAI_dt$codons))
tric_tAI <- cbind(tric_tAI, tAI = tAI_dt[i]$Scer)
tric_tAI[, adjusted := position - peak]
tric_tAI[, adjusted_start := position - peak_start]
tric_tAI[, adjusted_random := position - random]
setkeyv(tric_tAI, c("name"))
tric_tAI[, tAI_ma := movingAverage(tAI, n=3, center=T), by = name]

tric_tAI_max <- tric_fishers[temp, allow.cartesian = TRUE]
i <- cbind(match(tric_tAI_max$position1, tric_dt$position1))
tric_tAI_max <- cbind(tric_tAI_max, codon = tric_dt[i]$codon)
i <- cbind(match(tric_tAI_max$codon, codonusage_dt$codon))
tric_tAI_max <- cbind(tric_tAI_max, per1000 = codonusage_dt[i]$per1000)
tric_tAI_max <- cbind(tric_tAI_max, fraction = codonusage_dt[i]$fraction)
i <- cbind(match(tric_tAI_max$codon, cTE_dt$codons))
tric_tAI_max <- cbind(tric_tAI_max, cTE = cTE_dt[i]$Scer)
i <- cbind(match(tric_tAI_max$codon, nTE_dt$codons))
tric_tAI_max <- cbind(tric_tAI_max, nTE = nTE_dt[i]$Scer)
i <- cbind(match(tric_tAI_max$codon, tAI_dt$codons))
tric_tAI_max <- cbind(tric_tAI_max, tAI = tAI_dt[i]$Scer)
tric_tAI_max[, adjusted := position - peak]
tric_tAI_max[, adjusted_start := position - peak_start]
tric_tAI_max[, adjusted_random := position - random]

test <- tric_tAI[adjusted <= -20 & adjusted > -40]
setkeyv(test, c("name"))
setkeyv(tric_tAI, c("name"))
test[, tAI_avg := mean(tAI), by = name]
i <- cbind(match(tric_tAI$name, test$name))
tric_tAI <- cbind(tric_tAI, tAI_avg = test[i]$tAI_avg)
tric_tAI[, tAI_norm := tAI / tAI_avg]

ggplot(data = tric_tAI[, .(adjusted, tAI1 = movingAverage(tAI, n=5, center=T))]) +
  stat_summary(aes(adjusted, tAI1), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(adjusted, tAI1), fun.y = "median", geom = "line", size = 1, color = '#1F78B4') +
  stat_summary(data = tric_tAI[, .(adjusted_random, tAI2 = movingAverage(tAI, n=5, center=T))], aes(adjusted_random, tAI2), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = tric_tAI[, .(adjusted_random, tAI2 = movingAverage(tAI, n=5, center=T))], aes(adjusted_random, tAI2), fun.y = "median", geom = "line", size = 1, color = 'gray50') +
  xlim(-25, 25)

ggplot(data = tric_tAI[tric_tAI$name %in% tric_Early$position1 & stopdist < -25 & peak >= 70 & ribo_rpc >= 1, .(adjusted, tAI1 = movingAverage(tAI_norm, n=5, center=T))]) +
  stat_summary(aes(adjusted, tAI1, fill = '#1F78B4'), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5)) +
  stat_summary(aes(adjusted, tAI1, color = '#1F78B4'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(data = tric_tAI[tric_tAI$name %in% tric_Late$position1 & stopdist < -25 & peak >= 70 & ribo_rpc >= 1, .(adjusted, tAI2 = movingAverage(tAI_norm, n=5, center=T))], aes(adjusted, tAI2, fill = 'gray50'), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5)) +
  stat_summary(data = tric_tAI[tric_tAI$name %in% tric_Late$position1 & stopdist < -25 & peak >= 70 & ribo_rpc >= 1, .(adjusted, tAI2 = movingAverage(tAI_norm, n=5, center=T))], aes(adjusted, tAI2, color = 'gray50'), fun.y = "mean", geom = "line", size = 1.25) +
  xlim(-25, 25) + 
  scale_fill_manual(labels = c("High Z", "Low Z"), values = c("#1F78B4", "#A6CEE3"), name = "") +
  scale_color_manual(labels = c("High Z", "Low Z"), values = c("#1F78B4", "#A6CEE3"), name = "")
  scale_y_continuous(breaks = c(0.5,1,1.5)) +
  coord_cartesian(ylim = c(0.5,1.5))


# Ssb
temp <- data.table(orf = ssb_peaks5_max[,1],
                   name = ssb_peaks5_max[,2],
                   length = ssb_peaks5_max[,12],
                   stop_peak = ssb_peaks5_max[,13],
                   peak = ssb_peaks5_max[,9],
                   peak_start = ssb_peaks5_max[,7])
names(temp) <- c("orf", "name", "length", "stop_peak", "peak", "peak_start")

random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(30:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

setkeyv(temp, c("orf"))
setkeyv(ssb_fishers, c("orf"))
ssb_tAI <- ssb_fishers[temp, allow.cartesian = TRUE]
i <- cbind(match(ssb_tAI$position1, tric_dt$position1))
ssb_tAI <- cbind(ssb_tAI, codon = tric_dt[i]$codon)
i <- cbind(match(ssb_tAI$codon, codonusage_dt$codon))
ssb_tAI <- cbind(ssb_tAI, per1000 = codonusage_dt[i]$per1000)
ssb_tAI <- cbind(ssb_tAI, fraction = codonusage_dt[i]$fraction)
i <- cbind(match(ssb_tAI$codon, cTE_dt$codons))
ssb_tAI <- cbind(ssb_tAI, cTE = cTE_dt[i]$Scer)
i <- cbind(match(ssb_tAI$codon, nTE_dt$codons))
ssb_tAI <- cbind(ssb_tAI, nTE = nTE_dt[i]$Scer)
i <- cbind(match(ssb_tAI$codon, tAI_dt$codons))
ssb_tAI <- cbind(ssb_tAI, tAI = tAI_dt[i]$Scer)
ssb_tAI[, adjusted := position - peak]
ssb_tAI[, adjusted_start := position - peak_start]
ssb_tAI[, adjusted_random := position - random]
setkeyv(ssb_tAI, c("name"))
ssb_tAI[, tAI_ma := movingAverage(tAI, n=3, center=T), by = name]


ssb_tAI_max <- ssb_fishers[temp, allow.cartesian = TRUE]
i <- cbind(match(ssb_tAI_max$position1, tric_dt$position1))
ssb_tAI_max <- cbind(ssb_tAI_max, codon = tric_dt[i]$codon)
i <- cbind(match(ssb_tAI_max$codon, codonusage_dt$codon))
ssb_tAI_max <- cbind(ssb_tAI_max, per1000 = codonusage_dt[i]$per1000)
ssb_tAI_max <- cbind(ssb_tAI_max, fraction = codonusage_dt[i]$fraction)
i <- cbind(match(ssb_tAI_max$codon, cTE_dt$codons))
ssb_tAI_max <- cbind(ssb_tAI_max, cTE = cTE_dt[i]$Scer)
i <- cbind(match(ssb_tAI_max$codon, nTE_dt$codons))
ssb_tAI_max <- cbind(ssb_tAI_max, nTE = nTE_dt[i]$Scer)
i <- cbind(match(ssb_tAI_max$codon, tAI_dt$codons))
ssb_tAI_max <- cbind(ssb_tAI_max, tAI = tAI_dt[i]$Scer)
ssb_tAI_max[, adjusted := position - peak]
ssb_tAI_max[, adjusted_start := position - peak_start]
ssb_tAI_max[, adjusted_random := position - random]

test <- ssb_tAI[adjusted <= -20 & adjusted > -40]
setkeyv(test, c("name"))
setkeyv(ssb_tAI, c("name"))
test[, tAI_avg := mean(tAI), by = name]
i <- cbind(match(ssb_tAI$name, test$name))
ssb_tAI <- cbind(ssb_tAI, tAI_avg = test[i]$tAI_avg)
ssb_tAI[, tAI_norm := tAI / tAI_avg]

ggplot(data = ssb_tAI[, .(adjusted, tAI1 = movingAverage(tAI_norm, n=5, center=T))]) +
  stat_summary(aes(adjusted, tAI1), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(adjusted, tAI1), fun.y = "median", geom = "line", size = 1, color = '#1F78B4') +
  stat_summary(data = ssb_tAI[, .(adjusted_random, tAI2 = movingAverage(tAI_norm, n=5, center=T))], aes(adjusted_random, tAI2), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = ssb_tAI[, .(adjusted_random, tAI2 = movingAverage(tAI_norm, n=5, center=T))], aes(adjusted_random, tAI2), fun.y = "median", geom = "line", size = 1, color = 'gray50') +
  xlim(-25, 25)


# Bukau
temp <- data.table(orf = bukau_peaks5_max[,1],
                   name = bukau_peaks5_max[,2],
                   length = bukau_peaks5_max[,12],
                   stop_peak = bukau_peaks5_max[,13],
                   peak = bukau_peaks5_max[,9],
                   peak_start = bukau_peaks5_max[,7])
names(temp) <- c("orf", "name", "length", "stop_peak", "peak", "peak_start")

random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(30:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

setkeyv(temp, c("orf"))
setkeyv(bukau_fishers, c("orf"))
bukau_tAI <- bukau_fishers[temp, allow.cartesian = TRUE]
i <- cbind(match(bukau_tAI$position1, tric_dt$position1))
bukau_tAI <- cbind(bukau_tAI, codon = tric_dt[i]$codon)
i <- cbind(match(bukau_tAI$codon, codonusage_dt$codon))
bukau_tAI <- cbind(bukau_tAI, per1000 = codonusage_dt[i]$per1000)
bukau_tAI <- cbind(bukau_tAI, fraction = codonusage_dt[i]$fraction)
i <- cbind(match(bukau_tAI$codon, cTE_dt$codons))
bukau_tAI <- cbind(bukau_tAI, cTE = cTE_dt[i]$Scer)
i <- cbind(match(bukau_tAI$codon, nTE_dt$codons))
bukau_tAI <- cbind(bukau_tAI, nTE = nTE_dt[i]$Scer)
i <- cbind(match(bukau_tAI$codon, tAI_dt$codons))
bukau_tAI <- cbind(bukau_tAI, tAI = tAI_dt[i]$Scer)
bukau_tAI[, adjusted := position - peak]
bukau_tAI[, adjusted_start := position - peak_start]
bukau_tAI[, adjusted_random := position - random]
setkeyv(bukau_tAI, c("name"))
bukau_tAI[, tAI_ma := movingAverage(tAI, n=3, center=T), by = name]


bukau_tAI_max <- bukau_fishers[temp, allow.cartesian = TRUE]
i <- cbind(match(bukau_tAI_max$position1, tric_dt$position1))
bukau_tAI_max <- cbind(bukau_tAI_max, codon = tric_dt[i]$codon)
i <- cbind(match(bukau_tAI_max$codon, codonusage_dt$codon))
bukau_tAI_max <- cbind(bukau_tAI_max, per1000 = codonusage_dt[i]$per1000)
bukau_tAI_max <- cbind(bukau_tAI_max, fraction = codonusage_dt[i]$fraction)
i <- cbind(match(bukau_tAI_max$codon, cTE_dt$codons))
bukau_tAI_max <- cbind(bukau_tAI_max, cTE = cTE_dt[i]$Scer)
i <- cbind(match(bukau_tAI_max$codon, nTE_dt$codons))
bukau_tAI_max <- cbind(bukau_tAI_max, nTE = nTE_dt[i]$Scer)
i <- cbind(match(bukau_tAI_max$codon, tAI_dt$codons))
bukau_tAI_max <- cbind(bukau_tAI_max, tAI = tAI_dt[i]$Scer)
bukau_tAI_max[, adjusted := position - peak]
bukau_tAI_max[, adjusted_start := position - peak_start]
bukau_tAI_max[, adjusted_random := position - random]

test <- bukau_tAI[adjusted <= -20 & adjusted > -40]
setkeyv(test, c("name"))
setkeyv(bukau_tAI, c("name"))
test[, tAI_avg := mean(tAI), by = name]
i <- cbind(match(bukau_tAI$name, test$name))
bukau_tAI <- cbind(bukau_tAI, tAI_avg = test[i]$tAI_avg)
bukau_tAI[, tAI_norm := tAI / tAI_avg]

ggplot(data = bukau_tAI[, .(adjusted, tAI1 = movingAverage(tAI_norm, n=5, center=T))]) +
  stat_summary(aes(adjusted, tAI1), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(adjusted, tAI1), fun.y = "median", geom = "line", size = 1, color = '#1F78B4') +
  stat_summary(data = bukau_tAI[, .(adjusted_random, tAI2 = movingAverage(tAI_norm, n=5, center=T))], aes(adjusted_random, tAI2), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = bukau_tAI[, .(adjusted_random, tAI2 = movingAverage(tAI_norm, n=5, center=T))], aes(adjusted_random, tAI2), fun.y = "median", geom = "line", size = 1, color = 'gray50') +
  xlim(-25, 25)

ggplot(data = bukau_tAI[stop_peak < -25 & peak >= 70, .(adjusted, tAI1 = movingAverage(tAI, n=5, center=T))]) + xlim(-25, 25) +
  stat_summary(aes(adjusted, tAI1), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#33A02C') +
  stat_summary(aes(adjusted, tAI1), fun.y = "median", geom = "line", size = 1.25, color = '#33A02C') + 
  coord_cartesian(ylim = c(0.28,0.5))
  

# SRP
i <- cbind(match(srp_translatome_dt$codon, codonusage_dt$codon))
srp_translatome_dt <- cbind(srp_translatome_dt, per1000 = codonusage_dt[i]$per1000)
srp_translatome_dt <- cbind(srp_translatome_dt, fraction = codonusage_dt[i]$fraction)
i <- cbind(match(srp_translatome_dt$codon, cTE_dt$codons))
srp_translatome_dt <- cbind(srp_translatome_dt, cTE = cTE_dt[i]$Scer)
i <- cbind(match(srp_translatome_dt$codon, nTE_dt$codons))
srp_translatome_dt <- cbind(srp_translatome_dt, nTE = nTE_dt[i]$Scer)
i <- cbind(match(srp_translatome_dt$codon, tAI_dt$codons))
srp_translatome_dt <- cbind(srp_translatome_dt, tAI = tAI_dt[i]$Scer)
srp_translatome_dt[, tAI_ma := movingAverage(tAI, n=3, center=T), by = orf]

ggplot(data = srp_translatome_dt[, .(adjusted_start, tAI1 = movingAverage(tAI, n=5, center=T))]) +
  stat_summary(aes(adjusted, tAI1), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(adjusted, tAI1), fun.y = "median", geom = "line", size = 1, color = '#1F78B4') +
  stat_summary(data = srp_translatome_dt[, .(adjusted_random, tAI2 = movingAverage(tAI, n=5, center=T))], aes(adjusted_random, tAI2), fun.data = "median_hilow", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = srp_translatome_dt[, .(adjusted_random, tAI2 = movingAverage(tAI, n=5, center=T))], aes(adjusted_random, tAI2), fun.y = "median", geom = "line", size = 1, color = 'gray50') +
  xlim(-25, 25)

ggplot(data = srp_translatome_dt[WT_adjusted_rpc >= 0.5]) + xlim(-25, 25) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(aes(adjusted_random, per1000), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(aes(adjusted_random, per1000), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, per1000), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(adjusted_start, per1000), fun.y = "mean", geom = "line", size = 1.25, color = '#1F78B4')

