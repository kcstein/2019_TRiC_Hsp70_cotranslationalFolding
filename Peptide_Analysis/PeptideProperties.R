### Aggregation propensity
# TRiC
tango_tric <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/PeptideProperties/Tango/tric/tango_tric.txt", header = F)
colnames(tango_tric) <- c("orf", "position", "residue", "beta", "turn", "helix", "aggregation", "conc.stab.aggregation")
tango_tric$orf <- as.character(tango_tric$orf)
orf1 <- NULL
for (i in 1:nrow(tango_tric)) {
  str <- unlist(strsplit(tango_tric[i,]$orf, ".txt"))
  orf1 <- c(orf1, str)
}
tango_tric$orf <- orf1
tango_tric <- as.data.table(tango_tric)
setkeyv(tango_tric, c("orf"))
temp <- tric_peaks_all[, c(1,5,28,160)]

random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(60:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

tango_tric_dt <- tango_tric[temp, allow.cartesian = TRUE]
setkeyv(tango_tric_dt, c("orf"))
tango_tric_dt[, adjusted := position - peak]
tango_tric_dt[, adjusted_random := position - random]

# SSB
tango_ssb <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/PeptideProperties/Tango/ssb/tango_ssb.txt", header = F)
colnames(tango_ssb) <- c("orf", "position", "residue", "beta", "turn", "helix", "aggregation", "conc.stab.aggregation")
tango_ssb$orf <- as.character(tango_ssb$orf)
orf1 <- NULL
for (i in 1:nrow(tango_ssb)) {
  str <- unlist(strsplit(tango_ssb[i]$orf, ".txt"))
  orf1 <- c(orf1, str)
}
tango_ssb$orf <- orf1
tango_ssb <- as.data.table(tango_ssb)
setkeyv(tango_ssb, c("orf"))
temp <- ssb_peaks_all[, c(1,5,28,160)]

random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(60:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

tango_ssb_dt <- tango_ssb[temp, allow.cartesian = TRUE]
setkeyv(tango_ssb_dt, c("orf"))
tango_ssb_dt[, adjusted := position - peak]
tango_ssb_dt[, adjusted_random := position - random]


ggplot(data = tango_tric_dt[peak >= 70]) + xlim(-60, 5) +
  stat_summary(aes(adjusted, beta), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(adjusted, beta), fun.y = "mean", geom = "line", size = 1, color = '#1F78B4') + 
  stat_summary(aes(adjusted_random, beta), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'red') +
  stat_summary(aes(adjusted_random, beta), fun.y = "mean", geom = "line", size = 1, color = 'red')

ggplot(data = tango_ssb_dt[peak >= 70]) + xlim(-60, 5) +
  stat_summary(aes(adjusted, helix), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(adjusted, helix), fun.y = "mean", geom = "line", size = 1, color = '#1F78B4') + 
  stat_summary(aes(adjusted_random, helix), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'red') +
  stat_summary(aes(adjusted_random, helix), fun.y = "mean", geom = "line", size = 1, color = 'red')


### Psipred secondary structure with fishers
psipred_tric <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/PeptideProperties/Psipred/psipred_tric5.csv", header = F)
colnames(psipred_tric) <- c("orf", "position", "residue", "ss", "coil", "helix", "sheet")
psipred_tric <- as.data.table(psipred_tric)
setkeyv(psipred_tric, c("orf"))
temp <- tric_peaks5_max

random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(60:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

psipred_tric_dt <- psipred_tric[temp, allow.cartesian = TRUE]
setkeyv(psipred_tric_dt, c("orf"))
psipred_tric_dt[, adjusted := position - peak]
psipred_tric_dt[, adjusted_start := position - peak_start]
psipred_tric_dt[, adjusted_random := position - random]
setkeyv(psipred_tric_dt, c("position1"))
psipred_tric_dt[, coil_norm := coil / mean(coil), by = position1]
psipred_tric_dt[, helix_norm := helix / mean(helix), by = position1]
psipred_tric_dt[, sheet_norm := sheet / mean(sheet), by = position1]
setkeyv(psipred_tric_dt, c("orf"))

psipred_tric_max_dt <- psipred_tric[temp, allow.cartesian = TRUE]
setkeyv(psipred_tric_max_dt, c("orf"))
psipred_tric_max_dt[, adjusted := position - peak]
psipred_tric_max_dt[, adjusted_random := position - random]


ggplot(data = psipred_tric_dt[, .(adjusted, sheet1 = movingAverage(sheet, n=9, center=T))]) + xlim(-200, -15) + 
  stat_summary(aes(adjusted, sheet1, color = 'Sheet'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(data = psipred_tric_dt[, .(adjusted_random, helix1 = movingAverage(sheet, n=9, center=T))], aes(adjusted_random, helix1, color = 'Helix'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(data = psipred_tric_dt[, .(adjusted, coil1 = movingAverage(coil, n=9, center=T))], aes(adjusted, coil1, color = 'Coil'), fun.y = "mean", geom = "line", size = 1.25) +
  scale_color_manual(limits = c("Helix", "Sheet", "Coil"), labels = c("Helix", "Sheet", "Coil"), values = secondarystructure, name = "")
coord_cartesian(ylim = c(0.75,1.3)) +


# SSB
psipred_ssb <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/PeptideProperties/Psipred/psipred_ssb5.csv", header = F)
colnames(psipred_ssb) <- c("orf", "position", "residue", "ss", "coil", "helix", "sheet")
psipred_ssb <- as.data.table(psipred_ssb)
setkeyv(psipred_ssb, c("orf"))
temp <- ssb_peaks5_max

random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(60:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

psipred_ssb_dt <- psipred_ssb[temp, allow.cartesian = TRUE]
setkeyv(psipred_ssb_dt, c("orf"))
psipred_ssb_dt[, adjusted := position - peak]
psipred_ssb_dt[, adjusted_start := position - peak_start]
psipred_ssb_dt[, adjusted_random := position - random]
setkeyv(psipred_ssb_dt, c("position1"))
psipred_ssb_dt[, coil_norm := coil / mean(coil), by = position1]
psipred_ssb_dt[, helix_norm := helix / mean(helix), by = position1]
psipred_ssb_dt[, sheet_norm := sheet / mean(sheet), by = position1]
setkeyv(psipred_ssb_dt, c("orf"))

psipred_ssb_max_dt <- psipred_ssb[temp, allow.cartesian = TRUE]
setkeyv(psipred_ssb_max_dt, c("orf"))
psipred_ssb_max_dt[, adjusted := position - peak]
psipred_ssb_max_dt[, adjusted_random := position - random]


ggplot(data = psipred_ssb_dt[, .(adjusted, sheet1 = movingAverage(sheet, n=9, center=T))]) + xlim(-200, 5) +
  stat_summary(aes(adjusted, sheet1, color = 'Sheet'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(data = psipred_ssb_dt[, .(adjusted_random, helix1 = movingAverage(sheet, n=9, center=T))], aes(adjusted_random, helix1, color = 'Helix'), fun.y = "mean", geom = "line", size = 1.25) +
  stat_summary(data = psipred_ssb_dt[, .(adjusted, coil1 = movingAverage(coil, n=9, center=T))], aes(adjusted, coil1, color = 'Coil'), fun.y = "mean", geom = "line", size = 1.25) +
  scale_color_manual(limits = c("Helix", "Sheet", "Coil"), labels = c("Helix", "Sheet", "Coil"), values = secondarystructure, name = "")
coord_cartesian(ylim = c(0.75,1.3))

  
# Bukau
psipred_bukau <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/PeptideProperties/Psipred/psipred_bukau5.csv", header = F)
colnames(psipred_bukau) <- c("orf", "position", "residue", "ss", "coil", "helix", "sheet")
psipred_bukau <- as.data.table(psipred_bukau)
setkeyv(psipred_bukau, c("orf"))
temp <- bukau_peaks5

random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(60:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

psipred_bukau_dt <- psipred_bukau[temp, allow.cartesian = T]
setkeyv(psipred_bukau_dt, c("orf"))
psipred_bukau_dt[, adjusted := position - peak]
psipred_bukau_dt[, adjusted_start := position - peak_start]
psipred_bukau_dt[, adjusted_random := position - random]
setkeyv(psipred_bukau_dt, c("position1"))
psipred_bukau_dt[, coil_norm := coil / mean(coil), by = position1]
psipred_bukau_dt[, helix_norm := helix / mean(helix), by = position1]
psipred_bukau_dt[, sheet_norm := sheet / mean(sheet), by = position1]
setkeyv(psipred_bukau_dt, c("orf"))


### Hydrophobicity and charge
seq_tric_all <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/SequenceAnalysis/IndividualFiles/tric_peakstart_allfields.csv", header = F)
seq_tric_all <- as.data.table(seq_tric_all)
names(seq_tric_all) <- c("orf", "name", "position", "sequence")
seq_tric_all <- seq_tric_all[position >= 60]
seq_tric_all$sequence <- as.character(seq_tric_all$sequence)
seq_ssb_all <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/SequenceAnalysis/IndividualFiles/ssb_peakstart_allfields.csv", header = F)
seq_ssb_all <- as.data.table(seq_ssb_all)
names(seq_ssb_all) <- c("orf", "name", "position", "sequence")
seq_ssb_all <- seq_ssb_all[position >= 60]
seq_ssb_all$sequence <- as.character(seq_ssb_all$sequence)
seq_bukau_all <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/SequenceAnalysis/IndividualFiles/bukau_peakstart_allfields.csv", header = F)
seq_bukau_all <- as.data.table(seq_bukau_all)
names(seq_bukau_all) <- c("orf", "name", "position", "sequence")
seq_bukau_all <- seq_bukau_all[position >= 60]
seq_bukau_all$sequence <- as.character(seq_bukau_all$sequence)
seq_random_all <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/SequenceAnalysis/IndividualFiles/Random60.csv", header = F)
seq_random_all <- as.data.table(seq_random_all)
names(seq_random_all) <- c("sequence")
seq_random_all$sequence <- as.character(seq_random_all$sequence)

orfs <- NULL
names <- NULL
pos <- NULL
hydrophobic <- NULL
charge <- NULL
for (i in 1:nrow(seq_tric_all)) {
  orf <- seq_tric_all[i]$orf
  name <- seq_tric_all[i]$name
  sequence <- seq_tric_all[i]$sequence
  for (j in 1:nchar(sequence)) {
    seq <- as.character(substring(sequence,j,(j+6)))  
    hydrophobic1 <- hydrophobicity(seq, scale = "KyteDoolittle")
    hydrophobic <- c(hydrophobic, hydrophobic1)
    charge1 <- charge(seq, pH = 7, pKscale = "Lehninger")
    charge <- c(charge, charge1)
  }
  orfs1 <- as.character(rep(orf, each = 65))
  orfs <- c(orfs, orfs1)
  names1 <- as.character(rep(name, each = 65))
  names <- c(names, names1)
  pos1 <- c(1:65)
  pos <- c(pos, pos1)
}
hydro_tric <- data.table(orf = orfs, name = names, position = pos, hydro = hydrophobic,
                         charge = charge)

orfs <- NULL
names <- NULL
pos <- NULL
hydrophobic <- NULL
charge <- NULL
for (i in 1:nrow(seq_ssb_all)) {
  orf <- seq_ssb_all[i]$orf
  name <- seq_ssb_all[i]$name
  sequence <- seq_ssb_all[i]$sequence
  for (j in 1:nchar(sequence)) {
    seq <- as.character(substring(sequence,j,(j+6)))  
    hydrophobic1 <- hydrophobicity(seq, scale = "KyteDoolittle")
    hydrophobic <- c(hydrophobic, hydrophobic1)
    charge1 <- charge(seq, pH = 7, pKscale = "Lehninger")
    charge <- c(charge, charge1)
  }
  orfs1 <- as.character(rep(orf, each = 65))
  orfs <- c(orfs, orfs1)
  names1 <- as.character(rep(name, each = 65))
  names <- c(names, names1)
  pos1 <- c(1:65)
  pos <- c(pos, pos1)
}
hydro_ssb <- data.table(orf = orfs, name = names, position = pos, hydro = hydrophobic,
                        charge = charge)

orfs <- NULL
names <- NULL
pos <- NULL
hydrophobic <- NULL
charge <- NULL
for (i in 1:nrow(seq_bukau_all)) {
  orf <- seq_bukau_all[i]$orf
  name <- seq_bukau_all[i]$name
  sequence <- seq_bukau_all[i]$sequence
  for (j in 1:nchar(sequence)) {
    seq <- as.character(substring(sequence,j,(j+6)))  
    hydrophobic1 <- hydrophobicity(seq, scale = "KyteDoolittle")
    hydrophobic <- c(hydrophobic, hydrophobic1)
    charge1 <- charge(seq, pH = 7, pKscale = "Lehninger")
    charge <- c(charge, charge1)
  }
  orfs1 <- as.character(rep(orf, each = 65))
  orfs <- c(orfs, orfs1)
  names1 <- as.character(rep(name, each = 65))
  names <- c(names, names1)
  pos1 <- c(1:65)
  pos <- c(pos, pos1)
}
hydro_bukau <- data.table(orf = orfs, name = names, position = pos, hydro = hydrophobic,
                        charge = charge)

pos <- NULL
hydrophobic <- NULL
charge <- NULL
for (i in 1:nrow(seq_random_all)) {
  sequence <- seq_random_all[i]$sequence
  for (j in 1:nchar(sequence)) {
    seq <- as.character(substring(sequence,j,(j+6)))  
    hydrophobic1 <- hydrophobicity(seq, scale = "KyteDoolittle")
    hydrophobic <- c(hydrophobic, hydrophobic1)
    charge1 <- charge(seq, pH = 7, pKscale = "Lehninger")
    charge <- c(charge, charge1)
  }
  pos1 <- c(1:60)
  pos <- c(pos, pos1)
}
hydro_random <- data.table(position = pos, hydro = hydrophobic,
                         charge = charge)

ggplot(data = hydro_ssb) +
  stat_summary(aes(position, hydro), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(position, hydro), fun.y = "mean", geom = "line", size = 1, color = '#1F78B4')

ggplot(data = hydro_tric) + xlim(1,60) +
  stat_summary(aes(position, hydro), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#1F78B4') +
  stat_summary(aes(position, hydro), fun.y = "mean", geom = "line", size = 1, color = '#1F78B4') +
  stat_summary(data = hydro_ssb, aes(position, hydro), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
             fun.args=list(conf.int=0.5), fill = '#33A02C') +
  stat_summary(data = hydro_ssb, aes(position, hydro), fun.y = "mean", geom = "line", size = 1, color = '#33A02C') +
  stat_summary(data = hydro_bukau, aes(position, hydro), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#6A3D9A') +
  stat_summary(data = hydro_bukau, aes(position, hydro), fun.y = "mean", geom = "line", size = 1, color = '#6A3D9A') +
  stat_summary(data = hydro_random, aes(position, hydro), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = hydro_random, aes(position, hydro), fun.y = "mean", geom = "line", size = 1, color = 'gray50')


### Peptide properties
library(Peptides)

seq_tric_tunnel <- read.csv("/Users/KevinStein/Desktop/SequenceAnalysis/TRiCsubstrates_tunnel1.csv", header = F)
seq_tric_nascent <- read.csv("/Users/KevinStein/Desktop/SequenceAnalysis/TRiCsubstrates_nascent.csv", header = F)
seq_ssb_tunnel <- read.csv("/Users/KevinStein/Desktop/SequenceAnalysis/SSBsubstrates_tunnel1.csv", header = F)
seq_ssb_nascent <- read.csv("/Users/KevinStein/Desktop/SequenceAnalysis/SSBsubstrates_nascent.csv", header = F)
seq_random <- read.csv("/Users/KevinStein/Desktop/SequenceAnalysis/Random30mer_proteins.csv", header = F)

seq_tric_tunnel <- data.table(seq = seq_tric_tunnel$V1)
seq_tric_nascent <- data.table(seq = seq_tric_nascent$V1)
seq_ssb_tunnel <- data.table(seq = seq_ssb_tunnel$V1)
seq_ssb_nascent <- data.table(seq = seq_ssb_nascent$V1)
seq_random <- data.table(seq = seq_random$V1)

seq_tric_tunnel$seq <- as.character(seq_tric_tunnel$seq)
seq_tric_nascent$seq <- as.character(seq_tric_nascent$seq)
seq_ssb_tunnel$seq <- as.character(seq_ssb_tunnel$seq)
seq_ssb_nascent$seq <- as.character(seq_ssb_nascent$seq)
seq_random$seq <- as.character(seq_random$seq)

aliphatic1 <- NULL
boman1 <- NULL
charge1 <- NULL
hydrophobicity1 <- NULL
instability1 <- NULL
pI1 <- NULL
for (i in 1:nrow(seq_tric_tunnel)) {
  sequence <- seq_tric_tunnel[i]$seq
  aliphatic <- aIndex(sequence)
  aliphatic1 <- c(aliphatic1, aliphatic)
  boman <- boman(sequence)
  boman1 <- c(boman1, boman)
  charge <- charge(sequence, pH = 7, pKscale = "Lehninger")
  charge1 <- c(charge1, charge)
  hydrophobicity <- hydrophobicity(sequence, scale = "KyteDoolittle")
  hydrophobicity1 <- c(hydrophobicity1, hydrophobicity)
  instability <- instaIndex(sequence)
  instability1 <- c(instability1, instability)
  pI <- pI(sequence, pKscale = "EMBOSS")
  pI1 <- c(pI1, pI)
}
seq_tric_tunnel$aliphatic <- aliphatic1
seq_tric_tunnel$boman <- boman1
seq_tric_tunnel$pI <- pI1
seq_tric_tunnel$charge <- charge1
seq_tric_tunnel$instability <- instability1
seq_tric_tunnel$hydrophobicity <- hydrophobicity1


aliphatic1 <- NULL
boman1 <- NULL
charge1 <- NULL
hydrophobicity1 <- NULL
instability1 <- NULL
pI1 <- NULL
for (i in 1:nrow(seq_tric_nascent)) {
  sequence <- seq_tric_nascent[i]$seq
  aliphatic <- aIndex(sequence)
  aliphatic1 <- c(aliphatic1, aliphatic)
  boman <- boman(sequence)
  boman1 <- c(boman1, boman)
  charge <- charge(sequence, pH = 7, pKscale = "Lehninger")
  charge1 <- c(charge1, charge)
  hydrophobicity <- hydrophobicity(sequence, scale = "KyteDoolittle")
  hydrophobicity1 <- c(hydrophobicity1, hydrophobicity)
  instability <- instaIndex(sequence)
  instability1 <- c(instability1, instability)
  pI <- pI(sequence, pKscale = "EMBOSS")
  pI1 <- c(pI1, pI)
}
seq_tric_nascent$aliphatic <- aliphatic1
seq_tric_nascent$boman <- boman1
seq_tric_nascent$pI <- pI1
seq_tric_nascent$charge <- charge1
seq_tric_nascent$instability <- instability1
seq_tric_nascent$hydrophobicity <- hydrophobicity1


aliphatic1 <- NULL
boman1 <- NULL
charge1 <- NULL
hydrophobicity1 <- NULL
instability1 <- NULL
pI1 <- NULL
for (i in 1:nrow(seq_ssb_nascent)) {
  sequence <- seq_ssb_nascent[i]$seq
  aliphatic <- aIndex(sequence)
  aliphatic1 <- c(aliphatic1, aliphatic)
  boman <- boman(sequence)
  boman1 <- c(boman1, boman)
  charge <- charge(sequence, pH = 7, pKscale = "Lehninger")
  charge1 <- c(charge1, charge)
  hydrophobicity <- hydrophobicity(sequence, scale = "KyteDoolittle")
  hydrophobicity1 <- c(hydrophobicity1, hydrophobicity)
  instability <- instaIndex(sequence)
  instability1 <- c(instability1, instability)
  pI <- pI(sequence, pKscale = "EMBOSS")
  pI1 <- c(pI1, pI)
}
seq_ssb_nascent$aliphatic <- aliphatic1
seq_ssb_nascent$boman <- boman1
seq_ssb_nascent$pI <- pI1
seq_ssb_nascent$charge <- charge1
seq_ssb_nascent$instability <- instability1
seq_ssb_nascent$hydrophobicity <- hydrophobicity1


aliphatic1 <- NULL
boman1 <- NULL
charge1 <- NULL
hydrophobicity1 <- NULL
instability1 <- NULL
pI1 <- NULL
for (i in 1:nrow(seq_ssb_tunnel)) {
  sequence <- seq_ssb_tunnel[i]$seq
  aliphatic <- aIndex(sequence)
  aliphatic1 <- c(aliphatic1, aliphatic)
  boman <- boman(sequence)
  boman1 <- c(boman1, boman)
  charge <- charge(sequence, pH = 7, pKscale = "Lehninger")
  charge1 <- c(charge1, charge)
  hydrophobicity <- hydrophobicity(sequence, scale = "KyteDoolittle")
  hydrophobicity1 <- c(hydrophobicity1, hydrophobicity)
  instability <- instaIndex(sequence)
  instability1 <- c(instability1, instability)
  pI <- pI(sequence, pKscale = "EMBOSS")
  pI1 <- c(pI1, pI)
}
seq_ssb_tunnel$aliphatic <- aliphatic1
seq_ssb_tunnel$boman <- boman1
seq_ssb_tunnel$pI <- pI1
seq_ssb_tunnel$charge <- charge1
seq_ssb_tunnel$instability <- instability1
seq_ssb_tunnel$hydrophobicity <- hydrophobicity1


aliphatic1 <- NULL
boman1 <- NULL
charge1 <- NULL
hydrophobicity1 <- NULL
instability1 <- NULL
pI1 <- NULL
for (i in 1:nrow(seq_random)) {
  sequence <- seq_random[i]$seq
  aliphatic <- aIndex(sequence)
  aliphatic1 <- c(aliphatic1, aliphatic)
  boman <- boman(sequence)
  boman1 <- c(boman1, boman)
  charge <- charge(sequence, pH = 7, pKscale = "Lehninger")
  charge1 <- c(charge1, charge)
  hydrophobicity <- hydrophobicity(sequence, scale = "KyteDoolittle")
  hydrophobicity1 <- c(hydrophobicity1, hydrophobicity)
  instability <- instaIndex(sequence)
  instability1 <- c(instability1, instability)
  pI <- pI(sequence, pKscale = "EMBOSS")
  pI1 <- c(pI1, pI)
}
seq_random$aliphatic <- aliphatic1
seq_random$boman <- boman1
seq_random$pI <- pI1
seq_random$charge <- charge1
seq_random$instability <- instability1
seq_random$hydrophobicity <- hydrophobicity1



property1 <- NULL
tric_tunnel1 <- NULL
tric_nascent1 <- NULL
ssb_tunnel1 <- NULL
ssb_nascent1 <- NULL
nascent1 <- NULL
tunnel1 <- NULL
randommed1 <- NULL
tric_tunnelmed1 <- NULL
tric_nascentmed1 <- NULL
ssb_tunnelmed1 <- NULL
ssb_nascentmed1 <- NULL
for (i in 2:ncol(seq_tric_tunnel)) {
  property <- colnames(seq_tric_tunnel)[[i]] 
  filename <- paste("/Users/KevinStein/Desktop/PeptideProperties/", property, ".pdf", sep = "")  
  G <- ggplot(data = seq_random, aes("1_Random", seq_random[[i]])) + geom_violin() +
    geom_violin(data = seq_tric_nascent, aes("2_TRiCnascent", seq_tric_nascent[[i]])) +
    geom_violin(data = seq_tric_tunnel, aes("3_TRiCtunnel", seq_tric_tunnel[[i]])) +
    geom_violin(data = seq_ssb_nascent, aes("4_SSBnascent", seq_ssb_nascent[[i]])) +
    geom_violin(data = seq_ssb_tunnel, aes("5_SSBtunnel", seq_ssb_tunnel[[i]]))
  ggsave(filename, G, width = 6, height = 4, dpi = 300)
  property1 <- c(property1, property)
  tric_tunnel <- wilcox.test(na.omit(seq_tric_tunnel[[i]]), na.omit(seq_random[[i]]))$p.value
  tric_tunnel1 <- c(tric_tunnel1, tric_tunnel)
  tric_nascent <- wilcox.test(na.omit(seq_tric_nascent[[i]]), na.omit(seq_random[[i]]))$p.value
  tric_nascent1 <- c(tric_nascent1, tric_nascent)
  ssb_tunnel <- wilcox.test(na.omit(seq_ssb_tunnel[[i]]), na.omit(seq_random[[i]]))$p.value
  ssb_tunnel1 <- c(ssb_tunnel1, ssb_tunnel)
  ssb_nascent <- wilcox.test(na.omit(seq_ssb_nascent[[i]]), na.omit(seq_random[[i]]))$p.value
  ssb_nascent1 <- c(ssb_nascent1, ssb_nascent)
  nascent <- wilcox.test(na.omit(seq_tric_nascent[[i]]), na.omit(seq_ssb_nascent[[i]]))$p.value
  nascent1 <- c(nascent1, nascent)
  tunnel <- wilcox.test(na.omit(seq_tric_tunnel[[i]]), na.omit(seq_ssb_tunnel[[i]]))$p.value
  tunnel1 <- c(tunnel1, tunnel)
  randommed <- median(seq_random[[i]])
  randommed1 <- c(randommed1, randommed)
  tric_tunnelmed <- median(seq_tric_tunnel[[i]])
  tric_tunnelmed1 <- c(tric_tunnelmed1, tric_tunnelmed)
  tric_nascentmed <- median(seq_tric_nascent[[i]])
  tric_nascentmed1 <- c(tric_nascentmed1, tric_nascentmed)
  ssb_tunnelmed <- median(seq_ssb_tunnel[[i]])
  ssb_tunnelmed1 <- c(ssb_tunnelmed1, ssb_tunnelmed)
  ssb_nascentmed <- median(seq_ssb_nascent[[i]])
  ssb_nascentmed1 <- c(ssb_nascentmed1, ssb_nascentmed)
}
peptides_sig <- data.table(property = property1, tric_tunnel.p = tric_tunnel1,
                             tric_nascent.p = tric_nascent1, ssb_tunnel.p = ssb_tunnel1,
                             ssb_nascent.p = ssb_nascent1, nascent = nascent1,
                             tunnel = tunnel1, random.med = randommed1,
                             tric_tunnel.med = tric_tunnelmed1,
                             tric_nascent.med = tric_nascentmed1,
                             ssb_tunnel.med = ssb_tunnelmed1,
                             ssb_nascent.med = ssb_nascentmed1)

peptides_sig[, tric_tunnel_diff := tric_tunnel.med - random.med]
peptides_sig[, tric_nascent_diff := tric_nascent.med - random.med]
peptides_sig[, ssb_tunnel_diff := ssb_tunnel.med - random.med]
peptides_sig[, ssb_nascent_diff := ssb_nascent.med - random.med]
peptides_sig[, nascent_diff := tric_nascent.med - ssb_nascent.med]
peptides_sig[, tunnel_diff := tric_tunnel.med - ssb_tunnel.med]
View(peptides_sig)

