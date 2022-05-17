source('/Users/KevinStein/Desktop/Lab/Bioinformatics/R_scripts/MovingAverage.R')
library(data.table)

### Import raw reads of ribosome occupancy by codons
R1X <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/R1X.codon", header = F, stringsAsFactors = T)
R2X <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/R2X.codon", header = F, stringsAsFactors = T)
T1X <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/T1X.codon", header = F, stringsAsFactors = T)
T2X <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/T2X.codon", header = F, stringsAsFactors = T)
R1pur <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/R_PUR.codon", header = F, stringsAsFactors = T)
R2pur <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/Rpur.codon", header = F, stringsAsFactors = T)
T1pur <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/T_PUR.codon", header = F, stringsAsFactors = T)
T2pur <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/Tpur.codon", header = F, stringsAsFactors = T)
R1chx <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/R_CHX.codon", header = F, stringsAsFactors = T)
R2chx <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/Rchx.codon", header = F, stringsAsFactors = T)
T1chx <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/T_CHX.codon", header = F, stringsAsFactors = T)
T2chx <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/Tchx.codon", header = F, stringsAsFactors = T)
SSB_Rchx1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/SSB/Rchx1.codon", header = F, stringsAsFactors = T)
SSB_Rchx2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/SSB/Rchx2.codon", header = F, stringsAsFactors = T)
SSB_Schx1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/SSB/Schx1.codon", header = F, stringsAsFactors = T)
SSB_Schx2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/SSB/Schx2.codon", header = F, stringsAsFactors = T)
SSB_Rpuro1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/SSB/Rpuro1.codon", header = F, stringsAsFactors = T)
SSB_Rpuro2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/SSB/Rpuro2.codon", header = F, stringsAsFactors = T)
SSB_Spuro1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/SSB/Spuro1.codon", header = F, stringsAsFactors = T)
SSB_Spuro2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/SSB/Spuro2.codon", header = F, stringsAsFactors = T)

B_ssb1_I1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Doring_Bukau_2017/Codon_density/Ssb1interactome1.codon", header = F, stringsAsFactors = T)
B_ssb1_I2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Doring_Bukau_2017/Codon_density/Ssb1interactome2.codon", header = F, stringsAsFactors = T)
B_ssb1_T1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Doring_Bukau_2017/Codon_density/Ssb1translatome1.codon", header = F, stringsAsFactors = T)
B_ssb1_T2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Doring_Bukau_2017/Codon_density/Ssb1translatome2.codon", header = F, stringsAsFactors = T)
B_ssb2_I1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Doring_Bukau_2017/Codon_density/Ssb2interactome1.codon", header = F, stringsAsFactors = T)
B_ssb2_I2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Doring_Bukau_2017/Codon_density/Ssb2interactome2.codon", header = F, stringsAsFactors = T)
B_ssb2_T1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Doring_Bukau_2017/Codon_density/Ssb2translatome1.codon", header = F, stringsAsFactors = T)
B_ssb2_T2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Doring_Bukau_2017/Codon_density/Ssb2translatome2.codon", header = F, stringsAsFactors = T)

B_WT1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Doring_Bukau_2017/Codon_density/WTtranslatome1.codon", header = F, stringsAsFactors = T)
B_WT2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Doring_Bukau_2017/Codon_density/WTtranslatome2.codon", header = F, stringsAsFactors = T)
WT1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/KCS03_tRNA/Codon_density/WT1.codon", header = F, stringsAsFactors = T)
WT2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/KCS03_tRNA/Codon_density/WT2.codon",header = F,stringsAsFactors = T)
EV1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/KCS03_tRNA/Codon_density/EV1.codon", header = F, stringsAsFactors = T)
EV2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/KCS03_tRNA/Codon_density/EV2.codon",header = F,stringsAsFactors = T)
G_WT1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Young_Green_2015/Codon_density/WT1_A.codon", header = F, stringsAsFactors = T)
G_WT2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Young_Green_2015/Codon_density/WT2_A.codon", header = F, stringsAsFactors = T)
C_WT1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Radhakrishnan_Coller_2017/Rep1.codon", header = F, stringsAsFactors = T)
C_WT2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Radhakrishnan_Coller_2017/Rep2.codon", header = F, stringsAsFactors = T)
L_WT1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Nedialkova_Leidel_2015/Rep1.codon", header = F, stringsAsFactors = T)
L_WT2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Nedialkova_Leidel_2015/Rep2.codon", header = F, stringsAsFactors = T)
L_WT3 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/Published/Nedialkova_Leidel_2015/Rep3.codon", header = F, stringsAsFactors = T)

Atp2_R1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/Atp2R1.codon", header = F, stringsAsFactors = T)
Atp2_R2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/Atp2R2.codon", header = F, stringsAsFactors = T)
Atp2_T1 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/Atp2T1.codon", header = F, stringsAsFactors = T)
Atp2_T2 <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/Atp2T2.codon", header = F, stringsAsFactors = T)


### Create data tables
tric_dt <- data.table(orf = R1X[, 1], 
                    position = R1X[, 2],
                    codon = R1X[, 3],
                    R1X = R1X[, 4],
                    T1X = T1X[, 4],
                    R2X = R2X[, 4],
                    T2X = T2X[, 4],
                    R1pur = R1pur[, 4],
                    R2pur = R2pur[, 4],
                    T1pur = T1pur[, 4],
                    T2pur = T2pur[, 4],
                    R1chx = R1chx[, 4],
                    R2chx = R2chx[, 4],
                    T1chx = T1chx[, 4],
                    T2chx = T2chx[, 4],
                    SSB_Rchx1 = SSB_Rchx1[, 4],
                    SSB_Rchx2 = SSB_Rchx2[, 4],
                    SSB_Schx1 = SSB_Schx1[, 4],
                    SSB_Schx2 = SSB_Schx2[, 4],
                    SSB_Rpuro1 = SSB_Rpuro1[, 4],
                    SSB_Rpuro2 = SSB_Rpuro2[, 4],
                    SSB_Spuro1 = SSB_Spuro1[, 4],
                    SSB_Spuro2 = SSB_Spuro2[, 4])
setkeyv(tric_dt, c("orf"))
tric_dt[, length := length(position), by = orf]
tric_dt[, length := length - 15] # substract flanking 7 codons and stop codon
tric_dt[, stopdist := position - (length + 1)] # stop codon is 0
tric_dt <- tric_dt[orf != "YOR031W" & orf != "YFL057C" & orf != "YOL013W-A"] # overlapping or blocked orf
tric_dt[, position1 := as.character(base::paste(orf, position, sep = "_"))]

bukau_dt <- data.table(orf = B_ssb1_T1[, 1], 
                      position = B_ssb1_T1[, 2],
                      codon = B_ssb1_T1[, 3],
                      B_ssb1_T1 = B_ssb1_T1[, 4],
                      B_ssb1_T2 = B_ssb1_T2[, 4],
                      B_ssb1_I1 = B_ssb1_I1[, 4],
                      B_ssb1_I2 = B_ssb1_I2[, 4],
                      B_ssb2_T1 = B_ssb2_T1[, 4],
                      B_ssb2_T2 = B_ssb2_T2[, 4],
                      B_ssb2_I1 = B_ssb2_I1[, 4],
                      B_ssb2_I2 = B_ssb2_I2[, 4])
setkeyv(bukau_dt, c("orf"))
bukau_dt[, length := length(position), by = orf]
bukau_dt[, length := length - 15] # substract flanking 7 codons and stop codon
bukau_dt[, stopdist := position - (length + 1)] # stop codon is 0
bukau_dt <- bukau_dt[orf != "YOR031W" & orf != "YFL057C" & orf != "YOL013W-A"]
bukau_dt[, position1 := as.character(base::paste(orf, position, sep = "_"))]


translatome_dt <- data.table(orf = WT1[, 1], 
                             position = WT1[, 2],
                             codon = WT1[, 3],
                             WT1 = WT1[, 4],
                             WT2 = WT2[, 4],
                             EV1 = EV1[, 4],
                             EV2 = EV2[, 4],
                             B_WT1 = B_WT1[, 4],
                             B_WT2 = B_WT2[, 4],
                             G_WT1 = G_WT1[, 4],
                             G_WT2 = G_WT2[, 4],
                             C_WT1 = C_WT1[, 4],
                             C_WT2 = C_WT2[, 4],
                             L_WT1 = L_WT1[, 4],
                             L_WT2 = L_WT2[, 4],
                             L_WT3 = L_WT3[, 4])
setkeyv(translatome_dt, c("orf"))
translatome_dt[, length := length(position), by = orf]
translatome_dt[, length := length - 15] # substract flanking 7 codons and stop codon
translatome_dt[, stopdist := position - (length + 1)] # stop codon is 0
translatome_dt <- translatome_dt[orf != "YOR031W" & orf != "YFL057C" & orf != "YOL013W-A"]
translatome_dt[, position1 := as.character(base::paste(orf, position, sep = "_"))]

atp2_dt <- data.table(orf = Atp2_R1[, 1], 
                      position = Atp2_R1[, 2],
                      codon = Atp2_R1[, 3],
                      Atp2_R1 = Atp2_R1[, 4],
                      Atp2_R2 = Atp2_R2[, 4],
                      Atp2_T1 = Atp2_T1[, 4],
                      Atp2_T2 = Atp2_T2[, 4])
setkeyv(atp2_dt, c("orf"))
atp2_dt[, length := length(position), by = orf]
atp2_dt[, length := length - 15] # substract flanking 7 codons and stop codon
atp2_dt[, stopdist := position - (length + 1)] # stop codon is 0
atp2_dt <- atp2_dt[orf != "YOR031W" & orf != "YFL057C" & orf != "YOL013W-A"]
atp2_dt[, position1 := as.character(base::paste(orf, position, sep = "_"))]

### Add residue
i <- cbind(match(atp2_dt$codon, codon_table$codon))
atp2_dt <- cbind(atp2_dt, residue = codon_table[i]$residue)


### Calculate sum and rpm for each position, and rpc for each orf
tric_dt[, ribo_sum := R1X + R2X]
tric_dt[, tric_sum := T1X + T2X]
tric_dt[, Rpur_sum := R1pur + R2pur]
tric_dt[, Tpur_sum := T1pur + T2pur]
tric_dt[, Rchx_sum := R1chx + R2chx]
tric_dt[, Tchx_sum := T1chx + T2chx]
tric_dt[, ssb_Rchx_sum := SSB_Rchx1 + SSB_Rchx2]
tric_dt[, ssb_Schx_sum := SSB_Schx1 + SSB_Schx2]
tric_dt[, ssb_Rpur_sum := SSB_Rpuro1 + SSB_Rpuro2]
tric_dt[, ssb_Spur_sum := SSB_Spuro1 + SSB_Spuro2]

bukau_dt[, B_ssb1_T_sum := B_ssb1_T1 + B_ssb1_T2]
bukau_dt[, B_ssb1_I_sum := B_ssb1_I1 + B_ssb1_I2]
bukau_dt[, B_ssb2_T_sum := B_ssb2_T1 + B_ssb2_T2]
bukau_dt[, B_ssb2_I_sum := B_ssb2_I1 + B_ssb2_I2]
translatome_dt[, WT_sum := WT1 + WT2]
translatome_dt[, EV_sum := EV1 + EV2]
translatome_dt[, G_WT_sum := G_WT1 + G_WT2]
translatome_dt[, B_WT_sum := B_WT1 + B_WT2]
translatome_dt[, C_WT_sum := C_WT1 + C_WT2]
translatome_dt[, L_WT_sum := L_WT1 + L_WT2 + L_WT3]
atp2_dt[, ribo_sum := Atp2_R1 + Atp2_R2]
atp2_dt[, tric_sum := Atp2_T1 + Atp2_T2]

# Continue after expression analysis and calculation of rpc below
i <- cbind(match(tric_dt$orf, expression1$orf))
tric_dt <- cbind(tric_dt, ribo_total = expression1[i]$ribo_total)
tric_dt <- cbind(tric_dt, tric_total = expression1[i]$tric_total)
tric_dt <- cbind(tric_dt, Rpur_total = expression1[i]$Rpur_total)
tric_dt <- cbind(tric_dt, Tpur_total = expression1[i]$Tpur_total)
tric_dt <- cbind(tric_dt, ribo_rpc = expression1[i]$ribo_rpc)
tric_dt <- cbind(tric_dt, tric_rpc = expression1[i]$tric_rpc)
tric_dt <- cbind(tric_dt, Rpur_rpc = expression1[i]$Rpur_rpc)
tric_dt <- cbind(tric_dt, Tpur_rpc = expression1[i]$Tpur_rpc)
tric_dt <- cbind(tric_dt, Rchx_rpc = expression1[i]$Rchx_rpc)
tric_dt <- cbind(tric_dt, Tchx_rpc = expression1[i]$Tchx_rpc)
tric_dt <- cbind(tric_dt, R1X_rpc = expression1[i]$R1X_rpc)
tric_dt <- cbind(tric_dt, R2X_rpc = expression1[i]$R2X_rpc)
tric_dt <- cbind(tric_dt, T1X_rpc = expression1[i]$T1X_rpc)
tric_dt <- cbind(tric_dt, T2X_rpc = expression1[i]$T2X_rpc)
tric_dt <- cbind(tric_dt, R1pur_rpc = expression1[i]$R1pur_rpc)
tric_dt <- cbind(tric_dt, R2pur_rpc = expression1[i]$R2pur_rpc)
tric_dt <- cbind(tric_dt, T1pur_rpc = expression1[i]$T1pur_rpc)
tric_dt <- cbind(tric_dt, T2pur_rpc = expression1[i]$T2pur_rpc)
tric_dt <- cbind(tric_dt, R1chx_rpc = expression1[i]$R1chx_rpc)
tric_dt <- cbind(tric_dt, R2chx_rpc = expression1[i]$R2chx_rpc)
tric_dt <- cbind(tric_dt, T1chx_rpc = expression1[i]$T1chx_rpc)
tric_dt <- cbind(tric_dt, T2chx_rpc = expression1[i]$T2chx_rpc)
tric_dt <- cbind(tric_dt, ribo_tpm = expression1[i]$ribo_tpm)
tric_dt <- cbind(tric_dt, tric_tpm = expression1[i]$tric_tpm)
tric_dt <- cbind(tric_dt, Rpur_tpm = expression1[i]$Rpur_tpm)
tric_dt <- cbind(tric_dt, Tpur_tpm = expression1[i]$Tpur_tpm)
tric_dt <- cbind(tric_dt, R1X_tpm = expression1[i]$R1X_tpm)
tric_dt <- cbind(tric_dt, R2X_tpm = expression1[i]$R2X_tpm)
tric_dt <- cbind(tric_dt, T1X_tpm = expression1[i]$T1X_tpm)
tric_dt <- cbind(tric_dt, T2X_tpm = expression1[i]$T2X_tpm)
tric_dt <- cbind(tric_dt, R1pur_tpm = expression1[i]$R1pur_tpm)
tric_dt <- cbind(tric_dt, R2pur_tpm = expression1[i]$R2pur_tpm)
tric_dt <- cbind(tric_dt, T1pur_tpm = expression1[i]$T1pur_tpm)
tric_dt <- cbind(tric_dt, T2pur_tpm = expression1[i]$T2pur_tpm)
tric_dt <- cbind(tric_dt, ribo_rpkm = expression1[i]$ribo_rpkm)
tric_dt <- cbind(tric_dt, tric_rpkm = expression1[i]$tric_rpkm)
tric_dt <- cbind(tric_dt, R1X_rpkm = expression1[i]$R1X_rpkm)
tric_dt <- cbind(tric_dt, R2X_rpkm = expression1[i]$R2X_rpkm)
tric_dt <- cbind(tric_dt, T1X_rpkm = expression1[i]$T1X_rpkm)
tric_dt <- cbind(tric_dt, T2X_rpkm = expression1[i]$T2X_rpkm)
tric_dt <- cbind(tric_dt, ssb_Rchx_total = expression1[i]$ssb_Rchx_total)
tric_dt <- cbind(tric_dt, ssb_Schx_total = expression1[i]$ssb_Schx_total)
tric_dt <- cbind(tric_dt, ssb_Rpur_total = expression1[i]$ssb_Rpur_total)
tric_dt <- cbind(tric_dt, ssb_Spur_total = expression1[i]$ssb_Spur_total)
tric_dt <- cbind(tric_dt, ssb_Rchx_rpc = expression1[i]$ssb_Rchx_rpc)
tric_dt <- cbind(tric_dt, ssb_Schx_rpc = expression1[i]$ssb_Schx_rpc)
tric_dt <- cbind(tric_dt, ssb_Rpur_rpc = expression1[i]$ssb_Rpur_rpc)
tric_dt <- cbind(tric_dt, ssb_Spur_rpc = expression1[i]$ssb_Spur_rpc)
tric_dt <- cbind(tric_dt, ssb_Rchx1_rpc = expression1[i]$ssb_Rchx1_rpc)
tric_dt <- cbind(tric_dt, ssb_Rchx2_rpc = expression1[i]$ssb_Rchx2_rpc)
tric_dt <- cbind(tric_dt, ssb_Schx1_rpc = expression1[i]$ssb_Schx1_rpc)
tric_dt <- cbind(tric_dt, ssb_Schx2_rpc = expression1[i]$ssb_Schx2_rpc)
tric_dt <- cbind(tric_dt, ssb_Rpur1_rpc = expression1[i]$ssb_Rpur1_rpc)
tric_dt <- cbind(tric_dt, ssb_Rpur2_rpc = expression1[i]$ssb_Rpur2_rpc)
tric_dt <- cbind(tric_dt, ssb_Spur1_rpc = expression1[i]$ssb_Spur1_rpc)
tric_dt <- cbind(tric_dt, ssb_Spur2_rpc = expression1[i]$ssb_Spur2_rpc)
tric_dt <- cbind(tric_dt, ssb_Rchx_tpm = expression1[i]$ssb_Rchx_tpm)
tric_dt <- cbind(tric_dt, ssb_Schx_tpm = expression1[i]$ssb_Schx_tpm)
tric_dt <- cbind(tric_dt, ssb_Rpur_tpm = expression1[i]$ssb_Rpur_tpm)
tric_dt <- cbind(tric_dt, ssb_Spur_tpm = expression1[i]$ssb_Spur_tpm)
tric_dt <- cbind(tric_dt, ssb_Rchx1_tpm = expression1[i]$ssb_Rchx1_tpm)
tric_dt <- cbind(tric_dt, ssb_Rchx2_tpm = expression1[i]$ssb_Rchx2_tpm)
tric_dt <- cbind(tric_dt, ssb_Schx1_tpm = expression1[i]$ssb_Schx1_tpm)
tric_dt <- cbind(tric_dt, ssb_Schx2_tpm = expression1[i]$ssb_Schx2_tpm)
tric_dt <- cbind(tric_dt, ssb_Rpur1_tpm = expression1[i]$ssb_Rpur1_tpm)
tric_dt <- cbind(tric_dt, ssb_Rpur2_tpm = expression1[i]$ssb_Rpur2_tpm)
tric_dt <- cbind(tric_dt, ssb_Spur1_tpm = expression1[i]$ssb_Spur1_tpm)
tric_dt <- cbind(tric_dt, ssb_Spur2_tpm = expression1[i]$ssb_Spur2_tpm)
tric_dt <- cbind(tric_dt, ssb_Rchx_rpkm = expression1[i]$ssb_Rchx_rpkm)
tric_dt <- cbind(tric_dt, ssb_Schx_rpkm = expression1[i]$ssb_Schx_rpkm)
tric_dt <- cbind(tric_dt, ssb_Rchx1_rpkm = expression1[i]$ssb_Rchx1_rpkm)
tric_dt <- cbind(tric_dt, ssb_Rchx2_rpkm = expression1[i]$ssb_Rchx2_rpkm)
tric_dt <- cbind(tric_dt, ssb_Schx1_rpkm = expression1[i]$ssb_Schx1_rpkm)
tric_dt <- cbind(tric_dt, ssb_Schx2_rpkm = expression1[i]$ssb_Schx2_rpkm)
i <- cbind(match(tric_dt$orf, test$orf))
tric_dt <- cbind(tric_dt, ribo_cor = test[i]$ribo_cor)
tric_dt <- cbind(tric_dt, tric_cor = test[i]$tric_cor)
tric_dt <- cbind(tric_dt, ssb_Rchx_cor = test[i]$ssb_Rchx_cor)
tric_dt <- cbind(tric_dt, ssb_Schx_cor = test[i]$ssb_Schx_cor)


i <- cbind(match(bukau_dt$orf, expression1Bukau$orf))
bukau_dt <- cbind(bukau_dt, B_ssb1_T_total = expression1Bukau[i]$B_ssb1_T_total)
bukau_dt <- cbind(bukau_dt, B_ssb1_I_total = expression1Bukau[i]$B_ssb1_I_total)
bukau_dt <- cbind(bukau_dt, B_ssb1_T_rpc = expression1Bukau[i]$B_ssb1_T_rpc)
bukau_dt <- cbind(bukau_dt, B_ssb1_I_rpc = expression1Bukau[i]$B_ssb1_I_rpc)
bukau_dt <- cbind(bukau_dt, B_ssb1_T1_rpc = expression1Bukau[i]$B_ssb1_T1_rpc)
bukau_dt <- cbind(bukau_dt, B_ssb1_T2_rpc = expression1Bukau[i]$B_ssb1_T2_rpc)
bukau_dt <- cbind(bukau_dt, B_ssb1_I1_rpc = expression1Bukau[i]$B_ssb1_I1_rpc)
bukau_dt <- cbind(bukau_dt, B_ssb1_I2_rpc = expression1Bukau[i]$B_ssb1_I2_rpc)
bukau_dt <- cbind(bukau_dt, B_ssb1_T_tpm = expression1Bukau[i]$B_ssb1_T_tpm)
bukau_dt <- cbind(bukau_dt, B_ssb1_I_tpm = expression1Bukau[i]$B_ssb1_I_tpm)
bukau_dt <- cbind(bukau_dt, B_ssb1_T1_tpm = expression1Bukau[i]$B_ssb1_T1_tpm)
bukau_dt <- cbind(bukau_dt, B_ssb1_T2_tpm = expression1Bukau[i]$B_ssb1_T2_tpm)
bukau_dt <- cbind(bukau_dt, B_ssb1_I1_tpm = expression1Bukau[i]$B_ssb1_I1_tpm)
bukau_dt <- cbind(bukau_dt, B_ssb1_I2_tpm = expression1Bukau[i]$B_ssb1_I2_tpm)
bukau_dt <- cbind(bukau_dt, B_ssb1_T_rpkm = expression1Bukau[i]$B_ssb1_T_rpkm)
bukau_dt <- cbind(bukau_dt, B_ssb1_I_rpkm = expression1Bukau[i]$B_ssb1_I_rpkm)
bukau_dt <- cbind(bukau_dt, B_ssb1_T1_rpkm = expression1Bukau[i]$B_ssb1_T1_rpkm)
bukau_dt <- cbind(bukau_dt, B_ssb1_T2_rpkm = expression1Bukau[i]$B_ssb1_T2_rpkm)
bukau_dt <- cbind(bukau_dt, B_ssb1_I1_rpkm = expression1Bukau[i]$B_ssb1_I1_rpkm)
bukau_dt <- cbind(bukau_dt, B_ssb1_I2_rpkm = expression1Bukau[i]$B_ssb1_I2_rpkm)
i <- cbind(match(bukau_dt$orf, test$orf))
bukau_dt <- cbind(bukau_dt, B_ssb1_T_cor = test[i]$B_ssb1_T_cor)
bukau_dt <- cbind(bukau_dt, B_ssb1_I_cor = test[i]$B_ssb1_I_cor)


i <- cbind(match(translatome_dt$orf, expression1Translatome$orf))
translatome_dt <- cbind(translatome_dt, WT_rpc = expression1Translatome[i]$WT_rpc)
translatome_dt <- cbind(translatome_dt, EV_rpc = expression1Translatome[i]$EV_rpc)
translatome_dt <- cbind(translatome_dt, G_WT_rpc = expression1Translatome[i]$G_WT_rpc)
translatome_dt <- cbind(translatome_dt, B_WT_rpc = expression1Translatome[i]$B_WT_rpc)
translatome_dt <- cbind(translatome_dt, C_WT_rpc = expression1Translatome[i]$C_WT_rpc)
translatome_dt <- cbind(translatome_dt, L_WT_rpc = expression1Translatome[i]$L_WT_rpc)
translatome_dt <- cbind(translatome_dt, WT1_rpc = expression1Translatome[i]$WT1_rpc)
translatome_dt <- cbind(translatome_dt, WT2_rpc = expression1Translatome[i]$WT2_rpc)
translatome_dt <- cbind(translatome_dt, EV1_rpc = expression1Translatome[i]$EV1_rpc)
translatome_dt <- cbind(translatome_dt, EV2_rpc = expression1Translatome[i]$EV2_rpc)
translatome_dt <- cbind(translatome_dt, G_WT1_rpc = expression1Translatome[i]$G_WT1_rpc)
translatome_dt <- cbind(translatome_dt, G_WT2_rpc = expression1Translatome[i]$G_WT2_rpc)
translatome_dt <- cbind(translatome_dt, B_WT1_rpc = expression1Translatome[i]$B_WT1_rpc)
translatome_dt <- cbind(translatome_dt, B_WT2_rpc = expression1Translatome[i]$B_WT2_rpc)
translatome_dt <- cbind(translatome_dt, C_WT1_rpc = expression1Translatome[i]$C_WT1_rpc)
translatome_dt <- cbind(translatome_dt, C_WT2_rpc = expression1Translatome[i]$C_WT2_rpc)
translatome_dt <- cbind(translatome_dt, L_WT1_rpc = expression1Translatome[i]$L_WT1_rpc)
translatome_dt <- cbind(translatome_dt, L_WT2_rpc = expression1Translatome[i]$L_WT2_rpc)
translatome_dt <- cbind(translatome_dt, L_WT3_rpc = expression1Translatome[i]$L_WT3_rpc)


i <- cbind(match(atp2_dt$orf, expression1_atp2$orf))
atp2_dt <- cbind(atp2_dt, ribo_rpc = expression1_atp2[i]$ribo_rpc)
atp2_dt <- cbind(atp2_dt, ribo_total = expression1_atp2[i]$ribo_total)
atp2_dt <- cbind(atp2_dt, tric_rpc = expression1_atp2[i]$tric_rpc)
atp2_dt <- cbind(atp2_dt, tric_total = expression1_atp2[i]$tric_total)
atp2_dt <- cbind(atp2_dt, R1_rpc = expression1_atp2[i]$R1_rpc)
atp2_dt <- cbind(atp2_dt, R2_rpc = expression1_atp2[i]$R2_rpc)
atp2_dt <- cbind(atp2_dt, T1_rpc = expression1_atp2[i]$T1_rpc)
atp2_dt <- cbind(atp2_dt, T2_rpc = expression1_atp2[i]$T2_rpc)

# Calculate rpm
tric_dt[, ribo_rpm := (ribo_sum / sum(tric_dt$ribo_sum)) * 10^6]
tric_dt[, tric_rpm := (tric_sum / sum(tric_dt$tric_sum)) * 10^6]
tric_dt[, Rpur_rpm := (Rpur_sum / sum(tric_dt$Rpur_sum)) * 10^6]
tric_dt[, Tpur_rpm := (Tpur_sum / sum(tric_dt$Tpur_sum)) * 10^6]
tric_dt[, Rchx_rpm := (Rchx_sum / sum(tric_dt$Rchx_sum)) * 10^6]
tric_dt[, Tchx_rpm := (Tchx_sum / sum(tric_dt$Tchx_sum)) * 10^6]
tric_dt[, ssb_Rchx_rpm := (ssb_Rchx_sum / sum(tric_dt$ssb_Rchx_sum)) * 10^6]
tric_dt[, ssb_Schx_rpm := (ssb_Schx_sum / sum(tric_dt$ssb_Schx_sum)) * 10^6]
tric_dt[, ssb_Rpur_rpm := (ssb_Rpur_sum / sum(tric_dt$ssb_Rpur_sum)) * 10^6]
tric_dt[, ssb_Spur_rpm := (ssb_Spur_sum / sum(tric_dt$ssb_Spur_sum)) * 10^6]
tric_dt[, R1X_rpm := (R1X / sum(tric_dt$R1X)) * 10^6]  
tric_dt[, T1X_rpm := (T1X / sum(tric_dt$T1X)) * 10^6]  
tric_dt[, R2X_rpm := (R2X / sum(tric_dt$R2X)) * 10^6]  
tric_dt[, T2X_rpm := (T2X / sum(tric_dt$T2X)) * 10^6]
tric_dt[, R1pur_rpm := (R1pur / sum(tric_dt$R1pur)) * 10^6]  
tric_dt[, R2pur_rpm := (R2pur / sum(tric_dt$R2pur)) * 10^6]  
tric_dt[, T1pur_rpm := (T1pur / sum(tric_dt$T1pur)) * 10^6]  
tric_dt[, T2pur_rpm := (T2pur / sum(tric_dt$T2pur)) * 10^6]  
tric_dt[, R1chx_rpm := (R1chx / sum(tric_dt$R1chx)) * 10^6]  
tric_dt[, R2chx_rpm := (R2chx / sum(tric_dt$R2chx)) * 10^6]  
tric_dt[, T1chx_rpm := (T1chx / sum(tric_dt$T1chx)) * 10^6]  
tric_dt[, T2chx_rpm := (T2chx / sum(tric_dt$T2chx)) * 10^6]
tric_dt[, ssb_Rchx1_rpm := (SSB_Rchx1 / sum(tric_dt$SSB_Rchx1)) * 10^6]
tric_dt[, ssb_Rchx2_rpm := (SSB_Rchx2 / sum(tric_dt$SSB_Rchx2)) * 10^6]
tric_dt[, ssb_Schx1_rpm := (SSB_Schx1 / sum(tric_dt$SSB_Schx1)) * 10^6]
tric_dt[, ssb_Schx2_rpm := (SSB_Schx2 / sum(tric_dt$SSB_Schx2)) * 10^6]
tric_dt[, ssb_Rpur1_rpm := (SSB_Rpuro1 / sum(tric_dt$SSB_Rpuro1)) * 10^6]
tric_dt[, ssb_Rpur2_rpm := (SSB_Rpuro2 / sum(tric_dt$SSB_Rpuro2)) * 10^6]
tric_dt[, ssb_Spur1_rpm := (SSB_Spuro1 / sum(tric_dt$SSB_Spuro1)) * 10^6]
tric_dt[, ssb_Spur2_rpm := (SSB_Spuro2 / sum(tric_dt$SSB_Spuro2)) * 10^6]
bukau_dt[, B_ssb1_T_rpm := (B_ssb1_T_sum / sum(bukau_dt$B_ssb1_T_sum)) * 10^6]
bukau_dt[, B_ssb1_I_rpm := (B_ssb1_I_sum / sum(bukau_dt$B_ssb1_I_sum)) * 10^6]
bukau_dt[, B_ssb2_T_rpm := (B_ssb2_T_sum / sum(bukau_dt$B_ssb2_T_sum)) * 10^6]
bukau_dt[, B_ssb2_I_rpm := (B_ssb2_I_sum / sum(bukau_dt$B_ssb2_I_sum)) * 10^6]
bukau_dt[, B_ssb1_T1_rpm := (B_ssb1_T1 / sum(bukau_dt$B_ssb1_T1)) * 10^6]
bukau_dt[, B_ssb1_T2_rpm := (B_ssb1_T2 / sum(bukau_dt$B_ssb1_T2)) * 10^6]
bukau_dt[, B_ssb1_I1_rpm := (B_ssb1_I1 / sum(bukau_dt$B_ssb1_I1)) * 10^6]
bukau_dt[, B_ssb1_I2_rpm := (B_ssb1_I2 / sum(bukau_dt$B_ssb1_I2)) * 10^6]
bukau_dt[, B_ssb2_T1_rpm := (B_ssb2_T1 / sum(bukau_dt$B_ssb2_T1)) * 10^6]
bukau_dt[, B_ssb2_T2_rpm := (B_ssb2_T2 / sum(bukau_dt$B_ssb2_T2)) * 10^6]
bukau_dt[, B_ssb2_I1_rpm := (B_ssb2_I1 / sum(bukau_dt$B_ssb2_I1)) * 10^6]
bukau_dt[, B_ssb2_I2_rpm := (B_ssb2_I2 / sum(bukau_dt$B_ssb2_I2)) * 10^6]
translatome_dt[, WT_rpm := (WT_sum / sum(translatome_dt$WT_sum)) * 10^6]
translatome_dt[, EV_rpm := (EV_sum / sum(translatome_dt$EV_sum)) * 10^6]
translatome_dt[, G_WT_rpm := (G_WT_sum / sum(translatome_dt$G_WT_sum)) * 10^6]
translatome_dt[, B_WT_rpm := (B_WT_sum / sum(translatome_dt$B_WT_sum)) * 10^6]
translatome_dt[, C_WT_rpm := (C_WT_sum / sum(translatome_dt$C_WT_sum)) * 10^6]
translatome_dt[, L_WT_rpm := (L_WT_sum / sum(translatome_dt$L_WT_sum)) * 10^6]
translatome_dt[, WT1_rpm := (WT1 / sum(translatome_dt$WT1)) * 10^6]
translatome_dt[, WT2_rpm := (WT2 / sum(translatome_dt$WT2)) * 10^6]
translatome_dt[, EV1_rpm := (EV1 / sum(translatome_dt$EV1)) * 10^6]
translatome_dt[, EV2_rpm := (EV2 / sum(translatome_dt$EV2)) * 10^6]
translatome_dt[, G_WT1_rpm := (G_WT1 / sum(translatome_dt$G_WT1)) * 10^6]
translatome_dt[, G_WT2_rpm := (G_WT2 / sum(translatome_dt$G_WT2)) * 10^6]
translatome_dt[, B_WT1_rpm := (B_WT1 / sum(translatome_dt$B_WT1)) * 10^6]
translatome_dt[, B_WT2_rpm := (B_WT2 / sum(translatome_dt$B_WT2)) * 10^6]
translatome_dt[, C_WT1_rpm := (C_WT1 / sum(translatome_dt$C_WT1)) * 10^6]
translatome_dt[, C_WT2_rpm := (C_WT2 / sum(translatome_dt$C_WT2)) * 10^6]
translatome_dt[, L_WT1_rpm := (L_WT1 / sum(translatome_dt$L_WT1)) * 10^6]
translatome_dt[, L_WT2_rpm := (L_WT2 / sum(translatome_dt$L_WT2)) * 10^6]
translatome_dt[, L_WT3_rpm := (L_WT3 / sum(translatome_dt$L_WT3)) * 10^6]
atp2_dt[, ribo_rpm := (ribo_sum / sum(atp2_dt$ribo_sum)) * 10^6]
atp2_dt[, tric_rpm := (tric_sum / sum(atp2_dt$tric_sum)) * 10^6]
atp2_dt[, R1_rpm := (Atp2_R1 / sum(atp2_dt$Atp2_R1)) * 10^6]
atp2_dt[, R2_rpm := (Atp2_R2 / sum(atp2_dt$Atp2_R2)) * 10^6]
atp2_dt[, T1_rpm := (Atp2_T1 / sum(atp2_dt$Atp2_T1)) * 10^6]
atp2_dt[, T2_rpm := (Atp2_T2 / sum(atp2_dt$Atp2_T2)) * 10^6]

# Calculate normalized ribosome occupancy
tric_dt[, ribo_norm := ribo_sum / ribo_rpc]
tric_dt[, tric_norm := tric_sum / tric_rpc]
tric_dt[, Rpur_norm := Rpur_sum / Rpur_rpc]
tric_dt[, Tpur_norm := Tpur_sum / Tpur_rpc]
tric_dt[, Rchx_norm := Rchx_sum / Rchx_rpc]
tric_dt[, Tchx_norm := Tchx_sum / Tchx_rpc]
tric_dt[, ssb_Rchx_norm := ssb_Rchx_sum / ssb_Rchx_rpc]
tric_dt[, ssb_Schx_norm := ssb_Schx_sum / ssb_Schx_rpc]
tric_dt[, ssb_Rpur_norm := ssb_Rpur_sum / ssb_Rpur_rpc]
tric_dt[, ssb_Spur_norm := ssb_Spur_sum / ssb_Spur_rpc]
tric_dt[, R1X_norm := R1X / R1X_rpc]
tric_dt[, R2X_norm := R2X / R2X_rpc]
tric_dt[, T1X_norm := T1X / T1X_rpc]
tric_dt[, T2X_norm := T2X / T2X_rpc]
tric_dt[, R1pur_norm := R1pur / R1pur_rpc]
tric_dt[, R2pur_norm := R2pur / R2pur_rpc]
tric_dt[, T1pur_norm := T1pur / T1pur_rpc]
tric_dt[, T2pur_norm := T2pur / T2pur_rpc]
tric_dt[, ssb_Rchx1_norm := SSB_Rchx1 / ssb_Rchx1_rpc]
tric_dt[, ssb_Rchx2_norm := SSB_Rchx2 / ssb_Rchx2_rpc]
tric_dt[, ssb_Schx1_norm := SSB_Schx1 / ssb_Schx1_rpc]
tric_dt[, ssb_Schx2_norm := SSB_Schx2 / ssb_Schx2_rpc]
bukau_dt[, B_ssb1_T_norm := B_ssb1_T_sum / B_ssb1_T_rpc]
bukau_dt[, B_ssb1_I_norm := B_ssb1_I_sum / B_ssb1_I_rpc]
bukau_dt[, B_ssb1_T1_norm := B_ssb1_T1 / B_ssb1_T1_rpc]
bukau_dt[, B_ssb1_T2_norm := B_ssb1_T2 / B_ssb1_T2_rpc]
bukau_dt[, B_ssb1_I1_norm := B_ssb1_I1 / B_ssb1_I1_rpc]
bukau_dt[, B_ssb1_I2_norm := B_ssb1_I2 / B_ssb1_I2_rpc]
translatome_dt[, WT_norm := WT_sum / WT_rpc]
translatome_dt[, EV_norm := EV_sum / EV_rpc]
translatome_dt[, G_WT_norm := G_WT_sum / G_WT_rpc]
translatome_dt[, B_WT_norm := B_WT_sum / B_WT_rpc]
translatome_dt[, C_WT_norm := C_WT_sum / C_WT_rpc]
translatome_dt[, L_WT_norm := L_WT_sum / L_WT_rpc]
translatome_dt[, WT1_norm := WT1 / WT1_rpc]
translatome_dt[, WT2_norm := WT2 / WT2_rpc]
translatome_dt[, EV1_norm := EV1 / EV1_rpc]
translatome_dt[, EV2_norm := EV2 / EV2_rpc]
translatome_dt[, B_WT1_norm := B_WT1 / B_WT1_rpc]
translatome_dt[, B_WT2_norm := B_WT2 / B_WT2_rpc]
translatome_dt[, G_WT1_norm := G_WT1 / G_WT1_rpc]
translatome_dt[, G_WT2_norm := G_WT2 / G_WT2_rpc]
translatome_dt[, C_WT1_norm := C_WT1 / C_WT1_rpc]
translatome_dt[, C_WT2_norm := C_WT2 / C_WT2_rpc]
translatome_dt[, L_WT1_norm := L_WT1 / L_WT1_rpc]
translatome_dt[, L_WT2_norm := L_WT2 / L_WT2_rpc]
translatome_dt[, L_WT3_norm := L_WT3 / L_WT3_rpc]
atp2_dt[, ribo_norm := ribo_sum / ribo_rpc]
atp2_dt[, tric_norm := tric_sum / tric_rpc]
atp2_dt[, R1_norm := Atp2_R1 / R1_rpc]
atp2_dt[, R2_norm := Atp2_R2 / R2_rpc]
atp2_dt[, T1_norm := Atp2_T1 / T1_rpc]
atp2_dt[, T2_norm := Atp2_T2 / T2_rpc]


### Moving average for occupancy profiles
tric_dt[, tric_sum_ma := movingAverage(tric_sum, n=5, center=T)]
tric_dt[, ribo_sum_ma := movingAverage(ribo_sum, n=5, center=T)]
tric_dt[, ssb_Schx_sum_ma := movingAverage(ssb_Schx_sum, n=5, center=T)]
tric_dt[, ssb_Rchx_sum_ma := movingAverage(ssb_Rchx_sum, n=5, center=T)]


### Gene-level expression (tpm)
expression <- tric_dt[position > 5 & stopdist < -5]
expression[, ribo_total := sum(ribo_sum), by = orf]
expression[, tric_total := sum(tric_sum), by = orf]
expression[, Rpur_total := sum(Rpur_sum), by = orf]
expression[, Tpur_total := sum(Tpur_sum), by = orf]
expression[, R1X_total := sum(R1X), by = orf]
expression[, R2X_total := sum(R2X), by = orf]
expression[, T1X_total := sum(T1X), by = orf]
expression[, T2X_total := sum(T2X), by = orf]
expression[, R1pur_total := sum(R1pur), by = orf]
expression[, R2pur_total := sum(R2pur), by = orf]
expression[, T1pur_total := sum(T1pur), by = orf]
expression[, T2pur_total := sum(T2pur), by = orf]
expression[, R1X_rpm := (R1X / sum(expression$R1X)) * 10^6]
expression[, R2X_rpm := (R2X / sum(expression$R2X)) * 10^6]
expression[, T1X_rpm := (T1X / sum(expression$T1X)) * 10^6]
expression[, T2X_rpm := (T2X / sum(expression$T2X)) * 10^6]
expression[, ribo_rpm_gene := (ribo_total / sum(expression$ribo_sum)) * 10^6]
expression[, tric_rpm_gene := (tric_total / sum(expression$tric_sum)) * 10^6]
expression[, R1X_rpm_gene := (R1X_total / sum(expression$R1X)) * 10^6]
expression[, R2X_rpm_gene := (R2X_total / sum(expression$R2X)) * 10^6]
expression[, T1X_rpm_gene := (T1X_total / sum(expression$T1X)) * 10^6]
expression[, T2X_rpm_gene := (T2X_total / sum(expression$T2X)) * 10^6]
expression[, ribo_rpk := ribo_total / (length(ribo_sum) * 3 / 1000), by = orf]
expression[, tric_rpk := tric_total / (length(tric_sum) * 3 / 1000), by = orf]
expression[, Rpur_rpk := Rpur_total / (length(Rpur_sum) * 3 / 1000), by = orf]
expression[, Tpur_rpk := Tpur_total / (length(Tpur_sum) * 3 / 1000), by = orf]
expression[, R1X_rpk := R1X_total / (length(R1X) * 3 / 1000), by = orf]
expression[, R2X_rpk := R2X_total / (length(R2X) * 3 / 1000), by = orf]
expression[, T1X_rpk := T1X_total / (length(T1X) * 3 / 1000), by = orf]
expression[, T2X_rpk := T2X_total / (length(T2X) * 3 / 1000), by = orf]
expression[, R1pur_rpk := R1pur_total / (length(R1pur) * 3 / 1000), by = orf]
expression[, R2pur_rpk := R2pur_total / (length(R2pur) * 3 / 1000), by = orf]
expression[, T1pur_rpk := T1pur_total / (length(T1pur) * 3 / 1000), by = orf]
expression[, T2pur_rpk := T2pur_total / (length(T2pur) * 3 / 1000), by = orf]
expression[, ribo_rpkm := ribo_rpm_gene / (length(ribo_total) * 3 / 1000), by = orf]
expression[, tric_rpkm := tric_rpm_gene / (length(tric_total) * 3 / 1000), by = orf]
expression[, R1X_rpkm := R1X_rpm_gene / (length(R1X_total) * 3 / 1000), by = orf]
expression[, R2X_rpkm := R2X_rpm_gene / (length(R2X_total) * 3 / 1000), by = orf]
expression[, T1X_rpkm := T1X_rpm_gene / (length(T1X_total) * 3 / 1000), by = orf]
expression[, T2X_rpkm := T2X_rpm_gene / (length(T2X_total) * 3 / 1000), by = orf]
expression[, ribo_rpc := mean(ribo_sum), by = orf]
expression[, tric_rpc := mean(tric_sum), by = orf]
expression[, Rpur_rpc := mean(Rpur_sum), by = orf]
expression[, Tpur_rpc := mean(Tpur_sum), by = orf]
expression[, Rchx_rpc := mean(Rchx_sum), by = orf]
expression[, Tchx_rpc := mean(Tchx_sum), by = orf]
expression[, R1X_rpc := mean(R1X), by = orf]
expression[, R2X_rpc := mean(R2X), by = orf]
expression[, T1X_rpc := mean(T1X), by = orf]
expression[, T2X_rpc := mean(T2X), by = orf]
expression[, R1pur_rpc := mean(R1pur), by = orf]
expression[, T1pur_rpc := mean(T1pur), by = orf]
expression[, R2pur_rpc := mean(R2pur), by = orf]
expression[, T2pur_rpc := mean(T2pur), by = orf]
expression[, R1chx_rpc := mean(R1chx), by = orf]
expression[, T1chx_rpc := mean(T1chx), by = orf]
expression[, R2chx_rpc := mean(R2chx), by = orf]
expression[, T2chx_rpc := mean(T2chx), by = orf]
expression[, ssb_Rchx_total := sum(ssb_Rchx_sum), by = orf]
expression[, ssb_Schx_total := sum(ssb_Schx_sum), by = orf]
expression[, ssb_Rpur_total := sum(ssb_Rpur_sum), by = orf]
expression[, ssb_Spur_total := sum(ssb_Spur_sum), by = orf]
expression[, ssb_Rchx1_total := sum(SSB_Rchx1), by = orf]
expression[, ssb_Rchx2_total := sum(SSB_Rchx2), by = orf]
expression[, ssb_Schx1_total := sum(SSB_Schx1), by = orf]
expression[, ssb_Schx2_total := sum(SSB_Schx2), by = orf]
expression[, ssb_Rpur1_total := sum(SSB_Rpuro1), by = orf]
expression[, ssb_Rpur2_total := sum(SSB_Rpuro2), by = orf]
expression[, ssb_Spur1_total := sum(SSB_Spuro1), by = orf]
expression[, ssb_Spur2_total := sum(SSB_Spuro2), by = orf]
expression[, ssb_Rchx1_rpm := (SSB_Rchx1 / sum(expression$SSB_Rchx1)) * 10^6]
expression[, ssb_Rchx2_rpm := (SSB_Rchx2 / sum(expression$SSB_Rchx2)) * 10^6]
expression[, ssb_Schx1_rpm := (SSB_Schx1 / sum(expression$SSB_Schx1)) * 10^6]
expression[, ssb_Schx2_rpm := (SSB_Schx2 / sum(expression$SSB_Schx2)) * 10^6]
expression[, ssb_Rchx_rpm_gene := (ssb_Rchx_total / sum(expression$ssb_Rchx_sum)) * 10^6]
expression[, ssb_Schx_rpm_gene := (ssb_Schx_total / sum(expression$ssb_Schx_sum)) * 10^6]
expression[, ssb_Rchx1_rpm_gene := (ssb_Rchx1_total / sum(expression$SSB_Rchx1)) * 10^6]
expression[, ssb_Rchx2_rpm_gene := (ssb_Rchx2_total / sum(expression$SSB_Rchx2)) * 10^6]
expression[, ssb_Schx1_rpm_gene := (ssb_Schx1_total / sum(expression$SSB_Schx1)) * 10^6]
expression[, ssb_Schx2_rpm_gene := (ssb_Schx2_total / sum(expression$SSB_Schx2)) * 10^6]
expression[, ssb_Rchx_rpk := ssb_Rchx_total / (length(ssb_Rchx_sum) * 3 / 1000), by = orf]
expression[, ssb_Schx_rpk := ssb_Schx_total / (length(ssb_Schx_sum) * 3 / 1000), by = orf]
expression[, ssb_Rpur_rpk := ssb_Rpur_total / (length(ssb_Rpur_sum) * 3 / 1000), by = orf]
expression[, ssb_Spur_rpk := ssb_Spur_total / (length(ssb_Spur_sum) * 3 / 1000), by = orf]
expression[, ssb_Rchx1_rpk := ssb_Rchx1_total / (length(SSB_Rchx1) * 3 / 1000), by = orf]
expression[, ssb_Rchx2_rpk := ssb_Rchx2_total / (length(SSB_Rchx2) * 3 / 1000), by = orf]
expression[, ssb_Schx1_rpk := ssb_Schx1_total / (length(SSB_Schx1) * 3 / 1000), by = orf]
expression[, ssb_Schx2_rpk := ssb_Schx2_total / (length(SSB_Schx2) * 3 / 1000), by = orf]
expression[, ssb_Rpur1_rpk := ssb_Rpur1_total / (length(SSB_Rpuro1) * 3 / 1000), by = orf]
expression[, ssb_Rpur2_rpk := ssb_Rpur2_total / (length(SSB_Rpuro2) * 3 / 1000), by = orf]
expression[, ssb_Spur1_rpk := ssb_Spur1_total / (length(SSB_Spuro1) * 3 / 1000), by = orf]
expression[, ssb_Spur2_rpk := ssb_Spur2_total / (length(SSB_Spuro2) * 3 / 1000), by = orf]
expression[, ssb_Rchx_rpkm := ssb_Rchx_rpm_gene / (length(ssb_Rchx_total) * 3 / 1000), by = orf]
expression[, ssb_Schx_rpkm := ssb_Schx_rpm_gene / (length(ssb_Schx_total) * 3 / 1000), by = orf]
expression[, ssb_Rchx1_rpkm := ssb_Rchx1_rpm_gene / (length(ssb_Rchx1_total) * 3 / 1000), by = orf]
expression[, ssb_Rchx2_rpkm := ssb_Rchx2_rpm_gene / (length(ssb_Rchx2_total) * 3 / 1000), by = orf]
expression[, ssb_Schx1_rpkm := ssb_Schx1_rpm_gene / (length(ssb_Schx1_total) * 3 / 1000), by = orf]
expression[, ssb_Schx2_rpkm := ssb_Schx2_rpm_gene / (length(ssb_Schx2_total) * 3 / 1000), by = orf]
expression[, ssb_Rchx_rpc := mean(ssb_Rchx_sum), by = orf]
expression[, ssb_Schx_rpc := mean(ssb_Schx_sum), by = orf]
expression[, ssb_Rpur_rpc := mean(ssb_Rpur_sum), by = orf]
expression[, ssb_Spur_rpc := mean(ssb_Spur_sum), by = orf]
expression[, ssb_Rchx1_rpc := mean(SSB_Rchx1), by = orf]
expression[, ssb_Rchx2_rpc := mean(SSB_Rchx2), by = orf]
expression[, ssb_Schx1_rpc := mean(SSB_Schx1), by = orf]
expression[, ssb_Schx2_rpc := mean(SSB_Schx2), by = orf]
expression[, ssb_Rpur1_rpc := mean(SSB_Rpuro1), by = orf]
expression[, ssb_Rpur2_rpc := mean(SSB_Rpuro2), by = orf]
expression[, ssb_Spur1_rpc := mean(SSB_Spuro1), by = orf]
expression[, ssb_Spur2_rpc := mean(SSB_Spuro2), by = orf]

# Calculate correlation for each orf
orfs <- expression[, .SD[which.min(position)], by = orf]
ribo_cor <- NULL
tric_cor <- NULL
ssb_Rchx_cor <- NULL
ssb_Schx_cor <- NULL
gene <- NULL
for (g in orfs$orf) {
  ribo_cor1 <- cor(expression[g]$R1X_rpm, expression[g]$R2X_rpm)
  ribo_cor <- c(ribo_cor, ribo_cor1)
  tric_cor1 <- cor(expression[g]$T1X_rpm, expression[g]$T2X_rpm)
  tric_cor <- c(tric_cor, tric_cor1)
  ssb_Rchx_cor1 <- cor(expression[g]$ssb_Rchx1_rpm, expression[g]$ssb_Rchx2_rpm)
  ssb_Rchx_cor <- c(ssb_Rchx_cor, ssb_Rchx_cor1)
  ssb_Schx_cor1 <- cor(expression[g]$ssb_Schx1_rpm, expression[g]$ssb_Schx2_rpm)
  ssb_Schx_cor <- c(ssb_Schx_cor, ssb_Schx_cor1)
  gene <- c(gene, g)
}
test <- data.table(orf = gene, ribo_cor = ribo_cor, tric_cor = tric_cor,
                   ssb_Rchx_cor = ssb_Rchx_cor, ssb_Schx_cor = ssb_Schx_cor)
expression1 <- expression[, .SD[which.min(position)], by = orf]
expression1[, ribo_tpm := (ribo_rpk / sum(expression1$ribo_rpk)) * 10^6]
expression1[, tric_tpm := (tric_rpk / sum(expression1$tric_rpk)) * 10^6]
expression1[, Rpur_tpm := (Rpur_rpk / sum(expression1$Rpur_rpk)) * 10^6]
expression1[, Tpur_tpm := (Tpur_rpk / sum(expression1$Tpur_rpk)) * 10^6]
expression1[, R1X_tpm := (R1X_rpk / sum(expression1$R1X_rpk)) * 10^6]  
expression1[, R2X_tpm := (R2X_rpk / sum(expression1$R2X_rpk)) * 10^6]  
expression1[, T1X_tpm := (T1X_rpk / sum(expression1$T1X_rpk)) * 10^6]  
expression1[, T2X_tpm := (T2X_rpk / sum(expression1$T2X_rpk)) * 10^6]
expression1[, R1pur_tpm := (R1pur_rpk / sum(expression1$R1pur_rpk)) * 10^6]  
expression1[, R2pur_tpm := (R2pur_rpk / sum(expression1$R2pur_rpk)) * 10^6]  
expression1[, T1pur_tpm := (T1pur_rpk / sum(expression1$T1pur_rpk)) * 10^6]  
expression1[, T2pur_tpm := (T2pur_rpk / sum(expression1$T2pur_rpk)) * 10^6]
expression1[, ssb_Rchx_tpm := (ssb_Rchx_rpk / sum(expression1$ssb_Rchx_rpk)) * 10^6]
expression1[, ssb_Schx_tpm := (ssb_Schx_rpk / sum(expression1$ssb_Schx_rpk)) * 10^6]
expression1[, ssb_Rpur_tpm := (ssb_Rpur_rpk / sum(expression1$ssb_Rpur_rpk)) * 10^6]
expression1[, ssb_Spur_tpm := (ssb_Spur_rpk / sum(expression1$ssb_Spur_rpk)) * 10^6]
expression1[, ssb_Rchx1_tpm := (ssb_Rchx1_rpk / sum(expression1$ssb_Rchx1_rpk)) * 10^6]
expression1[, ssb_Rchx2_tpm := (ssb_Rchx2_rpk / sum(expression1$ssb_Rchx2_rpk)) * 10^6]
expression1[, ssb_Schx1_tpm := (ssb_Schx1_rpk / sum(expression1$ssb_Schx1_rpk)) * 10^6]
expression1[, ssb_Schx2_tpm := (ssb_Schx2_rpk / sum(expression1$ssb_Schx2_rpk)) * 10^6]
expression1[, ssb_Rpur1_tpm := (ssb_Rpur1_rpk / sum(expression1$ssb_Rpur1_rpk)) * 10^6]
expression1[, ssb_Rpur2_tpm := (ssb_Rpur2_rpk / sum(expression1$ssb_Rpur2_rpk)) * 10^6]
expression1[, ssb_Spur1_tpm := (ssb_Spur1_rpk / sum(expression1$ssb_Spur1_rpk)) * 10^6]
expression1[, ssb_Spur2_tpm := (ssb_Spur2_rpk / sum(expression1$ssb_Spur2_rpk)) * 10^6]


expressionBukau <- bukau_dt[position > 5 & stopdist < -5]
expressionBukau[, B_ssb1_T_total := sum(B_ssb1_T_sum), by = orf]
expressionBukau[, B_ssb1_I_total := sum(B_ssb1_I_sum), by = orf]
expressionBukau[, B_ssb1_T1_total := sum(B_ssb1_T1), by = orf]
expressionBukau[, B_ssb1_T2_total := sum(B_ssb1_T2), by = orf]
expressionBukau[, B_ssb1_I1_total := sum(B_ssb1_I1), by = orf]
expressionBukau[, B_ssb1_I2_total := sum(B_ssb1_I2), by = orf]
expressionBukau[, B_ssb1_T1_rpm := (B_ssb1_T1 / sum(expressionBukau$B_ssb1_T1)) * 10^6]
expressionBukau[, B_ssb1_T2_rpm := (B_ssb1_T2 / sum(expressionBukau$B_ssb1_T2)) * 10^6]
expressionBukau[, B_ssb1_I1_rpm := (B_ssb1_I1 / sum(expressionBukau$B_ssb1_I1)) * 10^6]
expressionBukau[, B_ssb1_I2_rpm := (B_ssb1_I2 / sum(expressionBukau$B_ssb1_I2)) * 10^6]
expressionBukau[, B_ssb1_T_rpm_gene := (B_ssb1_T_total / sum(expressionBukau$B_ssb1_T_sum)) * 10^6]
expressionBukau[, B_ssb1_I_rpm_gene := (B_ssb1_I_total / sum(expressionBukau$B_ssb1_I_sum)) * 10^6]
expressionBukau[, B_ssb1_T1_rpm_gene := (B_ssb1_T1_total / sum(expressionBukau$B_ssb1_T1)) * 10^6]
expressionBukau[, B_ssb1_T2_rpm_gene := (B_ssb1_T2_total / sum(expressionBukau$B_ssb1_T2)) * 10^6]
expressionBukau[, B_ssb1_I1_rpm_gene := (B_ssb1_I1_total / sum(expressionBukau$B_ssb1_I1)) * 10^6]
expressionBukau[, B_ssb1_I2_rpm_gene := (B_ssb1_I2_total / sum(expressionBukau$B_ssb1_I2)) * 10^6]
expressionBukau[, B_ssb1_T_rpk := B_ssb1_T_total / (length(B_ssb1_T_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_I_rpk := B_ssb1_I_total / (length(B_ssb1_I_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_T1_rpk := B_ssb1_T1_total / (length(B_ssb1_T1_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_T2_rpk := B_ssb1_T2_total / (length(B_ssb1_T2_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_I1_rpk := B_ssb1_I1_total / (length(B_ssb1_I1_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_I2_rpk := B_ssb1_I2_total / (length(B_ssb1_I2_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_T_rpkm := B_ssb1_T_rpm_gene / (length(B_ssb1_T_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_I_rpkm := B_ssb1_I_rpm_gene / (length(B_ssb1_I_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_T1_rpkm := B_ssb1_T1_rpm_gene / (length(B_ssb1_T1_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_T2_rpkm := B_ssb1_T2_rpm_gene / (length(B_ssb1_T2_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_I1_rpkm := B_ssb1_I1_rpm_gene / (length(B_ssb1_I1_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_I2_rpkm := B_ssb1_I2_rpm_gene / (length(B_ssb1_I2_total) * 3 / 1000), by = orf]
expressionBukau[, B_ssb1_T_rpc := mean(B_ssb1_T_sum), by = orf]
expressionBukau[, B_ssb1_I_rpc := mean(B_ssb1_I_sum), by = orf]
expressionBukau[, B_ssb1_T1_rpc := mean(B_ssb1_T1), by = orf]
expressionBukau[, B_ssb1_T2_rpc := mean(B_ssb1_T2), by = orf]
expressionBukau[, B_ssb1_I1_rpc := mean(B_ssb1_I1), by = orf]
expressionBukau[, B_ssb1_I2_rpc := mean(B_ssb1_I2), by = orf]

# Calculate correlation for each orf
orfs <- expressionBukau[, .SD[which.min(position)], by = orf]
B_ssb1_T_cor <- NULL
B_ssb1_I_cor <- NULL
gene <- NULL
for (g in orfs$orf) {
  B_ssb1_T_cor1 <- cor(expressionBukau[g]$B_ssb1_T1_rpm, expressionBukau[g]$B_ssb1_T2_rpm)
  B_ssb1_T_cor <- c(B_ssb1_T_cor, B_ssb1_T_cor1)
  B_ssb1_I_cor1 <- cor(expressionBukau[g]$B_ssb1_I1_rpm, expressionBukau[g]$B_ssb1_I2_rpm)
  B_ssb1_I_cor <- c(B_ssb1_I_cor, B_ssb1_I_cor1)
  gene <- c(gene, g)
}
test <- data.table(orf = gene, B_ssb1_T_cor = B_ssb1_T_cor, B_ssb1_I_cor = B_ssb1_I_cor)
expression1Bukau <- expressionBukau[, .SD[which.min(position)], by = orf]
expression1Bukau[, B_ssb1_T_tpm := (B_ssb1_T_rpk / sum(expression1Bukau$B_ssb1_T_rpk)) * 10^6]
expression1Bukau[, B_ssb1_I_tpm := (B_ssb1_I_rpk / sum(expression1Bukau$B_ssb1_I_rpk)) * 10^6]
expression1Bukau[, B_ssb1_T1_tpm := (B_ssb1_T1_rpk / sum(expression1Bukau$B_ssb1_T1_rpk)) * 10^6]
expression1Bukau[, B_ssb1_T2_tpm := (B_ssb1_T2_rpk / sum(expression1Bukau$B_ssb1_T2_rpk)) * 10^6]
expression1Bukau[, B_ssb1_I1_tpm := (B_ssb1_I1_rpk / sum(expression1Bukau$B_ssb1_I1_rpk)) * 10^6]
expression1Bukau[, B_ssb1_I2_tpm := (B_ssb1_I2_rpk / sum(expression1Bukau$B_ssb1_I2_rpk)) * 10^6]


expressionTranslatome <- translatome_dt[position > 5 & stopdist < -5]
expressionTranslatome[, WT_total := sum(WT_sum), by = orf]
expressionTranslatome[, EV_total := sum(EV_sum), by = orf]
expressionTranslatome[, G_WT_total := sum(G_WT_sum), by = orf]
expressionTranslatome[, B_WT_total := sum(B_WT_sum), by = orf]
expressionTranslatome[, C_WT_total := sum(C_WT_sum), by = orf]
expressionTranslatome[, L_WT_total := sum(L_WT_sum), by = orf]
expressionTranslatome[, WT1_total := sum(WT1), by = orf]
expressionTranslatome[, WT2_total := sum(WT2), by = orf]
expressionTranslatome[, EV1_total := sum(EV1), by = orf]
expressionTranslatome[, EV2_total := sum(EV2), by = orf]
expressionTranslatome[, B_WT1_total := sum(B_WT1), by = orf]
expressionTranslatome[, B_WT2_total := sum(B_WT2), by = orf]
expressionTranslatome[, G_WT1_total := sum(G_WT1), by = orf]
expressionTranslatome[, G_WT2_total := sum(G_WT2), by = orf]
expressionTranslatome[, C_WT1_total := sum(C_WT1), by = orf]
expressionTranslatome[, C_WT2_total := sum(C_WT2), by = orf]
expressionTranslatome[, L_WT1_total := sum(L_WT1), by = orf]
expressionTranslatome[, L_WT2_total := sum(L_WT2), by = orf]
expressionTranslatome[, L_WT3_total := sum(L_WT3), by = orf]
expressionTranslatome[, WT_rpk := WT_total / (length(WT_total) * 3 / 1000), by = orf]
expressionTranslatome[, EV_rpk := EV_total / (length(EV_total) * 3 / 1000), by = orf]
expressionTranslatome[, B_WT_rpk := B_WT_total / (length(B_WT_total) * 3 / 1000), by = orf]
expressionTranslatome[, G_WT_rpk := G_WT_total / (length(G_WT_total) * 3 / 1000), by = orf]
expressionTranslatome[, C_WT_rpk := C_WT_total / (length(C_WT_total) * 3 / 1000), by = orf]
expressionTranslatome[, L_WT_rpk := L_WT_total / (length(L_WT_total) * 3 / 1000), by = orf]
expressionTranslatome[, WT1_rpk := WT1_total / (length(WT1) * 3 / 1000), by = orf]
expressionTranslatome[, WT2_rpk := WT2_total / (length(WT2) * 3 / 1000), by = orf]
expressionTranslatome[, EV1_rpk := EV1_total / (length(EV1) * 3 / 1000), by = orf]
expressionTranslatome[, EV2_rpk := EV2_total / (length(EV2) * 3 / 1000), by = orf]
expressionTranslatome[, B_WT1_rpk := B_WT1_total / (length(B_WT1_total) * 3 / 1000), by = orf]
expressionTranslatome[, B_WT2_rpk := B_WT2_total / (length(B_WT2_total) * 3 / 1000), by = orf]
expressionTranslatome[, G_WT1_rpk := G_WT1_total / (length(G_WT1_total) * 3 / 1000), by = orf]
expressionTranslatome[, G_WT2_rpk := G_WT2_total / (length(G_WT2_total) * 3 / 1000), by = orf]
expressionTranslatome[, C_WT1_rpk := C_WT1_total / (length(C_WT1_total) * 3 / 1000), by = orf]
expressionTranslatome[, C_WT2_rpk := C_WT2_total / (length(C_WT2_total) * 3 / 1000), by = orf]
expressionTranslatome[, L_WT1_rpk := L_WT1_total / (length(L_WT1_total) * 3 / 1000), by = orf]
expressionTranslatome[, L_WT2_rpk := L_WT2_total / (length(L_WT2_total) * 3 / 1000), by = orf]
expressionTranslatome[, L_WT3_rpk := L_WT3_total / (length(L_WT3_total) * 3 / 1000), by = orf]
expressionTranslatome[, WT_rpc := mean(WT_sum), by = orf]
expressionTranslatome[, EV_rpc := mean(EV_sum), by = orf]
expressionTranslatome[, G_WT_rpc := mean(G_WT_sum), by = orf]
expressionTranslatome[, B_WT_rpc := mean(B_WT_sum), by = orf]
expressionTranslatome[, C_WT_rpc := mean(C_WT_sum), by = orf]
expressionTranslatome[, L_WT_rpc := mean(L_WT_sum), by = orf]
expressionTranslatome[, WT1_rpc := mean(WT1), by = orf]
expressionTranslatome[, WT2_rpc := mean(WT2), by = orf]
expressionTranslatome[, EV1_rpc := mean(EV1), by = orf]
expressionTranslatome[, EV2_rpc := mean(EV2), by = orf]
expressionTranslatome[, G_WT1_rpc := mean(G_WT1), by = orf]
expressionTranslatome[, G_WT2_rpc := mean(G_WT2), by = orf]
expressionTranslatome[, B_WT1_rpc := mean(B_WT1), by = orf]
expressionTranslatome[, B_WT2_rpc := mean(B_WT2), by = orf]
expressionTranslatome[, C_WT1_rpc := mean(C_WT1), by = orf]
expressionTranslatome[, C_WT2_rpc := mean(C_WT2), by = orf]
expressionTranslatome[, L_WT1_rpc := mean(L_WT1), by = orf]
expressionTranslatome[, L_WT2_rpc := mean(L_WT2), by = orf]
expressionTranslatome[, L_WT3_rpc := mean(L_WT3), by = orf]
expression1Translatome <- expressionTranslatome[, .SD[which.min(position)], by = orf]
expression1Translatome[, WT_tpm := (WT_rpk / sum(expression1Translatome$WT_rpk)) * 10^6]
expression1Translatome[, EV_tpm := (EV_rpk / sum(expression1Translatome$EV_rpk)) * 10^6]
expression1Translatome[, B_WT_tpm := (B_WT_rpk / sum(expression1Translatome$B_WT_rpk)) * 10^6]
expression1Translatome[, G_WT_tpm := (G_WT_rpk / sum(expression1Translatome$G_WT_rpk)) * 10^6]
expression1Translatome[, C_WT_tpm := (C_WT_rpk / sum(expression1Translatome$C_WT_rpk)) * 10^6]
expression1Translatome[, L_WT_tpm := (L_WT_rpk / sum(expression1Translatome$L_WT_rpk)) * 10^6]
expression1Translatome[, WT1_tpm := (WT1_rpk / sum(expression1Translatome$WT1_rpk)) * 10^6]
expression1Translatome[, WT2_tpm := (WT2_rpk / sum(expression1Translatome$WT2_rpk)) * 10^6]
expression1Translatome[, EV1_tpm := (EV1_rpk / sum(expression1Translatome$EV1_rpk)) * 10^6]
expression1Translatome[, EV2_tpm := (EV2_rpk / sum(expression1Translatome$EV2_rpk)) * 10^6]
expression1Translatome[, B_WT1_tpm := (B_WT1_rpk / sum(expression1Translatome$B_WT1_rpk)) * 10^6]
expression1Translatome[, B_WT2_tpm := (B_WT2_rpk / sum(expression1Translatome$B_WT2_rpk)) * 10^6]
expression1Translatome[, G_WT1_tpm := (G_WT1_rpk / sum(expression1Translatome$G_WT1_rpk)) * 10^6]
expression1Translatome[, G_WT2_tpm := (G_WT2_rpk / sum(expression1Translatome$G_WT2_rpk)) * 10^6]
expression1Translatome[, C_WT1_tpm := (C_WT1_rpk / sum(expression1Translatome$C_WT1_rpk)) * 10^6]
expression1Translatome[, C_WT2_tpm := (C_WT2_rpk / sum(expression1Translatome$C_WT2_rpk)) * 10^6]
expression1Translatome[, L_WT1_tpm := (L_WT1_rpk / sum(expression1Translatome$L_WT1_rpk)) * 10^6]
expression1Translatome[, L_WT2_tpm := (L_WT2_rpk / sum(expression1Translatome$L_WT2_rpk)) * 10^6]
expression1Translatome[, L_WT3_tpm := (L_WT3_rpk / sum(expression1Translatome$L_WT3_rpk)) * 10^6]


expression_atp2 <- atp2_dt[position > 5 & stopdist < -5]
expression_atp2[, ribo_total := sum(ribo_sum), by = orf]
expression_atp2[, tric_total := sum(tric_sum), by = orf]
expression_atp2[, R1_total := sum(Atp2_R1), by = orf]
expression_atp2[, R2_total := sum(Atp2_R2), by = orf]
expression_atp2[, T1_total := sum(Atp2_T1), by = orf]
expression_atp2[, T2_total := sum(Atp2_T2), by = orf]
expression_atp2[, ribo_rpk := ribo_total / (length(ribo_total) * 3 / 1000), by = orf]
expression_atp2[, tric_rpk := tric_total / (length(tric_total) * 3 / 1000), by = orf]
expression_atp2[, R1_rpk := R1_total / (length(R1_total) * 3 / 1000), by = orf]
expression_atp2[, R2_rpk := R2_total / (length(R2_total) * 3 / 1000), by = orf]
expression_atp2[, T1_rpk := T1_total / (length(T1_total) * 3 / 1000), by = orf]
expression_atp2[, T2_rpk := T2_total / (length(T2_total) * 3 / 1000), by = orf]
expression_atp2[, ribo_rpc := mean(ribo_sum), by = orf]
expression_atp2[, tric_rpc := mean(tric_sum), by = orf]
expression_atp2[, R1_rpc := mean(Atp2_R1), by = orf]
expression_atp2[, R2_rpc := mean(Atp2_R2), by = orf]
expression_atp2[, T1_rpc := mean(Atp2_T1), by = orf]
expression_atp2[, T2_rpc := mean(Atp2_T2), by = orf]
expression1_atp2 <- expression_atp2[, .SD[which.min(position)], by = orf]
expression1_atp2[, ribo_tpm := (ribo_rpk / sum(expression1_atp2$ribo_rpk)) * 10^6]
expression1_atp2[, tric_tpm := (tric_rpk / sum(expression1_atp2$tric_rpk)) * 10^6]
expression1_atp2[, R1_tpm := (R1_rpk / sum(expression1_atp2$R1_rpk)) * 10^6]  
expression1_atp2[, R2_tpm := (R2_rpk / sum(expression1_atp2$R2_rpk)) * 10^6]  
expression1_atp2[, T1_tpm := (T1_rpk / sum(expression1_atp2$T1_rpk)) * 10^6]  
expression1_atp2[, T2_tpm := (T2_rpk / sum(expression1_atp2$T2_rpk)) * 10^6]


### Processing of libraries from samples not crosslinked
R1O <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/R1O.codon", header = F, stringsAsFactors = T)
R2O <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/R2O.codon", header = F, stringsAsFactors = T)
T1O <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/T1O.codon", header = F, stringsAsFactors = T)
T2O <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Codon_density/T2O.codon", header = F, stringsAsFactors = T)

tric_noX_dt <- data.table(orf = R1O[, 1], 
                          position = R1O[, 2],
                          codon = R1O[, 3],
                          R1O = R1O[, 4],
                          T1O = T1O[, 4],
                          R2O = R2O[, 4],
                          T2O = T2O[, 4])
setkeyv(tric_noX_dt, c("orf"))
tric_noX_dt[, length := length(position), by = orf]
tric_noX_dt[, length := length - 15] # substract flanking 7 codons and stop codon
tric_noX_dt[, stopdist := position - (length + 1)] # stop codon is 0
tric_noX_dt <- tric_noX_dt[orf != "YOR031W" & orf != "YFL057C" & orf != "YOL013W-A"] # overlapping or blocked orf


i <- cbind(match(tric_noX_dt$codon, codon_table$codon))
tric_noX_dt <- cbind(tric_noX_dt, residue = codon_table[i]$residue)

tric_noX_dt[, ribo_sum := R1O + R2O]
tric_noX_dt[, tric_sum := T1O + T2O]

i <- cbind(match(tric_noX_dt$orf, expression1_noX$orf))
tric_noX_dt <- cbind(tric_noX_dt, ribo_rpc = expression1_noX[i]$ribo_rpc)
tric_noX_dt <- cbind(tric_noX_dt, ribo_total = expression1_noX[i]$ribo_total)
tric_noX_dt <- cbind(tric_noX_dt, tric_rpc = expression1_noX[i]$tric_rpc)
tric_noX_dt <- cbind(tric_noX_dt, tric_total = expression1_noX[i]$tric_total)
tric_noX_dt <- cbind(tric_noX_dt, R1O_rpc = expression1_noX[i]$R1O_rpc)
tric_noX_dt <- cbind(tric_noX_dt, R2O_rpc = expression1_noX[i]$R2O_rpc)
tric_noX_dt <- cbind(tric_noX_dt, T1O_rpc = expression1_noX[i]$T1O_rpc)
tric_noX_dt <- cbind(tric_noX_dt, T2O_rpc = expression1_noX[i]$T2O_rpc)
tric_noX_dt <- cbind(tric_noX_dt, R1O_total = expression1_noX[i]$R1O_total)
tric_noX_dt <- cbind(tric_noX_dt, R2O_total = expression1_noX[i]$R2O_total)
tric_noX_dt <- cbind(tric_noX_dt, T1O_total = expression1_noX[i]$T1O_total)
tric_noX_dt <- cbind(tric_noX_dt, T2O_total = expression1_noX[i]$T2O_total)


tric_noX_dt[, ribo_rpm := (ribo_sum / sum(tric_noX_dt$ribo_sum)) * 10^6]
tric_noX_dt[, tric_rpm := (tric_sum / sum(tric_noX_dt$tric_sum)) * 10^6]
tric_noX_dt[, R1O_rpm := (R1O / sum(tric_noX_dt$R1O)) * 10^6]  
tric_noX_dt[, T1O_rpm := (T1O / sum(tric_noX_dt$T1O)) * 10^6]  
tric_noX_dt[, R2O_rpm := (R2O / sum(tric_noX_dt$R2O)) * 10^6]  
tric_noX_dt[, T2O_rpm := (T2O / sum(tric_noX_dt$T2O)) * 10^6]
tric_noX_dt[, ribo_norm := ribo_sum / ribo_rpc]
tric_noX_dt[, tric_norm := tric_sum / tric_rpc]
tric_noX_dt[, R1O_norm := R1O / R1O_rpc]
tric_noX_dt[, R2O_norm := R2O / R2O_rpc]
tric_noX_dt[, T1O_norm := T1O / T1O_rpc]
tric_noX_dt[, T2O_norm := T2O / T2O_rpc]


expression_noX <- tric_noX_dt[position > 5 & stopdist < -5]
expression_noX[, ribo_total := sum(ribo_sum), by = orf]
expression_noX[, tric_total := sum(tric_sum), by = orf]
expression_noX[, R1O_total := sum(R1O), by = orf]
expression_noX[, R2O_total := sum(R2O), by = orf]
expression_noX[, T1O_total := sum(T1O), by = orf]
expression_noX[, T2O_total := sum(T2O), by = orf]
expression_noX[, ribo_rpk := ribo_total / (length(ribo_total) * 3 / 1000), by = orf]
expression_noX[, tric_rpk := tric_total / (length(tric_total) * 3 / 1000), by = orf]
expression_noX[, R1O_rpk := R1O_total / (length(R1O_total) * 3 / 1000), by = orf]
expression_noX[, R2O_rpk := R2O_total / (length(R2O_total) * 3 / 1000), by = orf]
expression_noX[, T1O_rpk := T1O_total / (length(T1O_total) * 3 / 1000), by = orf]
expression_noX[, T2O_rpk := T2O_total / (length(T2O_total) * 3 / 1000), by = orf]
expression_noX[, ribo_rpc := mean(ribo_sum), by = orf]
expression_noX[, tric_rpc := mean(tric_sum), by = orf]
expression_noX[, R1O_rpc := mean(R1O), by = orf]
expression_noX[, R2O_rpc := mean(R2O), by = orf]
expression_noX[, T1O_rpc := mean(T1O), by = orf]
expression_noX[, T2O_rpc := mean(T2O), by = orf]
expression1_noX <- expression_noX[, .SD[which.min(position)], by = orf]
expression1_noX[, ribo_tpm := (ribo_rpk / sum(expression1_noX$ribo_rpk)) * 10^6]
expression1_noX[, tric_tpm := (tric_rpk / sum(expression1_noX$tric_rpk)) * 10^6]
expression1_noX[, R1O_tpm := (R1O_rpk / sum(expression1_noX$R1O_rpk)) * 10^6]  
expression1_noX[, R2O_tpm := (R2O_rpk / sum(expression1_noX$R2O_rpk)) * 10^6]  
expression1_noX[, T1O_tpm := (T1O_rpk / sum(expression1_noX$T1O_rpk)) * 10^6]  
expression1_noX[, T2O_tpm := (T2O_rpk / sum(expression1_noX$T2O_rpk)) * 10^6]


### Create codon table
codon <- c("TTT","TTC", "TTA", "TTG",
           "CTT","CTC", "CTA", "CTG",
           "ATT","ATC", "ATA", "ATG",
           "GTT","GTC", "GTA", "GTG",
           "TCT","TCC", "TCA", "TCG",
           "CCT","CCC", "CCA", "CCG",
           "ACT","ACC", "ACA", "ACG",
           "GCT","GCC", "GCA", "GCG",
           "TAT","TAC", "TAA", "TAG",
           "CAT","CAC", "CAA", "CAG",
           "AAT","AAC", "AAA", "AAG",
           "GAT","GAC", "GAA", "GAG",
           "TGT","TGC", "TGA", "TGG",
           "CGT","CGC", "CGA", "CGG",
           "AGT","AGC", "AGA", "AGG",
           "GGT","GGC", "GGA", "GGG")
residue <- c("F", "F", "L", "L",
             "L", "L", "L", "L",
             "I", "I", "I", "M",
             "V", "V", "V", "V",
             "S", "S", "S", "S",
             "P", "P", "P", "P",
             "T", "T", "T", "T",
             "A", "A", "A", "A",
             "Y", "Y", "TAA", "TAG",
             "H", "H", "Q", "Q",
             "N", "N", "K", "K",
             "D", "D", "E", "E",
             "C", "C", "TGA", "W",
             "R", "R", "R", "R",
             "S", "S", "R", "R",
             "G", "G", "G", "G")
codon_table <- data.table(codon = codon,
                          residue = residue)
