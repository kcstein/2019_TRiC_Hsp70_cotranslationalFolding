### For ribosome occupancy profiles
temp <- tric_dt[orf == "YGR145W"]
temp <- bukau_dt[orf == "YGR145W"]
temp <- as.data.table(unique(wd_repeats[wd_repeats$orf %in% tric_substrates5$orf]$orf))
temp <- tric_dt[tric_dt$orf %in% temp$V1]
temp[, tric_sum1 := tric_sum + 1]
temp[, ribo_sum1 := ribo_sum + 1]
temp[, ssb_Rchx_sum1 := ssb_Rchx_sum + 1]
temp[, ssb_Schx_sum1 := ssb_Schx_sum + 1]
temp[, B_ssb1_T_sum1 := B_ssb1_T_sum + 1]
temp[, B_ssb1_I_sum1 := B_ssb1_I_sum + 1]
x <- temp[position==0]$length - 24
temp[, tric_total1 := tric_total + x]
temp[, ribo_total1 := ribo_total + x]
temp[, ssb_Rchx_total1 := ssb_Rchx_total + x]
temp[, ssb_Schx_total1 := ssb_Schx_total + x]
temp[, B_ssb1_T_total1 := B_ssb1_T_total + x]
temp[, B_ssb1_I_total1 := B_ssb1_I_total + x]
temp[, tric_sum1_ma := movingAverage(tric_sum, n=5, center=T)]
temp[, ribo_sum1_ma := movingAverage(ribo_sum, n=5, center=T)]
temp[, ssb_Rchx_sum1_ma := movingAverage(ssb_Rchx_sum, n=5, center=T)]
temp[, ssb_Schx_sum1_ma := movingAverage(ssb_Schx_sum, n=5, center=T)]
temp[, B_ssb1_T_sum1_ma := movingAverage(B_ssb1_T_sum1, n=5, center=T)]
temp[, B_ssb1_I_sum1_ma := movingAverage(B_ssb1_I_sum1, n=5, center=T)]
for (i in 1:nrow(temp)) {
  # tric <- matrix(c(temp[i]$tric_sum1_ma, temp[i]$ribo_sum1_ma, (temp[i]$tric_total - temp[i]$tric_sum1_ma), (temp[i]$ribo_total - temp[i]$ribo_sum1_ma)), nrow = 2)
  # tric_test <- fisher.test(tric)
  # temp[i, tric_odds := tric_test$estimate]
  # temp[i, tric_pvalue := tric_test$p.value]
  # ssb <- matrix(c(temp[i]$ssb_Schx_sum1_ma, temp[i]$ssb_Rchx_sum1_ma, (temp[i]$ssb_Schx_total - temp[i]$ssb_Schx_sum1_ma), (temp[i]$ssb_Rchx_total - temp[i]$ssb_Rchx_sum1_ma)), nrow = 2)
  # ssb_test <- fisher.test(ssb)
  # temp[i, ssb_odds := ssb_test$estimate]
  # temp[i, ssb_pvalue := ssb_test$p.value]
  ssb <- matrix(c(temp[i]$B_ssb1_I_sum1_ma, temp[i]$B_ssb1_T_sum1_ma, (temp[i]$B_ssb1_I_total1 - temp[i]$B_ssb1_I_sum1_ma), (temp[i]$B_ssb1_T_total1 - temp[i]$B_ssb1_T_sum1_ma)), nrow = 2)
  ssb_test <- fisher.test(ssb)
  temp[i, ssb_odds := ssb_test$estimate]
  temp[i, ssb_pvalue := ssb_test$p.value]
}
temp[, tric_padj := p.adjust(tric_pvalue, method = "BH"), by = orf]
temp[, ssb_padj := p.adjust(ssb_pvalue, method = "BH"), by = orf]
temp[, tric_odds_ma := movingAverage(tric_odds, n=5, center=T), by = orf]
temp[, ssb_odds_ma := movingAverage(ssb_odds, n=5, center=T), by = orf]


### Calculate fisher at each position for TRiC
tric_fishers <- tric_dt[ribo_cor >= 0.5 & tric_cor >= 0.5 & 
                          R1X_rpc >= 0.5 & R2X_rpc >= 0.5 & T1X_rpc >= 0.5 & T2X_rpc >= 0.5 & 
                          ribo_total >= 128 & tric_total >= 128, c(1:2,27:29,38:39,42:43)]
View(tric_fishers[(tric_total - tric_sum) < 0]) # Find orfs with position with aberrantly high number of reads
View(tric_fishers[(ribo_total - ribo_sum) < 0])
tric_fishers <- tric_fishers[orf != "YDL133C-A" & orf != "YPL220W" & orf != "YDL184C"]
for (i in 1:nrow(tric_fishers)) {
  tric <- matrix(c(tric_fishers[i]$tric_sum, tric_fishers[i]$ribo_sum, (tric_fishers[i]$tric_total - tric_fishers[i]$tric_sum), (tric_fishers[i]$ribo_total - tric_fishers[i]$ribo_sum)), nrow = 2)
  tric_test <- fisher.test(tric)
  tric_fishers[i, tric_odds := tric_test$estimate]
  tric_fishers[i, tric_pvalue := tric_test$p.value]
}
tric_fishers[, tric_padj := p.adjust(tric_pvalue, method = "BH"), by = orf]
tric_fishers[, tric_odds_ma := movingAverage(tric_odds, n=5, center=T), by = orf]
tric_peaks_all <- tric_fishers[tric_odds > 1 & tric_odds < Inf & tric_padj < 0.05 & position > 30 & tric_sum > tric_rpc]


### Calculate fisher at each position for SSB
ssb_fishers <- tric_dt[ssb_Schx_cor >= 0.5 & ssb_Rchx_cor >= 0.5 & 
                         ssb_Rchx1_rpc >= 0.5 & ssb_Schx1_rpc >= 0.5 & ssb_Rchx2_rpc >= 0.5 & ssb_Schx2_rpc >= 0.5 & 
                         ssb_Rchx_total >= 128 & ssb_Schx_total >= 128, c(1:2,27,34:35,78:79,82:83)]
View(ssb_fishers[(ssb_Rchx_total - ssb_Rchx_sum) < 0])
View(ssb_fishers[(ssb_Schx_total - ssb_Schx_sum) < 0])
ssb_fishers <- ssb_fishers[orf != "YDL184C" & orf != "YML129C" & orf != "YDR380W" & orf != "YEL048C" & orf != "YDL136W"]
for (i in 1:nrow(ssb_fishers)) {
  ssb <- matrix(c(ssb_fishers[i]$ssb_Schx_sum, ssb_fishers[i]$ssb_Rchx_sum, (ssb_fishers[i]$ssb_Schx_total - ssb_fishers[i]$ssb_Schx_sum), (ssb_fishers[i]$ssb_Rchx_total - ssb_fishers[i]$ssb_Rchx_sum)), nrow = 2)
  ssb_test <- fisher.test(ssb)
  ssb_fishers[i, ssb_odds := ssb_test$estimate]
  ssb_fishers[i, ssb_pvalue := ssb_test$p.value]
}
ssb_fishers[, ssb_padj := p.adjust(ssb_pvalue, method = "BH"), by = orf]
ssb_fishers[, ssb_odds_ma := movingAverage(ssb_odds, n=5, center=T), by = orf]
ssb_peaks_all <- ssb_fishers[ssb_odds > 1 & ssb_odds < Inf & ssb_padj < 0.05 & position > 30 & ssb_Schx_sum > ssb_Schx_rpc]


### Calculate fisher at each position for Bukau
bukau_fishers <- bukau_dt[B_ssb1_T_cor >= 0.5 & B_ssb1_I_cor >= 0.5 & 
                            B_ssb1_T1_rpc >= 0.5 & B_ssb1_T2_rpc >= 0.5 & B_ssb1_I1_rpc >= 0.5 & B_ssb1_I2_rpc >= 0.5 & 
                            B_ssb1_T_total >= 128 & B_ssb1_I_total >= 128, c(1:2,15:17,20:23)]
View(bukau_fishers[(B_ssb1_T_total - B_ssb1_T_sum) < 0])
View(bukau_fishers[(B_ssb1_I_total - B_ssb1_I_total) < 0])
bukau_fishers <- bukau_fishers[orf != "YJL047C-A"]
for (i in 1:nrow(bukau_fishers)) {
  ssb <- matrix(c(bukau_fishers[i]$B_ssb1_I_sum, bukau_fishers[i]$B_ssb1_T_sum, (bukau_fishers[i]$B_ssb1_I_total - bukau_fishers[i]$B_ssb1_I_sum), (bukau_fishers[i]$B_ssb1_T_total - bukau_fishers[i]$B_ssb1_T_sum)), nrow = 2)
  ssb_test <- fisher.test(ssb)
  bukau_fishers[i, ssb_odds := ssb_test$estimate]
  bukau_fishers[i, ssb_pvalue := ssb_test$p.value]
}
bukau_fishers[, ssb1_padj := p.adjust(ssb_pvalue, method = "BH"), by = orf]
bukau_fishers[, ssb1_odds_ma := movingAverage(ssb_odds, n=5, center=T), by = orf]
bukau_peaks_all <- bukau_fishers[ssb_odds > 1 & ssb_odds < Inf & ssb1_padj < 0.05 & position > 30 & B_ssb1_I_sum > B_ssb1_I_rpc]


### Bukau orfs
doring_fishers <- bukau_dt[bukau_dt$orf %in% bukau_dataset$orf]
unique(doring_fishers[(B_ssb1_T_total - B_ssb1_T_sum) < 0]$orf)
unique(doring_fishers[(B_ssb1_I_total - B_ssb1_I_total) < 0]$orf)
doring_fishers <- doring_fishers[orf != "YHR054C"]
unique(bukau_dataset[!bukau_dataset$orf %in% bukau_dt$orf])
for (i in 1:nrow(doring_fishers)) {
  ssb <- matrix(c(doring_fishers[i]$B_ssb1_I_sum, doring_fishers[i]$B_ssb1_T_sum, (doring_fishers[i]$B_ssb1_I_total - doring_fishers[i]$B_ssb1_I_sum), (doring_fishers[i]$B_ssb1_T_total - doring_fishers[i]$B_ssb1_T_sum)), nrow = 2)
  ssb_test <- fisher.test(ssb)
  doring_fishers[i, ssb_odds := ssb_test$estimate]
  doring_fishers[i, ssb_pvalue := ssb_test$p.value]
  print(i)
}
doring_fishers[, ssb1_padj := p.adjust(ssb_pvalue, method = "BH"), by = orf]
doring_fishers[, ssb1_odds_ma := movingAverage(ssb_odds, n=5, center=T), by = orf]


### Calculate fisher at each position in no crosslinking dataset
tric_noX_fishers <- tric_noX_dt[orf == "YFL039C"] # could try correlation plot of substrates identified in foundational dataset
tric_noX_fishers[, T2O_ma := movingAverage(T2O, n=5, center=T), by = orf]
tric_noX_fishers[, R2O_ma := movingAverage(R2O, n=5, center=T), by = orf]
for (i in 1:nrow(tric_noX_fishers_all)) {
  tric <- matrix(c(tric_noX_fishers[i]$T2O_ma, tric_noX_fishers[i]$R2O_ma, (tric_noX_fishers[i]$T2O_total - tric_noX_fishers[i]$T2O_ma), (tric_noX_fishers[i]$R2O_total - tric_noX_fishers[i]$R2O_ma)), nrow = 2)
  tric_test <- fisher.test(tric)
  tric_noX_fishers[i, tric_odds := tric_test$estimate]
  tric_noX_fishers[i, tric_pvalue := tric_test$p.value]
}
tric_noX_fishers[, tric_padj := p.adjust(tric_pvalue, method = "BH"), by = orf]
tric_noX_fishers[, tric_odds_ma := movingAverage(tric_odds, n=5, center=T), by = orf]

tric_X_fishers <- tric_dt[orf == "YFL039C"] # could try correlation plot of substrates identified in foundational dataset
tric_X_fishers[, tric_sum_ma := movingAverage(tric_sum, n=5, center=T), by = orf]
tric_X_fishers[, ribo_sum_ma := movingAverage(ribo_sum, n=5, center=T), by = orf]
for (i in 1:nrow(tric_X_fishers)) {
  tric <- matrix(c(tric_X_fishers[i]$tric_sum_ma, tric_X_fishers[i]$ribo_sum_ma, (tric_X_fishers[i]$tric_total - tric_X_fishers[i]$tric_sum_ma), (tric_X_fishers[i]$ribo_total - tric_X_fishers[i]$ribo_sum_ma)), nrow = 2)
  tric_test <- fisher.test(tric)
  tric_X_fishers[i, tric_odds := tric_test$estimate]
  tric_X_fishers[i, tric_pvalue := tric_test$p.value]
}
tric_X_fishers[, tric_padj := p.adjust(tric_pvalue, method = "BH"), by = orf]
tric_X_fishers[, tric_odds_ma := movingAverage(tric_odds, n=5, center=T), by = orf]


tric_noX_fishers_all <- tric_noX_dt[R2O_rpc >= 0.5 & T2O_rpc >= 0.5 & R2O_total >= 64 & T2O_total >= 64]
tric_noX_fishers_all <- tric_noX_fishers_all[tric_noX_fishers_all$orf %in% tric_substrates5$orf]
View(tric_noX_fishers_all[(T2O_total - T2O) < 0]) # Find orfs with position with aberrantly high number of reads
View(tric_noX_fishers_all[(R2O_total - R2O) < 0])
for (i in 1:nrow(tric_noX_fishers_all)) {
  print(i)
  tric <- matrix(c(tric_noX_fishers_all[i]$T2O, tric_noX_fishers_all[i]$R2O, (tric_noX_fishers_all[i]$T2O_total - tric_noX_fishers_all[i]$T2O), (tric_noX_fishers_all[i]$R2O_total - tric_noX_fishers_all[i]$R2O)), nrow = 2)
  tric_test <- fisher.test(tric)
  tric_noX_fishers_all[i, tric_odds := tric_test$estimate]
  tric_noX_fishers_all[i, tric_pvalue := tric_test$p.value]
}
tric_noX_fishers_all[, tric_padj := p.adjust(tric_pvalue, method = "BH"), by = orf]
tric_noX_fishers_all[, tric_odds_ma := movingAverage(tric_odds, n=5, center=T), by = orf]
tric_noX_peaks_all <- tric_noX_fishers_all[tric_odds > 1 & tric_odds < Inf & tric_padj < 0.05 & position > 30 & T2O > T2O_rpc]

orfs <- tric_noX_peaks_all[, .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
peak_start <- NULL
peak_end <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- tric_noX_peaks_all[g]$tric_odds
  peaks <- tric_noX_peaks_all[g]$position
  while (length(peaks) > 0) {
    y=peaks[1]
    while (y %in% peaks) {
      peaks_subset1 <- y
      peaks_subset <- c(peaks_subset,peaks_subset1)
      y <- y+1
    } # outputs vector of consecutive positions with odds_ma > 75%+1.5*IQR and padj < 0.05
    if (length(peaks_subset) >= 5) {
      odds_subset <- odds[which(peaks %in% peaks_subset)]
      final_peaks1 <- peaks_subset[which.max(odds_subset)] # determines position of max enrichment in that sub-vector
      final_peaks <- c(final_peaks, final_peaks1) # adds peak to new vector
      final_odds1 <- odds_subset[which.max(odds_subset)] # determines position of max enrichment in that sub-vector
      final_odds <- c(final_odds, final_odds1) # adds peak to new vector
      peak_start1 <- peaks_subset[1] 
      peak_start <- c(peak_start, peak_start1)
      peak_end1 <- peaks_subset[length(peaks_subset)] 
      peak_end <- c(peak_end, peak_end1)
      odds <- odds[which(!peaks %in% peaks_subset)]
      peaks <- peaks[!peaks %in% peaks_subset]
      peaks_subset <- NULL
      odds_subset <- NULL
      gene <- c(gene, g)
    } else {
      odds <- odds[which(!peaks %in% peaks_subset)]
      peaks <- peaks[!peaks %in% peaks_subset]
      peaks_subset <- NULL
    }
  }
}
tric_noX_peaks5 <- data.table(orf = gene, peak_start = peak_start, peak_end = peak_end,
                              peak = final_peaks,
                              odds = final_odds)
setkeyv(tric_noX_peaks5, c("orf"))

orfs <- tric_noX_peaks5[, .SD[which.min(peak)], by = orf]
peaks1 <- NULL # eliminate "peaks" within 5aa of each other that were simply due to a residue or two in a sequence not being significant
gene1 <- NULL
for (g in orfs$orf) {
  odds <- tric_noX_peaks5[g]$odds
  peaks <- tric_noX_peaks5[g]$peak
  for (i in 1:length(peaks)) {
    bottom <- peaks[i] - 5
    top <- peaks[i] + 5
    new_peaks <- peaks[which(peaks < top & peaks > bottom)]
    new_odds <- odds[which(peaks < top & peaks > bottom)]
    temp1 <- peaks[which(odds == new_odds[which.max(new_odds)])]
    temp2 <- odds[which(odds == new_odds[which.max(new_odds)])]
    new_peaks <- new_peaks[!new_peaks %in% temp1]
    new_odds <- new_odds[!new_odds %in% temp2]
    peaks <- peaks[!peaks %in% new_peaks]
    odds <- odds[!odds %in% new_odds]
  }
  gene <- rep(g, length(peaks))  
  peaks1 <- c(peaks1, peaks)
  gene1 <- c(gene1, gene)
}
tric_noX_peaks5_temp <- data.table(orf = gene1,
                                   peak = peaks1)
tric_noX_peaks5[, position1 := as.character(base::paste(orf, peak, sep = "_"))]
tric_noX_peaks5_temp[, position1 := as.character(base::paste(orf, peak, sep = "_"))]
tric_noX_peaks5 <- tric_noX_peaks5[tric_noX_peaks5$position1 %in% tric_noX_peaks5_temp$position1]
setkeyv(tric_noX_peaks5, c("orf"))


### Calculate fisher at all SSB and TRiC substrates
substrates_fishers <- tric_dt[tric_dt$orf %in% substrates$orf]
for (i in 1:nrow(substrates_fishers)) {
  tric <- matrix(c(substrates_fishers[i]$tric_sum, substrates_fishers[i]$ribo_sum, (substrates_fishers[i]$tric_total - substrates_fishers[i]$tric_sum), (substrates_fishers[i]$ribo_total - substrates_fishers[i]$ribo_sum)), nrow = 2)
  tric_test <- fisher.test(tric)
  substrates_fishers[i, tric_odds := tric_test$estimate]
  substrates_fishers[i, tric_pvalue := tric_test$p.value]
  ssb <- matrix(c(substrates_fishers[i]$ssb_Schx_sum, substrates_fishers[i]$ssb_Rchx_sum, (substrates_fishers[i]$ssb_Schx_total - substrates_fishers[i]$ssb_Schx_sum), (substrates_fishers[i]$ssb_Rchx_total - substrates_fishers[i]$ssb_Rchx_sum)), nrow = 2)
  ssb_test <- fisher.test(ssb)
  substrates_fishers[i, ssb_odds := ssb_test$estimate]
  substrates_fishers[i, ssb_pvalue := ssb_test$p.value]
}
substrates_fishers[, tric_padj := p.adjust(tric_pvalue, method = "BH"), by = orf]
substrates_fishers[, tric_odds_ma := movingAverage(tric_odds, n=5, center=T), by = orf]
substrates_fishers[, ssb_padj := p.adjust(ssb_pvalue, method = "BH"), by = orf]
substrates_fishers[, ssb_odds_ma := movingAverage(ssb_odds, n=5, center=T), by = orf]


### Calculate fisher at each position for Atp2 dataset
atp2_fishers <- atp2_dt[orf == "YJR121W" | orf == "YFL039C"] # could try correlation plot of substrates identified in foundational dataset
atp2_fishers <- atp2_dt[orf == "YGR281W"]
for (i in 1:nrow(atp2_fishers)) {
  atp2 <- matrix(c(atp2_fishers[i]$tric_sum, atp2_fishers[i]$ribo_sum, (atp2_fishers[i]$tric_total - atp2_fishers[i]$tric_sum), (atp2_fishers[i]$ribo_total - atp2_fishers[i]$ribo_sum)), nrow = 2)
  atp2_test <- fisher.test(atp2)
  atp2_fishers[i, atp2_odds := atp2_test$estimate]
  atp2_fishers[i, atp2_pvalue := atp2_test$p.value]
}
atp2_fishers[, atp2_padj := p.adjust(atp2_pvalue, method = "BH"), by = orf]
atp2_fishers[, atp2_odds_ma := movingAverage(atp2_odds, n=5, center=T), by = orf]


### Eliminate positions with fewer than 2, 3, or 5 consecutive codons
# TRiC
orfs <- tric_peaks_all[, .SD[which.min(position)], by = orf]
#orfs <- temp1[orf == "YML124C", .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
peak_start <- NULL
peak_end <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- tric_peaks_all[g]$tric_odds
  peaks <- tric_peaks_all[g]$position
  while (length(peaks) > 0) {
    y=peaks[1]
    while (y %in% peaks) {
      peaks_subset1 <- y
      peaks_subset <- c(peaks_subset,peaks_subset1)
      y <- y+1
    } # outputs vector of consecutive positions with odds_ma > 75%+1.5*IQR and padj < 0.05
    if (length(peaks_subset) >= 5) {
      odds_subset <- odds[which(peaks %in% peaks_subset)]
      final_peaks1 <- peaks_subset[which.max(odds_subset)] # determines position of max enrichment in that sub-vector
      final_peaks <- c(final_peaks, final_peaks1) # adds peak to new vector
      final_odds1 <- odds_subset[which.max(odds_subset)] # determines position of max enrichment in that sub-vector
      final_odds <- c(final_odds, final_odds1) # adds peak to new vector
      peak_start1 <- peaks_subset[1] 
      peak_start <- c(peak_start, peak_start1)
      peak_end1 <- peaks_subset[length(peaks_subset)] 
      peak_end <- c(peak_end, peak_end1)
      odds <- odds[which(!peaks %in% peaks_subset)]
      peaks <- peaks[!peaks %in% peaks_subset]
      peaks_subset <- NULL
      odds_subset <- NULL
      gene <- c(gene, g)
    } else {
        odds <- odds[which(!peaks %in% peaks_subset)]
        peaks <- peaks[!peaks %in% peaks_subset]
        peaks_subset <- NULL
    }
  }
}
tric_peaks5 <- data.table(orf = gene, peak_start = peak_start, peak_end = peak_end,
                    peak = final_peaks,
                    odds = final_odds)
setkeyv(tric_peaks5, c("orf"))

orfs <- tric_peaks5[, .SD[which.min(peak)], by = orf]
peaks1 <- NULL # eliminate "peaks" within 5aa of each other that were simply due to a residue or two in a sequence not being significant
gene1 <- NULL
for (g in orfs$orf) {
  odds <- tric_peaks5[g]$odds
  peaks <- tric_peaks5[g]$peak
  for (i in 1:length(peaks)) {
    bottom <- peaks[i] - 5
    top <- peaks[i] + 5
    new_peaks <- peaks[which(peaks < top & peaks > bottom)]
    new_odds <- odds[which(peaks < top & peaks > bottom)]
    temp1 <- peaks[which(odds == new_odds[which.max(new_odds)])]
    temp2 <- odds[which(odds == new_odds[which.max(new_odds)])]
    new_peaks <- new_peaks[!new_peaks %in% temp1]
    new_odds <- new_odds[!new_odds %in% temp2]
    peaks <- peaks[!peaks %in% new_peaks]
    odds <- odds[!odds %in% new_odds]
  }
  gene <- rep(g, length(peaks))  
  peaks1 <- c(peaks1, peaks)
  gene1 <- c(gene1, gene)
}
tric_peaks5_temp <- data.table(orf = gene1,
                             peak = peaks1)
tric_peaks5[, position1 := as.character(base::paste(orf, peak, sep = "_"))]
tric_peaks5_temp[, position1 := as.character(base::paste(orf, peak, sep = "_"))]
tric_peaks5 <- tric_peaks5[tric_peaks5$position1 %in% tric_peaks5_temp$position1]
setkeyv(tric_peaks5, c("orf"))


# SSB
orfs <- ssb_peaks_all[, .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
peak_start <- NULL
peak_end <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- ssb_peaks_all[g]$ssb_odds
  peaks <- ssb_peaks_all[g]$position
  while (length(peaks) > 0) {
    y=peaks[1]
    while (y %in% peaks) {
      peaks_subset1 <- y
      peaks_subset <- c(peaks_subset,peaks_subset1)
      y <- y+1
    } # outputs vector of consecutive positions with odds_ma > 75%+1.5*IQR and padj < 0.05
    if (length(peaks_subset) >= 5) {
      odds_subset <- odds[which(peaks %in% peaks_subset)]
      final_peaks1 <- peaks_subset[which.max(odds_subset)] # determines position of max enrichment in that sub-vector
      final_peaks <- c(final_peaks, final_peaks1) # adds peak to new vector
      final_odds1 <- odds_subset[which.max(odds_subset)] # determines position of max enrichment in that sub-vector
      final_odds <- c(final_odds, final_odds1) # adds peak to new vector
      peak_start1 <- peaks_subset[1] 
      peak_start <- c(peak_start, peak_start1)
      peak_end1 <- peaks_subset[length(peaks_subset)] 
      peak_end <- c(peak_end, peak_end1)
      odds <- odds[which(!peaks %in% peaks_subset)]
      peaks <- peaks[!peaks %in% peaks_subset]
      peaks_subset <- NULL
      odds_subset <- NULL
      gene <- c(gene, g)
    } else {
      odds <- odds[which(!peaks %in% peaks_subset)]
      peaks <- peaks[!peaks %in% peaks_subset]
      peaks_subset <- NULL
    }
  }
}
ssb_peaks5 <- data.table(orf = gene, peak_start = peak_start, peak_end = peak_end,
                         peak = final_peaks,
                         odds = final_odds)
setkeyv(ssb_peaks5, c("orf"))

orfs <- ssb_peaks5[, .SD[which.min(peak)], by = orf]
peaks1 <- NULL # eliminate "peaks" within 5aa of each other that were simply due to a residue or two in a sequence not being significant
gene1 <- NULL
for (g in orfs$orf) {
  odds <- ssb_peaks5[g]$odds
  peaks <- ssb_peaks5[g]$peak
  for (i in 1:length(peaks)) {
    bottom <- peaks[i] - 5
    top <- peaks[i] + 5
    new_peaks <- peaks[which(peaks < top & peaks > bottom)]
    new_odds <- odds[which(peaks < top & peaks > bottom)]
    temp1 <- peaks[which(odds == new_odds[which.max(new_odds)])]
    temp2 <- odds[which(odds == new_odds[which.max(new_odds)])]
    new_peaks <- new_peaks[!new_peaks %in% temp1]
    new_odds <- new_odds[!new_odds %in% temp2]
    peaks <- peaks[!peaks %in% new_peaks]
    odds <- odds[!odds %in% new_odds]
  }
  gene <- rep(g, length(peaks))  
  peaks1 <- c(peaks1, peaks)
  gene1 <- c(gene1, gene)
}
ssb_peaks5_temp <- data.table(orf = gene1,
                          peak = peaks1)
ssb_peaks5[, position1 := as.character(base::paste(orf, peak, sep = "_"))]
ssb_peaks5_temp[, position1 := as.character(base::paste(orf, peak, sep = "_"))]
ssb_peaks5 <- ssb_peaks5[ssb_peaks5$position1 %in% ssb_peaks5_temp$position1]
setkeyv(ssb_peaks5, c("orf"))


# Bukau
orfs <- bukau_peaks_all[, .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
peak_start <- NULL
peak_end <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- bukau_peaks_all[g]$ssb1_odds
  peaks <- bukau_peaks_all[g]$position
  while (length(peaks) > 0) {
    y=peaks[1]
    while (y %in% peaks) {
      peaks_subset1 <- y
      peaks_subset <- c(peaks_subset,peaks_subset1)
      y <- y+1
    } # outputs vector of consecutive positions with odds_ma > 75%+1.5*IQR and padj < 0.05
    if (length(peaks_subset) >= 5) {
      odds_subset <- odds[which(peaks %in% peaks_subset)]
      final_peaks1 <- peaks_subset[which.max(odds_subset)] # determines position of max enrichment in that sub-vector
      final_peaks <- c(final_peaks, final_peaks1) # adds peak to new vector
      final_odds1 <- odds_subset[which.max(odds_subset)] # determines position of max enrichment in that sub-vector
      final_odds <- c(final_odds, final_odds1) # adds peak to new vector
      peak_start1 <- peaks_subset[1] 
      peak_start <- c(peak_start, peak_start1)
      peak_end1 <- peaks_subset[length(peaks_subset)] 
      peak_end <- c(peak_end, peak_end1)
      odds <- odds[which(!peaks %in% peaks_subset)]
      peaks <- peaks[!peaks %in% peaks_subset]
      peaks_subset <- NULL
      odds_subset <- NULL
      gene <- c(gene, g)
    } else {
      odds <- odds[which(!peaks %in% peaks_subset)]
      peaks <- peaks[!peaks %in% peaks_subset]
      peaks_subset <- NULL
    }
  }
}
bukau_peaks5 <- data.table(orf = gene, peak_start = peak_start, peak_end = peak_end,
                           peak = final_peaks,
                           odds = final_odds)
setkeyv(bukau_peaks5, c("orf"))

orfs <- bukau_peaks5[, .SD[which.min(peak)], by = orf]
peaks1 <- NULL # eliminate "peaks" within 5aa of each other that were simply due to a residue or two in a sequence not being significant
gene1 <- NULL
for (g in orfs$orf) {
  odds <- bukau_peaks5[g]$odds
  peaks <- bukau_peaks5[g]$peak
  for (i in 1:length(peaks)) {
    bottom <- peaks[i] - 5
    top <- peaks[i] + 5
    new_peaks <- peaks[which(peaks < top & peaks > bottom)]
    new_odds <- odds[which(peaks < top & peaks > bottom)]
    temp1 <- peaks[which(odds == new_odds[which.max(new_odds)])]
    temp2 <- odds[which(odds == new_odds[which.max(new_odds)])]
    new_peaks <- new_peaks[!new_peaks %in% temp1]
    new_odds <- new_odds[!new_odds %in% temp2]
    peaks <- peaks[!peaks %in% new_peaks]
    odds <- odds[!odds %in% new_odds]
  }
  gene <- rep(g, length(peaks))  
  peaks1 <- c(peaks1, peaks)
  gene1 <- c(gene1, gene)
}
bukau_peaks5_temp <- data.table(orf = gene1,
                          peak = peaks1)
bukau_peaks5[, position1 := as.character(base::paste(orf, peak, sep = "_"))]
bukau_peaks5_temp[, position1 := as.character(base::paste(orf, peak, sep = "_"))]
bukau_peaks5 <- bukau_peaks5[bukau_peaks5$position1 %in% bukau_peaks5_temp$position1]
setkeyv(bukau_peaks5, c("orf"))


### Prepare data tables: 1) list of substrates, 2) just peaks with relevant data, 3) occupancy data for length of orf
tric_dt[, position_norm := position / length]
bukau_dt[, position_norm := position / length]

tric_peaks3 <- tric_peaks3[orf != "YJL014W" & orf != "YJR064W" & orf != "YDR188W" & orf != "YJL008C" & orf != "YDR212W" & orf != "YIL142W" & orf != "YDL143W" & orf != "YJL111W"]
tric_substrates3 <- as.data.table(unique(tric_peaks3[, 1])) # make data table of just the substrates
tric_peaks3[, PeaksPerOrf := length(which(peak > 0)), by = orf] # add number of peaks in each orf
ssb_peaks3 <- ssb_peaks3[orf != "YDL229W" & orf != "YNL209W"]
ssb_substrates3 <- as.data.table(unique(ssb_peaks3[, 1])) # make data table of just the substrates
ssb_peaks3[, PeaksPerOrf := length(which(peak > 0)), by = orf] # add number of peaks in each orf
bukau_peaks3 <- bukau_peaks3[orf != "YDL229W" & orf != "YNL209W"]
bukau_substrates3 <- as.data.table(unique(bukau_peaks3[, 1])) # make data table of just the substrates
bukau_peaks3[, PeaksPerOrf := length(which(peak > 0)), by = orf] # add number of peaks in each orf


## TRiC
# Data table of peaks with relevant data
setkeyv(tric_peaks3, c("position1")) 
setkeyv(tric_fishers, c("position1"))
temp <- tric_fishers[, c(1,3,10:13)] 
temp2 <- temp[tric_peaks3]
tric_peaks3 <- temp2[, c(1:6,8:10,12)]
setkeyv(tric_dt, c("position1"))
temp3 <- tric_dt[tric_dt$position1 %in% tric_peaks3$position1]
tric_peaks3 <- cbind(tric_peaks3, temp3[, c(3,24:26,42:43,48:51,168)])
setkeyv(tric_peaks3, c("orf","peak"))
setkeyv(tric_fishers, c("orf", "position"))
setkeyv(tric_dt, c("orf", "position"))
tric_peaks3_max <- tric_peaks3[, .SD[which.max(tric_odds)], by = orf]

# Data table of peaks with occupancy and padj for length of orf of substrates
temp4 <- tric_dt[tric_dt$orf %in% tric_substrates3$orf]
temp5 <- tric_fishers[tric_fishers$orf %in% tric_substrates3$orf, c(1,2,10:13)]
tric_substrates3_dt <- temp5[temp4]


## SSB
# Data table of peaks with relevant data
setkeyv(ssb_peaks3, c("position1")) 
setkeyv(ssb_fishers, c("position1"))
temp <- ssb_fishers[, c(1,3,10:13)] 
temp2 <- temp[ssb_peaks3]
ssb_peaks3 <- temp2[, c(1:6,8:10,12)]
setkeyv(tric_dt, c("position1"))
temp3 <- tric_dt[tric_dt$position1 %in% ssb_peaks3$position1]
ssb_peaks3 <- cbind(ssb_peaks3, temp3[, c(3,24:26,82:83,86:89,168)])
setkeyv(ssb_peaks3, c("orf","peak"))
setkeyv(ssb_fishers, c("orf", "position"))
setkeyv(tric_dt, c("orf", "position"))
ssb_peaks3_max <- ssb_peaks3[, .SD[which.max(ssb_odds)], by = orf]

# Data table of peaks with occupancy and padj for length of orf of substrates
temp4 <- tric_dt[tric_dt$orf %in% ssb_substrates5$orf]
temp5 <- ssb_fishers[ssb_fishers$orf %in% ssb_substrates5$orf, c(1,2,10:13)]
ssb_substrates5_dt <- temp5[temp4]


## Bukau
# Data table of peaks with relevant data
setkeyv(bukau_peaks3, c("position1")) 
setkeyv(bukau_fishers, c("position1"))
temp <- bukau_fishers[, c(1,3,10:13)] 
temp2 <- temp[bukau_peaks3]
bukau_peaks3 <- temp2[, c(1:6,8:10,12)]
setkeyv(bukau_dt, c("position1"))
temp3 <- bukau_dt[bukau_dt$position1 %in% bukau_peaks3$position1]
bukau_peaks3 <- cbind(bukau_peaks3, temp3[, c(3,12:14,22:27,60)])
setkeyv(bukau_peaks3, c("orf","peak"))
setkeyv(bukau_fishers, c("orf", "position"))
setkeyv(bukau_dt, c("orf", "position"))
bukau_peaks3_max <- bukau_peaks3[, .SD[which.max(ssb_odds)], by = orf]

# Data table of peaks with occupancy and padj for length of orf of substrates
temp4 <- bukau_dt[bukau_dt$orf %in% bukau_substrates3$orf]
temp5 <- bukau_fishers[bukau_fishers$orf %in% bukau_substrates3$orf, c(1,2,10:13)]
bukau_substrates3_dt <- temp5[temp4]



### Isolate all residues within identified peaks
# TRiC
orfs <- tric_peaks_all[, .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
final_padj <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- tric_peaks_all[g]$tric_odds
  peaks <- tric_peaks_all[g]$position
  padj <- tric_peaks_all[g]$tric_padj
  while (length(peaks) > 0) {
    y=peaks[1]
    while (y %in% peaks) {
      peaks_subset1 <- y
      peaks_subset <- c(peaks_subset,peaks_subset1)
      y <- y+1
    } # outputs vector of consecutive positions with odds_ma > 75%+1.5*IQR and padj < 0.05
    if (length(peaks_subset) >= 5) {
      odds_subset <- odds[which(peaks %in% peaks_subset)]
      padj_subset <- padj[which(peaks %in% peaks_subset)]
      gene1 <- rep(g, length(peaks_subset))
      final_peaks <- c(final_peaks, peaks_subset)
      final_odds <- c(final_odds, odds_subset)
      final_padj <- c(final_padj, padj_subset)
      gene <- c(gene, gene1)
      odds <- odds[which(!peaks %in% peaks_subset)]
      padj <- padj[which(!peaks %in% peaks_subset)]
      peaks <- peaks[!peaks %in% peaks_subset]
      peaks_subset <- NULL
      odds_subset <- NULL
    } else {
      odds <- odds[which(!peaks %in% peaks_subset)]
      peaks <- peaks[!peaks %in% peaks_subset]
      peaks_subset <- NULL
    }
  }
}
tric_peaks5all <- data.table(orf = gene, peak = final_peaks, 
                            tric_odds = final_odds, tric_padj = final_padj)
setkeyv(tric_peaks5all, c("orf"))
tric_peaks5all[, position1 := as.character(base::paste(orf, peak, sep = "_"))]


# SSB
orfs <- ssb_peaks_all[, .SD[which.min(position)], by = orf]
peaks_subset <- NULL 
odds_subset <- NULL
final_peaks <- NULL
final_odds <- NULL
final_padj <- NULL
gene <- NULL
for (g in orfs$orf) {
  odds <- ssb_peaks_all[g]$ssb_odds
  peaks <- ssb_peaks_all[g]$position
  padj <- ssb_peaks_all[g]$ssb_padj
  while (length(peaks) > 0) {
    y=peaks[1]
    while (y %in% peaks) {
      peaks_subset1 <- y
      peaks_subset <- c(peaks_subset,peaks_subset1)
      y <- y+1
    } # outputs vector of consecutive positions with odds_ma > 75%+1.5*IQR and padj < 0.05
    if (length(peaks_subset) >= 5) {
      odds_subset <- odds[which(peaks %in% peaks_subset)]
      padj_subset <- padj[which(peaks %in% peaks_subset)]
      gene1 <- rep(g, length(peaks_subset))
      final_peaks <- c(final_peaks, peaks_subset)
      final_odds <- c(final_odds, odds_subset)
      final_padj <- c(final_padj, padj_subset)
      gene <- c(gene, gene1)
      odds <- odds[which(!peaks %in% peaks_subset)]
      padj <- padj[which(!peaks %in% peaks_subset)]
      peaks <- peaks[!peaks %in% peaks_subset]
      peaks_subset <- NULL
      odds_subset <- NULL
    } else {
      odds <- odds[which(!peaks %in% peaks_subset)]
      peaks <- peaks[!peaks %in% peaks_subset]
      peaks_subset <- NULL
    }
  }
}
ssb_peaks5all <- data.table(orf = gene, peak = final_peaks, 
                            ssb_odds = final_odds, ssb_padj = final_padj)
setkeyv(ssb_peaks5all, c("orf"))
ssb_peaks5all[, position1 := as.character(base::paste(orf, peak, sep = "_"))]

