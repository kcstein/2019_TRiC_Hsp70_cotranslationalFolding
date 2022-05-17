### Isolate SSB substrates and binding sites from Bukau dataset
bukau_dataset <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/BukauDataset_Ssb1_ssb2KO.csv", header = T, stringsAsFactors = T)
bukau_dataset <- as.data.table(bukau_dataset)
setkeyv(bukau_dataset, c("orf"))

bukau_orfs1 <- NULL
bukau_cor1 <- NULL
bukau_peaks1 <- NULL
peak_start1 <- NULL
peak_width1 <- NULL
for (i in 1:nrow(bukau_dataset)) {
  orf <- bukau_dataset[i]$orf
  correlation <- bukau_dataset[i]$pearson.correlation
  peaks <- bukau_dataset[i]$number.of.peaks
  bukau_orfs <- as.character(rep(orf, each = peaks))
  bukau_cor <- rep(correlation, each = peaks)
  bukau_peaks <- rep(peaks, each = peaks)
  bukau_orfs1 <- c(bukau_orfs1, bukau_orfs)
  bukau_cor1 <- c(bukau_cor1, bukau_cor)
  bukau_peaks1 <- c(bukau_peaks1, bukau_peaks)
  for (j in 4:ncol(bukau_dataset)) {
    if (j %% 2 == 0) {
      if (!is.na(bukau_dataset[i,..j])) {
        peak_width <- bukau_dataset[i,..j]
        peak_width1 <- c(peak_width1, peak_width)  
      }
    } 
    if (j %% 2 == 1) {
      if (!is.na(bukau_dataset[i,..j])) {  
        peak_start <- bukau_dataset[i,..j]
        peak_start1 <- c(peak_start1, peak_start)
      }
    }
  }
}

bukau_orf <- as.character(bukau_orfs1)
bukau_cor <- as.numeric(bukau_cor1)
bukau_peaks <- as.integer(bukau_peaks1)
bukau_width <- as.integer(peak_width1)
bukau_start <- as.integer(peak_start1)
bukau_dataset_dt <- data.table(orf = bukau_orf,
                               correlation = bukau_cor,
                               peaks = bukau_peaks,
                               peak_start = bukau_start,
                               width = bukau_width)
names(bukau_dataset_dt) <- c("orf", "correlation", "peak_number", "peak_start", "width")
setkeyv(bukau_dataset_dt, c("orf"))
bukau_dataset_dt[, peak := floor(floor(peak_start/3) + floor((floor(width/3) / 2)))]
bukau_dataset_dt[, position1 := as.character(base::paste(orf, peak, sep = "_"))]
bukau_dataset_dt[, peak_startaa := ceiling((peak_start / 3))]
bukau_dataset_dt[, peak_endaa := floor(((peak_start + width) / 3))]

orfs <- translatome_dt[, .SD[which.min(length)], by = orf]
orfs <- data.table(orf = orfs$orf,
                   length = orfs$length)
setkeyv(orfs, c("orf"))
i <- cbind(match(bukau_dataset_dt$orf, orfs$orf))
bukau_dataset_dt <- cbind(bukau_dataset_dt, length=orfs[i]$length)

doring_substrates <- as.data.table(as.character(unique(bukau_dataset_dt$orf)))
names(doring_substrates) <- c("orf")
setkeyv(doring_substrates, c("orf"))


bukau_datasetWT <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Analysis/BukauDataset_Ssb1_WT.csv", header = T, stringsAsFactors = T)
bukau_datasetWT <- as.data.table(bukau_datasetWT)
setkeyv(bukau_datasetWT, c("orf"))

bukau_orfs1 <- NULL
bukau_cor1 <- NULL
bukau_peaks1 <- NULL
peak_start1 <- NULL
peak_width1 <- NULL
for (i in 1:nrow(bukau_datasetWT)) {
  orf <- bukau_datasetWT[i]$orf
  correlation <- bukau_datasetWT[i]$pearson.correlation
  peaks <- bukau_datasetWT[i]$number.of.peaks
  bukau_orfs <- as.character(rep(orf, each = peaks))
  bukau_cor <- rep(correlation, each = peaks)
  bukau_peaks <- rep(peaks, each = peaks)
  bukau_orfs1 <- c(bukau_orfs1, bukau_orfs)
  bukau_cor1 <- c(bukau_cor1, bukau_cor)
  bukau_peaks1 <- c(bukau_peaks1, bukau_peaks)
  for (j in 4:ncol(bukau_datasetWT)) {
    if (j %% 2 == 0) {
      if (!is.na(bukau_datasetWT[i,..j])) {
        peak_width <- bukau_datasetWT[i,..j]
        peak_width1 <- c(peak_width1, peak_width)  
      }
    } 
    if (j %% 2 == 1) {
      if (!is.na(bukau_datasetWT[i,..j])) {  
        peak_start <- bukau_datasetWT[i,..j]
        peak_start1 <- c(peak_start1, peak_start)
      }
    }
  }
}

bukau_orf <- as.character(bukau_orfs1)
bukau_cor <- as.numeric(bukau_cor1)
bukau_peaks <- as.integer(bukau_peaks1)
bukau_width <- as.integer(peak_width1)
bukau_start <- as.integer(peak_start1)
bukau_datasetWT_dt <- data.table(orf = bukau_orf,
                                 correlation = bukau_cor,
                                 peaks = bukau_peaks,
                                 peak_start = bukau_start,
                                 width = bukau_width)
names(bukau_datasetWT_dt) <- c("orf", "correlation", "peak_number", "peak_start", "width")
setkeyv(bukau_datasetWT_dt, c("orf"))
bukau_datasetWT_dt[, peak := floor(floor(peak_start/3) + floor((floor(width/3) / 2)))]
bukau_datasetWT_dt[, position1 := as.character(base::paste(orf, peak, sep = "_"))]
bukau_datasetWT_dt[, peak_startaa := ceiling((peak_start / 3))]
bukau_datasetWT_dt[, peak_endaa := floor(((peak_start + width) / 3))]

orfs <- translatome_dt[, .SD[which.min(length)], by = orf]
orfs <- data.table(orf = orfs$orf,
                   length = orfs$length)
setkeyv(orfs, c("orf"))
i <- cbind(match(bukau_datasetWT_dt$orf, orfs$orf))
bukau_datasetWT_dt <- cbind(bukau_datasetWT_dt, length=orfs[i]$length)

doringWT_substrates <- as.data.table(as.character(unique(bukau_datasetWT_dt$orf)))
names(doringWT_substrates) <- c("orf")
setkeyv(doringWT_substrates, c("orf"))


### Compare association sites between Bukau's dataset and ours
# Enrichment of shared substrates
temp <- ssb_substrates5[ssb_substrates5$orf %in% bukau_substrates5$orf]
temp1 <- ssb_substrates5_dt[ssb_substrates5_dt$orf %in% temp$orf]
i <- cbind(match(temp1$position1, bukau_substrates5_dt$position1))
temp1 <- cbind(temp1, ssb1_odds = bukau_substrates5_dt[i]$ssb_odds)
temp1 <- cbind(temp1, ssb1_odds_ma = bukau_substrates5_dt[i]$ssb_odds_ma)

temp <- ssb_peaks5[, .SD[which.min(peak)], by = orf]
temp1 <- bukau_peaks5[, .SD[which.min(peak)], by = orf]
temp2 <- temp[temp$orf %in% temp1$orf]
temp3 <- ssb_substrates5_dt[ssb_substrates5_dt$orf %in% temp2$orf]
temp4 <- bukau_substrates5_dt[bukau_substrates5_dt$orf %in% temp2$orf]


ggplot(temp3, aes(log2(ssb_odds_ma), log2(temp4$ssb1_odds_ma))) + 
  geom_hline(yintercept = 0, color = 'gray50', linetype = 'dashed', size = 1) + geom_vline(xintercept = 0, color = 'gray50', linetype = 'dashed', size = 1) +
  geom_point(color = "black", fill = "black", alpha = 0.5, size = 2) + ylim(-250,250) + xlim(-250,250)

