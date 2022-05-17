### Prepare table of WD proteins
wd_repeats <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Annotation/WDrepeats.csv", header = T)
wd_repeats <- as.data.table(wd_repeats)
setkeyv(wd_repeats, c("orf", "start"))
wd_repeats[, repeats := length(which(start > 0)), by = orf]
wd_repeats[, orf_wd := paste0(orf, "_", repeat.)]
wd_repeats[, start30 := start + 30]
wd_repeats[, end30 := end + 30]
temp <- wd_repeats[, .SD[which.min(start)], by = orf]

adjusted_end <- NULL
for (i in temp$orf) {
  temp1 <- wd_repeats[i]$end
  tempstart <- temp1[1]
  for (j in 1:length(temp1)) {
    adjusted_start1 <- temp1[j] - tempstart
    adjusted_end <- c(adjusted_end, adjusted_start1)
  }
}
wd_repeats$adjusted_end <- adjusted_end

wd_repeats[, start_avg := ifelse(repeat. == "WD2", round(mean(wd_repeats[repeat. == "WD2"]$adjusted_start)), adjusted_start)]
wd_repeats[, start_avg := ifelse(repeat. == "WD3", round(mean(wd_repeats[repeat. == "WD3"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD4", round(mean(wd_repeats[repeat. == "WD4"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD5", round(mean(wd_repeats[repeat. == "WD5"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD6", round(mean(wd_repeats[repeat. == "WD6"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD7", round(mean(wd_repeats[repeat. == "WD7"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD8", round(mean(wd_repeats[repeat. == "WD8"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD9", round(mean(wd_repeats[repeat. == "WD9"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD10", round(mean(wd_repeats[repeat. == "WD10"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD11", round(mean(wd_repeats[repeat. == "WD11"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD12", round(mean(wd_repeats[repeat. == "WD12"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD13", round(mean(wd_repeats[repeat. == "WD13"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD14", round(mean(wd_repeats[repeat. == "WD14"]$adjusted_start)), start_avg)]
wd_repeats[, start_avg := ifelse(repeat. == "WD15", round(mean(wd_repeats[repeat. == "WD15"]$adjusted_start)), start_avg)]

wd_repeats[, repeatlength := end - start]
wd_repeats[, length_avg := ifelse(repeat. == "WD1", round(mean(wd_repeats[repeat. == "WD1"]$repeatlength)), 0)]
wd_repeats[, length_avg := ifelse(repeat. == "WD2", round(mean(wd_repeats[repeat. == "WD2"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD3", round(mean(wd_repeats[repeat. == "WD3"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD4", round(mean(wd_repeats[repeat. == "WD4"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD5", round(mean(wd_repeats[repeat. == "WD5"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD6", round(mean(wd_repeats[repeat. == "WD6"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD7", round(mean(wd_repeats[repeat. == "WD7"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD8", round(mean(wd_repeats[repeat. == "WD8"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD9", round(mean(wd_repeats[repeat. == "WD9"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD10", round(mean(wd_repeats[repeat. == "WD10"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD11", round(mean(wd_repeats[repeat. == "WD11"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD12", round(mean(wd_repeats[repeat. == "WD12"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD13", round(mean(wd_repeats[repeat. == "WD13"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD14", round(mean(wd_repeats[repeat. == "WD14"]$repeatlength)), length_avg)]
wd_repeats[, length_avg := ifelse(repeat. == "WD15", round(mean(wd_repeats[repeat. == "WD15"]$repeatlength)), length_avg)]

