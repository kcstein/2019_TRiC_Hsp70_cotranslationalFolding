### Analyze ribosome occupancy of translatome around binding sites
# TRiC
tric_translatome_dt <- translatome_dt[tric_substrates5]

temp <- tric_peaks5
random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(30:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

tric_translatome_dt <- tric_translatome_dt[temp, allow.cartesian = TRUE]
tric_translatome_dt[, name := as.character(base::paste(orf, peak, sep = "_"))]
tric_translatome_dt[, adjusted := position - peak]
tric_translatome_dt[, adjusted_start := position - peak_start]
tric_translatome_dt[, adjusted_random := position - random]
tric_translatome_dt <- tric_translatome_dt[,c(1:108)]
test <- tric_translatome_dt[adjusted_start >= -50 & adjusted_start <= 50]
test2 <- tric_translatome_dt[adjusted_random >= -50 & adjusted_random <= 50]
setkeyv(test, c("name"))
setkeyv(test2, c("name"))
setkeyv(tric_translatome_dt, c("name"))
test[, WT_adjusted_rpc := mean(WT_sum), by = name]
test[, EV_adjusted_rpc := mean(EV_sum), by = name]
test[, B_WT_adjusted_rpc := mean(B_WT_sum), by = name]
test[, G_WT_adjusted_rpc := mean(G_WT_sum), by = name]
test[, C_WT_adjusted_rpc := mean(C_WT_sum), by = name]
test[, L_WT_adjusted_rpc := mean(L_WT_sum), by = name]
test2[, WT_adjusted_rpc_random := mean(WT_sum), by = name]
i <- cbind(match(tric_translatome_dt$name, test$name))
tric_translatome_dt <- cbind(tric_translatome_dt, WT_adjusted_rpc = test[i]$WT_adjusted_rpc)
tric_translatome_dt <- cbind(tric_translatome_dt, EV_adjusted_rpc = test[i]$EV_adjusted_rpc)
tric_translatome_dt <- cbind(tric_translatome_dt, B_WT_adjusted_rpc = test[i]$B_WT_adjusted_rpc)
tric_translatome_dt <- cbind(tric_translatome_dt, G_WT_adjusted_rpc = test[i]$G_WT_adjusted_rpc)
tric_translatome_dt <- cbind(tric_translatome_dt, C_WT_adjusted_rpc = test[i]$C_WT_adjusted_rpc)
tric_translatome_dt <- cbind(tric_translatome_dt, L_WT_adjusted_rpc = test[i]$L_WT_adjusted_rpc)
i <- cbind(match(tric_translatome_dt$name, test2$name))
tric_translatome_dt <- cbind(tric_translatome_dt, WT_adjusted_rpc_random = test2[i]$WT_adjusted_rpc_random)
tric_translatome_dt[, WT_adjusted_norm := WT_sum / WT_adjusted_rpc]
tric_translatome_dt[, EV_adjusted_norm := EV_sum / EV_adjusted_rpc]
tric_translatome_dt[, B_WT_adjusted_norm := B_WT_sum / B_WT_adjusted_rpc]
tric_translatome_dt[, G_WT_adjusted_norm := G_WT_sum / G_WT_adjusted_rpc]
tric_translatome_dt[, C_WT_adjusted_norm := C_WT_sum / C_WT_adjusted_rpc]
tric_translatome_dt[, L_WT_adjusted_norm := L_WT_sum / L_WT_adjusted_rpc]
tric_translatome_dt[, WT_adjusted_norm_random := WT_sum / WT_adjusted_rpc_random]

i <- cbind(match(tric_translatome_dt$position1, tric_substrates5_dt$position1))
tric_translatome_dt <- cbind(tric_translatome_dt, tric_odds1 = tric_substrates5_dt[i]$tric_odds)


tric_translatome_dt_max <- translatome_dt[tric_substrates5]

temp <- tric_peaks5_max
random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(30:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

tric_translatome_dt_max <- tric_translatome_dt_max[temp, allow.cartesian = TRUE]
tric_translatome_dt_max[, name := as.character(base::paste(orf, peak, sep = "_"))]
tric_translatome_dt_max[, adjusted := position - peak]
tric_translatome_dt_max[, adjusted_start := position - peak_start]
tric_translatome_dt_max[, adjusted_random := position - random]
tric_translatome_dt_max <- tric_translatome_dt_max[,c(1:108)]
test <- tric_translatome_dt_max[adjusted_start >= -50 & adjusted_start <= 50]
test2 <- tric_translatome_dt_max[adjusted_random >= -50 & adjusted_random <= 50]
setkeyv(test, c("name"))
setkeyv(tric_translatome_dt_max, c("name"))
setkeyv(test2, c("name"))
test[, WT_adjusted_rpc := mean(WT_sum), by = name]
test[, EV_adjusted_rpc := mean(EV_sum), by = name]
test[, B_WT_adjusted_rpc := mean(B_WT_sum), by = name]
test[, G_WT_adjusted_rpc := mean(G_WT_sum), by = name]
test[, C_WT_adjusted_rpc := mean(C_WT_sum), by = name]
test[, L_WT_adjusted_rpc := mean(L_WT_sum), by = name]
test2[, WT_adjusted_rpc_random := mean(WT_sum), by = name]
i <- cbind(match(tric_translatome_dt_max$name, test$name))
tric_translatome_dt_max <- cbind(tric_translatome_dt_max, WT_adjusted_rpc = test[i]$WT_adjusted_rpc)
tric_translatome_dt_max <- cbind(tric_translatome_dt_max, EV_adjusted_rpc = test[i]$EV_adjusted_rpc)
tric_translatome_dt_max <- cbind(tric_translatome_dt_max, B_WT_adjusted_rpc = test[i]$B_WT_adjusted_rpc)
tric_translatome_dt_max <- cbind(tric_translatome_dt_max, G_WT_adjusted_rpc = test[i]$G_WT_adjusted_rpc)
tric_translatome_dt_max <- cbind(tric_translatome_dt_max, C_WT_adjusted_rpc = test[i]$C_WT_adjusted_rpc)
tric_translatome_dt_max <- cbind(tric_translatome_dt_max, L_WT_adjusted_rpc = test[i]$L_WT_adjusted_rpc)
i <- cbind(match(tric_translatome_dt_max$name, test2$name))
tric_translatome_dt_max <- cbind(tric_translatome_dt_max, WT_adjusted_rpc_random = test2[i]$WT_adjusted_rpc_random)
tric_translatome_dt_max[, WT_adjusted_norm := WT_sum / WT_adjusted_rpc]
tric_translatome_dt_max[, EV_adjusted_norm := EV_sum / EV_adjusted_rpc]
tric_translatome_dt_max[, B_WT_adjusted_norm := B_WT_sum / B_WT_adjusted_rpc]
tric_translatome_dt_max[, G_WT_adjusted_norm := G_WT_sum / G_WT_adjusted_rpc]
tric_translatome_dt_max[, C_WT_adjusted_norm := C_WT_sum / C_WT_adjusted_rpc]
tric_translatome_dt_max[, L_WT_adjusted_norm := L_WT_sum / L_WT_adjusted_rpc]
tric_translatome_dt_max[, WT_adjusted_norm_random := WT_sum / WT_adjusted_rpc_random]

i <- cbind(match(tric_translatome_dt_max$position1, tric_substrates5_dt$position1))
tric_translatome_dt_max <- cbind(tric_translatome_dt_max, tric_odds1 = tric_substrates5_dt[i]$tric_odds)


# SSB
ssb_translatome_dt <- translatome_dt[ssb_substrates5]

temp <- ssb_peaks5
random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(30:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

ssb_translatome_dt <- ssb_translatome_dt[temp, allow.cartesian = TRUE]
ssb_translatome_dt[, name := as.character(base::paste(orf, peak, sep = "_"))]
ssb_translatome_dt[, adjusted := position - peak]
ssb_translatome_dt[, adjusted_start := position - peak_start]
ssb_translatome_dt[, adjusted_random := position - random]
ssb_translatome_dt <- ssb_translatome_dt[, c(1:107,120)]
test <- ssb_translatome_dt[adjusted_start >= -50 & adjusted_start <= 50]
test2 <- ssb_translatome_dt[adjusted_random >= -50 & adjusted_random <= 50]
setkeyv(test, c("name"))
setkeyv(test2, c("name"))
setkeyv(ssb_translatome_dt, c("name"))
test[, WT_adjusted_rpc := mean(WT_sum), by = name]
test[, EV_adjusted_rpc := mean(EV_sum), by = name]
test[, B_WT_adjusted_rpc := mean(B_WT_sum), by = name]
test[, G_WT_adjusted_rpc := mean(G_WT_sum), by = name]
test[, C_WT_adjusted_rpc := mean(C_WT_sum), by = name]
test[, L_WT_adjusted_rpc := mean(L_WT_sum), by = name]
test2[, WT_adjusted_rpc_random := mean(WT_sum), by = name]
i <- cbind(match(ssb_translatome_dt$name, test$name))
ssb_translatome_dt <- cbind(ssb_translatome_dt, WT_adjusted_rpc = test[i]$WT_adjusted_rpc)
ssb_translatome_dt <- cbind(ssb_translatome_dt, EV_adjusted_rpc = test[i]$EV_adjusted_rpc)
ssb_translatome_dt <- cbind(ssb_translatome_dt, B_WT_adjusted_rpc = test[i]$B_WT_adjusted_rpc)
ssb_translatome_dt <- cbind(ssb_translatome_dt, G_WT_adjusted_rpc = test[i]$G_WT_adjusted_rpc)
ssb_translatome_dt <- cbind(ssb_translatome_dt, C_WT_adjusted_rpc = test[i]$C_WT_adjusted_rpc)
ssb_translatome_dt <- cbind(ssb_translatome_dt, L_WT_adjusted_rpc = test[i]$L_WT_adjusted_rpc)
i <- cbind(match(ssb_translatome_dt$name, test2$name))
ssb_translatome_dt <- cbind(ssb_translatome_dt, WT_adjusted_rpc_random = test2[i]$WT_adjusted_rpc_random)
ssb_translatome_dt[, WT_adjusted_norm := WT_sum / WT_adjusted_rpc]
ssb_translatome_dt[, EV_adjusted_norm := EV_sum / EV_adjusted_rpc]
ssb_translatome_dt[, B_WT_adjusted_norm := B_WT_sum / B_WT_adjusted_rpc]
ssb_translatome_dt[, G_WT_adjusted_norm := G_WT_sum / G_WT_adjusted_rpc]
ssb_translatome_dt[, C_WT_adjusted_norm := C_WT_sum / C_WT_adjusted_rpc]
ssb_translatome_dt[, L_WT_adjusted_norm := L_WT_sum / L_WT_adjusted_rpc]
ssb_translatome_dt[, WT_adjusted_norm_random := WT_sum / WT_adjusted_rpc_random]

i <- cbind(match(ssb_translatome_dt$position1, ssb_substrates5_dt$position1))
ssb_translatome_dt <- cbind(ssb_translatome_dt, ssb_odds1 = ssb_substrates5_dt[i]$ssb_odds)


ssb_translatome_dt_max <- translatome_dt[ssb_substrates5]

temp <- ssb_peaks5_max
random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(30:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

ssb_translatome_dt_max <- ssb_translatome_dt_max[temp, allow.cartesian = TRUE]
ssb_translatome_dt_max[, name := as.character(base::paste(orf, peak, sep = "_"))]
ssb_translatome_dt_max[, adjusted := position - peak]
ssb_translatome_dt_max[, adjusted_start := position - peak_start]
ssb_translatome_dt_max[, adjusted_random := position - random]
ssb_translatome_dt_max <- ssb_translatome_dt_max[, c(1:108)]
test <- ssb_translatome_dt_max[adjusted_start >= -50 & adjusted_start <= 50]
test2 <- ssb_translatome_dt_max[adjusted_random >= -50 & adjusted_random <= 50]
setkeyv(test, c("name"))
setkeyv(test2, c("name"))
setkeyv(ssb_translatome_dt_max, c("name"))
test[, WT_adjusted_rpc := mean(WT_sum), by = name]
test[, EV_adjusted_rpc := mean(EV_sum), by = name]
test[, B_WT_adjusted_rpc := mean(B_WT_sum), by = name]
test[, G_WT_adjusted_rpc := mean(G_WT_sum), by = name]
test[, C_WT_adjusted_rpc := mean(C_WT_sum), by = name]
test[, L_WT_adjusted_rpc := mean(L_WT_sum), by = name]
test2[, WT_adjusted_rpc_random := mean(WT_sum), by = name]
i <- cbind(match(ssb_translatome_dt_max$name, test$name))
ssb_translatome_dt_max <- cbind(ssb_translatome_dt_max, WT_adjusted_rpc = test[i]$WT_adjusted_rpc)
ssb_translatome_dt_max <- cbind(ssb_translatome_dt_max, EV_adjusted_rpc = test[i]$EV_adjusted_rpc)
ssb_translatome_dt_max <- cbind(ssb_translatome_dt_max, B_WT_adjusted_rpc = test[i]$B_WT_adjusted_rpc)
ssb_translatome_dt_max <- cbind(ssb_translatome_dt_max, G_WT_adjusted_rpc = test[i]$G_WT_adjusted_rpc)
ssb_translatome_dt_max <- cbind(ssb_translatome_dt_max, C_WT_adjusted_rpc = test[i]$C_WT_adjusted_rpc)
ssb_translatome_dt_max <- cbind(ssb_translatome_dt_max, L_WT_adjusted_rpc = test[i]$L_WT_adjusted_rpc)
i <- cbind(match(ssb_translatome_dt_max$name, test2$name))
ssb_translatome_dt_max <- cbind(ssb_translatome_dt_max, WT_adjusted_rpc_random = test2[i]$WT_adjusted_rpc_random)
ssb_translatome_dt_max[, WT_adjusted_norm := WT_sum / WT_adjusted_rpc]
ssb_translatome_dt_max[, EV_adjusted_norm := EV_sum / EV_adjusted_rpc]
ssb_translatome_dt_max[, B_WT_adjusted_norm := B_WT_sum / B_WT_adjusted_rpc]
ssb_translatome_dt_max[, G_WT_adjusted_norm := G_WT_sum / G_WT_adjusted_rpc]
ssb_translatome_dt_max[, C_WT_adjusted_norm := C_WT_sum / C_WT_adjusted_rpc]
ssb_translatome_dt_max[, L_WT_adjusted_norm := L_WT_sum / L_WT_adjusted_rpc]
ssb_translatome_dt_max[, WT_adjusted_norm_random := WT_sum / WT_adjusted_rpc_random]

i <- cbind(match(ssb_translatome_dt_max$position1, ssb_substrates5_dt$position1))
ssb_translatome_dt_max <- cbind(ssb_translatome_dt_max, ssb_odds1 = ssb_substrates5_dt[i]$ssb_odds)


# Bukau
bukau_translatome_dt <- translatome_dt[bukau_substrates5]

temp <- bukau_peaks5
random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(30:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

bukau_translatome_dt <- bukau_translatome_dt[temp, allow.cartesian = TRUE]
bukau_translatome_dt[, name := as.character(base::paste(orf, peak, sep = "_"))]
bukau_translatome_dt[, adjusted := position - peak]
bukau_translatome_dt[, adjusted_start := position - peak_start]
bukau_translatome_dt[, adjusted_random := position - random]
bukau_translatome_dt <- bukau_translatome_dt[, c(1:107,120)]
test <- bukau_translatome_dt[adjusted_start >= -50 & adjusted_start <= 50]
test2 <- bukau_translatome_dt[adjusted_random >= -50 & adjusted_random <= 50]
setkeyv(test, c("name"))
setkeyv(test2, c("name"))
setkeyv(bukau_translatome_dt, c("name"))
test[, WT_adjusted_rpc := mean(WT_sum), by = name]
test[, EV_adjusted_rpc := mean(EV_sum), by = name]
test[, B_WT_adjusted_rpc := mean(B_WT_sum), by = name]
test[, G_WT_adjusted_rpc := mean(G_WT_sum), by = name]
test[, C_WT_adjusted_rpc := mean(C_WT_sum), by = name]
test[, L_WT_adjusted_rpc := mean(L_WT_sum), by = name]
test2[, WT_adjusted_rpc_random := mean(WT_sum), by = name]
i <- cbind(match(bukau_translatome_dt$name, test$name))
bukau_translatome_dt <- cbind(bukau_translatome_dt, WT_adjusted_rpc = test[i]$WT_adjusted_rpc)
bukau_translatome_dt <- cbind(bukau_translatome_dt, EV_adjusted_rpc = test[i]$EV_adjusted_rpc)
bukau_translatome_dt <- cbind(bukau_translatome_dt, B_WT_adjusted_rpc = test[i]$B_WT_adjusted_rpc)
bukau_translatome_dt <- cbind(bukau_translatome_dt, G_WT_adjusted_rpc = test[i]$G_WT_adjusted_rpc)
bukau_translatome_dt <- cbind(bukau_translatome_dt, C_WT_adjusted_rpc = test[i]$C_WT_adjusted_rpc)
bukau_translatome_dt <- cbind(bukau_translatome_dt, L_WT_adjusted_rpc = test[i]$L_WT_adjusted_rpc)
i <- cbind(match(bukau_translatome_dt$name, test2$name))
bukau_translatome_dt <- cbind(bukau_translatome_dt, WT_adjusted_rpc_random = test2[i]$WT_adjusted_rpc_random)
bukau_translatome_dt[, WT_adjusted_norm := WT_sum / WT_adjusted_rpc]
bukau_translatome_dt[, EV_adjusted_norm := EV_sum / EV_adjusted_rpc]
bukau_translatome_dt[, B_WT_adjusted_norm := B_WT_sum / B_WT_adjusted_rpc]
bukau_translatome_dt[, G_WT_adjusted_norm := G_WT_sum / G_WT_adjusted_rpc]
bukau_translatome_dt[, C_WT_adjusted_norm := C_WT_sum / C_WT_adjusted_rpc]
bukau_translatome_dt[, L_WT_adjusted_norm := L_WT_sum / L_WT_adjusted_rpc]
bukau_translatome_dt[, WT_adjusted_norm_random := WT_sum / WT_adjusted_rpc_random]

i <- cbind(match(bukau_translatome_dt$position1, bukau_substrates5_dt$position1))
bukau_translatome_dt <- cbind(bukau_translatome_dt, ssb_odds1 = bukau_substrates5_dt[i]$ssb_odds)


bukau_translatome_dt_max <- translatome_dt[bukau_substrates5]

temp <- bukau_peaks5_max
random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(30:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

bukau_translatome_dt_max <- bukau_translatome_dt_max[temp, allow.cartesian = TRUE]
bukau_translatome_dt_max[, name := as.character(base::paste(orf, peak, sep = "_"))]
bukau_translatome_dt_max[, adjusted := position - peak]
bukau_translatome_dt_max[, adjusted_start := position - peak_start]
bukau_translatome_dt_max[, adjusted_random := position - random]
bukau_translatome_dt_max <- bukau_translatome_dt_max[, c(1:108)]
test <- bukau_translatome_dt_max[adjusted_start >= -50 & adjusted_start <= 50]
test2 <- bukau_translatome_dt_max[adjusted_random >= -50 & adjusted_random <= 50]
setkeyv(test, c("name"))
setkeyv(test2, c("name"))
setkeyv(bukau_translatome_dt_max, c("name"))
test[, WT_adjusted_rpc := mean(WT_sum), by = name]
test[, EV_adjusted_rpc := mean(EV_sum), by = name]
test[, B_WT_adjusted_rpc := mean(B_WT_sum), by = name]
test[, G_WT_adjusted_rpc := mean(G_WT_sum), by = name]
test[, C_WT_adjusted_rpc := mean(C_WT_sum), by = name]
test[, L_WT_adjusted_rpc := mean(L_WT_sum), by = name]
test2[, WT_adjusted_rpc_random := mean(WT_sum), by = name]
i <- cbind(match(bukau_translatome_dt_max$name, test$name))
bukau_translatome_dt_max <- cbind(bukau_translatome_dt_max, WT_adjusted_rpc = test[i]$WT_adjusted_rpc)
bukau_translatome_dt_max <- cbind(bukau_translatome_dt_max, EV_adjusted_rpc = test[i]$EV_adjusted_rpc)
bukau_translatome_dt_max <- cbind(bukau_translatome_dt_max, B_WT_adjusted_rpc = test[i]$B_WT_adjusted_rpc)
bukau_translatome_dt_max <- cbind(bukau_translatome_dt_max, G_WT_adjusted_rpc = test[i]$G_WT_adjusted_rpc)
bukau_translatome_dt_max <- cbind(bukau_translatome_dt_max, C_WT_adjusted_rpc = test[i]$C_WT_adjusted_rpc)
bukau_translatome_dt_max <- cbind(bukau_translatome_dt_max, L_WT_adjusted_rpc = test[i]$L_WT_adjusted_rpc)
i <- cbind(match(bukau_translatome_dt_max$name, test2$name))
bukau_translatome_dt_max <- cbind(bukau_translatome_dt_max, WT_adjusted_rpc_random = test2[i]$WT_adjusted_rpc_random)
bukau_translatome_dt_max[, WT_adjusted_norm := WT_sum / WT_adjusted_rpc]
bukau_translatome_dt_max[, EV_adjusted_norm := EV_sum / EV_adjusted_rpc]
bukau_translatome_dt_max[, B_WT_adjusted_norm := B_WT_sum / B_WT_adjusted_rpc]
bukau_translatome_dt_max[, G_WT_adjusted_norm := G_WT_sum / G_WT_adjusted_rpc]
bukau_translatome_dt_max[, C_WT_adjusted_norm := C_WT_sum / C_WT_adjusted_rpc]
bukau_translatome_dt_max[, L_WT_adjusted_norm := L_WT_sum / L_WT_adjusted_rpc]
bukau_translatome_dt_max[, WT_adjusted_norm_random := WT_sum / WT_adjusted_rpc_random]

i <- cbind(match(bukau_translatome_dt_max$position1, bukau_substrates5_dt$position1))
bukau_translatome_dt_max <- cbind(bukau_translatome_dt_max, ssb_odds1 = bukau_substrates5_dt[i]$ssb_odds)


# Doring
doring_translatome_dt <- translatome_dt[doring_substrates]
setkeyv(doring_translatome_dt, c("orf"))
doring_translatome_dt <- doring_translatome_dt[orf != "YLR262C" & orf != "YNL139C" & orf != "YAR075W"]
doring_translatome_dt1 <- doring_translatome_dt[position > 40 & stopdist < -20]
doring_translatome_dt_B1 <- doring_translatome_dt1[, .SD[all(B_WT_sum > 0)], by = orf]
orfs2 <- doring_translatome_dt_B1[, .SD[which.min(position)], by = orf]
doring_translatome_dt_B1 <- doring_translatome_dt[doring_translatome_dt$orf %in% orfs2$orf]

temp <- bukau_dataset_dt[bukau_dataset_dt$orf %in% orfs2$orf]

random <- NULL
for (i in 1:nrow(temp)) {
  length1 <- temp[i]$length
  random1 <- sample(30:length1, size = 1, replace = T)
  random <- c(random, random1)
}
temp$random <- random

doring_translatome_dt_B1 <- doring_translatome_dt_B1[temp, allow.cartesian = TRUE]
doring_translatome_dt_B1[, name := as.character(base::paste(orf, peak, sep = "_"))]
doring_translatome_dt_B1[, adjusted := position - peak]
doring_translatome_dt_B1[, adjusted_start := position - peak_startaa]
doring_translatome_dt_B1[, adjusted_random := position - random]
doring_translatome_dt_B1 <- doring_translatome_dt_B1[,c(1:97)]
test <- doring_translatome_dt_B1[adjusted_start >= -50 & adjusted_start <= 50]
test <- doring_translatome_dt_B1[position >= 30 & stopdist <= -30]
test2 <- doring_translatome_dt_B1[adjusted_random >= -50 & adjusted_random <= 50]
test2 <- doring_translatome_dt_B1[position >= 30 & stopdist <= -30]
setkeyv(test, c("name"))
setkeyv(test2, c("name"))
setkeyv(doring_translatome_dt_B1, c("name"))
test[, WT_adjusted_rpc := mean(WT_sum), by = name]
test[, EV_adjusted_rpc := mean(EV_sum), by = name]
test[, B_WT_adjusted_rpc := mean(B_WT_sum), by = name]
test[, G_WT_adjusted_rpc := mean(G_WT_sum), by = name]
test[, C_WT_adjusted_rpc := mean(C_WT_sum), by = name]
test[, L_WT_adjusted_rpc := mean(L_WT_sum), by = name]
test2[, WT_adjusted_rpc_random := mean(WT_sum), by = name]
test2[, EV_adjusted_rpc_random := mean(EV_sum), by = name]
test2[, B_WT_adjusted_rpc_random := mean(B_WT_sum), by = name]
test2[, G_WT_adjusted_rpc_random := mean(G_WT_sum), by = name]
test2[, C_WT_adjusted_rpc_random := mean(C_WT_sum), by = name]
test2[, L_WT_adjusted_rpc_random := mean(L_WT_sum), by = name]
i <- cbind(match(doring_translatome_dt_B1$name, test$name))
doring_translatome_dt_B1 <- cbind(doring_translatome_dt_B1, WT_adjusted_rpc = test[i]$WT_adjusted_rpc)
doring_translatome_dt_B1 <- cbind(doring_translatome_dt_B1, EV_adjusted_rpc = test[i]$EV_adjusted_rpc)
doring_translatome_dt_B1 <- cbind(doring_translatome_dt_B1, B_WT_adjusted_rpc = test[i]$B_WT_adjusted_rpc)
doring_translatome_dt_B1 <- cbind(doring_translatome_dt_B1, G_WT_adjusted_rpc = test[i]$G_WT_adjusted_rpc)
doring_translatome_dt_B1 <- cbind(doring_translatome_dt_B1, C_WT_adjusted_rpc = test[i]$C_WT_adjusted_rpc)
doring_translatome_dt_B1 <- cbind(doring_translatome_dt_B1, L_WT_adjusted_rpc = test[i]$L_WT_adjusted_rpc)
i <- cbind(match(doring_translatome_dt_B1$name, test2$name))
doring_translatome_dt_B1 <- cbind(doring_translatome_dt_B1, WT_adjusted_rpc_random = test2[i]$WT_adjusted_rpc_random)
doring_translatome_dt_B1 <- cbind(doring_translatome_dt_B1, L_WT_adjusted_rpc_random = test2[i]$L_WT_adjusted_rpc_random)
doring_translatome_dt_B1 <- cbind(doring_translatome_dt_B1, C_WT_adjusted_rpc_random = test2[i]$C_WT_adjusted_rpc_random)
doring_translatome_dt_B1 <- cbind(doring_translatome_dt_B1, G_WT_adjusted_rpc_random = test2[i]$G_WT_adjusted_rpc_random)
doring_translatome_dt_B1 <- cbind(doring_translatome_dt_B1, B_WT_adjusted_rpc_random = test2[i]$B_WT_adjusted_rpc_random)
doring_translatome_dt_B1[, WT_adjusted_norm := WT_sum / WT_adjusted_rpc]
doring_translatome_dt_B1[, EV_adjusted_norm := EV_sum / EV_adjusted_rpc]
doring_translatome_dt_B1[, B_WT_adjusted_norm := B_WT_sum / B_WT_adjusted_rpc]
doring_translatome_dt_B1[, G_WT_adjusted_norm := G_WT_sum / G_WT_adjusted_rpc]
doring_translatome_dt_B1[, C_WT_adjusted_norm := C_WT_sum / C_WT_adjusted_rpc]
doring_translatome_dt_B1[, L_WT_adjusted_norm := L_WT_sum / L_WT_adjusted_rpc]
doring_translatome_dt_B1[, WT_adjusted_norm_random := WT_sum / WT_adjusted_rpc_random]
doring_translatome_dt_B1[, B_WT_adjusted_norm_random := B_WT_sum / B_WT_adjusted_rpc_random]
doring_translatome_dt_B1[, L_WT_adjusted_norm_random := L_WT_sum / L_WT_adjusted_rpc_random]
doring_translatome_dt_B1[, C_WT_adjusted_norm_random := C_WT_sum / C_WT_adjusted_rpc_random]
doring_translatome_dt_B1[, G_WT_adjusted_norm_random := G_WT_sum / G_WT_adjusted_rpc_random]

ggplot(data = doring_translatome_dt_B1[B_WT_adjusted_rpc >= 0.5 & peak > 75, .(adjusted_start, occupancy = movingAverage(B_WT_adjusted_norm, n=5, center=T))]) +
  geom_vline(xintercept = 0, color = "gray70", linetype = 'dashed', size = 1.25) +
  stat_summary(data = doring_translatome_dt_B1[B_WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(B_WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = 'gray50') +
  stat_summary(data = doring_translatome_dt_B1[B_WT_adjusted_rpc_random >= 0.5 & random > 75, .(adjusted_random, occupancy2 = movingAverage(B_WT_adjusted_norm_random, n=5, center=T))],
               aes(adjusted_random, occupancy2), fun.y = "mean", geom = "line", size = 1.25, color = 'gray50') + 
  stat_summary(aes(adjusted_start, occupancy), fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.3,
               fun.args=list(conf.int=0.5), fill = '#6A3D9A') +
  stat_summary(aes(adjusted_start, occupancy), fun.y = "mean", geom = "line", size = 1.25, color = '#6A3D9A') +
  xlim(-100, 200)
coord_cartesian(ylim = c(0.8,1.4))


### SRP
srp_peaks <- read.csv("/Users/KevinStein/Desktop/SRP_pronouncedBinding.csv", header = TRUE, stringsAsFactors = TRUE)

srp_peaks <- as.data.table(srp_peaks)
i <- cbind(match(srp_peaks$orf, translatome_dt$orf))
srp_peaks <- cbind(srp_peaks, length = translatome_dt[i]$length)

random <- NULL
for (i in 1:nrow(srp_peaks)) {
  length1 <- temp[i]$length
  random1 <- sample(30:length1, size = 1, replace = T)
  random <- c(random, random1)
}
srp_peaks$random <- random

srp_translatome_dt <- translatome_dt[srp_peaks, allow.cartesian = TRUE]
srp_translatome_dt[, name := as.character(base::paste(orf, binding.position, sep = "_"))]
srp_translatome_dt[, adjusted_start := position - binding.position]
srp_translatome_dt[, adjusted_random := position - random]
srp_translatome_dt <- srp_translatome_dt[, c(1:95)]
test <- srp_translatome_dt[adjusted_start >= -40 & adjusted_start <= -20]
test2 <- srp_translatome_dt[adjusted_random >= -40 & adjusted_random <= -20]
test <- srp_translatome_dt[position >= 30 & stopdist <= -30]
test2 <- srp_translatome_dt[position >= 30 & stopdist <= -30]
setkeyv(test, c("name"))
setkeyv(test2, c("name"))
setkeyv(srp_translatome_dt, c("name"))
test[, WT_adjusted_rpc := mean(WT_sum), by = name]
test[, EV_adjusted_rpc := mean(EV_sum), by = name]
test[, B_WT_adjusted_rpc := mean(B_WT_sum), by = name]
test[, G_WT_adjusted_rpc := mean(G_WT_sum), by = name]
test[, C_WT_adjusted_rpc := mean(C_WT_sum), by = name]
test[, L_WT_adjusted_rpc := mean(L_WT_sum), by = name]
test2[, WT_adjusted_rpc_random := mean(WT_sum), by = name]
i <- cbind(match(srp_translatome_dt$name, test$name))
srp_translatome_dt <- cbind(srp_translatome_dt, WT_adjusted_rpc = test[i]$WT_adjusted_rpc)
srp_translatome_dt <- cbind(srp_translatome_dt, EV_adjusted_rpc = test[i]$EV_adjusted_rpc)
srp_translatome_dt <- cbind(srp_translatome_dt, B_WT_adjusted_rpc = test[i]$B_WT_adjusted_rpc)
srp_translatome_dt <- cbind(srp_translatome_dt, G_WT_adjusted_rpc = test[i]$G_WT_adjusted_rpc)
srp_translatome_dt <- cbind(srp_translatome_dt, C_WT_adjusted_rpc = test[i]$C_WT_adjusted_rpc)
srp_translatome_dt <- cbind(srp_translatome_dt, L_WT_adjusted_rpc = test[i]$L_WT_adjusted_rpc)
i <- cbind(match(srp_translatome_dt$name, test2$name))
srp_translatome_dt <- cbind(srp_translatome_dt, WT_adjusted_rpc_random = test2[i]$WT_adjusted_rpc_random)
srp_translatome_dt[, WT_adjusted_norm := WT_sum / WT_adjusted_rpc]
srp_translatome_dt[, EV_adjusted_norm := EV_sum / EV_adjusted_rpc]
srp_translatome_dt[, B_WT_adjusted_norm := B_WT_sum / B_WT_adjusted_rpc]
srp_translatome_dt[, G_WT_adjusted_norm := G_WT_sum / G_WT_adjusted_rpc]
srp_translatome_dt[, C_WT_adjusted_norm := C_WT_sum / C_WT_adjusted_rpc]
srp_translatome_dt[, L_WT_adjusted_norm := L_WT_sum / L_WT_adjusted_rpc]
srp_translatome_dt[, WT_adjusted_norm_random := WT_sum / WT_adjusted_rpc_random]

