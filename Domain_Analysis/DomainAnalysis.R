### Read lists from SGD of domains in putative TRiC substrates. 
# Domains from SGD in curated data
domains <- read.delim("/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/Yeast/YeastProteomeCategorization/domains_SGDcurated.txt", header = T, stringsAsFactors = T)
domains <- as.data.table(domains)
domains[, DomainLength := End - Start + 1]

# Domains from SGD using column selection
# domains1 <- read.delim("/Users/KevinStein/Desktop/domains_GeneListDownload.tsv", header=F)
# colnames(domains1) <- c("SystematicGene", "name", "Start", "End", "Method", "Domain","Description")
# domains1 <- as.data.table(domains1)
# domains1[, DomainLength := End - Start + 1]

# Subset only Gene3D and Superfamily methods and add description from Interpro
setnames(domains, "SystematicGene","orf")
domains_subset <- domains[Method == "Gene3D" | Method == "SUPERFAMILY"]
setkeyv(domains_subset, c("orf", "Start"))
domains_cath <- domains_subset[Method == "Gene3D"]
domains_scop <- domains_subset[Method == "SUPERFAMILY"]
domains_cathgenes <- as.data.frame(unique(domains_cath[,1]))
domains_scopgenes <- as.data.frame(unique(domains_scop[,1]))


### Add column of downstream domain start position
ds <- NULL
for (g in domains_cathgenes$orf) {
  cath <- domains_cath[g]$Start
  for (i in 1:length(cath)) {
    ds1 <- ifelse(!is.na(cath[i+1]), cath[i+1], NA)
    ds <- c(ds, ds1)
  }
}
domains_cath$ds_domain <- ds

ds <- NULL
for (g in domains_scopgenes$orf) {
  scop <- domains_scop[g]$Start
  for (i in 1:length(scop)) {
    ds1 <- ifelse(!is.na(scop[i+1]), scop[i+1], NA)
    ds <- c(ds, ds1)
  }
}
domains_scop$ds_domain <- ds


### Add Gene3D hierarchical description of domain
domains_cath[, Gene3D_domain := substring(Domain,7,200)]
Gene3D_topo <- NULL
Gene3D_arch <- NULL
for (i in 1:nrow(domains_cath)) {
  str <- unlist(strsplit(domains_cath[i]$Gene3D_domain, "[.]"))
  new_str1 <- base::paste(str[1], str[2], str[3], sep = ".")
  Gene3D_topo <- c(Gene3D_topo, new_str1)
  new_str2 <- base::paste(str[1], str[2], sep = ".")
  Gene3D_arch <- c(Gene3D_arch, new_str2)
}
domains_cath$Gene3D_topo <- Gene3D_topo
domains_cath$Gene3D_arch <- Gene3D_arch


### Idenitfy domains enriched in TRiC and SSB substrates
# Do this for CATH and SCOP
domains_cath[, domain_alias := paste0(orf, "_", Domain, "_", Start, "_", End)]
domains_scop[, domain_alias := paste0(orf, "_", Domain, "_", Start, "_", End)]

setkeyv(tric_peaks5, c("orf"))
domains_cath_tric3 <- domains_cath[tric_peaks3]
domains_cath_tric3 <- domains_cath_tric3[!is.na(Start)]
domains_cath_tric5 <- domains_cath[tric_peaks5]
domains_cath_tric5 <- domains_cath_tric5[!is.na(Start)]
domains_cath_tric3_max <- domains_cath[tric_peaks3_max]
domains_cath_tric3_max <- domains_cath_tric3_max[!is.na(Start)]
domains_cath_tric5_max <- domains_cath[tric_peaks5_max]
domains_cath_tric5_max <- domains_cath_tric5_max[!is.na(Start)]
domains_scop_tric3 <- domains_scop[tric_peaks3]
domains_scop_tric3 <- domains_scop_tric3[!is.na(Start)]
domains_scop_tric5 <- domains_scop[tric_peaks5]
domains_scop_tric5 <- domains_scop_tric5[!is.na(Start)]
domains_scop_tric3_max <- domains_scop[tric_peaks3_max]
domains_scop_tric3_max <- domains_scop_tric3_max[!is.na(Start)]
domains_scop_tric5_max <- domains_scop[tric_peaks5_max]
domains_scop_tric5_max <- domains_scop_tric5_max[!is.na(Start)]

setkeyv(ssb_peaks5, c("orf"))
domains_cath_ssb3 <- domains_cath[ssb_peaks3]
domains_cath_ssb3 <- domains_cath_ssb3[!is.na(Start)]
domains_cath_ssb5 <- domains_cath[ssb_peaks5]
domains_cath_ssb5 <- domains_cath_ssb5[!is.na(Start)]
domains_cath_ssb3_max <- domains_cath[ssb_peaks3_max]
domains_cath_ssb3_max <- domains_cath_ssb3_max[!is.na(Start)]
domains_cath_ssb5_max <- domains_cath[ssb_peaks5_max]
domains_cath_ssb5_max <- domains_cath_ssb5_max[!is.na(Start)]
domains_scop_ssb3 <- domains_scop[ssb_peaks3]
domains_scop_ssb3 <- domains_scop_ssb3[!is.na(Start)]
domains_scop_ssb5 <- domains_scop[ssb_peaks5]
domains_scop_ssb5 <- domains_scop_ssb5[!is.na(Start)]
domains_scop_ssb3_max <- domains_scop[ssb_peaks3_max]
domains_scop_ssb3_max <- domains_scop_ssb3_max[!is.na(Start)]
domains_scop_ssb5_max <- domains_scop[ssb_peaks5_max]
domains_scop_ssb5_max <- domains_scop_ssb5_max[!is.na(Start)]
setkeyv(bukau_peaks5, c("orf"))
domains_cath_bukau3 <- domains_cath[bukau_peaks3]
domains_cath_bukau3 <- domains_cath_bukau3[!is.na(Start)]
domains_cath_bukau5 <- domains_cath[bukau_peaks5]
domains_cath_bukau5 <- domains_cath_bukau5[!is.na(Start)]
domains_cath_bukau3_max <- domains_cath[bukau_peaks3_max]
domains_cath_bukau3_max <- domains_cath_bukau3_max[!is.na(Start)]
domains_cath_bukau5_max <- domains_cath[bukau_peaks5_max]
domains_cath_bukau5_max <- domains_cath_bukau5_max[!is.na(Start)]
domains_scop_bukau5 <- domains_scop[bukau_peaks5]
domains_scop_bukau5 <- domains_scop_bukau5[!is.na(Start)]


# Add domain count 
temp <- count(domains_cath_tric5_max$Gene3D_domain) # counting each orf once, then add to both max peaks and all peaks datasets
temp <- as.data.table(temp)
i <- cbind(match(domains_cath_tric5_max$Gene3D_domain, temp$x))
domains_cath_tric5_max <- cbind(domains_cath_tric5_max, DomainCount = temp[i]$freq)
i <- cbind(match(domains_cath_tric5$Gene3D_domain, temp$x))
domains_cath_tric5 <- cbind(domains_cath_tric5, DomainCount = temp[i]$freq)
temp <- count(domains_scop_tric5_max$Domain)
temp <- as.data.table(temp)
i <- cbind(match(domains_scop_tric5_max$Domain, temp$x))
domains_scop_tric5_max <- cbind(domains_scop_tric5_max, DomainCount = temp[i]$freq)
i <- cbind(match(domains_scop_tric5$Domain, temp$x))
domains_scop_tric5 <- cbind(domains_scop_tric5, DomainCount = temp[i]$freq)

temp <- count(domains_cath_ssb5_max$Gene3D_domain)
temp <- as.data.table(temp)
i <- cbind(match(domains_cath_ssb5_max$Gene3D_domain, temp$x))
domains_cath_ssb5_max <- cbind(domains_cath_ssb5_max, DomainCount = temp[i]$freq)
i <- cbind(match(domains_cath_ssb5$Gene3D_domain, temp$x))
domains_cath_ssb5 <- cbind(domains_cath_ssb5, DomainCount = temp[i]$freq)
temp <- count(domains_scop_ssb5_max$Domain)
temp <- as.data.table(temp)
i <- cbind(match(domains_scop_ssb5_max$Domain, temp$x))
domains_scop_ssb5_max <- cbind(domains_scop_ssb5_max, DomainCount = temp[i]$freq)
i <- cbind(match(domains_scop_ssb5$Domain, temp$x))
domains_scop_ssb5 <- cbind(domains_scop_ssb5, DomainCount = temp[i]$freq)
temp <- count(domains_cath_bukau5_max$Gene3D_domain)
temp <- as.data.table(temp)
i <- cbind(match(domains_cath_bukau5_max$Gene3D_domain, temp$x))
domains_cath_bukau5_max <- cbind(domains_cath_bukau5_max, DomainCount = temp[i]$freq)
i <- cbind(match(domains_cath_bukau5$Gene3D_domain, temp$x))
domains_cath_bukau5 <- cbind(domains_cath_bukau5, DomainCount = temp[i]$freq)


### Add domain frequency
freq <- NULL
y <- NULL
for (i in 1:nrow(domains_cath_tric5_max)) {
  x <- as.character(domains_cath_tric5_max[i]$Gene3D_topo)
  y1 <- length(unique(domains_cath_tric5_max[Gene3D_topo == x]$orf))
  z1 <- length(unique(domains_cath_tric5_max$orf))
  freq1 <- y1 / z1
  freq1 <- freq1*100
  freq <- c(freq, freq1)
  y <- c(y, y1)
}
domains_cath_tric5_max$freq_topo <- freq
domains_cath_tric5_max$topocount <- y

freq <- NULL
y <- NULL
for (i in 1:nrow(domains_scop_tric5_max)) {
  x <- as.character(domains_scop_tric5_max[i]$Domain)
  y1 <- length(unique(domains_scop_tric5_max[Domain == x]$orf))
  z1 <- length(unique(domains_scop_tric5_max$orf))
  freq1 <- y1 / z1
  freq1 <- freq1*100
  freq <- c(freq, freq1)
  y <- c(y, y1)
}
domains_scop_tric5_max$freq <- freq
domains_scop_tric5_max$domaincount <- y

freq <- NULL
y <- NULL
for (i in 1:nrow(domains_cath_ssb5_max)) {
  x <- as.character(domains_cath_ssb5_max[i]$Gene3D_topo)
  y1 <- length(unique(domains_cath_ssb5_max[Gene3D_topo == x]$orf))
  z1 <- length(unique(domains_cath_ssb5_max$orf))
  freq1 <- y1 / z1
  freq1 <- freq1*100
  freq <- c(freq, freq1)
  y <- c(y, y1)
}
domains_cath_ssb5_max$freq_topo <- freq
domains_cath_ssb5_max$topocount <- y

freq <- NULL
y <- NULL
for (i in 1:nrow(domains_scop_ssb5_max)) {
  x <- as.character(domains_scop_ssb5_max[i]$Domain)
  y1 <- length(unique(domains_scop_ssb5_max[Domain == x]$orf))
  z1 <- length(unique(domains_scop_ssb5_max$orf))
  freq1 <- y1 / z1
  freq1 <- freq1*100
  freq <- c(freq, freq1)
  y <- c(y, y1)
}
domains_scop_ssb5_max$freq <- freq
domains_scop_ssb5_max$domaincount <- y

freq <- NULL
y <- NULL
for (i in 1:nrow(domains_cath_bukau5_max)) {
  x <- as.character(domains_cath_bukau5_max[i]$Gene3D_topo)
  y1 <- length(unique(domains_cath_bukau5_max[Gene3D_topo == x]$orf))
  z1 <- length(unique(domains_cath_bukau5_max$orf))
  freq1 <- y1 / z1
  freq1 <- freq1*100
  freq <- c(freq, freq1)
  y <- c(y, y1)
}
domains_cath_bukau5_max$freq_topo <- freq
domains_cath_bukau5_max$topocount <- y


### Significance of particular CATH topology
length(unique(domains_cath_tric5$orf)) # 449
length(unique(domains_cath_tric5[Gene3D_topo == "2.130.10"]$orf)) # 31
length(unique(domains_cath_tric5[Gene3D_topo == "3.40.50"]$orf)) # 139
length(unique(domains_cath_tric5[Gene3D_topo == "1.25.10"]$orf)) # 14
length(unique(domains_cath_tric5[Gene3D_topo == "3.30.70"]$orf)) # 22
length(unique(domains_cath_tric5[Gene3D_topo == "3.20.20"]$orf)) # 32
length(unique(domains_cath_ssb5$orf)) # 1010
length(unique(domains_cath_ssb5[Gene3D_topo == "2.130.10"]$orf)) # 39
length(unique(domains_cath_ssb5[Gene3D_topo == "3.40.50"]$orf)) # 282
length(unique(domains_cath_ssb5[Gene3D_topo == "1.25.10"]$orf)) # 29
length(unique(domains_cath_ssb5[Gene3D_topo == "3.30.70"]$orf)) # 46
length(unique(domains_cath_ssb5[Gene3D_topo == "3.20.20"]$orf)) # 47
length(unique(domains_cath_bukau5$orf)) # 530
length(unique(domains_cath_bukau5[Gene3D_topo == "2.130.10"]$orf)) # 11
length(unique(domains_cath_bukau5[Gene3D_topo == "3.40.50"]$orf)) # 143
length(unique(domains_cath_bukau5[Gene3D_topo == "1.25.10"]$orf)) # 9
length(unique(domains_cath_bukau5[Gene3D_topo == "3.30.70"]$orf)) # 29
length(unique(domains_cath_bukau5[Gene3D_topo == "3.20.20"]$orf)) # 30
temp <- domains_cath[domains_cath$orf %in% properties_background$orf]
length(unique(temp$orf)) # 3381
length(unique(temp[Gene3D_topo == "2.130.10"]$orf)) # 134
length(unique(temp[Gene3D_topo == "3.40.50"]$orf)) # 829
length(unique(temp[Gene3D_topo == "1.25.10"]$orf)) # 73
length(unique(temp[Gene3D_topo == "3.30.70"]$orf)) # 124
length(unique(temp[Gene3D_topo == "3.20.20"]$orf)) # 125

length(unique(domains_scop_tric5$orf)) # 466
length(unique(domains_scop_tric5[Domain == "SSF50978"]$orf)) # 26 ; 0.001751485762415849
length(unique(domains_scop_tric5[Domain == "SSF48371"]$orf)) # 23 ; 0.01347074313269091
length(unique(domains_scop_tric5[Domain == "SSF51445" | Domain == "SSF51569"]$orf)) # 16 ; 8.045684761905192e-06
length(unique(domains_scop_tric5[Domain == "SSF54928"]$orf)) # 11 ; 0.1660725553868298
length(unique(domains_scop_ssb5$orf)) # 1049
length(unique(domains_scop_ssb5[Domain == "SSF50978"]$orf)) # 33 ; 0.5112903694979751
length(unique(domains_scop_ssb5[Domain == "SSF48371"]$orf)) # 40 ; 0.06774655455990033
length(unique(domains_scop_ssb5[Domain == "SSF51445" | Domain == "SSF51569"]$orf)) # 21 ; 0.0007398319155130851
length(unique(domains_scop_ssb5[Domain == "SSF54928"]$orf)) # 24 ; 0.06120642394989462
length(unique(domains_scop_bukau5$orf)) # 547
length(unique(domains_scop_bukau5[Domain == "SSF50978"]$orf)) # 9 ; 0.01634871518114181 under
length(unique(domains_scop_bukau5[Domain == "SSF48371"]$orf)) # 13 ; 0.1819077684581061
length(unique(domains_scop_bukau5[Domain == "SSF51445" | Domain == "SSF51569"]$orf)) # 13 ; 0.003053931699104894
length(unique(domains_scop_bukau5[Domain == "SSF54928"]$orf)) # 12 ; 0.2176905448078913
temp <- domains_scop[domains_scop$orf %in% properties_background$orf]
length(unique(temp$orf)) # 3561
length(unique(temp[Domain == "SSF50978"]$orf)) # 111 
length(unique(temp[Domain == "SSF48371"]$orf)) # 110 
length(unique(temp[Domain == "SSF51445" | Domain == "SSF51569"]$orf)) # 38 
length(unique(temp[Domain == "SSF54928"]$orf)) # 61 



# Significance of early vs late TRiC binding
tric_Early <- tric_peaks5[position_norm < 0.75]
tric_Late <- tric_peaks5[position_norm >= 0.75]
tric_Early_max <- tric_peaks5_max[position_norm < 0.75]
tric_Late_max <- tric_peaks5_max[position_norm >= 0.75]
tric_Low_max <- tric_peaks5_max[tric_odds < 4.821]
tric_High_max <- tric_peaks5_max[tric_odds >= 4.821]
tric_Low <- tric_peaks5[tric_odds < 4.141]
tric_High <- tric_peaks5[tric_odds >= 4.141]
temp <- as.data.table(unique(tric_Early_max$orf))
temp1 <- domains_cath_tric5_max[(domains_cath_tric5_max$orf %in% temp$V1)]
temp2 <- as.data.table(unique(tric_Late_max$orf))
temp3 <- domains_cath_tric5_max[(domains_cath_tric5_max$orf %in% temp2$V1)]
temp <- as.data.table(unique(tric_Low_max$orf))
temp1 <- domains_cath_tric5_max[(domains_cath_tric5_max$orf %in% temp$V1)]
temp2 <- as.data.table(unique(tric_High_max$orf))
temp3 <- domains_cath_tric5_max[(domains_cath_tric5_max$orf %in% temp2$V1)]
length(unique(temp1$orf)) # 301 for early TRiC binding: (all peaks)
length(unique(temp1[Gene3D_topo == "2.130.10"]$orf)) # 13
length(unique(temp1[Gene3D_topo == "3.40.50"]$orf)) # 99
length(unique(temp1[Gene3D_topo == "1.25.10"]$orf)) # 7
length(unique(temp1[Gene3D_topo == "3.20.20"]$orf)) # 23
length(unique(temp1[Gene3D_topo == "3.30.70"]$orf)) # 16
length(unique(temp3$orf)) # 245 for late TRiC binding: (all peaks)
length(unique(temp3[Gene3D_topo == "2.130.10"]$orf)) # 22
length(unique(temp3[Gene3D_topo == "3.40.50"]$orf)) # 76
length(unique(temp3[Gene3D_topo == "1.25.10"]$orf)) # 10
length(unique(temp3[Gene3D_topo == "3.20.20"]$orf)) # 18
length(unique(temp3[Gene3D_topo == "3.30.70"]$orf)) # 8
length(unique(domains_cath_tric5$orf)) # 449
length(unique(domains_cath_tric5[Gene3D_topo == "2.130.10"]$orf)) # 31
length(unique(domains_cath_tric5[Gene3D_topo == "3.40.50"]$orf)) # 139
length(unique(domains_cath_tric5[Gene3D_topo == "1.25.10"]$orf)) # 14
length(unique(domains_cath_tric5[Gene3D_topo == "3.20.20"]$orf)) # 32
length(unique(domains_cath_tric5[Gene3D_topo == "3.30.70"]$orf)) # 22


### Assign enrichment sites to domains
domains_cath_tric5[, peak.end := peak - End]
domains_cath_tric5[, peak.start := peak - Start]
domains_cath_tric5[, peak.ds := peak - ds_domain]
domains_cath_tric5[, length.ratio := peak.start / DomainLength]
setkeyv(domains_cath_tric5, c("position1"))
domains_cath_tric5_max[, peak.end := peak - End]
domains_cath_tric5_max[, peak.start := peak - Start]
domains_cath_tric5_max[, peak.ds := peak - ds_domain]
domains_cath_tric5_max[, length.ratio := peak.start / DomainLength]
setkeyv(domains_cath_tric5_max, c("position1"))
domains_cath_ssb5[, peak.end := peak - End]
domains_cath_ssb5[, peak.start := peak - Start]
domains_cath_ssb5[, peak.ds := peak - ds_domain]
domains_cath_ssb5[, length.ratio := peak.start / DomainLength]
setkeyv(domains_cath_ssb5, c("position1"))
domains_cath_ssb5_max[, peak.end := peak - End]
domains_cath_ssb5_max[, peak.start := peak - Start]
domains_cath_ssb5_max[, peak.ds := peak - ds_domain]
domains_cath_ssb5_max[, length.ratio := peak.start / DomainLength]
setkeyv(domains_cath_ssb5_max, c("position1"))
domains_scop_tric5[, peak.end := peak - End]
domains_scop_tric5[, peak.start := peak - Start]
domains_scop_tric5[, peak.ds := peak - ds_domain]
domains_scop_tric5[, length.ratio := peak.start / DomainLength]
setkeyv(domains_scop_tric5, c("position1"))
domains_scop_tric5_max[, peak.end := peak - End]
domains_scop_tric5_max[, peak.start := peak - Start]
domains_scop_tric5_max[, peak.ds := peak - ds_domain]
domains_scop_tric5_max[, length.ratio := peak.start / DomainLength]
setkeyv(domains_scop_tric5_max, c("position1"))
domains_scop_ssb5[, peak.end := peak - End]
domains_scop_ssb5[, peak.start := peak - Start]
domains_scop_ssb5[, peak.ds := peak - ds_domain]
domains_scop_ssb5[, length.ratio := peak.start / DomainLength]
setkeyv(domains_scop_ssb5, c("position1"))
domains_scop_ssb5_max[, peak.end := peak - End]
domains_scop_ssb5_max[, peak.start := peak - Start]
domains_scop_ssb5_max[, peak.ds := peak - ds_domain]
domains_scop_ssb5_max[, length.ratio := peak.start / DomainLength]
setkeyv(domains_scop_ssb5_max, c("position1"))
domains_cath_bukau5[, peak.end := peak - End]
domains_cath_bukau5[, peak.start := peak - Start]
domains_cath_bukau5[, peak.ds := peak - ds_domain]
domains_cath_bukau5[, length.ratio := peak.start / DomainLength]
setkeyv(domains_cath_bukau5, c("position1"))
domains_cath_bukau5_max[, peak.end := peak - End]
domains_cath_bukau5_max[, peak.start := peak - Start]
domains_cath_bukau5_max[, peak.ds := peak - ds_domain]
domains_cath_bukau5_max[, length.ratio := peak.start / DomainLength]
setkeyv(domains_cath_bukau5_max, c("position1"))


domains_cath_tric5_start <- domains_cath_tric5[peak.start > 30 & (peak.ds < 30 | is.na(peak.ds)), .SD[which.min(peak.start)], by = position1]
domains_cath_tric5_max_start <- domains_cath_tric5_max[peak.start > 30 & (peak.ds < 30 | is.na(peak.ds)), .SD[which.min(peak.start)], by = position1]
domains_scop_tric5_start <- domains_scop_tric5[peak.start > 30 & (peak.ds < 30 | is.na(peak.ds)), .SD[which.min(peak.start)], by = position1]
domains_scop_tric5_max_start <- domains_scop_tric5_max[peak.start > 30 & (peak.ds < 30 | is.na(peak.ds)), .SD[which.min(peak.start)], by = position1]
domains_cath_ssb5_start <- domains_cath_ssb5[peak.start > 30 & (peak.ds < 30 | is.na(peak.ds)), .SD[which.min(peak.start)], by = position1]
domains_cath_ssb5_max_start <- domains_cath_ssb5_max[peak.start > 30 & (peak.ds < 30 | is.na(peak.ds)), .SD[which.min(peak.start)], by = position1]
domains_scop_ssb5_start <- domains_scop_ssb5[peak.start > 30 & (peak.ds < 30 | is.na(peak.ds)), .SD[which.min(peak.start)], by = position1]
domains_scop_ssb5_max_start <- domains_scop_ssb5_max[peak.start > 30 & (peak.ds < 30 | is.na(peak.ds)), .SD[which.min(peak.start)], by = position1]
domains_cath_tric5_max_start[, domainend_norm := End / length]
domains_cath_tric5_start[, domainend_norm := End / length]
domains_scop_tric5_max_start[, domainend_norm := End / length]
domains_scop_tric5_start[, domainend_norm := End / length]
domains_cath_ssb5_max_start[, domainend_norm := End / length]
domains_cath_ssb5_start[, domainend_norm := End / length]
domains_scop_ssb5_max_start[, domainend_norm := End / length]
domains_scop_ssb5_start[, domainend_norm := End / length]
domains_cath_tric5_max[, domainend_norm := End / length]
domains_cath_tric5[, domainend_norm := End / length]
domains_scop_tric5_max[, domainend_norm := End / length]
domains_scop_tric5[, domainend_norm := End / length]
domains_cath_ssb5_max[, domainend_norm := End / length]
domains_cath_ssb5[, domainend_norm := End / length]
domains_scop_ssb5_max[, domainend_norm := End / length]
domains_scop_ssb5[, domainend_norm := End / length]
domains_cath_bukau5_start <- domains_cath_bukau5[peak.start > 30 & (peak.ds < 30 | is.na(peak.ds)), .SD[which.min(peak.start)], by = position1]
domains_cath_bukau5_start[, domainend_norm := End / length]
domains_cath_bukau5[, domainend_norm := End / length]
domains_cath_bukau5_max_start <- domains_cath_bukau5_max[peak.start > 30 & (peak.ds < 30 | is.na(peak.ds)), .SD[which.min(peak.start)], by = position1]
domains_cath_bukau5_max_start[, domainend_norm := End / length]
domains_cath_bukau5_max[, domainend_norm := End / length]


### Significance compared to random position
random <- NULL
for (i in 1:nrow(domains_cath_tric5_max_start)) {
  start <- domains_cath_tric5_max_start[i]$Start + 30
  if (is.na(domains_cath_tric5_max_start[i]$ds_domain)) { 
    end <- domains_cath_tric5_max_start[i]$length
  } else { end <- domains_cath_tric5_max_start[i]$ds_domain + 30 }
  random1 <- sample(start:end, size = 1, replace = T)
  random <- c(random, random1)
}
domains_cath_tric5_max_start$random <- random

random <- NULL
for (i in 1:nrow(domains_scop_tric5_start)) {
  start <- domains_scop_tric5_start[i]$Start + 30
  if (is.na(domains_scop_tric5_start[i]$ds_domain)) { 
    end <- domains_scop_tric5_start[i]$length
  } else { end <- domains_scop_tric5_start[i]$ds_domain + 30 }
  random1 <- sample(start:end, size = 1, replace = T)
  random <- c(random, random1)
}
domains_scop_tric5_start$random <- random

random <- NULL
for (i in 1:nrow(domains_cath_ssb5_max_start)) {
  start <- domains_cath_ssb5_max_start[i]$Start + 30
  if (is.na(domains_cath_ssb5_max_start[i]$ds_domain)) { 
    end <- domains_cath_ssb5_max_start[i]$length
  } else { end <- domains_cath_ssb5_max_start[i]$ds_domain + 30 }
  random1 <- sample(start:end, size = 1, replace = T)
  random <- c(random, random1)
}
domains_cath_ssb5_max_start$random <- random

random <- NULL
for (i in 1:nrow(domains_scop_ssb5_max_start)) {
  start <- domains_scop_ssb5_max_start[i]$Start + 30
  if (is.na(domains_scop_ssb5_max_start[i]$ds_domain)) { 
    end <- domains_scop_ssb5_max_start[i]$length
  } else { end <- domains_scop_ssb5_max_start[i]$ds_domain + 30 }
  random1 <- sample(start:end, size = 1, replace = T)
  random <- c(random, random1)
}
domains_scop_ssb5_max_start$random <- random

domains_cath_tric5_start[, random.start := random - Start]
domains_cath_tric5_start[, random.length.ratio := random.start / DomainLength]
domains_cath_tric5_max_start[, random.start := random - Start]
domains_cath_tric5_max_start[, random.length.ratio := random.start / DomainLength]
domains_scop_tric5_start[, random.start := random - Start]
domains_scop_tric5_start[, random.length.ratio := random.start / DomainLength]
domains_scop_tric5_max_start[, random.start := random - Start]
domains_scop_tric5_max_start[, random.length.ratio := random.start / DomainLength]
domains_cath_ssb5_start[, random.start := random - Start]
domains_cath_ssb5_start[, random.length.ratio := random.start / DomainLength]
domains_cath_ssb5_max_start[, random.start := random - Start]
domains_cath_ssb5_max_start[, random.length.ratio := random.start / DomainLength]
domains_scop_ssb5_start[, random.start := random - Start]
domains_scop_ssb5_start[, random.length.ratio := random.start / DomainLength]
domains_scop_ssb5_max_start[, random.start := random - Start]
domains_scop_ssb5_max_start[, random.length.ratio := random.start / DomainLength]


### Enrichment centered around domain end
temp <- as.data.table(unique(domains_cath_tric5_max_start$domain_alias))
temp2 <- domains_cath_tric5_max[domains_cath_tric5_max$domain_alias %in% temp$V1, c(1:17,22:25,37:46)]
setkeyv(temp2, c("orf"))
domains_cath_tric5_max_align <- tric_substrates5_dt[temp2]
domains_cath_tric5_max_align[, adjusted_start := position - Start]
domains_cath_tric5_max_align[, adjusted_end := position - End]
domains_cath_tric5_max_align[, adjusted_norm := (position - Start) / DomainLength]
temp <- as.data.table(unique(domains_cath_tric5_start$domain_alias))
temp2 <- domains_cath_tric5[domains_cath_tric5$domain_alias %in% temp$V1, c(1:17,22:25,37:44)]
setkeyv(temp2, c("orf"))
domains_cath_tric5_align <- tric_substrates5_dt[temp2, allow.cartesian = TRUE]
domains_cath_tric5_align[, adjusted_start := position - Start]
domains_cath_tric5_align[, adjusted_end := position - End]
domains_cath_tric5_align[, adjusted_norm := (position - Start) / DomainLength]

temp <- as.data.table(unique(domains_cath_ssb5_max_start$domain_alias))
temp2 <- domains_cath_ssb5_max[domains_cath_ssb5_max$domain_alias %in% temp$V1, c(1:17,22:25,37:46)]
setkeyv(temp2, c("orf"))
domains_cath_ssb5_max_align <- ssb_substrates5_dt[temp2]
domains_cath_ssb5_max_align[, adjusted_start := position - Start]
domains_cath_ssb5_max_align[, adjusted_end := position - End]
domains_cath_ssb5_max_align[, adjusted_norm := (position - Start) / DomainLength]
temp <- as.data.table(unique(domains_cath_ssb5_start$domain_alias))
temp2 <- domains_cath_ssb5[domains_cath_ssb5$domain_alias %in% temp$V1, c(1:17,22:25,37:44)]
setkeyv(temp2, c("orf"))
domains_cath_ssb5_align <- ssb_substrates5_dt[temp2, allow.cartesian = TRUE]
domains_cath_ssb5_align[, adjusted_start := position - Start]
domains_cath_ssb5_align[, adjusted_end := position - End]
domains_cath_ssb5_align[, adjusted_norm := (position - Start) / DomainLength]

