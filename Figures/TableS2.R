### Table S2 ###

### Binding peaks

# TRiC
data_tric <- tric_peaks5[, c(1,2,4:10,12)]

temp1 <- read.csv("/Users/KevinStein/Desktop/Lab/Bioinformatics/ProfilingData/AKK/Annotation/Localization_Uniprot.csv", header = T)
temp1 <- as.data.table(temp1)
i <- cbind(match(data_tric$orf, temp1$orf))
data_tric <- cbind(data_tric, name = temp1[i]$Names)

i <- cbind(match(data_tric$position1, domains_cath_tric5_start$position1))
data_tric <- cbind(data_tric, tric_domain = domains_cath_tric5_start[i]$Gene3D_domain)
data_tric <- cbind(data_tric, tric_interpro = domains_cath_tric5_start[i]$InterproDomainName)
data_tric <- cbind(data_tric, tric_start = domains_cath_tric5_start[i]$Start)
data_tric <- cbind(data_tric, tric_end = domains_cath_tric5_start[i]$End)

data_tric[is.na(data_tric)] <- ""

write.csv(data_tric, "/Users/KevinStein/Desktop/data_tric.csv")


# Ssb
data_ssb <- ssb_peaks5[, c(1,2,4:10,12)]

i <- cbind(match(data_ssb$orf, temp1$orf))
data_ssb <- cbind(data_ssb, name = temp1[i]$Names)

i <- cbind(match(data_ssb$position1, domains_cath_ssb5_start$position1))
data_ssb <- cbind(data_ssb, ssb_domain = domains_cath_ssb5_start[i]$Gene3D_domain)
data_ssb <- cbind(data_ssb, ssb_interpro = domains_cath_ssb5_start[i]$InterproDomainName)
data_ssb <- cbind(data_ssb, ssb_start = domains_cath_ssb5_start[i]$Start)
data_ssb <- cbind(data_ssb, ssb_end = domains_cath_ssb5_start[i]$End)

data_ssb[is.na(data_ssb)] <- ""

write.csv(data_ssb, "/Users/KevinStein/Desktop/data_ssb.csv")

