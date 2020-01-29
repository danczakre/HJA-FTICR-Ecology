### Script to correlate significantly different transformations

library(Hmisc)

setwd("/path/to/ICR_bNTI_results/")
bNTI = read.csv("HJ_Andrews_TWCD_bNTI_999.csv", row.names = 1)
trans = read.csv("HJ_Andrews_Trans_Profiles.csv", row.names = 1)
factors = read.csv("HJ_Andrews_Metadata.csv", row.names = 1)


# ####################### #
#### Preliminary setup ####
# ####################### #

# Removing outlier sample
factors = factors[-which(row.names(factors) %in% "PP48.000012"),]

# Clean up transformation sheet
trans = trans[,-which(colnames(trans) %in% "Mass")]
colnames(trans) = gsub("Sample_", "", colnames(trans))

# Checking row names consistency between molecular info and data
if(identical(x = colnames(bNTI), y = colnames(trans)) == FALSE){
  stop("Something is incorrect: the data and peak counts don't match")
}


# Converting trans data to rel abund
trans = as.data.frame(apply(trans, 2, function(x) x/sum(x)))


# ############################ #
#### Correlate average bNTI ####
# ############################ #

avg.bNTI = bNTI
avg.bNTI[upper.tri(avg.bNTI)] = t(avg.bNTI)[upper.tri(avg.bNTI)]
avg.bNTI = apply(avg.bNTI, 2, mean, na.rm = T)

trans = t(trans)
avg.cor = data.frame(Correlation = rep(NA, length(trans[1,])*3), r.value = NA, p.value = NA)

for(i in 1:length(trans[1,])){
  
  # Bulk correlation
  temp.cor = rcorr(trans[,i], avg.bNTI, type = "spearman")
  
  avg.cor$Correlation[i] = paste0("Bulk ", colnames(trans)[i], " Correlation")
  avg.cor$r.value[i] = temp.cor$r[1,2]
  avg.cor$p.value[i] = temp.cor$P[1,2]
  
  # PP48 correlation
  temp.bNTI = avg.bNTI[grep("PP48", names(avg.bNTI))]
  temp.trans = trans[grep("PP48", row.names(trans)),]
  temp.cor = rcorr(temp.trans[,i], temp.bNTI, type = "spearman")
  
  avg.cor$Correlation[i+length(trans[1,])] = paste0("PP48 ", colnames(trans)[i], " Correlation")
  avg.cor$r.value[i+length(trans[1,])] = temp.cor$r[1,2]
  avg.cor$p.value[i+length(trans[1,])] = temp.cor$P[1,2]
  
  # SW48 correlation
  temp.bNTI = avg.bNTI[grep("SW48", names(avg.bNTI))]
  temp.trans = trans[grep("SW48", row.names(trans)),]
  temp.cor = rcorr(temp.trans[,i], temp.bNTI, type = "spearman")
  
  avg.cor$Correlation[i+(length(trans[1,])*2)] = paste0("SW48 ", colnames(trans)[i], " Correlation")
  avg.cor$r.value[i+(length(trans[1,])*2)] = temp.cor$r[1,2]
  avg.cor$p.value[i+(length(trans[1,])*2)] = temp.cor$P[1,2]
  
  rm("temp.bNTI", "temp.trans", "temp.cor")
}


# ######################################### #
#### Correlating by transformation group ####
# ######################################### #

CHO = trans[,-grep("P|S|N|ine|phan|Aspartic|Glutamic|amination|uracil|urea|co-enzyme", colnames(trans))]

N = trans[,grep("N|ine|phan|urea|biotinyl|co-enzyme", colnames(trans))]
N = N[,-grep("N/A|glucose-N-|Na_", colnames(N))]

S = trans[,grep("S|Cysteine|Methionine|co-enzyme", colnames(trans))]
S = S[,-grep("Serine", colnames(S))]

P = trans[,grep("P|co-enzyme", colnames(trans))]
P = P[,-grep("2ndIP|Phenylalanine|Proline", colnames(P))]

freq.by.group = data.frame(CHO = rowSums(CHO), N = rowSums(N), S = rowSums(S), P = rowSums(P))
cor.by.group = data.frame(Correlation = rep(NA, length(freq.by.group[1,])*3), r.value = NA, p.value = NA)

for(i in 1:ncol(freq.by.group)){
  # Bulk correlation
  temp.cor = rcorr(freq.by.group[,i], avg.bNTI, type = "spearman")
  
  cor.by.group$Correlation[i] = paste0("Bulk ", colnames(freq.by.group)[i], " Correlation")
  cor.by.group$r.value[i] = temp.cor$r[1,2]
  cor.by.group$p.value[i] = temp.cor$P[1,2]
  
  # PP48 correlation
  temp.bNTI = avg.bNTI[grep("PP48", names(avg.bNTI))]
  temp.freq = freq.by.group[grep("PP48", row.names(freq.by.group)),]
  temp.cor = rcorr(temp.freq[,i], temp.bNTI, type = "spearman")
  
  cor.by.group$Correlation[i+length(freq.by.group[1,])] = paste0("PP48 ", colnames(freq.by.group)[i], " Correlation")
  cor.by.group$r.value[i+length(freq.by.group[1,])] = temp.cor$r[1,2]
  cor.by.group$p.value[i+length(freq.by.group[1,])] = temp.cor$P[1,2]
  
  # SW48 correlation
  temp.bNTI = avg.bNTI[grep("SW48", names(avg.bNTI))]
  temp.freq = freq.by.group[grep("SW48", row.names(freq.by.group)),]
  temp.cor = rcorr(temp.freq[,i], temp.bNTI, type = "spearman")
  
  cor.by.group$Correlation[i+(length(freq.by.group[1,])*2)] = paste0("SW48 ", colnames(freq.by.group)[i], " Correlation")
  cor.by.group$r.value[i+(length(freq.by.group[1,])*2)] = temp.cor$r[1,2]
  cor.by.group$p.value[i+(length(freq.by.group[1,])*2)] = temp.cor$P[1,2]
  
  rm("temp.bNTI", "temp.freq", "temp.cor")
}


# ##################################### #
#### Cleaning up and writing results ####
# ##################################### #

avg.cor$p.value = round(avg.cor$p.value, digits = 6)
cor.by.group$p.value = round(cor.by.group$p.value, digits = 5)

avg.cor = avg.cor[which(avg.cor$p.value < 0.05),] # Given the sheer number of transformation correlations, filtering the results to significant values

write.csv(avg.cor, "HJA_bNTI-Transf_Correlations.csv", quote = F, row.names = F)
write.csv(cor.by.group, "HJA_bNTI-TransfGroup_Correlations.csv", quote = F, row.names = F)
