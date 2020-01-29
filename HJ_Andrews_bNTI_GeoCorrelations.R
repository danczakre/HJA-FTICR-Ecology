#### Script to correlate bNTI to external parameters

library(Hmisc) # Necessary for "rcorr" function

# ################## #
#### Load in data ####
# ################## #

# Set directory
setwd("~/Documents/HJ Andrews Diel Sampling/FT-ICR Analyses (NoRI - No Outlier)/Null Models/")

# Load in bNTI and geochem
bNTI = read.csv("HJ_Andrews_NoRI_NoOut_Weighted_bNTI_999.csv", row.names = 1)
ext = read.csv("~/Documents/HJ Andrews Diel Sampling/Geochem/HJ_Andrews_Geochem.csv", row.names = 1)
factors = read.csv("~/Documents/HJ Andrews Diel Sampling/Factor Sheets/HJ_Andrews_NoRep_Sample_Sheet.csv", row.names = 1)

# ####################### #
#### Preliminary setup ####
# ####################### #

# Removing the outlier from other data
ext = ext[-which(row.names(ext) %in% "PP48.000012"),]
factors = factors[-which(row.names(factors) %in% "PP48.000012"),]

# Checking row names consistency between molecular info and data
if(identical(x = row.names(bNTI), y = row.names(ext)) == FALSE){
  stop("Something is incorrect: the data and peak counts don't match")
}

# Removing problematic or redundant variables
variables = c("P_mg_per_L_as_PO4", "F_mg_per_L", "K_mg_per_L")
ext = ext[,-which(colnames(ext) %in% variables)]
rm("variables")

# Adding in approx. DO measurements
ext$DO = factors$Approx_DO_in_WC
ext$Stage = factors$Water_Stage


# ############################# #
#### Correlate with averages ####
# ############################# #

avg.bNTI = bNTI # Keeping the original bNTI information unaltered
avg.bNTI[upper.tri(avg.bNTI)] = t(avg.bNTI)[upper.tri(avg.bNTI)] # Mirroring the bNTI results
avg.bNTI = apply(avg.bNTI, 2, mean, na.rm = T) # Averaging by sample

avg.cor = data.frame(Correlation = rep(NA, length(ext[1,])*3), r.value = NA, p.value = NA) # Creating object to store corrrelations

for(i in 1:length(ext[1,])){
  
  # Bulk correlation
  temp.cor = rcorr(ext[,i], avg.bNTI, type = "spearman")
  
  avg.cor$Correlation[i] = paste0("Bulk ", colnames(ext)[i], " Correlation")
  avg.cor$r.value[i] = temp.cor$r[1,2]
  avg.cor$p.value[i] = temp.cor$P[1,2]
  
  # PP48 correlation
  temp.bNTI = avg.bNTI[grep("PP48", names(avg.bNTI))]
  temp.ext = ext[grep("PP48", row.names(ext)),]
  temp.cor = rcorr(temp.ext[,i], temp.bNTI, type = "spearman")
  
  avg.cor$Correlation[i+length(ext[1,])] = paste0("PP48 ", colnames(ext)[i], " Correlation")
  avg.cor$r.value[i+length(ext[1,])] = temp.cor$r[1,2]
  avg.cor$p.value[i+length(ext[1,])] = temp.cor$P[1,2]
  
  # SW48 correlation
  temp.bNTI = avg.bNTI[grep("SW48", names(avg.bNTI))]
  temp.ext = ext[grep("SW48", row.names(ext)),]
  temp.cor = rcorr(temp.ext[,i], temp.bNTI, type = "spearman")
  
  avg.cor$Correlation[i+(length(ext[1,])*2)] = paste0("SW48 ", colnames(ext)[i], " Correlation")
  avg.cor$r.value[i+(length(ext[1,])*2)] = temp.cor$r[1,2]
  avg.cor$p.value[i+(length(ext[1,])*2)] = temp.cor$P[1,2]
  
  rm("temp.bNTI", "temp.ext", "temp.cor")
}

avg.cor$p.value = round(avg.cor$p.value, digits = 6)

avg.cor = avg.cor[which(avg.cor$p.value < 0.05),]

write.csv(avg.cor, "HJA_Geochem_bNTI_Correlations.csv", quote = F, row.names = F)