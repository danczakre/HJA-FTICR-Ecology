### Merging null bMNTD matrices into an array and calculating bNTI

# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r; Stegen et al, 2013
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r; Stegen et al, 2013
# RED 2018; danczak.6@osu.edu

# File Name Labels
Sample_Name = "Altamaha_Alt_Conf"
tree_type = "Weighted"

# Load in libraries/functions
library(picante)
library(reshape2)
library(abind)

acomb = function(...) abind(..., along = 3)


# ################## #
#### Data Loading ####
# ################## #

setwd("/Users/danc783/Null Outputs/Altamaha_Alt_Conf/") # Working directory for FT-ICR data
data = read.csv("~/Documents/48 Hour Meta-analysis/Individual Datasets/Altamaha River/Alternate_Confident_Peaks/Processed_Altamaha_48Hour_Alt_Data_clean_noRepComp.csv", row.names = 1) # Importing the organismal data
tree = read.tree("~/Documents/48 Hour Meta-analysis/Individual Datasets/Altamaha River/Alternate_Confident_Peaks/Altamaha_Alt_Conf_48Hour_Weighted_All-Trans_UPGMA.tre") # Importing the tree


# ###################### #
#### bNTI Calculation ####
# ###################### #

# Converting data to presence/absence
data[data>1] = 1

# Matching the tree to peak data
phylo = match.phylo.data(tree, data)

# Calculating bMNTD for my samples
coph = cophenetic(phylo$phy)

bMNTD = as.matrix(comdistnt(t(phylo$data), coph, abundance.weighted = F, exclude.conspecifics = F))

# Merging the separate bMNTD files
files = list.files(path = paste(tree_type, "_Null_Output/", sep = "")
                   , pattern = "bMNTD_rep", full.names = T) # Listing files

rand.bMNTD = NULL # Dummy object

for(curr.file in files){
  temp = as.data.frame(read.csv(curr.file, row.names = 1))
  rand.bMNTD = c(rand.bMNTD, list(temp))
} # Merging individual 

rand.bMNTD = do.call(acomb, rand.bMNTD)
rm("curr.file", "files")


# Calculate bNTI
bNTI = matrix(c(NA), nrow = ncol(rand.bMNTD[,,1]), ncol = ncol(rand.bMNTD[,,1]))

for(i in 1:(ncol(rand.bMNTD[,,1])-1)){
  for(j in (i+1):ncol(rand.bMNTD[,,1])){
    m = rand.bMNTD[j,i,] # Just setting all the randomizations for a given comparison to a matrix
    bNTI[j,i] = ((bMNTD[j,i]-mean(m))/sd(m)) # The bNTI calculation
  }
}

dimnames(bNTI) = dimnames(rand.bMNTD[,,1])
rm("m", "j")

write.csv(bNTI, paste(Sample_Name, "_", tree_type, "_bNTI_999.csv", sep = ""), quote = F)
