### HJ Andrews data processing

options(digits = 10)

require(vegan); require(picante); require(phytools); require(pez) # Loading in packages for ecology/dendrogram analyses
require(reshape2); require(ggplot2); require(ggthemes) # Packages for managing plots
require(plot3D) # Needed for 3D PCoA plots


# ################## #
#### Load in data ####
# ################## #

# Set directory
setwd("/path/to/ICR_data")

# Standard data
data = read.csv("Processed_HJ_Andrews_cleaned_Data.csv", row.names = 1) # Load in peak counts
mol = read.csv("Processed_HJ_Andrews_cleaned_Mol.csv", row.names = 1) # Load in mol. info
meta = read.csv("HJ_Andrews_Metadata.csv", row.names = 1) # Load in metadata

# Load in geochem.
geo = read.csv("HJ_Andrews_Geochem.csv", row.names = 1)

# Transformation profile
trans.pro = read.csv("HJ_Andrews_Trans_Profiles.csv", row.names = 1) # Load in trans. profile

# Load in tree
tree = read.tree("HJ_Andrews_TWCD_UPGMA.tre")


# ###################### #
#### Define functions ####
# ###################### #

beta.thru.time = function(beta.div){
  # Converting input to matrix
  beta.div = as.matrix(beta.div)
  
  # Creating dummy vector
  thru.time = NULL
  first.time = which(factors$Date_Time %in% unique(factors$Date_Time)[1]) # Finding first time points in dataset
  
  for(curr.samp in unique(factors$Sample_Type)){
    
    # Find necessary data
    w = which(factors$Sample_Type %in% curr.samp)
    
    # Removing first time point in order to track dissimilarity through time
    v = first.time[(first.time %in% w)] 
    w = w[!(w %in% first.time)]
    
    # Store data in temp data frame
    temp = data.frame(Date = factors$Date_Time[w], Sample = factors$Sample_Type[w], beta = beta.div[w,v])
    
    # Add corresponding DO and temp. data for correlations
    temp$DO = factors$DO[w]; temp$Temp = factors$Temp[w] 
    
    # Merge data
    thru.time = rbind(thru.time, temp)
    
    # Clean-up
    rm("temp", "w", "v")
  }
  
  rm("curr.samp")
  
  # Ensuring the date is in the correct order
  thru.time$Date = factor(as.character(thru.time$Date), levels = unique(thru.time$Date)[order(as.character(unique(thru.time$Date)))])
  
  # Plotting thru time data
  print(
    ggplot(data = thru.time, aes(x = as.POSIXct(Date), y = beta, group = Sample))+
      geom_point(aes(color = Sample))+
      geom_line(aes(color = Sample))+
      scale_color_stata()+
      xlab(NULL)+
      hori_x_theme
  )
  
  return(thru.time)
}

plot3D = function(pcoa){
  scores.pcoa = as.data.frame(pcoa$vectors) # Getting scores
  
  axis.font = list(family = "Arial, sans-serif", size = 17, color = "black")
  tick.font = list(family = "Arial, sans-serif", size = 14, color = "black")
  
  xaxis = list(title = paste0("PCoA1 (", round((pcoa$values$Relative_eig[1]*100), 2), "%)"),
               titlefont = axis.font, tickfont = tick.font, gridcolor = "black")
  yaxis = list(title = paste0("PCoA2 (", round((pcoa$values$Relative_eig[2]*100), 2), "%)"),
               titlefont = axis.font, tickfont = tick.font, gridcolor = "black")
  zaxis = list(title = paste0("PCoA3 (", round((pcoa$values$Relative_eig[3]*100), 2), "%)"),
               titlefont = axis.font, tickfont = tick.font, gridcolor = "black")
  legend = list(font = axis.font, x = 0, y = 0)
  
  plot_ly(x = scores.pcoa$Axis.1, y = scores.pcoa$Axis.2, z = scores.pcoa$Axis.3, type = "scatter3d", 
          mode = "markers", symbol = plot.factors, color = plot.factors, symbols = c('circle', 'diamond'),
          colors = stata_pal("s2color")(2))%>%
    layout(scene = list(xaxis = xaxis, yaxis = yaxis, zaxis = zaxis), legend = legend, 
           paper_bgcolor = 'rgb(243, 243, 243', plot_bgcolor = 'rgb(243, 243, 243')
}

# ############ #
#### Errors ####
# ############ #

# Checking row names consistency between molecular info and data
if(identical(x = row.names(data), y = row.names(mol)) == FALSE){
  stop("Something is incorrect: the mol. info and peak counts don't match")
}

# Determining whether isotopic peaks are present
if(length(which(mol$C13 == 1)) > 0){
  stop("Isotopic signatures weren't removed, please process using the ftmsRanalysis R package")
}

# Data should be presence/absence
if(max(data) > 1){
  stop("Data is not presence/absence, data is being convert to presence/absence")
  data[data > 1] = 1
}

# Miscellaneous fixes
if(length(grep("QC_SRFAII", colnames(data))) > 0){
  stop("Suwannee River standards are still in the data")
}

if(length(grep("rep1|rep2", colnames(data))) > 0){
  stop("Technical replicates are still in the data")
}



# ######################## #
#### Preprocessing data ####
# ######################## #

# Removing the outlier (PP48.000012) from other datasets
meta = meta[-which(row.names(meta) %in% "PP48.000012"),]
geo = geo[-which(row.names(geo) %in% "PP48.000012"),]

# Creating a factors sheet
factors = colnames(data) # The factors sheet will be used to help plot
factors = as.data.frame(factors)

# Adding in date and time, then merging them
factors$Date = meta[gsub("_.*$", "", factors$factors), "Date", drop = T]
factors$Time = meta[gsub("_.*$", "", factors$factors), "Time", drop = T]

meta$Date_Time = as.POSIXct(paste(meta$Date, meta$Time), format = "%m/%d/%Y %H:%M")
factors$Date_Time = meta[gsub("_.*$", "", factors$factors), "Date_Time", drop = T]
factors$Date_Time = gsub("0018", "2018", factors$Date_Time)
factors$Date_Time = gsub(":..$", "", factors$Date_Time)

# Adding in depths
factors$Depth = meta[gsub("_.*$", "", factors$factors), "Depth", drop = T]

# Providing a column name for samples
colnames(factors)[1] = c("Sample_Name")

# Adding in sample type information
factors$Sample_Type = meta[gsub("_.*$", "", factors$Sample_Name), "Sample", drop = T]
factors$WC_Height = meta[gsub("_.*$", "", factors$Sample_Name), "Water_Column_Height", drop = T]
factors$Stage = meta[gsub("_.*$", "", factors$Sample_Name), "Water_Stage", drop = T]
factors$DO = meta[gsub("_.*$", "", factors$Sample_Name), "Approx_DO_in_WC", drop = T]
factors$Temp = meta[gsub("_.*$", "", factors$Sample_Name), "Water_Temp_50per", drop = T]

# Merging DO and water stage into geochem.
geo$DO = factors$DO
geo$Stage = factors$Stage

# Removing problematic or redundant variables
variables = c("P_mg_per_L_as_PO4", "F_mg_per_L", "K_mg_per_L")
geo = geo[,-which(colnames(geo) %in% variables)]
rm("variables")

# Cleaning up transformation profile
trans.pro = trans.pro[,-which(colnames(trans.pro) %in% "Mass")] # Removing mass column
trans.pro = as.data.frame(apply(trans.pro, 2, function(x) x/sum(x))) # Convert to relative abundance

# ggplot2 theme objects
hori_x_theme = theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())

vert_x_theme =  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.6),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())


# ################################### #
#### Analyzing chemical properties ####
# ################################### #
# Looping through each sample to obtain some summary stats of the peaks
chem.prop = data.frame(GFE = rep(NA, length(data[1,])), AI_Mod = NA, DBE = NA, Peaks = NA, 
                       row.names = colnames(data), stringsAsFactors = F) # Object to store properties

for(i in 1:ncol(data)){
  temp = data[which(data[,i] > 0), i, drop = F] # Need to keep names, looking at columns
  temp = mol[row.names(temp),]
  
  chem.prop$GFE[i] = mean(temp$GFE, na.rm = T)
  chem.prop$AI_Mod[i] = mean(temp$AI_Mod, na.rm = T)
  chem.prop$DBE[i] = mean(temp$DBE, na.rm = T)
  chem.prop$Peaks[i] = length(temp$P)
} # I'm not sure how to do this without the for-loop, but I'm simply just finding the mean for peak stats

chem.melt$Date_Time = factors$Date_Time
chem.melt$Type = factors$Sample_Type


# Plotting the data through time
chem.melt = melt(chem.melt, id.vars = c("Date_Time", "Type")) # Melting this for ggplot

print(
  ggplot(data = chem.melt, aes(x = as.POSIXct(Date_Time), y = value, group = Type))+
    geom_point(aes(color = Type))+
    geom_line(aes(color = Type))+
    facet_grid(variable~., scales = "free_y")+
    scale_color_stata()+
    xlab(NULL)+
    hori_x_theme+
    ggtitle(label = paste("Variables through time"))
)

rm("chem.melt")

# Measuring statistical differences for chemical properties
chem.stats = data.frame(MWU = rep(NA, ncol(chem.prop)), p.value = NA, row.names = colnames(chem.prop))

for(i in 1:ncol(chem.prop)){
  chem.stats[i,] = c(wilcox.test(chem.prop[,i]~factors$Sample_Type)$statistic, wilcox.test(chem.prop[,i]~factors$Sample_Type)$p.value)
}


# ############################### #
#### Analyzing alpha diversity ####
# ############################### #
tree = midpoint.root(tree) # Rooting the tree for consistent results
div = data.frame(row.names = colnames(data), Date_Time = factors$Date_Time, Type = factors$Sample_Type, 
                 PD = NA, SR = NA, MPD = NA, MNTD = NA, VPD = NA, VNTD = NA) # Creaating empty data frame to store data

div[,c("PD", "SR")] = pd(t(data), tree) # Faith's PD and species richness

comp.comm = comparative.comm(tree, t(data)) # Creating a "pez" compatible object
div[,c("MPD", "MNTD", "VPD", "VNTD")] = generic.metrics(comp.comm, metrics = c(.mpd, .mntd, .vpd, .vntd)) # Processing other phylogenetic a-diversity metrics 
rm(comp.comm) # Removing this object as it is no longer necessary

div.melt = div
div.melt = melt(div.melt, id.vars = c("Date_Time", "Type")) # Converting to long format

# By location boxplot
print(
  ggplot(data = div.melt, aes(x = Type, y = value))+
    geom_boxplot(aes(color = Type))+
    facet_grid(variable~., scales = "free_y")+
    scale_color_stata()+
    hori_x_theme+
    ggtitle("Alpha diversity")
) # Plotting the data

rm("div.melt")

# Measuring statistical differences for diversity
div.stats = data.frame(MWU = rep(NA, ncol(div[,3:8])), p.value = NA, row.names = colnames(div[,3:8]))

for(i in 1:(ncol(div)-2)){
  div.stats[i,] = c(wilcox.test(div[,i+2]~factors$Sample_Type)$statistic, wilcox.test(div[,i+2]~factors$Sample_Type)$p.value)
}


# #################### #
#### Beta-diversity ####
# #################### #

# Creating distance matrix
dist = vegdist(x = t(data), method = "jaccard") # Using Jaccard for historical reasons (ICR data is often analyzed using it)

# 3D Plotting
jac.pcoa = ape::pcoa(dist) # Generating a Jaccard-based principal coordinate analysis

plot3D(jac.pcoa)

# Alternative 3D plot
scores.pcoa = as.data.frame(jac.pcoa$vectors)
with(scores.pcoa, points3D(x = Axis.1, y = Axis.2, z = Axis.3, phi = 10, theta = 35,
                          colvar = as.integer(factor(plot.factors)),
                          col = stata_pal("s2color")(2), pch = c(rep(16, 17), rep(17, 17)), cex = 2,  
                          clim = c(1, 2), ticktype = "detailed",
                          xlab = paste0("PCoA1 (", round((jac.pcoa$values$Relative_eig[1]*100), 2), "%)"), 
                          ylab = paste0("PCoA2 (", round((jac.pcoa$values$Relative_eig[2]*100), 2), "%)"),
                          zlab = paste0("PCoA3 (", round((jac.pcoa$values$Relative_eig[3]*100), 2), "%)"),
                          colkey = list(at = c(1.25,1.75), side = 1, 
                                        addlines = TRUE, length = 0.5, width = 0.5,
                                        labels = unique(plot.factors))))

# Determining whether differences are significant or not
jac.sig = adonis(dist~factors$Sample_Type, permutations = 9999)

rm("jac.pcoa", "scores.pcoa")

### Phylogenetic beta-diversity
## Matching data to the provided tree
phylo = match.phylo.data(tree, data) # Matching ICR dataset to the tree

## NMDS on tree-subsetted data
# Creating distance matrix
dist = vegdist(x = t(phylo$data), method = "jaccard") # Using Jaccard for historical reasons (ICR data is often analyzed using it)

# Plotting Jaccard NMDS
nms = metaMDS(dist, try = 40) # Determining NMDS
nms = as.data.frame(scores(nms)) # Conveting to scores
nms$Type = factors$Sample_Type # Adding meta-data

print(
  ggplot(data = nms, aes(x = NMDS1, y = NMDS2))+
    geom_point(aes(color = Type, shape = Type), size = 2.5)+
    scale_color_stata()+
    ggtitle(paste0("Subset NMDS"))+
    hori_x_theme
) # Plotting NMS graph

nms.time = beta.thru.time(dist)

## bMNTD Analysis
coph = cophenetic(phylo$phy)
bMNTD = comdistnt(t(phylo$data), coph, abundance.weighted = F, exclude.conspecifics = F)

bMNTD.pcoa = ape::pcoa(bMNTD)
bMNTD.scores = as.data.frame(bMNTD.pcoa$vectors)
bMNTD.scores$Type = factors$Sample_Type

print(
  ggplot(data = bMNTD.scores, aes(x = Axis.1, y = Axis.2))+
    geom_point(aes(color = Type, shape = Type), size = 2.5)+
    xlab(label = paste0("PCoA1 (", round((bMNTD.pcoa$values$Relative_eig[1]*100), 2), "%)"))+
    ylab(label = paste0("PCoA2 (", round((bMNTD.pcoa$values$Relative_eig[2]*100), 2), "%)"))+
    scale_color_stata()+
    ggtitle(paste0("βMNTD"))+
    hori_x_theme
) # Plotting the UniFrac PCoA

bMNTD.time = beta.thru.time(bMNTD)

# 3D Plot
plot3D(bMNTD.pcoa)

# Alternative 3D plot
scores.pcoa = as.data.frame(bMNTD.pcoa$vectors)
with(scores.pcoa, points3D(x = Axis.1, y = Axis.2, z = Axis.3, phi = 10, theta = 30,
                           colvar = as.integer(factor(plot.factors)),
                           col = stata_pal("s2color")(2), pch = c(rep(16, 17), rep(17, 17)), cex = 2,  
                           clim = c(1, 2), ticktype = "detailed",
                           xlab = paste0("PCoA1 (", round((bMNTD.pcoa$values$Relative_eig[1]*100), 2), "%)"), 
                           ylab = paste0("PCoA2 (", round((bMNTD.pcoa$values$Relative_eig[2]*100), 2), "%)"),
                           zlab = paste0("PCoA3 (", round((bMNTD.pcoa$values$Relative_eig[3]*100), 2), "%)"),
                           colkey = list(at = c(1.25,1.75), side = 1, 
                                         addlines = TRUE, length = 0.5, width = 0.5,
                                         labels = unique(plot.factors))))

# Stats
bmntd.sig = adonis(bMNTD~factors$Sample_Type, ppermutations = 9999)

rm("bMNTD.pcoa", "bMNTD.scores", "coph")


# ############################### #
#### Examining transformations ####
# ############################### #

# Jaccard-based NMDS
dist = vegdist(t(trans.pro), method = "jaccard")
nms = metaMDS(dist, trymax = 40)
nms = as.data.frame(scores(nms)) # Conveting to scores
nms$Type = factors$Sample_Type # Adding meta-data

print(
  ggplot(data = nms, aes(x = NMDS1, y = NMDS2))+
    geom_point(aes(color = Type, shape = Type), size = 2.5)+
    scale_color_stata()+
    ggtitle(paste0("Transformation NMDS"))+
    hori_x_theme
) # Plotting NMS graph

# Transformation dissimilarity through time
trans.time = beta.thru.time(dist)

# Stats
trans.sig = adonis(dist~factors$Sample_Type, permutations = 9999)