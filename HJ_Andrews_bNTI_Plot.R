### Script to test individual null models
library(reshape2); library(ggplot2); library(ggthemes)

setwd("~/Documents/HJ Andrews Diel Sampling/FT-ICR Analyses (NoRI - No Outlier)/Null Models/")
bNTI = read.csv("HJ_Andrews_NoRI_NoOut_Weighted_bNTI_999.csv", row.names = 1)
factors = read.csv("~/Documents/HJ Andrews Diel Sampling/Factor Sheets/HJ_Andrews_NoRep_Sample_Sheet.csv", row.names = 1)

### Preprocessing and clean up
# Dropping PP48-000012 (outlier)
factors = factors[-which(row.names(factors) %in% "PP48.000012"),]

# Fixing date and time
factors$Date_Time = as.POSIXct(paste(factors$Date, factors$Time), format = "%m/%d/%Y %H:%M")
factors$Date_Time = gsub("0018", "2018", factors$Date_Time)
factors$Date_Time = gsub(":..$", "", factors$Date_Time)

# Creating ggplot themes
hori_x_theme = theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())

vert_x_theme = theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.6),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())



### Actual analysis
# Reflecting null matrices
bNTI[upper.tri(bNTI)] = t(bNTI)[upper.tri(bNTI)]

# Plotting mean thru time
avg.bNTI = apply(bNTI, 2, function(x) mean(x, na.rm = T))
avg.bNTI = data.frame(Samples = names(avg.bNTI), bNTI = avg.bNTI, Date = factors$Date_Time, Type = factors$Sample)

ggplot(data = avg.bNTI, aes(x = as.POSIXct(Date), y = bNTI, group = Type))+
  geom_point(aes(color = Type))+
  geom_line(aes(color = Type))+
  xlab(NULL)+
  scale_color_stata()+
  hori_x_theme

# Melting data
bNTI = melt(as.matrix(bNTI))

# Cleaning data
bNTI = bNTI[!is.na(bNTI$value),]

# Adding Location Information
bNTI$Type = "Pushpoint"
bNTI$Type[grep("SW48", bNTI$Var2)] = "Surface Water"

# Counting the proportion of processes
proc = data.frame(hom.sel = NA, stoc = NA, var.sel = NA)
proc$hom.sel = round(length(which(bNTI$value < -2))/length(bNTI$value), digits = 3)
proc$stoc = round(length(which(abs(bNTI$value) < 2))/length(bNTI$value), digits = 3)
proc$var.sel = round(length(which(bNTI$value > 2))/length(bNTI$value), digits = 3)

### All Samples
# Boxplots by group
ggplot(data = bNTI, aes(x = Type, y = value))+
  geom_boxplot(aes(group = Type, color = Type))+
  geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
  ggtitle("All Samples")+
  scale_color_stata()+
  hori_x_theme

### Within-group comparisons only
temp = bNTI[-c(grep("SW48", bNTI$Var1), grep("SW48", bNTI$Var2)),]
temp = rbind(temp, bNTI[intersect(grep("SW48", bNTI$Var1), grep("SW48", bNTI$Var2)),])

# Boxplots by group
ggplot(data = temp, aes(x = Type, y = value))+
  geom_boxplot(aes(group = Type, color = Type))+
  geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
  ggtitle("Without Cross-comparisons")+
  xlab(NULL)+
  scale_color_stata()+
  hori_x_theme

# Within-group stats
nocross.mwu = wilcox.test(temp$value~temp$Type)
