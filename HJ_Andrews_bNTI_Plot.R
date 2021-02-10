### Script to test individual null models
library(reshape2); library(ggplot2); library(ggthemes); library(ggpubr)

setwd("/path/to/Data Folder")
bNTI = read.csv("HJ_Andrews_TWCD_bNTI_999.csv", row.names = 1)
factors = read.csv("HJ_Andrews_Metadata.csv", row.names = 1)

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

# Selecting within group dynamics
avg.bNTI = bNTI
avg.bNTI[grep("SW48", row.names(avg.bNTI)),grep("PP48", colnames(avg.bNTI))] = NA
avg.bNTI[grep("PP48", row.names(avg.bNTI)),grep("SW48", colnames(avg.bNTI))] = NA

# Plotting mean thru time
std.bNTI = apply(avg.bNTI, 2, function(x) sd(x, na.rm = T))
avg.bNTI = apply(avg.bNTI, 2, function(x) mean(x, na.rm = T))
avg.bNTI = data.frame(Samples = names(avg.bNTI), bNTI = avg.bNTI, StdDev = std.bNTI, Date = factors$Date_Time, Type = factors$Sample)

ggplot(data = avg.bNTI, aes(x = as.POSIXct(Date), y = bNTI, group = Type))+
  geom_point(aes(color = Type)) + geom_line(aes(color = Type))+
  geom_errorbar(aes(color = Type, ymin = bNTI-std.bNTI, ymax = bNTI+std.bNTI), width=2500,
                position=position_dodge(500))+
  xlab(NULL) + scale_color_stata()+
  hori_x_theme

# Plotting means as a boxplot
ggplot(data = avg.bNTI, aes(x = Type, y = bNTI))+
  geom_boxplot(aes(group = Type, color = Type))+
  geom_jitter(aes(color = Type))+
  geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
  stat_compare_means(method = "wilcox.test")+
  xlab(NULL) + scale_color_stata()+
  hori_x_theme + theme(legend.position = "none")

# Paired Wilcoxon test
avg.bNTI = avg.bNTI[-grep("12", avg.bNTI$Samples),]
stat = wilcox.test(bNTI~Type, data = avg.bNTI, paired = T)
