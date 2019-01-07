# ============================================================
# Tutorial on drawing a heatmap using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================
setwd("~/Documents/NWHI_SeamountRecovery/RDocs")
library(reshape2)
library(ggplot2)
library(plyr)
library(scales)
library(ggpubr)

spe_moku<-read.csv("Moku_simper_heatmap2.csv",row.names=1,check.names=FALSE)
spe_moku<-t(spe_moku)
# Convert to relative frequencies
spe_moku2 <- spe_moku/rowSums(spe_moku)

#abund_order2 = abund_order2[,colSums(abund_order2) > 0.05] #removing low values

# S362_simper<-read.csv("NorthSimperdepth.csv",row.names=1,check.names=FALSE)
# S362_simper<-t(S362_simper)
# # Convert to relative frequencies
# S362_heat <- S362_simper/rowSums(S362_simper)

df<-melt(spe_moku2)
colnames(df)<-c("Samples","Species","Value")

#get depthbin as a factor to order y-axis by depth rather than location
groups = data.frame(do.call('rbind', strsplit(as.character(df$Samples),'m',fixed=TRUE)))
df$DepthBin = as.factor(groups$X1)
levorder= c("250m-North","250m-South","250m-West",
            "350m-North","350m-South","350m-West",
            "450m-North","450m-South","450m-West",
            "550m-North","550m-South","550m-West",
            "650m-North","650m-South","650m-West")
df$Samples2 <- factor(df$Samples, levels=levorder)

#df$DepthBin2 = factor(df$DepthBin, levels=c("250","350","450","550","650") )
#df$Samples3 <- factor(df$Samples, levels=df$Samples[order(df$DepthBin2)])
#df$Samples <- factor(df$Samples, levels=unique(df$DepthBin))
#df<-melt(S362_heat)
#colnames(df)<-c("Samples","Species","Value")


# We are going to apply transformation to our data to make it easier on eyes 
#df<-ddply(df,.(Samples),transform,rescale=scale(Value))
#df<-ddply(df,.(Samples),transform,rescale=sqrt(Value))

# Plot heatmap
m <- ggplot(df, aes(x=Species,y=Samples2)) + 
  geom_tile(aes(fill = Value) , colour = "gray") + 
  coord_equal() +
  scale_fill_gradient(low = "white",
                      high = "black",
                      guide = "colourbar")+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), limits = rev(levels(df$Samples2))) + 
  ylab(label = " ") +
  xlab(label = " ") +
  guides(fill=guide_legend(title="Relative\nAbundance")) +
  theme_pubr(margin = FALSE, legend = "right", base_size = 5) +
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust =.5,size=14, color="black"),
        axis.text.y = element_text(size=14, hjust=0, vjust = 0.5, color="black"))
print(m)
