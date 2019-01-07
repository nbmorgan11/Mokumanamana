library(vegan)
library(cluster)
#library(mclust)
#library(clustsig)
library(ggplot2)
library(ggdendro)
library(ggpubr)
#library(tidyverse)  # data manipulation
#library(factoextra) # clustering visualization
#library(dendextend) # for comparing two dendrograms

setwd("~/Documents/NWHI_SeamountRecovery/RDocs")

abund_table<-read.csv("SPE_moku_stnd.csv",row.names=1,check.names=FALSE)
abund_table<-t(abund_table)

bray = vegdist(abund_table)
clus = agnes(bray)
dend = as.dendrogram(clus)
dend_data <- dendro_data(dend, type = "rectangle")
labs <- label(dend_data)

#export labs, reorder with grouping info to get correct order of transects from dendro
#otherwise depth bins don't match up when added to labs
#make sure to sort by x before closing
grouping_info_dend<-read.csv("grouping_info_dend2.csv",row.names=1,check.names=FALSE)
labs$depthbin=grouping_info_dend[,1]
labs$dive = grouping_info_dend[,3]

p <- ggplot(dend_data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data=label(dend_data), aes(label=labs$depthbin, x=x, y=0, colour=labs$dive),
            hjust = 1, angle = 90, size = 6)+ 
  ylim(-0.25, 1) +
  scale_color_manual(values= c("#00C791","#115089", "#59ECFF")) +
  ylab(label = "Dissimilarity") +
  xlab(label = " ") +
  theme_pubr(border = TRUE, legend = "right") +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(colour=guide_legend(title= "Dive",override.aes = list(shape = 19)))
print(p)


## Remove Outliers##

outlier= c("S360-04","S360-06","S362-23","S362-24")
abund_table = abund_table[!rownames(abund_table) %in% outlier, ]
#now re-run above

#### modded abundance dataset for urchin barrens
S362<-read.csv("S362_SPE_std.csv",row.names=1,check.names=FALSE)
S362<-t(S362)

bray2 = vegdist(S362)
clus2 = agnes(bray2)
dend2 = as.dendrogram(clus2)
dend_data2 <- dendro_data(dend2, type = "rectangle")
labs2 <- label(dend_data2)

#export labs, reorder with grouping info to get correct order of transects from dendro
#otherwise depth bins don't match up when added to labs
S362grouping_dend <- read.csv("S362_grouping_dend.csv",row.names=1,check.names=FALSE)
labs2$depthbin = S362grouping_dend[,1]

p <- ggplot(dend_data2$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data=label(dend_data2), aes(label=labs2$depthbin, x=x, y=0),
            hjust = 1, angle = 90, size = 8)+
  #scale_color_viridis(direction= -1,discrete = TRUE) +
  ylim(-0.25, 1) +
  ylab(label = "Dissimilarity") +
  xlab(label = " ") +
  theme_pubr(border = TRUE, legend = "right")
  #theme(legend.text = element_text(size=20),legend.title = element_text(size=20))
  #guides(colour=guide_legend(title= "Depth Bin",override.aes = list(shape=19)))
print(p)



###Removing outlier samples####
outlier= c("S360-04","S360-06","S362-23","S362-24")
abund_table2 = abund_table[!rownames(abund_table) %in% outlier, ]
grouping_info2 = grouping_info[!rownames(grouping_info) %in% outlier, ]

bray = vegdist(abund_table2)
clus = agnes(bray)
dend = as.dendrogram(clus)
dend_data <- dendro_data(dend, type = "rectangle")
labs <- label(dend_data)

#export labs, reorder with grouping info to get correct order of transects from dendro
#otherwise depth bins don't match up when added to labs
#make sure to sort by x before closing
grouping_info_dend2<-read.csv("grouping_info_dend2.csv",row.names=1,check.names=FALSE)
labs$depthbin=grouping_info_dend2[,1]
labs$dive = grouping_info_dend2[,3]

p <- ggplot(dend_data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data=label(dend_data), aes(label=label, x=x, y=0, colour=labs$dive),
            hjust = 1, angle = 90, size = 6)+ 
  ylim(-0.25, 1) +
  scale_color_manual(values= c("firebrick4","goldenrod1", "royalblue3")) +
  ylab(label = "Dissimilarity") +
  xlab(label = " ") +
  theme_pubr(border = TRUE, legend = "right") +
  theme(legend.text = element_text(size=15),legend.title = element_text(size=15)) +
  guides(colour=guide_legend(title= "Dive"))
print(p)

##########################################
cluster = Mclust(abund_table2)

plot.Mclust(cluster, what="BIC") 
summary.Mclust(cluster)
print.Mclust(cluster)

grouping_info2$mclust = cluster$classification

bray2 = vegdist(S362)
k1 = kmeans(bray2, centers = 3, nstart = 25)
print(k1)
S362grouping$kmeans = k1$cluster

bray = vegdist(abund_table2)
k2 = kmeans(bray, centers =8 , nstart = 25)
grouping_info2$kmeans = k2$cluster

write.csv(S362grouping,file = "S362_kmeanscluster.csv")
write.csv(grouping_info2,file = "Moku_kcluster.csv")
