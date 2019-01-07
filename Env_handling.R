library(corrplot)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(scales)
library(grid)
library(plyr)
library(gridExtra)
library(viridis)
library(ggpubr)

setwd("~/Documents/NWHI_SeamountRecovery/RDocs")

##NorthSite##
S362ENV = read.csv("S362_ENV.csv",row.names=1,check.names=FALSE)
S362Scale = as.data.frame(scale(S362ENV))
S362Cor = cor(S362Scale, use="complete.obs")
#corrplot(S362Cor, type = "lower", method = "number",number.cex=0.7, bg="gray48", tl.col = "black")
corrplot(S362Cor, type = "lower", method = "color", bg="gray48", tl.col = "black", cl.pos = "r",tl.cex = 0.7)

col1= colorRampPalette(c("#440154FF","#39568CFF","#1F968BFF","#FDE725FF"))

S362env = read.csv("S362_ENV_short.csv",row.names=1,check.names=FALSE)
S362env=S362env[,-(5)] #removing backscatter b/c PCA just adds in average value across samples
S362Scale2 = as.data.frame(scale(S362env))
S362Cor2 = cor(S362Scale2, use="complete.obs")
corrplot(S362Cor2, type = "lower", method = "number",number.cex=0.9,col = col1(100), tl.col = "black",tl.cex = 0.9,cl.pos = "r")


### All Sites ###
meta_table = read.csv("ENV_moku.csv",row.names=1,check.names=FALSE)
envScale = as.data.frame(scale(meta_table))
envCor = cor(envScale, use="complete.obs")
corrplot(envCor, type = "lower", method = "color", tl.col = "black", cl.pos = "r",tl.cex = 0.7)

env = read.csv("ENV_moku_short.csv",row.names=1,check.names=FALSE)
env=env[,-(5)] #removing backscatter b/c PCA just adds in average value across samples
envScale2 = as.data.frame(scale(env))
envCor2 = cor(envScale2, use="complete.obs")
corrplot(envCor2, type = "lower", method = "number",number.cex=0.9,col = col1(100), tl.col = "black",tl.cex = 0.9,cl.pos = "r")

grouping_info<-read.csv("grouping_info.csv",row.names=1,check.names=FALSE)
S362grouping<-read.csv("S362_grouping.csv",row.names=1,check.names=FALSE)


### PCA ###

env.pca = PCA(S362Scale2)
# S362_grouping<-read.csv("S362_grouping.csv",row.names=1,check.names=FALSE)
# fviz_pca_biplot(env.pca, repel = TRUE, habillage = S362grouping$DepthBin, geom.ind = "point", pointsize = 7,
#                 col.var = "black", alpha.var = "contrib", title = "") +
#   guides(color = guide_legend(override.aes = list(size = 9))) +
#   scale_shape_manual(values=shape_values2) +
#   scale_color_viridis(direction= -1,discrete = TRUE) +
#   #geom_text(aes(label=S362_grouping$DepthBin, size=12),hjust=0,vjust=2) +
#   xlab(label= "PC1 (51.3%)") +
#   ylab(label= "PC2 (14.1%)") +
#   theme_pubr(border=TRUE, legend = "right") +
#   theme(legend.text = element_text(size=10),legend.title = element_text(size=12),axis.title = element_text(size=10))
# 
# env.pca$var
# env.pca$ind$coord
 write.csv(env.pca$var, file="North_PCA_values.csv")

PC1 <- env.pca$ind$coord[,1]
PC2 <- env.pca$ind$coord[,2]
labs <- rownames(env.pca$ind$coord)
PCs <- data.frame(cbind(PC1,PC2))
rownames(PCs) <- labs
vPC1 <- env.pca$var$coord[,1]
vPC2 <- env.pca$var$coord[,2]
vlabs <- rownames(env.pca$var$coord)
vPCs <- data.frame(cbind(vPC1,vPC2))
rownames(vPCs) <- vlabs
colnames(vPCs) <- c("x","y")
PCs$depthbin = S362grouping[,1]

angle <- seq(-pi, pi, length = 50) 
df <- data.frame(x = sin(angle), y = cos(angle)) 

p <- ggplot() +  
  geom_point(data=PCs, aes(x=PC1,y=PC2, color=depthbin), size=8) +
  geom_text(data=vPCs, aes(x=vPC1*3.2,y=vPC2*3.2,label=rownames(vPCs)), size=4) + 
  geom_path(aes(x*3, y*3), data = df, colour="grey70") +
  geom_segment(data=vPCs, aes(x = 0, y = 0, xend = x*3, yend = y*3), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  scale_color_viridis(direction = -1,discrete = TRUE) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_pubr(border = TRUE, legend = "right") +
  ylab(label = "PC2 14.29%") +
  xlab( label = "PC1 54.61%") +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=12))
print(p)

#### MOKU OVERALL PCA ###

env.pca2 = PCA(envScale2)
# grouping_info<-read.csv("grouping_info.csv",row.names=1,check.names=FALSE)
# envshape = grouping_info$Dive
# 
# fviz_pca_biplot(env.pca2, repel = TRUE, habillage = grouping_info$Dive, geom.ind = "point", pointsize = 8,
#                 col.var = "black", title = "") +
#   scale_color_manual(values= c("firebrick4","goldenrod1", "royalblue3")) +
#   guides(color = guide_legend(override.aes = list(size = 9))) +
#   geom_text(aes(label=grouping_info$DepthBin, color=grouping_info$Dive),hjust=0,vjust=2) +
#   #xlab(label= "PC1 (36.9%)") +
#   #ylab(label= "PC2 (15.4%)") +
#   theme_pubr(border=TRUE, legend = "right")
#   #theme(legend.text = element_text(size=15),legend.title = element_text(size=18),axis.title = element_text(size=15))
# 
write.csv(env.pca2$var, file="Moku_PCA_values.csv")

PC1 <- env.pca2$ind$coord[,1]
PC2 <- env.pca2$ind$coord[,2]
labs <- rownames(env.pca2$ind$coord)
PCs <- data.frame(cbind(PC1,PC2))
rownames(PCs) <- labs
vPC1 <- env.pca2$var$coord[,1]
vPC2 <- env.pca2$var$coord[,2]
vlabs <- rownames(env.pca2$var$coord)
vPCs <- data.frame(cbind(vPC1,vPC2))
rownames(vPCs) <- vlabs
colnames(vPCs) <- c("x","y")
PCs$DepthBin = grouping_info[,1]
PCs$Dive = grouping_info[,3]

angle <- seq(-pi, pi, length = 50) 
df <- data.frame(x = sin(angle), y = cos(angle)) 

p <- ggplot() +  
  geom_point(data=PCs, aes(x=PC1,y=PC2, color=DepthBin, shape=Dive), size=11) +
  geom_text(data=vPCs, aes(x=vPC1*3.5,y=vPC2*3.5,label=rownames(vPCs)), size=8) + 
  geom_path(aes(x*3, y*3), data = df, colour="grey70") +
  geom_segment(data=vPCs, aes(x = 0, y = 0, xend = x*3, yend = y*3), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  scale_shape_manual(values=c(15,17,19))+
  scale_color_viridis(direction = -1,discrete = TRUE) +
  guides(color = guide_legend(override.aes = list(size = 6)),
         shape = guide_legend(override.aes = list(size = 6))) +
  theme_pubr(border = TRUE, legend = "right") +
  ylab(label = "PC2 14.72%") +
  xlab( label = "PC1 38.32%") +
  theme(legend.text = element_text(size=22),
        legend.title = element_text(size=22),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))
print(p)
