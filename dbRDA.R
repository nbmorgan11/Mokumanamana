library(vegan)
library(ggplot2)
library(viridis)
library(ggpubr)
setwd("~/Documents/NWHI_SeamountRecovery/RDocs")

#RDA analysis with metric scaling, ignores neg. eigenvalues
#vare.cap <- capscale(S362 ~ Oxygen + MeanPOC + Rugosity, S362env,dist="bray")

### NORTH ONLY ###
S362 <-read.csv("S362_SPE_std.csv",row.names=1,check.names=FALSE)
S362 = t(S362) #Transpose the data to have sample names on rows
S362env = read.csv("S362_ENV_short.csv",row.names=1,check.names=FALSE)
S362Scale2 = as.data.frame(scale(S362env))
S362grouping<-read.csv("S362_grouping.csv",row.names=1,check.names=FALSE)

#dbRBDA - using environmental variables selectd for by distLM
#bioenv in PRIMER and R give same results with fewer variables.
S362db = dbrda(S362 ~ Oxygen + MeanPOC + Rugosity, S362Scale2,dist="bray")

summary(S362db)
anova(S362db) ## overall test of the significance of the analysis
anova(S362db, by="axis", perm.max=500) ## test axes for significance
anova(S362db, by="terms", permu=200) ## test for sig. environ. variables

scrs<-scores(S362db,display=c("sp","wa","lc","bp","cn"))
df_sites<-data.frame(scrs$sites)
df_sites = cbind(df_sites,S362grouping)
#colnames(df_sites)<-c("x","y","Dive","DepthBin")

# Reference: http://www.inside-r.org/packages/cran/vegan/docs/envfit
# The printed output of continuous variables (vectors) gives the direction cosines 
# which are the coordinates of the heads of unit length vectors. In plot these are 
# scaled by their correlation (square root of the column r2) so that "weak" predictors 
# have shorter arrows than "strong" predictors. You can see the scaled relative lengths 
# using command scores. The plotted (and scaled) arrows are further adjusted to the 
# current graph using a constant multiplier: this will keep the relative r2-scaled 
# lengths of the arrows but tries to fill the current plot. You can see the multiplier 
# using vegan:::ordiArrowMul(result_of_envfit), and set it with the argument arrow.mul. 

#Draw biplots
#multiplier <- vegan:::ordiArrowMul(scrs$biplot) #because my species values are so small this multiplier makes the arrows unreadable
#df_arrows<- scrs$biplot*multiplier
df_arrows<- scrs$biplot
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

#plotting

p<-ggplot()+
  geom_point(data=df_sites,aes(dbRDA1,dbRDA2,color=DepthBin), size=8) + 
  scale_color_viridis(direction = -1, discrete = TRUE)+
  geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(df_arrows*1.1),
            aes(x, y, label = rownames(df_arrows)),color="black", size = 5)+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ylab(label = "dbRDA2 39.83%") +
  xlab( label = "dbRDA1 49.12%") +
  theme_pubr(border = TRUE, legend = "right") +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))
print(p)


######ALL Sites#####
abund_table<-read.csv("SPE_moku_stnd.csv",row.names=1,check.names=FALSE)
abund_table<-t(abund_table)
grouping_info<-read.csv("grouping_info.csv",row.names=1,check.names=FALSE)
env = read.csv("ENV_moku_short.csv",row.names=1,check.names=FALSE)
env = env[,-(5)]

outlier= c("S360-04","S360-06","S362-23","S362-24")
abund_table2 = abund_table[!rownames(abund_table) %in% outlier, ]
grouping_info2 = grouping_info[!rownames(grouping_info) %in% outlier, ]
env2 = env[!rownames(env) %in% outlier, ]
envScale2 = as.data.frame(scale(env2))

mokudb2 = dbrda(abund_table2 ~  MeanDirection + Oxygen + Length + percWS +
        StdDevSubstrateSize + MeanPOC + Slope + StdDevSlope + Rugosity + 
        StdDevRugosity + AvgCurrentEW +  AvgCurrentNS, 
        envScale2,dist="bray")

summary(mokudb2)

anova(mokudb2) ## overall test of the significance of the analysis
anova(mokudb2, by="axis", perm.max=500) ## test axes for significance
anova(mokudb2, by="terms", permu=200) ## test for sig. environ. variables


scrs<-scores(mokudb2,display=c("sp","wa","lc","bp","cn"))
df_sites<-data.frame(scrs$sites)
df_sites = cbind(df_sites,grouping_info2)
#colnames(df_sites)<-c("x","y","Dive","DepthBin")


# Reference: http://www.inside-r.org/packages/cran/vegan/docs/envfit
# The printed output of continuous variables (vectors) gives the direction cosines 
# which are the coordinates of the heads of unit length vectors. In plot these are 
# scaled by their correlation (square root of the column r2) so that "weak" predictors 
# have shorter arrows than "strong" predictors. You can see the scaled relative lengths 
# using command scores. The plotted (and scaled) arrows are further adjusted to the 
# current graph using a constant multiplier: this will keep the relative r2-scaled 
# lengths of the arrows but tries to fill the current plot. You can see the multiplier 
# using vegan:::ordiArrowMul(result_of_envfit), and set it with the argument arrow.mul. 

#Draw biplots
#multiplier <- vegan:::ordiArrowMul(scrs$biplot)
#df_arrows<- scrs$biplot*multiplier
df_arrows<- scrs$biplot
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

#plotting

p<-ggplot()+
  geom_point(data=df_sites,aes(dbRDA1,dbRDA2,colour=DepthBin, shape=Dive), size=8) + 
  scale_shape_manual(values=c(15,17,19))+ 
  scale_color_viridis(direction = -1, discrete = TRUE)+
  geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(df_arrows*1.1),
            aes(x, y, label = rownames(df_arrows)),color="black", size=5)+
  guides(color = guide_legend(override.aes = list(size = 6)),
         shape = guide_legend(override.aes = list(size = 6))) +
  ylab(label = "dbRDA2 17.13%") +
  xlab( label = "dbRDA1 27.55%") +
  theme_pubr(border = TRUE, legend = "right") +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))
print(p)

df_sites$sand = env2[,5]
df_sites$rugosity = env2[,9]
df_sites$oxygen = env2[,1]
df_sites$poc = env2[,11]
df_sites$currentN = env2[,12]
df_sites$currentE = env2[,13] #less important
df_sites$mdir = env2[,3]
df_sites$len = env2[,4]
df_sites$sdSub = env2[,6]
df_sites$sdSlp = env2[,8]


ps<-ggplot(data=df_sites, aes(dbRDA1,dbRDA2, fill=DepthBin, size=sand)) +
  geom_point(shape=21)+
  #scale_size_area(max_size = 20) +
  scale_size(range = c(3,15)) +
  scale_fill_viridis(direction = -1,discrete = TRUE) +
  guides(fill = guide_legend(override.aes = list(size = 4)),
         size=guide_legend(title= "Perc. Sand", override.aes = list(shape=19))) +
  #geom_text(aes(label=Dive),size=3, color="black",hjust=0.5,vjust=2.5) +
  theme_pubr(border = TRUE, legend = "right") +
  ylab(label = "dbRDA2") +
  xlab( label = "dbRDA1") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank()) # remove axis ticks
print(ps)

pr<-ggplot(data=df_sites,aes(dbRDA1,dbRDA2, fill=DepthBin, size=rugosity)) +
  geom_point(shape=21)+
  #scale_size_area(max_size = 20) +
  scale_size(range = c(3,15)) +
  scale_fill_viridis(direction = -1,discrete = TRUE) +
  guides(fill = guide_legend(override.aes = list(size = 4), order = 2),
         size=guide_legend(title= "Rugosity", override.aes = list(shape=19), order=1)) +
  #geom_text(aes(label=Dive),size=3, color="black",hjust=0.5,vjust=2.5) +
  theme_pubr(border = TRUE, legend = "right") +
  ylab(label = "dbRDA2") +
  xlab( label = "dbRDA1") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank()) # remove axis ticks
print(pr)

po<-ggplot(data=df_sites,aes(dbRDA1,dbRDA2, fill=DepthBin, size=oxygen)) +
  geom_point(shape=21)+
  #scale_size_area(max_size = 20) +
  scale_size(range = c(3,15)) +
  scale_fill_viridis(direction = -1,discrete = TRUE) +
  guides(fill = guide_legend(override.aes = list(size = 4), order = 2),
         size=guide_legend(title= "Oxygen", override.aes = list(shape=19), order = 1)) +
  #geom_text(aes(label=Dive),size=3, color="black",hjust=0.5,vjust=2.5) +
  theme_pubr(border = TRUE, legend = "right") +
  ylab(label = "dbRDA2") +
  xlab( label = "dbRDA1") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank()) # remove axis ticks

pmc<-ggplot(data=df_sites,aes(dbRDA1,dbRDA2, fill=DepthBin, size=poc)) +
  geom_point(shape=21)+
  #scale_size_area(max_size = 20) +
  scale_size(range = c(3,15)) +
  scale_fill_viridis(direction = -1,discrete = TRUE) +
  guides(fill = guide_legend(override.aes = list(size = 4), order = 2),
         size=guide_legend(title= "POC", override.aes = list(shape=19), order = 1)) +
  theme_pubr(border = TRUE, legend = "right") +
  ylab(label = "dbRDA2") +
  xlab( label = "dbRDA1") +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank()) # remove axis ticks
print(pmc)

pdepth=ggarrange(ps, pr, po, pmc, ncol = 2, nrow = 2, labels = c("a","b","c","d"))
print(pdepth)



# psl<-ggplot(data=df_sites,aes(dbRDA1,dbRDA2, fill=DepthBin, size=sdSlp)) +
#   geom_point(shape=21)+
#   scale_size(range = c(3,15)) +
#   scale_fill_viridis(direction = -1,discrete = TRUE) +
#   guides(fill = guide_legend(override.aes = list(size = 4)),
#          size=guide_legend(title= "Slope SD", override.aes = list(shape=19))) +
#   #geom_text(aes(label=Dive),size=3, color="black",hjust=0.5,vjust=2.5) +
#   theme_pubr(border = TRUE, legend = "right") +
#   ylab(label = "dbRDA2") +
#   xlab( label = "dbRDA1") +
#   theme(legend.text = element_text(size=10),
#         legend.title = element_text(size=10),
#         axis.text.x = element_blank(),  # remove x-axis text
#         axis.text.y = element_blank(), # remove y-axis text
#         axis.ticks = element_blank()) # remove axis ticks


# psz<-ggplot(data=df_sites,aes(dbRDA1,dbRDA2, fill=DepthBin, size=sdSub)) +
#   geom_point(shape=21)+
#   scale_size(range = c(3,15)) +
#   scale_fill_viridis(direction = -1,discrete = TRUE) +
#   guides(fill = guide_legend(override.aes = list(size = 4)),
#          size=guide_legend(title= "Substrate SD", override.aes = list(shape=19))) +
#   #geom_text(aes(label=Dive),size=3, color="black",hjust=0.5,vjust=2.5) +
#   theme_pubr(border = TRUE, legend = "right") +
#   ylab(label = "dbRDA2") +
#   xlab( label = "dbRDA1") +
#   theme(legend.text = element_text(size=10),
#         legend.title = element_text(size=10),
#         axis.text.x = element_blank(),  # remove x-axis text
#         axis.text.y = element_blank(), # remove y-axis text
#         axis.ticks = element_blank()) # remove axis ticks

# pdir<-ggplot(data=df_sites,aes(dbRDA1,dbRDA2, fill=Dive, size=mdir)) +
#   geom_point(shape=21)+
#   scale_size(range = c(3,15)) +
#   scale_fill_viridis(direction = -1,discrete = TRUE) +
#   guides(fill = guide_legend(override.aes = list(size = 4)),
#          size=guide_legend(title= "Mean Dir", override.aes = list(shape=19))) +
#   #geom_text(aes(label=Dive),size=3, color="black",hjust=0.5,vjust=2.5) +
#   theme_pubr(border = TRUE, legend = "right") +
#   ylab(label = "dbRDA2") +
#   xlab( label = "dbRDA1") +
#   theme(legend.text = element_text(size=10),
#         legend.title = element_text(size=10),
#         axis.text.x = element_blank(),  # remove x-axis text
#         axis.text.y = element_blank(), # remove y-axis text
#         axis.ticks = element_blank()) # remove axis ticks

## BY DIVE ##
pcn<-ggplot(data=df_sites,aes(dbRDA1,dbRDA2, fill=Dive, size=currentN)) +
  geom_point(shape=21)+
  #scale_size_area(max_size = 20) +
  scale_size(range = c(3,15)) +
  scale_fill_manual(values= c("#00C791","#115089", "#59ECFF")) +
  guides(fill = guide_legend(override.aes = list(size = 4)),
         size=guide_legend(title= "N-S Current", override.aes = list(shape=19))) +
  #geom_text(aes(label=DepthBin),size=3, color="black",hjust=0.5,vjust=2.5) +
  theme_pubr(border = TRUE, legend = "right") +
  ylab(label = "dbRDA2") +
  xlab( label = "dbRDA1 ") +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank()) # remove axis ticks

pce<-ggplot(data=df_sites,aes(dbRDA1,dbRDA2, fill=Dive, size=currentE)) +
  geom_point(shape=21)+
  #scale_size_area(max_size = 20) +
  scale_size(range = c(3,15)) +
  scale_fill_manual(values= c("#00C791","#115089", "#59ECFF")) +
  guides(fill = guide_legend(override.aes = list(size = 4)),
         size=guide_legend(title= "E-W Current", override.aes = list(shape=19))) +
  #geom_text(aes(label=DepthBin),size=3, color="black",hjust=0.5,vjust=2.5) +
  theme_pubr(border = TRUE, legend = "right") +
  ylab(label = "dbRDA2") +
  xlab( label = "dbRDA1") +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank()) # remove axis ticks


pdive=ggarrange(pcn,pce, ncol = 1, nrow = 2, labels = c("a","b"))
print(pdive)
