setwd("~/Documents/NWHI_SeamountRecovery/RDocs")
library(vegan)
library(betapart)
library(ggplot2)
library(viridis)
library(ggpubr)

### NORTH ONLY ###

#community data
S362 <-read.csv("S362_SPE_std.csv",row.names=1,check.names=FALSE)
S362 = as.data.frame(t(S362)) #Transpose the data to have sample names on rows
#grouping information
S362grouping<-read.csv("S362_grouping.csv",row.names=1,check.names=FALSE)

#To look at beta diversity along depth gradient -> does not seem to work for sampling/density plot
S362$DepthBin = S362grouping[,1]
depth.s362 = do.call(data.frame, aggregate(.~ DepthBin, S362, function(x) sum(x)))
row.names(depth.s362) = depth.s362$DepthBin
depth.s362 = depth.s362[,-(1)]


#all samples
# multi.S362 = beta.multi.abund(S362)
# S362.samp <- beta.sample.abund(S362,
#                                 sites=20, samples=100)
# d.362 <- S362.samp$sampled.values
# df.S362 = stack(d.362)
# is.factor(df.S362[,2])
# colnames(df.S362) = c("Value","Diversity")
# 
# p = ggplot(data=df.S362, aes(x = Value)) +
#   geom_density(aes(color = Diversity)) +
#   theme_pubr(border = TRUE, legend = "right") +
#   ylab(label = "Density") +
#   xlab( label = "Beta Diversity") 
#         
# print(p)


#one site per depth bin
s362.multi = beta.multi.abund(depth.s362)
ds362.samp <- beta.sample.abund(depth.s362,
                            sites=7, samples=100) #make sure using fewer sites than full number
dist.362 <- ds362.samp$sampled.values

df.ds362 = stack(dist.362)
is.factor(df.ds362[,2])
colnames(df.ds362) = c("Value","Diversity")

# Color plot for standalone graph
p = ggplot(data=df.ds362, aes(x = Value, fill=Diversity)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(name = expression(paste(beta, "-North")),
                    values = c("#00523B","#00916A","#46EBBE"),
                    labels = c(expression(paste(beta, "-Overall")),expression(paste(beta, "-Substitution")),expression(paste(beta, "-Nested"))))+
  scale_x_continuous(breaks = c(0.25,0.50, 0.75)) +
  theme_pubr(border = TRUE, legend = "right") +
  theme(legend.text = element_text(size=22, hjust=0),
        legend.title = element_text(size=22),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22)) +
  ylab(label = "Density") +
  xlab( label = expression(paste(beta, "-Diversity"))) 
print(p)

#black and gray for publication combination graph
p = ggplot(data=df.ds362, aes(x = Value, fill=Diversity)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("gray15","gray48","gray83"),
                    labels = c(expression(beta[Overall]),expression(beta[Substitution]),expression(beta[Nested])))+
  scale_x_continuous(breaks = c(0.25,0.50, 0.75)) +
  theme_pubr(border = TRUE, legend = "right") +
  theme(legend.text = element_text(size=15),
        legend.title = element_blank()) +
  ylab(label = "Density - North") +
  xlab( label = " ") 
print(p)

ggsave(filename = "NorthDepthBeta.jpg", plot = p)

#plotting from Baselga 2017
# plot(density(dist.362$beta.BRAY),
#      xlim=c(0,2), ylim=c(0, 15), xlab='Beta
#      diversity', main='', lwd=3)
# lines(density(dist.362$beta.BRAY.BAL), lty=1, lwd=2)
# lines(density(dist.362$beta.BRAY.GRA), lty=2, lwd=2)


### All Moku ###

abund_table<-read.csv("SPE_moku_stnd.csv",row.names=1,check.names=FALSE)
abund_table<-as.data.frame(t(abund_table)) #Transpose the data to have sample names on rows

grouping_info<-read.csv("grouping_info.csv",row.names=1,check.names=FALSE)

abund_table$Area = grouping_info[,5]
depth.moku = do.call(data.frame, aggregate(.~ Area, abund_table, function(x) sum(x)))
row.names(depth.moku) = depth.moku$Area
depth.moku = depth.moku[,-(1)]

depth.west = depth.moku[(11:15),]
depth.south = depth.moku[(6:10),]

west.multi = beta.multi.abund(depth.west)
south.multi = beta.multi.abund(depth.south)
west.samp <- beta.sample.abund(depth.west,
                                sites=3, samples=100) #make sure using fewer sites than full number
south.samp <- beta.sample.abund(depth.south,
                                sites=3, samples=100) #make sure using fewer sites than full number
dist.west <- west.samp$sampled.values
dist.south <- south.samp$sampled.values

dfw = stack(dist.west)
is.factor(dfw[,2])
colnames(dfw) = c("Value","Diversity")

dfs = stack(dist.south)
is.factor(dfs[,2])
colnames(dfs) = c("Value","Diversity")

pw = ggplot(data=dfw, aes(x = Value, fill=Diversity)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(name = expression(paste(beta, "-West")),
                    values = c("#004A4C","#00C1C9","#59ECFF"),
                    labels = c(expression(paste(beta, "-Overall")),expression(paste(beta, "-Substitution")),expression(paste(beta, "-Nested"))))+
  scale_x_continuous(breaks = c(0.25,0.50, 0.75)) +
  theme_pubr(border = TRUE, legend = "right") +
  theme(legend.text = element_text(size=22, hjust=0),
        legend.title = element_text(size=22),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22)) +
  ylab(label = "Density") +
  xlab( label = expression(paste(beta, "-Diversity"))) 
print(pw)

ps = ggplot(data=dfs, aes(x = Value, fill=Diversity)) +
  geom_density(alpha = 0.8) +
  scale_fill_manual(name = expression(paste(beta, "-South")),
                    values = c("#00064C","blue","dodgerblue"),
                    labels = c(expression(paste(beta, "-Overall")),expression(paste(beta, "-Substitution")),expression(paste(beta, "-Nested"))))+
  scale_x_continuous(breaks = c(0.25,0.50, 0.75)) +
  theme_pubr(border = TRUE, legend = "right") +
  theme(legend.text = element_text(size=22, hjust=0),
        legend.title = element_text(size=22),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22)) +
  ylab(label = "Density") +
  xlab(label = expression(paste(beta, "-Diversity"))) 
print(ps)

pd=ggarrange(p, ps, pw, ncol = 2, nrow = 2)
print(pd)


# plotting the distributions of components (Baselga 2017)
# dist.s <- ceram.s.samp$sampled.values
# dist.n <- ceram.n.samp$sampled.values
# plot(density(dist.s$beta.SOR),
#      xlim=c(0,0.8), ylim=c(0, 19), xlab='Beta
#      diversity', main='', lwd=3)
# lines(density(dist.s$beta.SNE), lty=1, lwd=2)
# lines(density(dist.s$beta.SIM), lty=2, lwd=2)
# lines(density(dist.n$beta.SOR), col='grey60',
#       lwd=3)
# lines(density(dist.n$beta.SNE), col='grey60',
#       lty=1, lwd=2)
# lines(density(dist.n$beta.SIM), col='grey60',
#       lty=2, lwd=2)





