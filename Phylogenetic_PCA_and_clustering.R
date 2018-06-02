#stochastic mapping of flower color of the Polemoniaceae phylogeny
library(phytools)
library(geiger)
library(ggfortify)
library(ggplot2)

setwd("/Phylogenet_PCA")

#trait data is in a csv file, with column one being species names, columns 2 and 3 the continuous traits of interest.
#Species missing data for either or both traits was removed prior to analysis
#Use the next two lines for flower color and Linanthus nuclear
mydata <- read.csv("Trait_data_working.csv",row.names=1)
mytree <- read.nexus("MCC.fixed.nex")

Pollinator <- mydata[,4]
names(Pollinator) <- row.names(mydata)

Flower <- mydata[,1]
names(Flower) <- row.names(mydata)

#keep only taxa with data
comparison <- name.check(phy=mytree,data=mydata)

# prune taxa that don't have data but are present in the tree
mytree <- drop.tip(mytree,comparison$tree_not_data)
comparison <- name.check(phy=mytree,data=mydata)
comparison

#perform phylogenetic PCA
pca <- phyl.pca(mytree, mydata[,2:3], method="lambda", mode="cov")
summary(pca)

#pull out cordinates
tbl<-pca$S

#scatter plot with PC1 and PC2
ggplot(tbl)+ geom_point(aes(x=PC1, y=PC2))

ggplot(tbl,aes(x=PC1, y=PC2),colour=Pollinator) + geom_point() + scale_color_brewer(levels(Pollinator))

#standard PCA of length and width
pca2 <- prcomp(mydata[,2:3])
#plot of PC1 and PC2
pdf("PC1_vs_PC2.pdf")
autoplot(pca2, data=mydata, colour='Pollinator', size=3, alpha=0.9, variance_percentage=FALSE, loadings=FALSE, scale=0) + theme_classic() + scale_color_brewer(palette="Paired")
dev.off()

## Plot PC1 vs. PC2 with ellipses around functional groups
pdf("PC1_vs_PC2_ellipses.pdf")
autoplot(pca2, data=mydata, colour='Pollinator', size=3, alpha=0.9, variance_percentage=FALSE, loadings=FALSE, scale=0, frame=TRUE, frame.type='norm') + theme_classic() + scale_color_brewer(palette="Paired")
dev.off()


## Plot PC1 vs. PC2 with ellipses around functional groups
pdf("PC1_vs_PC2_guild.pdf")
autoplot(pca2, data=mydata, colour='Functional', size=3, alpha=0.9, variance_percentage=FALSE, loadings=FALSE, scale=0, frame=TRUE, frame.type='norm') + theme_classic() + scale_color_brewer(palette="Paired")
dev.off()


pdf("biplots.pdf", 30,30)
biplot(pca, cex=1, var.axes=TRUE, arrow.len=0.5,xlim=c(-1.3,2.25))
dev.off()


pdf("biplots2.pdf", 30,30)
biplot(results, choices=c(2,3), cex=1, var.axes=TRUE, arrow.len=0.5,ylim=c(-0.5,0.5))
dev.off()

#########################################
#K-means clustering

library(plyr)

mydata <- read.csv("Trait_data_working.csv",row.names=1)

pollinator <- as.factor(setNames(mydata[,4], rownames(mydata)))
guilds <- as.factor(setNames(mydata[,5], rownames(mydata)))
clust.res <- hclust(dist(mydata[,2:3]), method="average")

pdf("Clustering_poll_label.pdf", 10,10)
plot(clust.res, hang=-1, labels=pollinator, cex=0.4)
dev.off()

pdf("Clustering_guild.pdf", 10,10)
plot(clust.res, hang=-1, labels=guilds, cex=0.4)
dev.off()
