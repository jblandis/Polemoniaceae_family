#Reconstruct corolla lenght on the tree
setwd("/continuous")
library(phytools)
library(phytools)
library(ape)
library(geiger)
library(nlme)

mytree <- read.nexus("MCC.fixed.nex")
mydata <- read.csv("Corolla_length.csv",row.names=1)

comparison <- name.check(phy=mytree,data=mydata)

# prune taxa that don't have data but are present in the tree
pruned.tree <- drop.tip(mytree,comparison$tree_not_data)

#create vector for traits
# if the vector names are changed, then commands following this will need to be altered to maintain the procedure
flower_length <- mydata[,1]
names(flower_length) <- row.names(mydata)

fit <- fastAnc(pruned.tree,flower_length, vars=TRUE, CI=TRUE)
fit

fit$CI[1,]

range(flower_length)

#plot the tree
pdf("Length_continuous_ASR.pdf")
obj <- contMap(pruned.tree,flower_length, plot=FALSE)
plot(obj, type="fan", fsize=0.25, lwd=1.5, legend=0.7*max(nodeHeights(pruned.tree)))
dev.off()

#project of tree into phenotype space
pdf("Length_Phenotype_space.pdf")
phenogram(pruned.tree,flower_length,fsize=0.25, spread.costs=c(1,0))
dev.off()

##########################################################################################
#Reconstruct corolla width on the tree
#Reconstruct corolla lenght on the tree
setwd("/continuous")
library(phytools)
library(phytools)
library(ape)
library(geiger)
library(nlme)

mytree <- read.nexus("MCC.fixed.nex")
mydata <- read.csv("Corolla_width.csv",row.names=1)

comparison <- name.check(phy=mytree,data=mydata)

# prune taxa that don't have data but are present in the tree
pruned.tree <- drop.tip(mytree,comparison$tree_not_data)
comparison <- name.check(phy=pruned.tree,data=mydata)

#create vector for traits
# if the vector names are changed, then commands following this will need to be altered to maintain the procedure
flower_width <- mydata[,1]
names(flower_width) <- row.names(mydata)

fit <- fastAnc(pruned.tree,flower_width, vars=TRUE, CI=TRUE)
fit

fit$CI[1,]

range(flower_width)

#plot the tree
pdf("Width_continuous_ASR.pdf")
obj <- contMap(pruned.tree,flower_width, plot=FALSE)
plot(obj, type="fan", fsize=0.25, lwd=1.5, legend=0.7*max(nodeHeights(pruned.tree)))
dev.off()

#project of tree into phenotype space
pdf("Width_Phenotype_space.pdf")
phenogram(pruned.tree,flower_width,fsize=0.25, spread.costs=c(1,0))
dev.off()
