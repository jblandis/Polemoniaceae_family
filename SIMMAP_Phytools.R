#stochastic mapping of flower color of the Polemoniaceae phylogeny
setwd("~/Stochastic_mapping")
library(phytools)
library(ape)
library(geiger)
library(nlme)

#trait data is in a csv file, with column one being species names, columns 2 the trait of interest.
mydata <- read.csv("Color.csv",row.names=1)
mytree <- read.nexus("MCC.fixed.nex")

comparison <- name.check(phy=mytree,data=mydata)

# prune taxa that don't have data but are present in the tree
mytree <- drop.tip(mytree,comparison$tree_not_data)

#create vector for traits
# if the vector names are changed, then commands following this will need to be altered to maintain the procedure
flower_color <- mydata[,1]
names(flower_color) <- row.names(mydata)

## simulate single stochastic character map using empirical Bayes method
#ER= equal rates
#ARD = all rates different
#SYM = symmetrical rates
mtree <- make.simmap(mytree, flower_color, model = "ARD")

## simulate 1000 stochastic character maps over the best tree
mtrees <- make.simmap(mytree, flower_color, model = "ARD", nsim = 1000)

cols <- setNames(c("snow2","blue", "yellow", "green"), sort(unique(flower_color)))
pdf("Summary_stochastic_flower_color.pdf",10,10)
XX <- describe.simmap(mtrees)
plot(XX, fsize=0.45, cex=c(0.3,0.1), type="fan", ftype="i",colors=cols)
tip_states <-sort(unique(getStates(mtrees[[1]], "tips")))
legend("bottomleft", legend=c("Non pigmented","Anthocyanin","Carotenoid","Chlorophyll"), fill=c("snow2","blue", "yellow", "green"))
h <- max(nodeHeights(mytree))
#max age of Pole phylogeny is 92.5 MYA
g <- signif(h, digits=2)
d <- (g-85)
e <- (g-70)
f <- (g-45)
m <- (g-30)
q <- (g-15)
s <- (g-0)
obj <- axis(1,pos=-0.015*g, at=seq(0,92.5,by=15), labels=c(g,q,m,f,e,d,""), hadj=0.6, mgp=c(3,0.25,1), cex.axis=0.75, lwd=0.85)
text(x=0.4*g, y=-0.14*g, "Time (Million Years Ago)", cex=0.75)
#for(j in 1:(length(obj))){
#  draw.arc(0,0,radius=obj[j], deg1=0,deg2=360,lwd=1, col=make.transparent("gray",0.35))
#}	
arc.cladelabels(text="Polemonium", node=830, mark.node=FALSE, ln.offset=1.43, lab.offset=1.45, cex=0.8)
arc.cladelabels(text="Linanthus", node=811, mark.node=FALSE, ln.offset=1.36, lab.offset=1.38, cex=0.8)
arc.cladelabels(text="Phlox", node=741, mark.node=FALSE, ln.offset=1.36, lab.offset=1.38, cex=0.8)
arc.cladelabels(text="Leptosiphon", node=711, mark.node=FALSE, ln.offset=1.43, lab.offset=1.45, cex=0.8)
arc.cladelabels(text="Gilia", node=640, mark.node=FALSE, ln.offset=1.36, lab.offset=1.38, cex=0.8)
arc.cladelabels(text="Navarretia", node=594, mark.node=FALSE, ln.offset=1.41, lab.offset=1.43, cex=0.8)
arc.cladelabels(text="Eriastrum", node=548, mark.node=FALSE, ln.offset=1.41, lab.offset=1.43, cex=0.8)
arc.cladelabels(text="Ipomopsis", node=501, mark.node=FALSE, ln.offset=1.41, lab.offset=1.43, cex=0.8)
arc.cladelabels(text="Loeselia", node=486, mark.node=FALSE, ln.offset=1.41, lab.offset=1.43, cex=0.8)
arc.cladelabels(text="Aliciella", node=461, mark.node=FALSE, ln.offset=1.3, lab.offset=1.32, cex=0.8)
arc.cladelabels(text="Cantua", node=451, mark.node=FALSE, ln.offset=1.3, lab.offset=1.32, cex=0.8)
arc.cladelabels(text="Cobaea", node=436, mark.node=FALSE, ln.offset=1.3, lab.offset=1.32, cex=0.8)
arc.cladelabels(text="Allophyllum", node=576, mark.node=FALSE, ln.offset=1.36, lab.offset=1.38, cex=0.8)
dev.off()


#do this over the distribution of trees
mydata <- read.csv("Color.csv",row.names=1)
mytrees <- read.nexus("1000_distribution.nex")

flower_color <- mydata[,1]
names(flower_color) <- row.names(mydata)


# simulate 100 stochastic character maps
mtrees <- make.simmap(mytrees, flower_color, model = "ARD", nsim = 100)

XX <- describe.simmap(mtrees, plot = FALSE)

XX
sink()
##########################################################################################################
#stochastic mapping of pollinators over the Polemoniaceae phylogey after removing taxa with missing data
setwd("~/Stochastic_mapping")

library(phytools)
library(ape)
library(geiger)
library(nlme)

#trait data is in a csv file, with column one being species names, columns 2 the trait of interest.
mydata <- read.csv("Autogamy_and_out.csv",row.names=1)
mytree <- read.nexus("MCC.fixed.nex")

#compares names between the tree and the data to list any discrepancies
comparison <- name.check(phy=mytree,data=mydata)

# prune taxa that don't have data but are present in the tree
mytree <- drop.tip(mytree,comparison$tree_not_data)

#double check to make sure that taxa all match with tree and data
name.check(phy=mytree,data=mydata)
comparison <- name.check(phy=mytree,data=mydata)
comparison

#create vector for traits
# if the vector names are changed, then commands following this will need to be altered to maintain the procedure
pollinator <- mydata[,1]
names(pollinator) <- row.names(mydata)

## simulate single stochastic character map using empirical Bayes method
mtree <- make.simmap(mytree, pollinator, model = "ARD")

## simulate 1000 stochastic character maps over the best tree
mtrees <- make.simmap(mytree, pollinator, model = "ARD", nsim = 1000)

#settings for plotting the whole tree in a circular phylogeny
cols <- setNames(c("snow2","black"), sort(unique(pollinator)))
pdf("Summary_stochastic_selfing_outcrossing.pdf",12,12)
XX <- describe.simmap(mtrees)
plot(XX, fsize=0.60, cex=c(0.4,0.2), type="fan", ftype="i",colors=cols)
tip_states <-sort(unique(getStates(mtrees[[1]], "tips")))
legend("bottomleft", legend=c("Autogamy", "Outcrossing"), fill=c("snow2","black"))
h <- max(nodeHeights(mytree))
#max age of Pole phylogeny is 92.5 MYA
g <- signif(h, digits=2)
d <- (g-85)
e <- (g-70)
f <- (g-45)
m <- (g-30)
q <- (g-15)
s <- (g-0)
obj <- axis(1,pos=-0.015*g, at=seq(0,92.5,by=15), labels=c(g,q,m,f,e,d,""), hadj=0.6, mgp=c(3,0.25,1), cex.axis=0.75, lwd=0.85)
text(x=0.4*g, y=-0.12*g, "Time (Million Years Ago)", cex=0.75)
#for(j in 1:(length(obj))){
#  draw.arc(0,0,radius=obj[j], deg1=0,deg2=360,lwd=1, col=make.transparent("gray",0.35))
#}	
arc.cladelabels(text="Polemonium", node=370, mark.node=FALSE, ln.offset=1.36, lab.offset=1.38, cex=0.8)
arc.cladelabels(text="Linanthus", node=367, mark.node=FALSE, ln.offset=1.36, lab.offset=1.38, cex=0.8)
arc.cladelabels(text="Phlox", node=342, mark.node=FALSE, ln.offset=1.44, lab.offset=1.46, cex=0.8)
arc.cladelabels(text="Leptosiphon", node=330, mark.node=FALSE, ln.offset=1.53, lab.offset=1.55, cex=0.8)
arc.cladelabels(text="Gilia", node=280, mark.node=FALSE, ln.offset=1.43, lab.offset=1.45, cex=0.8)
arc.cladelabels(text="Navarretia", node=271, mark.node=FALSE, ln.offset=1.41, lab.offset=1.43, cex=0.8)
arc.cladelabels(text="Eriastrum", node=245, mark.node=FALSE, ln.offset=1.41, lab.offset=1.43, cex=0.8)
arc.cladelabels(text="Ipomopsis", node=221, mark.node=FALSE, ln.offset=1.51, lab.offset=1.53, cex=0.8)
arc.cladelabels(text="Aliciella", node=212, mark.node=FALSE, ln.offset=1.3, lab.offset=1.32, cex=0.8)
arc.cladelabels(text="Cantua", node=208, mark.node=FALSE, ln.offset=1.3, lab.offset=1.32, cex=0.8)
arc.cladelabels(text="Cobaea", node=199, mark.node=FALSE, ln.offset=1.3, lab.offset=1.32, cex=0.8)
arc.cladelabels(text="Allophyllum", node=258, mark.node=FALSE, ln.offset=1.36, lab.offset=1.38, cex=0.8)
dev.off()



#do this over the distribution of trees
mydata <- read.csv("Autogamy_and_out.csv",row.names=1)
mytree <- read.nexus("MCC.fixed.nex")
mytrees <- read.nexus("1000_distribution.nex")

#compares names between the tree and the data to list any discrepancies
comparison <- name.check(phy=mytree,data=mydata)
to_prune <- comparison$tree_not_data
to_prune

# prune taxa that don't have data but are present in the tree iteratively through all trees
pruned.trees <- lapply(mytrees,drop.tip,to_prune)

#make sure the pruned trees are multiphylo object
class(pruned.trees)<-"multiPhylo"

#write the pruned trees to a new file, and open said file again
write.nexus(pruned.trees, file="Pruned_trees.nex")
pruned_trees <- read.nexus("Pruned_trees.nex")

#create vector for traits
# if the vector names are changed, then commands following this will need to be altered to maintain the procedure
pollinator <- mydata[,1]
names(pollinator) <- row.names(mydata)

# simulate 100 stochastic character maps
mtrees <- make.simmap(pruned_trees, pollinator, model="ARD", nsim = 100)


XX <- describe.simmap(mtrees, plot = FALSE)

XX

##################################################################
#modified Q matrix based on MuSSE runs
#trait data is in a csv file, with column one being species names, columns 2 the trait of interest.
mydata <- read.csv("Color.csv",row.names=1)
mytree <- read.nexus("MCC.fixed.nex")

comparison <- name.check(phy=mytree,data=mydata)
# prune taxa that don't have data but are present in the tree
mytree <- drop.tip(mytree,comparison$tree_not_data)

#create vector for traits
# if the vector names are changed, then commands following this will need to be altered to maintain the procedure
flower_color <- mydata[,1]
names(flower_color) <- row.names(mydata)

#create matrix based on transition rates from MuSSE (or HiSSE); order is column by column from table, with rows equaling 0
#last two entries specify the number of columns and rows in matrix
mod <- matrix(c(-0.102, 0.054, 0.073, 0.000, 0.046, -0.066, 0.019, 0.000, 0.056, 0.012, -0.092, 0.078, 0.000, 0.000, 0.000, -0.078),4,4)
#data is coded as 0-3, so rows and columns need to be labeled the same to distinguis transition rates, may changed based on data coding
rownames(mod) <- colnames(mod) <- c(0,1,2,3)
#make simmap function with user specified Q-matrix
mtree <- make.simmap(mytree, flower_color, Q=mod)


## simulate 1000 stochastic character maps over the best tree
mtrees <- make.simmap(mytree, flower_color, Q=mod, nsim = 1000)

cols <- setNames(c("snow2","blue", "yellow", "green"), sort(unique(flower_color)))
pdf("Summary_stochastic_flower_color2.pdf",10,10)
XX <- describe.simmap(mtrees)
plot(XX, fsize=0.45, cex=c(0.3,0.1), type="fan", ftype="i",colors=cols)
tip_states <-sort(unique(getStates(mtrees[[1]], "tips")))
legend("bottomleft", legend=c("Non pigmented","Anthocyanin","Carotenoid","Chlorophyll"), fill=c("snow2","blue", "yellow", "green"))
h <- max(nodeHeights(mytree))
#max age of Pole phylogeny is 92.5 MYA
g <- signif(h, digits=2)
d <- (g-85)
e <- (g-70)
f <- (g-45)
m <- (g-30)
q <- (g-15)
s <- (g-0)
obj <- axis(1,pos=-0.015*g, at=seq(0,92.5,by=15), labels=c(g,q,m,f,e,d,""), hadj=0.6, mgp=c(3,0.25,1), cex.axis=0.75, lwd=0.85)
text(x=0.4*g, y=-0.14*g, "Time (Million Years Ago)", cex=0.75)
#for(j in 1:(length(obj))){
#  draw.arc(0,0,radius=obj[j], deg1=0,deg2=360,lwd=1, col=make.transparent("gray",0.35))
#}	
arc.cladelabels(text="Polemonium", node=830, mark.node=FALSE, ln.offset=1.43, lab.offset=1.45, cex=0.8)
arc.cladelabels(text="Linanthus", node=811, mark.node=FALSE, ln.offset=1.36, lab.offset=1.38, cex=0.8)
arc.cladelabels(text="Phlox", node=741, mark.node=FALSE, ln.offset=1.36, lab.offset=1.38, cex=0.8)
arc.cladelabels(text="Leptosiphon", node=711, mark.node=FALSE, ln.offset=1.43, lab.offset=1.45, cex=0.8)
arc.cladelabels(text="Gilia", node=640, mark.node=FALSE, ln.offset=1.36, lab.offset=1.38, cex=0.8)
arc.cladelabels(text="Navarretia", node=594, mark.node=FALSE, ln.offset=1.41, lab.offset=1.43, cex=0.8)
arc.cladelabels(text="Eriastrum", node=548, mark.node=FALSE, ln.offset=1.41, lab.offset=1.43, cex=0.8)
arc.cladelabels(text="Ipomopsis", node=501, mark.node=FALSE, ln.offset=1.41, lab.offset=1.43, cex=0.8)
arc.cladelabels(text="Loeselia", node=486, mark.node=FALSE, ln.offset=1.41, lab.offset=1.43, cex=0.8)
arc.cladelabels(text="Aliciella", node=461, mark.node=FALSE, ln.offset=1.3, lab.offset=1.32, cex=0.8)
arc.cladelabels(text="Cantua", node=451, mark.node=FALSE, ln.offset=1.3, lab.offset=1.32, cex=0.8)
arc.cladelabels(text="Cobaea", node=436, mark.node=FALSE, ln.offset=1.3, lab.offset=1.32, cex=0.8)
arc.cladelabels(text="Allophyllum", node=576, mark.node=FALSE, ln.offset=1.36, lab.offset=1.38, cex=0.8)
dev.off()

#summarize trees with the MuSSE transition rates
XX <- describe.simmap(mtrees, plot = FALSE)

XX