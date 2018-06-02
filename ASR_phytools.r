library(phytools)
library(geiger)
library(nlme)

#starting files of the tree and the data file
mydata <- read.csv("Color_discrete.csv",row.names=1)
mytree <- read.nexus("MCC.fixed.nex")


#create vector for traits
# if the vector names are changed, then commands following this will need to be altered to maintain the procedure
states <- mydata[,1]
names(states) <- row.names(mydata)

#marginal TRUE
fit.ace.true <- ace(states, mytree, type="discrete", model="ARD", marginal=TRUE)
pdf("ASR_marginal_true.pdf")
cols <- setNames(c("blue","green","snow2","yellow"), sort(unique(states)))
h<-max(nodeHeights(mytree))
offset.factor<-1.01 ## increase this for greater offset
plotTree(rescale(mytree,model="depth",depth=offset.factor*h),
    color="transparent",ftype="i",type="fan",fsize=0.25,lwd=2)
par(fg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
plotTree(mytree, colors=cols, type="fan",fsize=0.25,ftype="i",add=TRUE,
	lwd=2,cex=0.3, label.offset=1)
nodelabels(pie=fit.ace.true$lik.anc,cex=0.3, piecol=c("blue","green","snow2","yellow"))
dev.off()

#marginal FALSE
fit.ace.false <- ace(states, mytree, type="discrete", model="ARD", marginal=FALSE)
pdf("ASR_marginal_false.pdf")
cols <- setNames(c("blue","green","snow2","yellow"), sort(unique(states)))
h<-max(nodeHeights(mytree))
offset.factor<-1.01 ## increase this for greater offset
plotTree(rescale(mytree,model="depth",depth=offset.factor*h),
    color="transparent",ftype="i",type="fan",fsize=0.25,lwd=2)
par(fg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
plotTree(mytree, colors=cols, type="fan",fsize=0.25,ftype="i",add=TRUE,
	lwd=2,cex=0.3, label.offset=1)
nodelabels(pie=fit.ace.false$lik.anc,cex=0.3, piecol=c("blue","green","snow2","yellow"))
dev.off()


