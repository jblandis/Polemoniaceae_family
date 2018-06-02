library(phytools)
library(geiger)
library(nlme)
library(diversitree)

#starting files of the tree and the data file
mydata <- read.csv("Autogamy_and_out.csv",row.names=1)
mytree <- read.nexus("MCC.fixed.nex")

#compares names between the tree and the data to list any discrepancies
comparison <- name.check(phy=mytree,data=mydata)

# prune taxa that don't have data but are present in the tree
mytree <- drop.tip(mytree,comparison$tree_not_data)

#double check to make sure that taxa all match with tree and data
name.check(phy=mytree,data=mydata)
comparison <- name.check(phy=mytree,data=mydata)

#create vector for traits
# if the vector names are changed, then commands following this will need to be altered to maintain the procedure
states <- mydata[,1]
names(states) <- row.names(mydata)

#ancestral state reconstruction
sampling.f <- 192 / 471
lik <- make.bisse(mytree,states,sampling=sampling.f)
p <- starting.point.bisse(mytree)
fit <- find.mle(lik,p)
st <- asr.marginal(lik, coef(fit))
#this method produces a circle tree, easy to read
#need to use this code for doing stochastic character mapping
pdf("Pollinator_ASR.pdf")
state.colors <- c("gray", "black")
plot(mytree, type = "fan", label.offset=1.5, cex = 0.3, no.margin=TRUE)
nodelabels(pie = t(st), frame ="circle", piecol = state.colors, cex =  0.35)
dev.off()


#MK2 ancestral state reconstructions ignoring shifts in diversification
pdf("Pollination_MK2.pdf")
lik.mk2 <- make.mk2(mytree, states)
p <- c(.1,.1)
fit.mk2 <- find.mle(lik.mk2, p)
coef(fit.mk2)
logLik(fit.mk2)
st.mk2 <- asr.marginal(lik.mk2, coef(fit.mk2))
plot(mytree, type = "fan", label.offset=1.5, cex = 0.3, no.margin=TRUE)
state.colors <- c("gray", "black")
nodelabels(pie = t(st.mk2),  frame ="circle", piecol = state.colors, cex = 0.35, adj = -0.5)
dev.off()
