#perform a MuSSE analysis for scored traits in Polemoniaceae
setwd("/Users/jacoblandis/Desktop/Pole_MuSSE")
library(diversitree)
library(phytools)
library(ape)
library(geiger)
library(nlme)
library(optparse)
library(picante)

#starting files of the tree and the data file
mydata <- read.csv("Color.csv",row.names=1)
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

#proportion of tips sampled.  For the Pole phylogeny its a total of 429 tips, 471 for total sampling, for color all 429 in the phylogeny
#have data.  Ran with and without this proportion change, and it didn't affect p-value much
sampling.f <- 429 / 471

#make.musse needs two arguments, a tree and character states
lik <- make.musse(mytree,states,4,sampling.f=sampling.f)

#argument names
argnames(lik)

#heuristic guess as to sensible starting point
p <- starting.point.musse(mytree,4)

fit.base <- find.mle(lik, p[argnames(lik)])
fit.base

#round to three decimal places
round(coef(fit.base),3)

save(fit.base, file="fit.base.RData")
save.image()
#load(file="fit.base.RData")

#test hypothesis of speciation rates are different for different states
lik.lambda <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1, lambda4 ~ lambda1)

#then start search again
fit.lambda <- find.mle(lik.lambda, p[argnames(lik.lambda)])
fit.lambda

#round to three decimal places
round(coef(fit.lambda),3)

save(fit.lambda, file="fit.lambda.RData")
save.image()


#perform statistical test to get a p-value
anova(fit.base, free.lambda=fit.lambda)

#then start search again
lik.mu <- constrain(lik, mu2 ~ mu1, mu3 ~ mu1, mu4 ~ mu1)

#then start search again
fit.mu <- find.mle(lik.mu, p[argnames(lik.mu)])
fit.mu

round(coef(fit.mu),3)

save(fit.mu, file="fit.mu.RData")
save.image()


#perform statistical test to get a p-value
anova(fit.base, free.lambda=fit.mu)
anova

#mcmc analysis
prior.exp <- make.prior.exponential(2)
#prior <- make.prior.exponential(1 / (2 * (p[1] - p[4])))
set.seed(1)
samples <- mcmc(lik, coef(fit.base), fail.value=NULL, nsteps=10000, prior=prior.exp, w=1, print.every=50)

#save distribution
save(samples, file="Samples_distribution.RData")
#if need be, can reload file
#load("Samples_distribution.RData")
#############################################################
#plotting results from MuSSE run
#to visualize speciation rate
pdf("Speciation.pdf")
col <- c("snow3","blue", "yellow", "green2")
profiles.plot(samples[c("lambda1", "lambda2", "lambda3", "lambda4")], col.line=col, las=1, 
	xlab="Speciation rate", legend=NULL, font=2, cex.lab=1, font.lab=1, cex.legend=1, margin=1/4, opacity=0.4, main="Speciation Rate")
legend("topright", legend=c("Non pigmented", "Anthocyanin", "Carotenoid", "Chlorophyll"), fill=c("snow3","blue", "yellow", "green2"))
dev.off()

#to visualize extinction rate
pdf("Extinction_rate.pdf")
col <- c("snow3","blue", "yellow", "green2")
profiles.plot(samples[c("mu1", "mu2", "mu3", "mu4")], col.line=col, las=1, 
    xlab="Extinction rate", legend=NULL, font=2, cex.lab=1, font.lab=1, cex.legend=1, margin=1/4, opacity=0.4, main="Extinction Rate")
legend("topright", legend=c("Non pigmented", "Anthocyanin", "Carotenoid", "Chlorophyll"), fill=c("snow3","blue", "yellow", "green2"))
dev.off()

#visualize net diversification
pdf("Net_diversification.pdf")
samples$net1 <- samples$lambda1 - samples$mu1
samples$net2 <- samples$lambda2 - samples$mu2
samples$net3 <- samples$lambda3 - samples$mu3
samples$net4 <- samples$lambda4 - samples$mu4
col <- c("snow3","blue", "yellow", "green2")
profiles.plot(samples[c("net1", "net2", "net3", "net4")], col.line=col, las=1, 
  xlab="Net diversification rate", legend=NULL, font=2, cex.lab=1.0, font.lab=1, cex.legend=1, margin=1/4, opacity=0.4, main="Net Diversification")
legend("topleft", legend=c("Non pigmented", "Anthocyanin", "Carotenoid", "Chlorophyll"), fill=c("snow3","blue", "yellow", "green2"))
dev.off()



##########################################################################
#perform ancestral state reconstructions incorporating diversification results
#ancestral state reconstruction
sampling.f <- 429 / 471
lik <- make.musse(mytree,states,4,sampling.f=sampling.f)
p <- starting.point.musse(mytree,4)
fit <- find.mle(lik, p[argnames(lik)])
st <- asr.marginal(lik, coef(fit))
#this method produces a circle tree, easy to read
#need to use this code for doing stochastic character mapping
pdf("Flower_color_ASR.pdf")
state.colors <- c("snow2","blue","yellow","green")
plot(mytree, type = "fan", label.offset=1.5, cex = 0.25, no.margin=TRUE)
nodelabels(pie = t(st), piecol = state.colors, cex = 0.35)
dev.off()


#MKn ancestral state reconstructions ignoring shifts in diversification
pdf("Flower_color_MKN.pdf")
lik.m <- make.mkn(mytree, states,4)
p <- c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1)
fit.m <- find.mle(lik.m,p, method = "subplex")
st.m <- asr.marginal(lik.m, coef(fit.m))
state.colors <- c("snow2","blue","yellow","green")
plot(mytree, type = "fan", label.offset=1.5, cex = 0.25, no.margin=TRUE)
nodelabels(pie = t(st.m), piecol = state.colors, cex = 0.34, adj = -0.5)
dev.off()
