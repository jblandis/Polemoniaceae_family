#perform a BiSSe analysis with pollinators coded as binary
#since there is missing trait data, the tree has to be trimmed
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

#proportion of tips sampled.  For the Pole phylogeny its a total of 192 tips, 471 for total sampling
sampling.f <- 192 / 471

#make.bisse needs two arguments, a tree and character states
lik <- make.bisse(mytree,states, sampling.f=sampling.f)

#heuristic guess as to sensible starting point
p <- starting.point.bisse(mytree)
p

#ML search (make take some time)
fit <- find.mle(lik,p)

#round to three decimal places
round(coef(fit),3)

#test hypothesis of speciation rates are different for different states
lik.l <- constrain(lik, lambda1 ~ lambda0)

#then start search again
fit.l <- find.mle(lik.l, p[argnames(lik.l)])

#round to thre decimal places
round(rbind(full=coef(fit), equal.l=coef(fit.l, TRUE)), 3)

#perform statistical test to get a p-value
anova(fit, equal.l=fit.l)

##test hypothesis of equal extinction rates
#test hypothesis of extinction rates are different for different states
lik.l <- constrain(lik, mu1 ~ mu0)

#then start search again
fit.l <- find.mle(lik.l, p[argnames(lik.l)])

#round to thre decimal places
round(rbind(full=coef(fit), equal.l=coef(fit.l, TRUE)), 3)

#perform statistical test to get a p-value
anova(fit, equal.l=fit.l)

##test hypothesis of equal transition rates
#test hypothesis of transition rates are different for different states
lik.l <- constrain(lik, q01 ~ q10)

#then start search again
fit.l <- find.mle(lik.l, p[argnames(lik.l)])

#round to thre decimal places
round(rbind(full=coef(fit), equal.l=coef(fit.l, TRUE)), 3)

#perform statistical test to get a p-value
anova(fit, equal.l=fit.l)

#run the mcmc chain
prior <- make.prior.exponential(1/(2*(p[1]-p[3])))
set.seed(1)
tmp <- mcmc(lik, fit$par, nsteps=100, prior=prior,lower=0, w=rep(1, 6), print.every=0)
w <- diff(sapply(tmp[2:7], range))
samples <- mcmc(lik, fit$par, nsteps=10000, w=w, lower=0, prior=prior,print.every=100)

#to visualize speciation rate
pdf("Pollinator_Speciation.pdf")
col <- c("gray", "black")
profiles.plot(samples[c("lambda0", "lambda1")], col.line=col, las=1,
              xlab="Speciation rate", legend="topright")
#abline(v=c(.1, .2), col=col)
dev.off()
   
# visualize extiction rate
pdf("Pollinator_Extinction_rate.pdf")
col <- c("gray", "black")
profiles.plot(samples[c("mu0", "mu1")], col.line=col, las=1,
    xlab="Extinction rate", legend="topright")
#abline(v=c(.1, .2), col=col)
dev.off()
  

#visualize r
pdf("Net_diversification.pdf")
samples$net0 <- samples$lambda0 - samples$mu0
samples$net1 <- samples$lambda1 - samples$mu1
col <- c("gray", "black")
profiles.plot(samples[c("net0", "net1")], col.line=col, las=1,
	xlab="Net diversification", legend.pos="topright")
#abline(v=c(.1, .2), col=col)
dev.off()


#difference in r
pdf("Pollinator_r_diff.pdf")
samples$tot <- samples$net0 - samples$net1
col <- c("gray", "black")
profiles.plot(samples[c("tot")], col.line=col, las=1,
	xlab="Net", legend.pos="topright")
#abline(v=c(.1, .2), col=col)
dev.off()


