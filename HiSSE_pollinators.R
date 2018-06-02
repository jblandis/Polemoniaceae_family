#perform a HiSSE analysis with color coded as binary for presence or absence of anthocyanins
library(phytools)
library(geiger)
library(nlme)
library(diversitree)
library(hisse)
library(devtools)
#devtools::install_github("thej022214/hisse")


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

#output types can be any of the following: "turnover", "net.div", or "raw"
###################################################################################
#####run the full model with different turnover.anc, eps.anc, and transition rates run 1
#full hisse model
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual,output.type="net.div")

pp
###################################################################################
###Bissee model all free, run 2
turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)

#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
trans.rates.bisse

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse)

pp
###################################################################################
###Bissee model, extinction 0=1, run 3
turnover.anc = c(1,2,0,0)
eps.anc = c(1,1,0,0)

#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
trans.rates.bisse

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse)

pp
###################################################################################
###Bissee model with equal q's, run 4
turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)

#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
trans.rates.bisse
trans.rates.bisse.equal = ParEqual(trans.rates.bisse, c(1,2))
trans.rates.bisse.equal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse.equal)

pp
###################################################################################
###Bissee model with equal q's and e0=e1, run 5
turnover.anc = c(1,2,0,0)
eps.anc = c(1,1,0,0)

#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
trans.rates.bisse
trans.rates.bisse.equal = ParEqual(trans.rates.bisse, c(1,2))
trans.rates.bisse.equal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse.equal)

pp
###################################################################################

#2-state character independent CID-2 model, equal q's but different transition rates, run 6
#null two model, traits A and B have a diversification rate
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)

#full 8 transition model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))

#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#2-state character independent CID-2 model, equal q's and e's but different transition rates, run 7
#null two model, traits A and B have a diversification rate
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,1,1)

#full 8 transition model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))

#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
######CD-4 model, q's equal, run 8

#full 8 transition model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))

#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp.hisse.null4 <- hisse(mytree, mydata, f=c(0.40,0.40), turnover.anc=rep(c(1,2,3,4),2),eps.anc=rep(c(1,2,3,4),2), trans.rate=trans.rates.nodual.allequal)

pp.hisse.null4
###################################################################################
######CD-4 model, q's and e's equal, run9
eps.anc = c(1,1,1,1)
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))

#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp.hisse.null4 = hisse(mytree, mydata, f=c(0.40,0.40), turnover.anc=rep(c(1,2,3,4),2),eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal)

pp.hisse.null4
###################################################################################
#####run the full model with transition rates equal, run 10
#full hisse model
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model with transition rates and e's equal, run 11
turnover.anc = c(1,2,3,4)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a=0b, and extinction 0a=1a=0b run 12
turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,2)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a=0b, and extinction equal run 13
turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=0b, and e0a=e0b, run 14
turnover.anc = c(1,2,1,3)
eps.anc = c(1,2,1,3)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=0b, and extinction equal, run 15
turnover.anc = c(1,2,1,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a, and e0a=e1a, run 16
turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,2,3)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a, and extinction equal, run 17
turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model with different turnover.anc, eps.anc, and but rates q0b1b=0,q1b0b=0, all other equals, run 18
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model with different e's equal and but rates q0b1b=0,q1b0b=0, all other equals, run 19
turnover.anc = c(1,2,3,4)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model with r0a=r1a=r0b, e0a=e0b=e0b and rates q0b1b=0,q1b0b=0, all other equals, run 20
turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,2)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model with r0a=r1a=r0b, extinction equals and rates q0b1b=0,q1b0b=0, all other equals, run 21
turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 22
turnover.anc = c(1,2,1,3)
eps.anc = c(1,2,1,3)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 23
turnover.anc = c(1,2,1,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 24
turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,2,3)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################
#####run the full model withr0a=r0b, extinction equal, and rates q0b1b=0,q1b0b=0, all other equals, run 25
turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp
###################################################################################

#Run MarginRecon after rerunning the best fit model, may need to rename the variable for pars to best model name
pp.recon <- MarginRecon(mytree, mydata, f=c(0.40,0.40), pars=pp.hisse.null4$solution, hidden.states=TRUE)

#can plot turnover,net diversification ("net.div"), speciation, extinction, or extinction fraction ("extinction.fraction")
#to make the figure more easliy understood, pass a range of min and max of all models by rate.range=c(min,max), if you leave this blank it will select automatically
plot.hisse.states(pp.recon, rate.param="net.div", show.tip.label=FALSE, rate.range=c( 0.00123711868605535,0.0902830996788747))

##adaptive sappling of the MLE surface for support
SupportRegion(pp.hisse.null4, n.points=1000, scale.int=0.1, desired.delta=2, output.type="raw", hidden.states=TRUE, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, verbose=TRUE)
  

