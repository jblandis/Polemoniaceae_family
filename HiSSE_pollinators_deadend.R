#perform a HiSSE analysis with selfing or outcrossing coded as binary 
setwd("~/Hisse")

library(ape)
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
###################################################################################
#####run the full model with different turnover.anc, eps.anc, and transition rates run 1
#HiSSE model testing selfing as an evolutionary dead end
#contrained no transitions from selfing to outcrossing
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits, and transition from observed selfing to observed outcrossing
trans.rates.nodual = ParDrop(trans.rates, c(3,4,5,8,10))
trans.rates.nodual

pp = hisse(mytree, mydata, f=c(0.40,0.40), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual,output.type="net.div")

pp
###################################################################################

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

