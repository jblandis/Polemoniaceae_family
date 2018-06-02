#perform a BiSSe analysis with pollinators coded as binary
#since there is missing trait data, the tree has to be trimmed
library(phytools)
library(ape)
library(geiger)
library(nlme)
library(diversitree)

setwd("~/Bisse")
source("~/Bisse/traitDependent_functions.R")

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

#proportion of tips sampled.  For the Pole phylogeny its a total of 429 tips, 471 for total sampling
#but for pollination it is only 191 out of 471
sampling.f <- 191 / 471


res <- FISSE.binary(mytree, states)
res
