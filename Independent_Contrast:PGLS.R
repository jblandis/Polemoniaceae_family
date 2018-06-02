#This R code will allow for testing correlations between two continuous traits using Indendent contrast, in this case flower length and flower width
#Towards the end of the code you can also look for phylogentic signal in each trait by calculating Bloomberg's K value and Pagel's lambda for each trait individaully

library(phytools)
library(ape)
library(geiger)
library(nlme)

#trait data is in a csv file, with column one being species names, columns 2 and 3 the continuous traits of interest.
#Species missing data for either or both traits was removed prior to analysis
mydata <- read.csv("Length_Width.csv",row.names=1)
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
flower_length <- mydata[,1]
flower_width <- mydata[,2]

names(flower_width) <- row.names(mydata)
names(flower_length) <- row.names(mydata)

#calculate contrasts
#testing whether traits are correlated with independent contrast.  The null hypothesis is that there is no correlation and the PICs should be low
#a signficant value from the correlataion analysis later on will determine if there is signifcant correlation or not
Contrastwidth <- pic(flower_width, mytree, scaled=FALSE)
Contrastlength <- pic(flower_length, mytree, scaled=FALSE)

#displaying contrasts of one variable
Contrastwidth

#plotting contrasts with one trait on the x and the other on the y axis
plot (mytree)
nodelabels(round(Contrastwidth, 3), adj = c(0, -0.5), frame="n")
nodelabels(round(Contrastlength, 3), adj = c(0, 1), frame="n")

#linear regression with the intercept going through zero
RegressWidthLength <- lm(Contrastwidth~Contrastlength -1)

#regression stats
summary.lm(RegressWidthLength)

plot(Contrastwidth, Contrastlength)

#add in line of best fit
abline(RegressWidthLength)

#PGLS to look at correlation between traits
cor.bm <- corBrownian(phy=mytree)
pgls <- gls(flower_length ~ flower_width, correlation = cor.bm)
summary(pgls)

#Bloombergs K
# K=1 indicates the observed trait is predicted by the structure of the phylogney under a Brownian motion model of trait evolution
# K>1 suggest more phylogenetic signal than expected from brownian motion
# K<1 suggests less phylogenetic signal than expected from brownian motion

phylosig(tree = mytree, x = flower_length, method = "K", test = T)
phylosig(tree = mytree, x = flower_width, method = "K", test = T)

#Pagel's lambda, ranging from 0 to 1
# A low lambda value indicates very little phylogenetic signal in the trait data given the original tree
# A high lamda indicates relatively more phylognetic signal in the trait data given the original tree
# Null hypothesis is that lamda=0

phylosig(tree=mytree, x=flower_length, method = "lambda", test = T)

phylosig(tree=mytree, x=flower_width, method = "lambda", test = T)