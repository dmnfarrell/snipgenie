# Create phylogeny from outbreak simulation
# D. Farrell May 2021

library(ape)
library(phangorn)
library('TransPhylo')
library('TreeDist')
setwd('~/gitprojects/snipgenie/notebooks/')

set.seed(2)
neg=80/365
off.r=5
w.shape=10
w.scale=0.1
pi=0.5

simu <- simulateOutbreak(neg=neg, nSampled=100, pi=pi, off.r=off.r, w.shape=w.shape,
                         w.scale=w.scale, dateStartOutbreak=2010, dateT=2017)

print(simu)

ttree<-extractTTree(simu)
plot(ttree)
ptree<-extractPTree(simu)
plot(ptree)
p<-phyloFromPTree(ptree)
plot(p)
axisPhylo(backward = F)
write.tree(p,'sim.newick')
