#create M.bovis phylogeny from outbreak simulation

library(ape)
library('TransPhylo')

set.seed(0)
neg=100/365
off.r=5
w.shape=10
w.scale=0.1
pi=0.25

simu <- simulateOutbreak(neg=neg,pi=pi,off.r=off.r,w.shape=w.shape,
                         w.scale=w.scale,dateStartOutbreak=2017,dateT=2020)

plot(simu)
ttree<-extractTTree(simu)
plot(ttree)
ptree<-extractPTree(simu)
plot(ptree)
p<-phyloFromPTree(ptree)
plot(p)
axisPhylo(backward = F)
write.tree(p,'simtree.nwk')
