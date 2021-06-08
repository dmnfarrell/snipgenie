# Create M.bovis phylogeny from outbreak simulation
# D. Farrell May 2021

library(ape)
library(phangorn)
library('TransPhylo')

setwd('~/gitprojects/snipgenie/notebooks/')

set.seed(0)
neg=100/365
off.r=5
w.shape=10
w.scale=0.1
pi=0.25

#simulate outbreak
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
write.tree(p,'sim.newick')

seq <- read.dna(file="Mbovis_AF212297.fa",format='fasta')#,as.character=T)

seqdata <- simSeq(p, Q = NULL, bf = NULL,
           rootseq = c(seq), type = "DNA", model = NULL, levels = NULL,
           rate = 1, code = 1)

write.phyDat(seqdata, file="temp.dat", format="sequential", nbcol = -1,
             colsep = "")

