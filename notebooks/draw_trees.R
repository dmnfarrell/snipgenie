# draw trees

library("ape")
library("phytools")
library(seqinr)
library(RColorBrewer)

setwd('/storage/btbgenie/')
samples <- read.table('/storage/btbgenie/mbovis_ireland/all_ireland_samples.csv',sep=',',header=TRUE,row.names=1)
samples[is.na(samples)] <- "-"

mltree <- read.tree('/storage/btbgenie/all_ireland_results/tree.newick')
mltree <- root(mltree,"TB20-002800")
mltree <- drop.tip(mltree, c('ref','182-MBovis'))

# Getnumber of positions in alignment
alignment <- read.alignment("/storage/btbgenie/all_ireland_results/core.fa", format = "fasta")
#nsites <- nchar(alignment$seq[[1]])
# Rescale the branch lengths
#mltree$edge.length <- mltree$edge.length * nsites

labels <- samples[mltree$tip.label,]$county
labels[is.na(labels)] <- "-"
counties <- levels(as.factor(labels))
n <- length(counties)

colors <- c('blue','green','red','blueviolet','orange','cadetblue','chartreuse','chocolate',
           'coral','gold','cornflowerblue','cornsilk','khaki','orange','pink','burlywood',
           'gold','cyan','lightblue','purple','salmon','maroon','beige',
           'yellow','gray','lightgreen')

cols<-setNames(colors,levels(as.factor(labels)))
cols

#png('all_ireland_tree.png',width=1200,height=1200)
plot(mltree, type='fan',cex=1,label.offset=1, edge.width=1,show.tip.label=FALSE)
tiplabels(pie=to.matrix(labels, levels(as.factor(labels))),cex=0.2,piecol=cols)
add.scale.bar(lwd=3, cex=1, xpd=TRUE)
counties<-c("Clare","Monaghan","NI","Wicklow","Cork","Cavan","Meath")
legcolors <- cols[counties]
legend("topleft", legend=counties, pch=22, pt.bg=legcolors, pt.cex=2.5)
#dev.off()


