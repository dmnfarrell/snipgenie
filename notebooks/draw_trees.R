# draw trees

library("ape")
library("phytools")
library(seqinr)
library(RColorBrewer)

setwd('/storage/btbgenie/')
samples <- read.table('/storage/btbgenie/mbovis_ireland/all_ireland_samples.csv',sep=',',header=TRUE,row.names=1)
#samples <- samples[samples$county == "Wicklow",]

mltree <- read.tree('/storage/btbgenie/all_ireland_results/tree.newick')
mltree <- root(mltree,"TB20-002800")
mltree <- drop.tip(mltree, c('ref','182-MBovis'))
#mltree$edge.length[mltree$edge.length<0] <- 0

# Getnumber of positions in alignment
alignment <- read.alignment("/storage/btbgenie/all_ireland_results/core.fa", format = "fasta")
#nsites <- nchar(alignment$seq[[1]])
# Rescale the branch lengths
#mltree$edge.length <- mltree$edge.length * nsites

labels <- samples[mltree$tip.label,]$county
counties <- levels(as.factor(labels))
n <- length(counties)
#colors <- c('red','blue','green','orange','yellow','purple','gold','brown')

colors <- c('blue','green','red','blueviolet','orange','cadetblue','chartreuse','chocolate',
           'coral','gold','cornflowerblue','cornsilk','khaki','orange','pink','burlywood',
           'mediumvioletred','navy','lightblue','purple','salmon','maroon','beige',
           'yellow','gray','lightgreen')

cols<-setNames(colors,levels(as.factor(labels)))
cols
mltree$tip.label <- labels

#png('all_ireland_tree.png',width=1200,height=1200)
plot(mltree, type='fan',cex=1,label.offset=1, edge.width=1,show.tip.label=TRUE)
tiplabels(pie=to.matrix(labels, levels(as.factor(labels))),cex=0.14,piecol=cols)
add.scale.bar(lwd=3, cex=2, xpd=TRUE)
legcolors <- c("green","burlywood","mediumvioletred","gray","red","blue","pink")
legend("topleft", legend=c("Clare","Monaghan","NI","Wicklow","Cork","Cavan","Meath"), 
       pch=22, pt.bg=legcolors, pt.cex=2.5)
#dev.off()
to.matrix(labels, levels(as.factor(labels)))
#plotTree(mltree,ftype="off",ylim=c(-0.2*Ntip(mltree),Ntip(mltree)),lwd=1,
#         xlim=c(max(nodeHeights(mltree)),0),direction="leftwards")

