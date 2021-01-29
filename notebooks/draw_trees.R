# draw trees

library("ape")
library("phytools")
library(seqinr)

setwd('~/gitprojects/snpgenie/')
samples <- read.table('/storage/btbgenie/mbovis_ireland/all_ireland_samples.csv',sep=',',header=TRUE,row.names=1)
#samples <- samples[samples$county == "Wicklow",]

mltree <- read.tree('/storage/btbgenie/wicklow_results/RAxML_bestTree.variants')
#rooted <- root(mltree,"ref")
mltree <- drop.tip(mltree, c('ref','182-MBovis'))
#mltree$edge.length[mltree$edge.length<0] <- 0

# Getnumber of positions in alignment
alignment <- read.alignment("/storage/btbgenie/wicklow_results/core.fa", format = "fasta")
nsites <- nchar(alignment$seq[[1]])
# Rescale the branch lengths
mltree$edge.length <- mltree$edge.length * nsites

labels <- samples[mltree$tip.label,]$species
colors <- c('red','blue','green','orange')
cols<-setNames(colors,levels(as.factor(labels)))

# Add coloured shapes at tips
#tiplabels(pch=19, 
#           col=ifelse(grepl(treeBS$tip.label, pattern="lwof"), "red", "blue"), 
#           cex=2)

#png('tree.png',width=800,height=800)
plot(mltree, type='fan',cex=1,label.offset=1, edge.width=2,show.tip.label=FALSE)
tiplabels(pie=to.matrix(labels, levels(as.factor(labels))),cex=0.3,piecol=cols)
add.scale.bar(lwd=3, cex=2, xpd=TRUE)
legend("topleft", legend=c("Cow","Deer","Badger","Dog"), pch=22, pt.bg=colors, pt.cex=2.5)

#dev.off()

plotTree(mltree,ftype="off",ylim=c(-0.2*Ntip(mltree),Ntip(mltree)),lwd=1,
         xlim=c(max(nodeHeights(mltree)),0),direction="leftwards")
