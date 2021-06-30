# draw trees

library("ape")
library("phytools")
library(seqinr)
library(RColorBrewer)

setwd('/storage/btbgenie/')
samples <- read.table('/storage/btbgenie/mbovis_ireland/all_ireland_samples.csv',sep=',',header=TRUE,row.names=1)
#samples <- samples[samples$county == "Wicklow",]

mltree <- read.tree('/storage/btbgenie/all_ireland_results/tree.newick')
#rooted <- root(mltree,"ref")
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
           'red','mediumvioletred','navy','lightblue','purple','salmon','maroon')

#library(randomcoloR)
#colors <- distinctColorPalette(n)
#colors<-brewer.pal(11, "rainbow")
cols<-setNames(colors,levels(as.factor(labels)))
cols
# Add coloured shapes at tips
#tiplabels(pch=19, 
#           col=ifelse(grepl(mltree$tip.label, pattern="lwof"), "red", "blue"), 
#           cex=2)


#png('all_ireland_tree.png',width=1200,height=1200)
plot(mltree, type='fan',cex=1,label.offset=1, edge.width=1,show.tip.label=FALSE)
tiplabels(pie=to.matrix(labels, levels(as.factor(labels))),cex=0.18,piecol=cols)
add.scale.bar(lwd=3, cex=2, xpd=TRUE)
legcolors <- c("red","khaki","orange","lightblue","blueviolet","green","cornsilk")
legend("topleft", legend=c("Clare","Monaghan","NI","Wicklow","Cork","Cavan","Meath"), pch=22, pt.bg=legcolors, pt.cex=2.5)
#dev.off()

plotTree(mltree,ftype="off",ylim=c(-0.2*Ntip(mltree),Ntip(mltree)),lwd=1,
         xlim=c(max(nodeHeights(mltree)),0),direction="leftwards")

