# draw all ireland trees
# btbgenie project

library("ape")
library("phytools")
library(seqinr)
library(RColorBrewer)

setwd('/storage/btbgenie/')
samples <- read.table('/storage/btbgenie/mbovis_ireland/ireland_metadata.csv',sep=',',
                      header=TRUE,row.names=1)
samples[is.na(samples)] <- "-"

mltree <- read.tree('/storage/btbgenie/all_ireland_results/tree.newick')
#mltree <- root(mltree,"TB20-002800")
#mltree <- drop.tip(mltree, c('ref','182-MBovis'))

plot_tree <- function(mltree,type='phylogram',title='',column='clade') {
  labels <- samples[mltree$tip.label,][[column]]
  print(labels)
  labels[is.na(labels)] <- "Other"
  
  leglabels <- levels(as.factor(labels))
  n<-length(leglabels)
  colors <- brewer.pal(n = n, name = "Set1")
  cols<-setNames(colors[1:length(leglabels)],leglabels)
  
  plot(mltree, type=type,cex=.6,label.offset=5, edge.width=.5,show.tip.label=FALSE)
  title(title,cex.main= 2)
  tiplabels(pie=to.matrix(labels, levels(as.factor(labels))),cex=.16,piecol=cols)
  add.scale.bar(x=100,lwd=2, cex=1)
  #counties<-c("Clare","Monaghan","NI","Wicklow","Cork","Cavan","Meath")
  legcolors <- cols[leglabels]
  legend("topright", legend=names(cols), pch=22, pt.bg=legcolors, pt.cex=2.5, bty='n')
  #text<-samples[mltree$tip.label,]$county
  #tiplabels(text=text,bg=rgb(0,0,0,0),offset=2,adj=0,frame='n',cex=.5)
}

#label clades
clades_label <- function(){
  nodes = list(clade1=50,clade2=300,clade3=400)
  for (n in names(nodes)){
    cladelabels(tree=NULL, n, node=nodes[n], cex=1)
  }
}

#png('all_ireland_tree_county.png',width=1200,height=1800,pointsize=26)
plot_tree(mltree,title='',column='county1')
#dev.off()

#sub-sample
#png('all_ireland_tree_subsampled.png',width=1800,height=1200,pointsize=26)
randtips<-sample(mltree$tip.label,420)
randtree <- drop.tip(mltree, randtips)
plot_tree(randtree,column='clade')#,type='fan')
#dev.off()

#common counties only
#png('all_ireland_tree_common.png',width=1200,height=1200,pointsize=26)
sub<-samples[!samples$county %in% c('Monaghan','Clare'),]
subtree <- drop.tip(mltree,row.names(sub))
plot_tree(subtree, column='county')
#dev.off()

#combined-results


#----------------
#phytools tree
cols<-setNames(colors[1:length(leglabels)],leglabels)
labels <- samples[mltree$tip.label,]$clade
plotTree(subtree,type="fan",lwd=.5,ftype="off",fsize=1,outline=TRUE)

