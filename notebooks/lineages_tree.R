library("ape")
library("phytools")
library(seqinr)
library(RColorBrewer)

setwd('/storage/btbgenie/combined_results/')
meta <- read.table('/storage/btbgenie/combined_results/metadata.csv',sep=',',
                      header=TRUE,row.names=2)
meta[is.na(meta)] <- "-"
mltree <- read.tree('/storage/btbgenie/combined_results/tree.newick')

meta1 <- read.table('/storage/btbgenie/global_results/metadata.csv',sep=',',
                   header=TRUE,row.names=1)
meta1[is.na(meta1)] <- "-"
mltree1 <- read.tree('/storage/btbgenie/global_results/tree.newick')

plot_tree <- function(mltree,samples,type='phylogram',title='',column='SB',cmap="Set1") {
  labels <- samples[mltree$tip.label,][[column]]
  #print(labels)
  labels[is.na(labels)] <- "Other"
  
  leglabels <- levels(as.factor(labels))
  n<-length(leglabels)
  colors <- brewer.pal(n = n, name = cmap)
  cols<-setNames(colors[1:length(leglabels)],leglabels)
  
  plot(mltree, type=type,cex=.2,label.offset=5, edge.width=.5,show.tip.label=FALSE)
  title(title,cex.main= 2)
  tiplabels(pie=to.matrix(labels, levels(as.factor(labels))),cex=.2,piecol=cols)
  add.scale.bar(x=100,lwd=2, cex=1)
  legcolors <- cols[leglabels]
  legend("topleft", legend=names(cols), pch=22, pt.bg=legcolors, pt.cex=2.5, bty='n')
}

#png('tree.png',width=1200,height=1800,pointsize=26)
plot_tree(mltree,meta,'fan',column='Country',cmap='Paired')

#sub-sample

randtips<-sample(mltree$tip.label,400)
randtree <- drop.tip(mltree, randtips)
plot_tree(randtree,meta,column='SB1',type='fan')

#ireland
tips<-row.names(meta[meta$Country=='Ireland',])
itree <- drop.tip(mltree, tips)
plot_tree(itree,meta,column='county',type='fan')

#non-Ireland
plot_tree(mltree1,meta1,'fan',column='geo_loc_name_country',cmap='Paired')

row.names(tips)
