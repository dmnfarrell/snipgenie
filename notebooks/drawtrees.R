library("ape")
#library(seqinr)
library(RColorBrewer)
library(dplyr)
library('ggplot2')
library('ggtree')
library(tidytree)
#library(ggnewscale)
#library('ggtreeExtra')

plot_tree <- function(tree,samples,type='phylogram',title='',column='SB',cmap="Set1") {
  labels <- samples[tree$tip.label,][[column]]
  #print(labels)
  #print (samples[mltree$tip.label,])
  labels[is.na(labels)] <- "Other"  
  leglabels <- levels(as.factor(labels))
  n<-length(leglabels)
  colors <- brewer.pal(n = n, name = cmap)
  cols<-setNames(colors[1:length(leglabels)],leglabels)
  print (cols)
  l<-length(labels)
  w<- max(dist.nodes(tree))*.8
    
  if (l<70){
      showtip=TRUE
      }
  else {
      showtip=FALSE
  }
  plot(tree,type=type,cex=.5,label.offset=10, edge.width=.6,show.tip.label=showtip)
  title(title,cex.main= 2)
  cex<-(.3/l*100)
  tiplabels(pie=to.matrix(labels, levels(as.factor(labels))),cex=cex,size=2,piecol=cols)
  add.scale.bar(x=100,lwd=2, cex=1)
  legcolors <- cols[leglabels]
  legend("topright", legend=names(cols), pch=22, pt.bg=cols, pt.cex=2.0, cex=1.2, 
         bty="n",ncol=1,x.intersp=.3)
}

gettreedata <- function(tree, meta){
    d<-meta[row.names(meta) %in% tree$tip.label,]
    d$label <- row.names(d)
    y <- full_join(as_tibble(tree), d, by='label')
    y <- as.treedata(y)
    return(y)
}

ggplottree <- function(tree, meta, col1, col2=NULL, col3=NULL, layout="rectangular", cmap='Set1', tiplabel=FALSE, title='') {

    y <- gettreedata(tree, meta)    
    p <- ggtree(y, layout=layout)
    p1 <- p + scale_color_brewer(palette=cmap) + 
              #scale_fill_manual(values=colors) +
              #geom_text2(aes(subset=!isTip, label=node)) +
              geom_tippoint(mapping=aes( shape=NULL, color=.data[[col1]]),size=2) 
    if (tiplabel==TRUE){  
          p1 <- p1 + geom_tiplab(size=4, color="black", offset=1)
    }
    if (is.null(col2)){        
        p2<-p1
    }
    else {
        df<-meta[tree$tip.label,][col2]
        p2 <- gheatmap(p1, df, offset=70, width=.05,
                  colnames_angle=0, colnames_offset_y = .05) 
 }
         #geom_cladelab(node=10, label="test label") +
         #geom_hilight(mapping=aes(subset = node %in% c(1, 2,4,6), fill = snp50)) + 
    
    if (is.null(col3)){        
        p2<-p2
    }
    else {
        df<-meta[tree$tip.label,][col3]
        p2 <- gheatmap(p2, df, offset=140, width=.05,
                  colnames_angle=0, colnames_offset_y = .05) 
    }  
    
    p2 <- p2 + theme_tree2(legend.text = element_text(size=20), legend.key.size = unit(1, 'cm'), 
                        legend.position="left", plot.title = element_text(size=30)) +
            guides(color = guide_legend(override.aes = list(size=10))) +
            ggtitle(title) 
    
    return(p2)
}

labelclades <- function(p, labels){    
    for(l in names(labels)){        
        p <- p + geom_cladelab(node=labels[[l]], label=l, angle=0, offset=30, fontsize=9)
    }
    return(p)
}

highlightclades <- function(p){
    p <- p +
     geom_hilight(node=20, fill="lightblue", alpha=.2)
    return(p)
}

ggtreefruit <- function(tree, meta, layout='c', col1=NULL){
    y <- gettreedata(tree, meta)    
    p <- ggtree(y, layout=layout) +
            geom_tippoint( mapping=aes( shape=NULL, color=.data[[col1]]),size=3) 
    
    p <- p +
         geom_fruit(             
             geom = geom_tile,
             mapping = aes( x=node, 
                         group=RDWicklow,
                         fill=SB1),size=2,
            axis.params=list(
                         axis       = "x",
                         text.size  = 1.8,
                         hjust      = 1,
                         vjust      = 0.5,
                         nbreak     = 3,
                     ),
         grid.params=list())
    
    return(p)
}