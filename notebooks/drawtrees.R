library("ape")
#library(seqinr)
library(RColorBrewer)
library(dplyr)
library('ggplot2')
library('ggtree')
library(tidytree)
library(ggnewscale)
#library(ggstar)

plot_tree <- function(tree,samples,type='phylogram',title='',colorcol=NULL,
						tiplabelcol=NULL,showtip=TRUE,cex=1.2,
						cmap="Set1") {

    if (!is.null(colorcol)){
        labels <- samples[tree$tip.label,][[colorcol]]
        #print(labels)
        #labels[is.na(labels)] <- "Other"
        leglabels <- levels(as.factor(labels))
        n<-length(leglabels)
        colors <- brewer.pal(n = n, name = cmap)
        clrs <- setNames(colors[1:length(leglabels)],leglabels)
        tiplabels(pie=to.matrix(labels, levels(as.factor(labels))),cex=cex,size=2,piecol=clrs)

        legcolors <- cols[leglabels]
        legend("topright", legend=names(cols), pch=22, pt.bg=cols, pt.cex=2.0, cex=1.2,
             bty="n",ncol=1,x.intersp=.3)
    }
    l<-length(labels)
    w<- max(dist.nodes(tree))*.8

    #plot tree
    plot(tree,type=type,cex=.8,label.offset=.2, edge.width=.6,show.tip.label=showtip)
    title(title,cex.main= 2)
    cex<-(.3/l*100)
    if (!is.null(tiplabelcol)){
        tl <- samples[tree$tip.label,][[tiplabelcol]]
        #print(tl)
        #tiplabels(tl,type=2,cex=.8,size=1)
        #sub.taxa.label(tree, tl)
        dotTree(tree,labels,colors=cols)
    }
    add.scale.bar(x=100,lwd=2, cex=1)
}

gettreedata <- function(tree, meta){
    d<-meta[row.names(meta) %in% tree$tip.label,]
    d$label <- row.names(d)
    y <- full_join(as_tibble(tree), d, by='label')
    y <- as.treedata(y)
    return(y)
}

get_color_mapping <- function(data, col, cmap){
    labels <- (data[[col]])
    names <- levels(as.factor(labels))
    n <- length(names)
    if (n<10){
        colors <- suppressWarnings(c(brewer.pal(n, cmap)))[1:n]
    }
    else {
        colors <- colorRampPalette(brewer.pal(8, cmap))(n)
    }
    names(colors) = names
    return (colors)
}

ggplottree <- function(tree, meta, cols=NULL, colors=NULL, cmaps=NULL, layout="rectangular",
                       offset=10, tiplabel=FALSE, tipsize=3, tiplabelsize=5, tiplabelcol=NULL,
					   align=FALSE) {

    y <- gettreedata(tree, meta)
    p <- ggtree(y, layout=layout)

    if (is.null(cols)) {
        if (tiplabel){
            p <- p + geom_tiplab(size=tiplabelsize)
        }
        return (p)
    }
    col <- cols[1]
    if (!is.null(colors)) {
        #use predefined colors
        clrs <- colors
    }
    else {
        #calculate colors from cmap
        cmap <- cmaps[1]
        df <- meta[tree$tip.label,][col]
        clrs <- get_color_mapping(df, col, cmap)
    }
    #print (clrs)
    p <- p + new_scale_fill() +
            geom_tippoint(mapping=aes(fill=.data[[col]]),size=tipsize,shape=21,stroke=0) +
            scale_fill_manual(values=clrs, na.value="black")

    p2 <- p
    if (length(cols)>1){
        for (i in 2:length(cols)){
            col <- cols[i]
            if (length(cmaps)>=i){
                cmap <- cmaps[i]
            }
            else {
                cmap = 'Greys'
            }
            df <- meta[tree$tip.label,][col]
            type <- class(df[col,])
            p2 <- p2 + new_scale_fill()
            p2 <- gheatmap(p2, df, offset=i*offset, width=.08,
                      colnames_angle=0, colnames_offset_y = .05)
            if (type == 'numeric'){
				p2 <- p2 + scale_fill_gradient(low='#F8F699',high='#06A958')
            }
            else {
                colors <- get_color_mapping(df, col, cmap)
                p2 <- p2 + scale_fill_manual(values=colors, name=col)
            }
        }
    }

    p2 <- p2 + theme_tree2(legend.text = element_text(size=20), legend.key.size = unit(1, 'cm'),
                        legend.position="left", plot.title = element_text(size=40))
            guides(color = guide_legend(override.aes = list(size=10)))
    if (tiplabel) {
		if (!is.null(tiplabelcol)) {
			p2 <- p2 + geom_tiplab(mapping=aes(label=.data[[tiplabelcol]]),
								size=tiplabelsize, align=align)
		}
		else {
        	p2 <- p2 + geom_tiplab(size=tiplabelsize, align=align)
		}
    }
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
