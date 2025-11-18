#https://www.randigriffin.com/2017/05/11/primate-phylogeny-ggtree.html
#https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html

library("ape")
#library(seqinr)
library(RColorBrewer)
library(dplyr)
library('ggplot2')
library('ggtree')
library(tidytree)
library(ggnewscale)

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

prune_tree <- function(tree, k=10){

    dist_matrix <- cophenetic(tree)
    # Perform hierarchical clustering
    hc <- hclust(as.dist(dist_matrix))
    # Cut the tree into clusters (set k to the desired number of clusters)
    clusters <- cutree(hc, k = k)
    # Select one tip per cluster
    representative_tips <- sapply(unique(clusters), function(cluster) {
      tip <- names(clusters)[clusters == cluster][1]
      return(tip)
    })
    prunedtree <- drop.tip(tree, setdiff(tree$tip.label, representative_tips))
    return(prunedtree)
}

prune_tree_by_distance <- function(tree, threshold = 0.1) {

  # Compute the cophenetic distance matrix
  dist_matrix <- cophenetic(tree)
  # Perform hierarchical clustering
  hc <- hclust(as.dist(dist_matrix))
  # Cut the tree based on a distance threshold
  clusters <- cutree(hc, h = threshold)
  # Select one representative tip per cluster
  representative_tips <- sapply(unique(clusters), function(cluster) {
    tip <- names(clusters)[clusters == cluster][1]
    return(tip)
  })
  # Prune the tree to keep only the representative tips
  pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, representative_tips))
  return(pruned_tree)
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

#' @title ggplottree
#' @description  convenvience function for ggtree drawing
#' @param tree R ape tree
#' @param meta metadata for coloring and labels, an R data.frame
#' @param cols name of color column(s) as list
#' @return ggtree
ggplottree <- function(tree, meta, cols=NULL, colors=NULL, cmaps=NULL, legends=NULL,
                       layout="rectangular", offset=10,
                       tiplabel=FALSE, tipsize=3, tiplabelsize=5, tiplabelcol=NULL, tipalpha=NULL,
                       legend.fontsize=18, legend.position='left', legend.pointsize=8, size=.5,
					   align=FALSE, tipoffset=0, scalebar=FALSE, scalebar.width=10,
                       gradient=FALSE) {

    #make legends list if not provided
    if (is.null(legends)) {
        legends <- vector("list", length(cols)+1)
        for (i in seq_along(cols)) {
          legends[[i]] <- TRUE
        }
    }
    y <- gettreedata(tree, meta)
    if (layout == 'cladogram'){
        p <- ggtree(y, layout='c', branch.length='none')
    }
    else {
        p <- ggtree(y, layout=layout, size=size)
    }

    if (is.null(cols)) {
        if (tiplabel){
            p <- p + geom_tiplab(size=tiplabelsize,offset=tipoffset)
        }
        return (p)
    }
    col <- cols[1]
    df <- meta[tree$tip.label,][col]
    if (!is.null(colors) && !is.null(colors[[1]])) {
        #use predefined colors
        clrs <- colors[[1]]
    }
    else {
        #calculate colors from cmap
        cmap <- cmaps[1]
        clrs <- get_color_mapping(df, col, cmap)
    }

    p <- p + new_scale_fill() +
            geom_tippoint(mapping=aes(fill=.data[[col]]),
                          size=tipsize,shape=21,stroke=0,show.legend=legends[[1]]) +
            scale_fill_manual(values=clrs, na.value="black",
                guide = guide_legend(override.aes = list(size=legend.pointsize)))
    p2 <- p
    current_offset <- max(p2$data$x)/100
    heatmap_width <- .1
    constant_gap <- 50
    if (length(cols)>1){
        for (i in 2:length(cols)){
            col <- cols[i]
            df <- meta[tree$tip.label,][col]
            #type <- class(df[col,])
            type <- sapply(meta,class)[col]
            if (length(colors)>=i) {
                clrs <- colors[[i]]
            }
            else if (length(cmaps)>=i){
                #use predefined colors
                cmap <- cmaps[i]
                clrs <- get_color_mapping(df, col, cmap)
            }
            else {
                cmap = 'Greys'
                clrs <- get_color_mapping(df, col, cmap)
            }
            p2 <- p2 + new_scale_fill()
            if (i==1){
                current_offset <- current_offset + constant_gap
            }
            print(current_offset)
            p2 <- gheatmap(p2, df, offset=i*offset, width=heatmap_width,
                      colnames_angle=90, colnames_offset_y=.001, color=NA)
            current_offset <- current_offset + heatmap_width + constant_gap
            if (type == 'integer' && gradient == TRUE) {
				#p2 <- p2 + scale_fill_gradient(low='#F8F699',high='#06A958', na.value="white")
                clrs = brewer.pal(9, cmap)
                p2 <- p2 + scale_fill_gradientn(colors=clrs, na.value="white")
            }
            else {
                #hide legend if required
                p2 <- p2 + scale_fill_manual(values=clrs, name=col, na.value="white")
            }
            if (legends[[i]] == FALSE){
                    p2 <- p2 + guides(colour='none', fill='none')
            }
        }
    }

    p2 <- p2 + theme_tree(legend.text = element_text(size=legend.fontsize), legend.key.size = unit(.7, 'cm'),
                        legend.position=legend.position, plot.title = element_text(size=22),
                        legend.title=element_text(size=legend.fontsize+2) )
    #p2 <- p2 +guides(color = guide_legend(override.aes = list(size=20)))

    if (tiplabel) {
		if (!is.null(tiplabelcol)) {
			p2 <- p2 + geom_tiplab(mapping=aes(label=.data[[tiplabelcol]]),
								size=tiplabelsize, align=align,offset=tipoffset)
		}
		else {
        	p2 <- p2 + geom_tiplab(size=tiplabelsize, align=align, offset=tipoffset)
		}
    }
    if (scalebar) {
        p2 <- p2 + geom_treescale(x=NULL, y=-10, width=scalebar.width, offset = NULL,
                       label='SNPs', color="black", linesize = 1.2, fontsize = 7)
    }
    return(p2)
}

plotcladelabels <- function(tree, meta, col, p, align=FALSE) {
    data <- meta[rownames(meta) %in% tree$tip.label,]
    for (cat in unique(data[[col]])) {
        tips <- rownames(data[data[[col]]==cat,])
        mrca_node <- MRCA(tree, tips)
        #print (mrca_node)
        p <- p + geom_cladelab(node=mrca_node, label=cat, fontsize=6, horizontal=FALSE, barsize=0,
                               align=align, offset = 0.001)
    }
    return (p)
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
