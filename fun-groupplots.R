
gene_pheatmap <- function(exprdat,sampleid,annotation_row=NULL) {
  sampleDists <- dist(t(exprdat))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- sampleid
  if(!is.null(annotation_row)) rownames(annotation_row) <- sampleid
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows=sampleDists,
                     clustering_distance_cols=sampleDists,
                     annotation_row=annotation_row,
                     col=colors)
}

gene_pcaplot <- function(exprdat,sampleid,groupdat=NULL,colorfactor=NULL,shapefactor=NULL,
                         plot_sampleids=TRUE, pcnum=1:2, plottitle = "PCA Plot") {
  #adapted from DESeq2:::plotPCA.DESeqTransform
  pca <- prcomp(t(exprdat))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if(is.null(groupdat)) groupdat = data.frame("group"=rep(1,ncol(exprdat)))
  intgroup = colnames(groupdat)
  allgroup <- if (length(intgroup) > 1) {
    factor(apply(groupdat, 1, paste, collapse = ":"))
  }else{allgroup <- intgroup}
  d <- data.frame(PC1 = pca$x[, pcnum[1]], PC2 = pca$x[, pcnum[2]], uniquegroup = allgroup, 
                  groupdat, name = sampleid)
  percentVar <- round(100 * percentVar)
  if(is.null(colorfactor)) {d$color=as.factor(1)}else{
    colnames(d)[colnames(d)==colorfactor] <- "color"}
  if(is.null(shapefactor)) {d$shape=as.factor(1)}else{
    colnames(d)[colnames(d)==shapefactor] <- "shape"
  }
  if(identical(shapefactor,colorfactor)) {d$shape = d$color}
  p <- ggplot(d, aes(PC1, PC2, color=color, shape=shape, size=3)) 
  if(plot_sampleids) {
    p <- p + geom_text(aes(label=name,size=10))
  }else{
    p <- p + geom_point()
  }
  
  if(!is.null(colorfactor)) {
    p <- p + guides(color=guide_legend(title=colorfactor))
  }else {
    p <- p + guides(color = "none")
  }
  if(!is.null(shapefactor)) {
    p <- p + guides(shape=guide_legend(title=shapefactor))
  }else{
    p <- p + guides(shape = "none")
  }
  p <- p + guides(size= "none") + theme_bw() + 
    xlab(paste0("PC",pcnum[1],": ",percentVar[pcnum[1]],"% variance")) +
    ylab(paste0("PC",pcnum[2],": ",percentVar[pcnum[2]],"% variance")) + ggtitle(plottitle)
  
  return(p)
}

