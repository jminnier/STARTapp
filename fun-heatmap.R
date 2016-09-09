## ==================================================================================== ##
# START Shiny App for analysis and visualization of transcriptome data.
# Copyright (C) 2016  Jessica Minnier
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# You may contact the author of this code, Jessica Minnier, at <minnier@ohsu.edu>
## ==================================================================================== ##

##================================================================================##
## HEATMAP
##================================================================================##
heatmap_colors = c("#2D004B","white","gold3")
heatcols = colorRampPalette(heatmap_colors)(500)


#Returns data for a heatmap
heatmap_subdat <- function(data_analyzed,
                           yname="log2cpm",
                           orderby="significance",
                           usesubset=FALSE,
                           subsetids=NULL,
                           FDRcut=0.05,maxgenes=NULL,
                           view_group=NULL,
                           sel_test=NULL,
                           filter_by_go=FALSE,
                           filter_fdr=FALSE,
                           filter_maxgene=TRUE,
                           #filter_maxgeneN="genesN",
                           filter_cpm=FALSE,
                           filter_fc=FALSE,
                           fold_change_range=NULL,
                           fold_change_groups=NULL,
                           group_filter_range =NULL,
                           fixed_genes_list=NULL) {
  
  data_results = data_analyzed$results
  data_results$unique_id = as.character(data_results$unique_id)
  
  if(is.null(view_group)) view_group=data_analyzed$group_names
  
  if((usesubset)&(!is.null(subsetids))) {thesegenes = subsetids}else{
    
    if(orderby=="significance") {
      res = data_results%>%filter(test==sel_test)
      #Order by FDR
      if(is.null(res$adj.P.Val)) res$adj.P.Val = res$P.Value
      #tmpout = res[order(res$adj.P.Val),]
      tmpout = res%>%arrange(adj.P.Val,P.Value,abs(logFC))
      
      thesegenes = tmpout$unique_id
      print(paste("start",length(thesegenes)))
      
      if(length(thesegenes)==0) {return(NULL)}
      
      ## FILTER GENES by
      
      #FDR cut
      if((filter_fdr)&&(!is.null(FDRcut))) {
        #print("yes"); print(FDRcut)
        tmpg = tmpout$unique_id[which(tmpout$adj.P.Val<FDRcut)]
        thesegenes = intersect(thesegenes,tmpg)
        print(paste("fdrcut",length(thesegenes)))
      }
      
      if(length(thesegenes)==0) {return(NULL)}
      
      #Range of fold change, fold_change_range[1:2], based on two selected groups
      if((filter_fc)&&(!is.null(fold_change_groups))&&(!is.null(fold_change_range))) {
        
        fcgroup1 = fold_change_groups[1]; fcgroup2 = fold_change_groups[2]
        res_fc = data_results%>%filter(test==paste0(fcgroup1,"/",fcgroup2))
        fc <- res_fc$logFC
        tmpg = res_fc$unique_id[which((fc >= fold_change_range[1])*(fc <= fold_change_range[2])==1)]
        thesegenes = intersect(thesegenes,tmpg)
        print(paste("filterfc",length(thesegenes)))
      }
      if(length(thesegenes)==0) {return(NULL)}
      
      tmp = tmpout%>%filter(unique_id%in%thesegenes)%>%arrange(adj.P.Val)
      thesegenes = tmp$unique_id
      
    }else {
      
      #filter by standard deviation (SD of log2cpm is coefficient of variation of unlogged data)
      tmpdat = data_analyzed$data_long
      tmpdat = tmpdat%>%filter(group%in%view_group)
      tmpdat = reshape2::dcast(tmpdat,unique_id~sampleid,value.var=yname)
      tmpsd = apply(tmpdat[,-1],1,sd)
      
      thesegenes = tmpdat$unique_id[order(tmpsd,decreasing=TRUE)]
    }
    
    if(length(thesegenes)==0) {return(NULL)}
    
    # if((!filter_maxgene)&(filter_maxgeneN=="genesN")) {
    #   maxgenes=10000; filter_maxgene=TRUE
    #   }
    
    if((filter_maxgene)&&(!is.null(maxgenes))) {
      tmpg = thesegenes[1:min(maxgenes,length(thesegenes))]
      thesegenes = tmpg
      print(paste("maxgenes",length(thesegenes)))
    }
    
  }#end subsetids
  
  
  if(length(thesegenes)==0) {return(NULL)}
  
  subdat = filter(data_analyzed$data_long,unique_id%in%thesegenes,group%in%view_group)
  data = reshape2::dcast(subdat,unique_id~sampleid,value.var=yname)
  data_unique_id = data[,1]
  data_geneids = data_analyzed$geneids[match(data_unique_id,data_analyzed$geneids$unique_id),]
  data = data[,-1]
  rown <- as.character(data_unique_id)
  rown[is.na(rown)] <- ""
  rownames(data) = rown
  return(list("data"=data,"data_geneids"=data_geneids))
  
}


heatmap_data <- function(...) {
  #create the data
  tmpdat = heatmap_subdat(...)
  if(is.null(tmpdat)) {frame()}else{	
    heatdat = tmpdat$data
    heatdat1 = as.matrix(heatdat)
    
    #order the rows same way heatmap does
    
    Rowv <- rowMeans(heatdat1, na.rm = TRUE)
    hcr <- hclust(dist(heatdat1))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    
    #hmap1 = heatmap(heatdat1,col=heatcols,margins=c(8,11),Colv=NA,scale="none",keep.dendro=TRUE,show=FALSE)
    
    heatdat = data.frame(tmpdat$data_geneids, heatdat)
    heatdat <- heatdat[rowInd,]
    heatdat
  }
}

heatmap_render <- function(yname,interactive=FALSE,heatmap_rowlabels=TRUE,...)  {
  #possible inputs
  tmpdat <- heatmap_subdat(yname,...)
  
  if(is.null(tmpdat)) {frame()}else{
    heatdat = as.matrix(tmpdat$data)
    myheatmap_colors = heatmap_colors
    heatdat_rowmean = sweep(heatdat,1,rowMeans(heatdat))
    if(min(heatdat_rowmean)>0) {myheatmap_colors = heatmap_colors[2:3]}
    if(max(heatdat_rowmean)<0) {myheatmap_colors = heatmap_colors[1:2]}
    color.palette  <- colorRampPalette(myheatmap_colors)
    
    #http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
    lmat = rbind(c(0,3),c(2,1),c(0,4))
    lwid = c(1,4)
    lhei = c(0.2,4,1)
    
    # heatmap.2(heatdat,col=color.palette,
    #           margins=c(8,11),labRow= rownames(heatdat),
    #           Colv=FALSE,
    #           lmat = lmat, lwid = lwid, lhei = lhei,
    #           key.xlab=yname,
    #           scale="row",tracecol=NULL,trace='none',density.info='none',
    #           #changed to scale="row"
    #           symkey=TRUE)
    
    #par(oma=c(0,0,0,5))
    tcols = as.numeric(as.factor(do.call(rbind,strsplit(colnames(heatdat),"_"))[,1]))+1
    

    # heatmap(heatdat_rowmean,labRow = rownames(heatdat),scale="none",col=color.palette(100),
    #         ColSideColors = rainbow(max(tcols))[tcols],
    #         margins=c(5,10))
    
    acexRow = min(0.2 + 1/log10(nrow(heatdat_rowmean)), 1.2) # default aheatmap cexRow
    lmargin = 300
    if(!heatmap_rowlabels) {acexRow = 0; lmargin = 0}
    
    if(interactive) {
      heatmaply(heatdat_rowmean,colors=color.palette(100))%>% layout(margin = list(l = lmargin, b = 100))
    }else{
      aheatmap(heatdat_rowmean,col=color.palette(100),scale = "none",
               cexRow=acexRow,
               revC=TRUE,
             annCol=data.frame("group"=as.factor(do.call(rbind,strsplit(colnames(heatdat),"_"))[,1])))
    }
  }
}


#heatmap_ggvis_render(data_analyzed,sel_group=c("group1","group2"),maxgenes=20)
#tmpdat = heatmap_subdat(data_analyzed,sel_group=c("BE","HP"),maxgenes=200)
#Not actually used in app, rendering is done in server.R
# heatmap_ggvis_render <- function(...) {
#   
#   tmpdat = heatmap_ggvis_data(...)
#   if(is.null(tmpdat)) {return(NULL)}else{	
#     heatdatmelt = tmpdat
#     
#     
#     all_values <- function(x) {
#       if(is.null(x)) return(NULL)
#       thisrow <- heatdatmelt[heatdatmelt$id == x$id,]
#       thisrow <- thisrow[,c("unique_id","sample",yname)]
#       thisrow[,yname] <- round(thisrow[,yname],3)
#       paste0(names(thisrow), ": ", format(thisrow), collapse = "<br />")
#     }
#     
#     lb <- linked_brush(keys = 1:nrow(heatdatmelt),fill = "red")
#     
#     p <- heatdatmelt%>%
#       ggvis(~sample, ~unique_id, fill:=~defcolor,stroke:=~defcolor, key:=~id)%>%
#       layer_rects(width = band(), height = band()) %>%
#       scale_nominal("x", padding = 0, points = FALSE) %>%
#       scale_nominal("y", padding = 0, points = FALSE) %>%
#       add_tooltip(html=all_values,on = "hover") %>%
#       lb$input()%>%set_options(width=600,height=1000)
#     return(p)
#     
#   }}

#Make a dataset for ggvis, needs to be melted first
heatmap_ggvis_data <- function(yname,...) {
  tmpdat = heatmap_subdat(yname,...)
  if(is.null(tmpdat)) {return(NULL)}else{  
    heatdat = tmpdat$data
    heatdat <- sweep(heatdat, 1L, rowMeans(heatdat, na.rm = TRUE), check.margin = FALSE)
    
    heatdat1 = as.matrix(heatdat)
    #browser()
    
    
    heatdat_scale = heatdat
    
    # x <- heatdat
    # x <- sweep(x, 1L, rowMeans(x, na.rm = TRUE), check.margin = FALSE)
    # sx <- apply(x, 1L, sd, na.rm = TRUE)
    # x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    # heatdat_scale = x
    
    #order the rows same way heatmap does
    
    Rowv <- rowMeans(heatdat1, na.rm = TRUE)
    hcr <- hclust(dist(heatdat1))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    myorder = rowInd
    #myorder = myorder[length(myorder):1]
    heatdat <- heatdat[myorder,]
    heatdat_scale <- heatdat_scale[myorder,]
    
    Colv <- colMeans(heatdat1,na.rm=TRUE)
    hcc <- hclust(dist(t(heatdat1)))
    ddc <- as.dendrogram(hcc)
    colInd = order.dendrogram(ddc)
    ddc <- reorder(ddc, Colv)
    myorder = order.dendrogram(ddc)
    # myorder = myorder[length(myorder):1]
    heatdat <- heatdat[,myorder]
    heatdat_scale <- heatdat_scale[,myorder]
    
    
    heatdat$unique_id = rownames(heatdat)
    heatdat$order = 1:nrow(heatdat)
    print(dim(heatdat))
    
    heatdatmelt = melt(heatdat,value.name=yname,id.vars=c("unique_id","order"),
                       variable.name="sample")
    
    
    heatdat_scale$unique_id = rownames(heatdat_scale)
    heatdat_scale$order = 1:nrow(heatdat_scale)
    print(dim(heatdat_scale))
    
    heatdatmelt_scale = melt(heatdat_scale,value.name=paste0(yname,"_scaled"),id.vars=c("unique_id","order"),variable.name="sample")
    
    heatdatmelt[,paste0(yname,"_scaled")]= heatdatmelt_scale[,paste0(yname,"_scaled")]
    #can make custom fill color as a column in heatdatmelt, then need to pass in key to get tooltip to show log2cpm and not color
    #http://stackoverflow.com/questions/24810470/custom-fill-color-in-ggvis-and-other-options
    
    heatdatmelt$defcolor = heatcols[1]
    heatdatmelt$id = 1:nrow(heatdatmelt)
    lenhc = length(heatcols)
    
    ll = heatdatmelt[,paste0(yname,"_scaled")]
    minll = -1*ceiling(max(abs(ll)))
    maxll = ceiling(max(abs(ll)))
    dd =  heatdatmelt$defcolor 
    #print(summary(ll))
    #dd[ll< 0] = as.character(cut(ll[ll<0],seq(floor(min(ll)),0,length.out=(lenhc/2+1)),labels=heatcols[1:floor(lenhc/2)]))
    #dd[ll>= 0] = as.character(cut(ll[ll>=0],seq(0,ceiling(max(ll)),length.out=(lenhc/2+1)),labels=heatcols[(floor(lenhc/2)+1):lenhc]))
    dd[ll< 0] = as.character(cut(ll[ll<0],seq(minll,0,length.out=(lenhc/2+1)),labels=heatcols[1:floor(lenhc/2)]))
    dd[ll>= 0] = as.character(cut(ll[ll>=0],seq(0,maxll,length.out=(lenhc/2+1)),labels=heatcols[(floor(lenhc/2)+1):lenhc]))
    heatdatmelt$defcolor = dd
    #to do: add scale
    
    heatdatmelt
    
    
  }}