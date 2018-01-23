
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



## ==================================================================================== ##
## Volcano Plot
## ==================================================================================== ##  	


# change rna_volcanoplot to have an input function that depends on type of variables and then a plotting function
# need to be more careful about what p-value (adjusted or raw) is used for colors
# add option to label a set of genes (top 5, or name them)

# profvis::profvis(
#   rna_volcanoplot(results,test_sel="group1/group2",absFCcut=0,pvalcut=0.05,fdrcut=0.05)
#   )

rna_volcanoplot <- function(data_results, geneids=NULL, 
                            test_sel=NULL,absFCcut=0,pvalcut=0.05,fdrcut=0.05,
                            sel_genes=NULL) {
  
  validate(need(mean(is.na(data_results$P.Value))<1,message = "All p-values are NA. 
                Check to make sure you have replicates or >1 groups for statistical analysis."))
  
  
  validate(need(test_sel%in%data_results$test,message = "Incompatable test selection. Check group names of file."))
  
  #res = data_results%>%filter(test==paste0(group1,"/",group2))
  #if(test_sel%in%data_results$test) {
  res = data_results%>%filter(test==test_sel)
  #}else{res = data_results}
  
  res = res%>%filter(!is.na(res$P.Value))
  
  validate(need(mean(is.na(res$P.Value))<1,message = "All p-values for this test are NA. 
                Check to make sure you have replicates or >1 groups for statistical analysis."))
  
  usepadj=TRUE
  pvalname = "adj-pval"
  if(is.null(res$adj.P.Val)) {
    res$adj.P.Val = res$P.Value
    usepadj = FALSE 
    pvalname = "pval"
    print("no adjusted p-value found for volcano plot")
  }
  
  res$color="None"
  res$color[which((abs(res$logFC)>absFCcut)*(res$P.Value<pvalcut)==1)] = 
    paste0("pval","<",pvalcut," & abs(logfc)>",absFCcut)
  res$color[which((abs(res$logFC)<absFCcut)*(res$P.Value<pvalcut)==1)] =  
    paste0("pval","<",pvalcut, " & abs(logfc)<",absFCcut)
  res$color[which((abs(res$logFC)>absFCcut)*(res$adj.P.Val<fdrcut)==1)] = 
    paste0(pvalname,"<",fdrcut," & abs(logfc)>",absFCcut)
  res$color[which((abs(res$logFC)<absFCcut)*(res$adj.P.Val<fdrcut)==1)] = 
    paste0(pvalname,"<",fdrcut, " & abs(logfc)<",absFCcut)
  # if pvalcut is high and only have genes < fdrcut, fdrcut dominates, is this the best way, or should it
  # be intersection?
  
  # levels of color will be a subset of all_levels, but we want the color to match the all_levels
  tmplevels = levels(as.factor(res$color))
  all_levels = c("None",
                 paste0("pval","<",pvalcut, " & abs(logfc)<",absFCcut), # only exists if pval != pvalname 
                 paste0(pvalname,"<",fdrcut," & abs(logfc)<",absFCcut), 
                 paste0("pval","<",pvalcut, " & abs(logfc)>",absFCcut), # only exists if pval != pvalname 
                 paste0(pvalname,"<",fdrcut, " & abs(logfc)>",absFCcut))
  res$color = factor(
    res$color,
    levels = intersect(all_levels,
                       tmplevels
    ))
  tmplevels = levels(res$color)
  
  # add selected genes
  shapedata = data.frame()
  if(!is.null(sel_genes)) {
    tmpind = sapply(unlist(sel_genes),function(k) grep(k,res$unique_id,fixed=TRUE))
    tmpind = unique(unlist(tmpind))
    shapedata <- res[tmpind,]
  }
  
  
  p <- ggplot(res,aes(x=logFC,y=-log10(P.Value),color=color,text=unique_id))+
    geom_point(shape=19,fill="black")
  p <- p + scale_color_manual(values=
                                c("grey40","grey60","green3","grey70","red2")[match(tmplevels,all_levels)],
                              limits=levels(res$color),
                              name="Significance")
  
  if(nrow(shapedata)>0) {
    p <- p + geom_point(data=shapedata,fill="orange",shape=23,size=3,color="grey40") +
      #scale_size_manual(values = c(1,3))+
      guides(size=FALSE,shape=FALSE,fill=FALSE)
  }
  
  
  p <- p + theme_base() + theme(plot.margin = unit(c(2,2,2,2), "cm"))
  
  g <- plotly_build(p)
  
  
  #Match order of text to proper gene order
  newtext =  paste("Gene ID:",res$unique_id,"<br>",
                   "Comparison",res$test,"<br>",
                   "logFC",signif(res$logFC,3),"<br>",
                   "P.Value",signif(res$P.Value,3),"<br>",
                   "adj.P.Val",signif(res$adj.P.Val,3))
  
  for(ii in 1:length(g$x$data)) {
    tmpid = do.call(rbind,strsplit(g$x$data[[ii]]$text,"<br />"))[,4]
    g$x$data[[ii]]$text <- newtext[match(tmpid,res$unique_id)]
  }
  
  g
  
}

# switched from ggvis to plotly, this function is not currently used
rna_volcanoplot_ggvis <- function(data_results, geneids=NULL, 
                                  test_sel=NULL,absFCcut=0,fdrcut=0.05) {
  print(dim(data_results))
  res = data_results%>%filter(test==test_sel)
  
  usepadj=TRUE
  if(is.null(res$adj.P.Val)) {
    res$adj.P.Val = res$P.Value
    usepadj = FALSE 
  }
  
  res$color="None"
  res$color[which(res$adj.P.Val<fdrcut)] = paste0("adj-pval<",fdrcut)
  res$color[which(abs(res$logFC)>absFCcut)] = paste0("abs(logfc)>",absFCcut)
  res$color[which((abs(res$logFC)>absFCcut)*(res$adj.P.Val<.05)==1)] = paste0("adj-pval<",fdrcut," & abs(logfc)>",absFCcut)
  res$color = factor(res$color,
                     levels = c("None",
                                paste0("adj-pval<",fdrcut),
                                paste0("abs(logfc)>",absFCcut),
                                paste0("adj-pval<",fdrcut," & abs(logfc)>",absFCcut)))
  res$id = 1:nrow(res)
  
  
  all_values <- function(x){
    if(is.null(x)) return(NULL)
    row <- res[res$id==x$id,]
    if(usepadj) {
      show <- c("unique_id","test","logFC","P.Value","adj.P.Val")
      showname <- c("Gene ID","Comparison","logFC","raw p-value","BH FDR adjusted p-value")
    }else{
      show <- c("unique_id","test","logFC","P.Value")
      showname <- c("Gene ID","Comparison","logFC","p-value")
    }
    tmpout = paste0(showname,": ",format(row[,show],digits=3),collapse="<br />")
    tmpout
    #     paste0(tmpout,"<br/>",paste0("proteomics MGI Name: ",tmpgenename,collapse="<br/>"),
    #            "<br/>",paste0("proteomics Accession: ",tmpaccession,collapse="<br/>"))
  }
  
  res%>%ggvis(~logFC,~ -log10(P.Value),fill=~color,key := ~id)%>%
    layer_points()%>%add_axis("x",title="Log2(FC)")%>%
    add_axis("x",orient = "top",title=paste0("Comparison: ",unique(res$test)),
             ticks=0)%>%
    add_axis("y",title="-log10(p-value)")%>%
    add_tooltip(all_values, "hover")%>%add_legend("fill",title="Significance")
  
}


## ==================================================================================== ##
## Scatter plot of log2 fold changes
## ==================================================================================== ##  	
# 
# profvis::profvis(
#   rna_scatterplot(data_long,results,results_test_name="group1/group2",
#                   color_result_name="-log10(p-value)",
#                   group_sel=c('group1','group2'),
#                   sel_genes=c("Itpkb","ENSMUSG00000051977_Prdm9"))
# )


rna_scatterplot <- function(data_long, results, 
                            results_test_name = NULL,
                            color_result_name=NULL,
                            color_low="blue",
                            color_hi="orange",
                            geneids=NULL, group_sel=NULL,
                            valuename="log2cpm",
                            sel_genes=NULL) {
  group1 = group_sel[1]; group2 = group_sel[2]
  
  data_long$value = data_long[,valuename]
  
  pp = data_long%>%filter(group%in%group_sel)
  pp_sum = pp%>%group_by(unique_id,group)%>%summarise("Ave_value"=mean(value))
  
  pp_wide = pp_sum%>%spread(key = group,Ave_value)
  pp_wide$id = 1:nrow(pp_wide)
  
  colnames(pp_wide)[c(match(group1,colnames(pp_wide)),match(group2,colnames(pp_wide)))] = c("g1","g2")
  #pp_wide = pp_wide%>%mutate(diff = g1-g2,color=1*(g1>=g2)) # mutate is too slow
  pp_wide$diff = pp_wide$g1 - pp_wide$g2
  pp_wide$color = 1*(pp_wide$g1>=pp_wide$g2)
  
  results = results%>%filter(test==results_test_name)
  pp_wide = left_join(pp_wide,results)
  
  
  # Choose variable for colors
  colorlabels = c("logFC","p-value","adjusted p-value (q-value)",
                  "-log10(p-value)","-log10(q-value)",
                  "p-value < .1","q-value < .1")
  colorvars = c("logFC","P.Value","adj.P.Val","log10.P.Value","log10.adj.P.Val",
                "P.Value.1","adj.P.Val.1")
  pp_wide$log10.P.Value = -log10(pp_wide$P.Value)
  pp_wide$log10.adj.P.Val = -log10(pp_wide$adj.P.Val)
  pp_wide$P.Value.1 = pp_wide$P.Value
  pp_wide$P.Value.1[pp_wide$P.Value>.1] = .11
  pp_wide$adj.P.Val.1 = pp_wide$adj.P.Val
  pp_wide$adj.P.Val.1[pp_wide$adj.P.Val>.1] = .11
  
  len_nacolor = 0
  colorname = NULL
  if(color_result_name=="Sign of FC") color_result_name = NULL
  color_is_factor = TRUE
  if(!is.null(color_result_name)) {
    tmpcolorvar = colorvars[match(color_result_name,colorlabels)]
    tmpcolor = get(tmpcolorvar,pp_wide)
    len_nacolor = sum(is.na(tmpcolor))
    colorname = color_result_name
    color_is_factor = FALSE
    if(len_nacolor>0) {
      warning(paste0("Color factor has ",len_nacolor, "missing values, these genes will not appear on graph."))}
    pp_wide$color = tmpcolor
  }
  if(length(unique(pp_wide$color))<5) {
    color_is_factor = TRUE
    pp_wide$color = factor(pp_wide$color)
  }
  
  # add selected genes
  shapedata = data.frame()
  if(!is.null(sel_genes)) {
    tmpgenes = stringr::str_split(pp_wide$unique_id,"_",simplify = TRUE)
    tmpind = sapply(sel_genes,function(k) unique(which(tmpgenes==k,arr.ind=T)[,1]))
    tmpind = unique(unlist(tmpind))
    shapedata <- pp_wide[tmpind,]
  }
  
  #  pp_wide = pp_wide%>%filter(value>=valuecut[1],value<=valuecut[2])
  
  # all_values <- function(x){
  #   if(is.null(x)) return(NULL)
  #   row <- pp_wide[pp_wide$id==x$id,]
  #   show <- c("unique_id","g1","g2","diff")
  #   showname <- c("Gene ID",
  #                 paste0(group1,"_Ave",valuename),paste0(group2,"_Ave",valuename),
  #                 "difference")
  #   tmpout = paste0(showname,": ",format(row[,show],digits=3),collapse="<br />")
  #   tmpout
  # }
  # pp_wide%>%ggvis(~g1,~g2,fill=~factor(color),key := ~id)%>%
  #   layer_points()%>%add_axis("x",title=paste0(group1,"_Ave",valuename))%>%
  #   add_axis("x",orient = "top",title=paste0("Number of genes: ",nrow(pp_wide)),
  #            ticks=0)%>%
  #   add_axis("y",title=paste0(group2,"_Ave",valuename))%>%
  #   add_tooltip(all_values, "hover")%>%hide_legend("fill")
  
  # switch to ggplotly since ggvis was slow
  p <-   ggplot(pp_wide,aes(x=g1,y=g2,
                            color=color,
                            text=unique_id),fill=1)+geom_point(shape=19,size=1)+guides(fill=FALSE)
  p <- p + xlab(paste0(group1,"_Ave",valuename)) + ylab(paste0(group2,"_Ave",valuename))
  
  
  if(is.null(colorname)) {
    p <- p + guides(color=FALSE)
  }else {
    if(color_is_factor){
      mycolors = colorRampPalette(c(color_low,color_hi))(nlevels(pp_wide$color))
      p <- p + scale_color_manual(name=colorname,values = mycolors)
    }else{
      mycolors = colorRampPalette(c(color_low,color_hi))(nlevels(pp_wide$color))
      p <- p + scale_color_gradient(name=colorname,low=color_low,high=color_hi)
    }
  }
  
  
  if(nrow(shapedata)>0) {
    p <- p + geom_point(data=shapedata,fill=2,shape=23,size=4) +
      scale_size_manual(values = c(1,3))+
      scale_fill_manual(values = c("black","red"))+
      guides(size=FALSE,shape=FALSE,fill=FALSE)
  }
  
  p <- p + theme_base() + #ggtitle(paste0("Number of genes: ",nrow(pp_wide))) + 
    theme(plot.margin = unit(c(2,2,2,2), "cm"))
  
  g <- plotly_build(p)
  
  # just in case we don't have adj.p.val, don't error newtext
  if(is.null(pp_wide$adj.P.Val)) pp_wide$adj.P.Val = NA
  
  #Match order of text to proper gene order
  newtext =  paste("Gene ID:",pp_wide$unique_id,"<br>",
                   paste0(group1,"_Ave",valuename,":"),round(pp_wide$g1,3),"<br>",
                   paste0(group2,"_Ave",valuename,":"),round(pp_wide$g2,3),"<br>",
                   "Difference:",round(pp_wide$diff,3),"<br>",
                   "logFC:",round(pp_wide$logFC,3),"<br>",
                   "P.Value:",signif(pp_wide$P.Value,3),"<br>",
                   "adj.P.Val:",signif(pp_wide$adj.P.Val,3),"<br>"
  )
  
  
  
  for(ii in 1:length(g$x$data)) {
    if(!is.null(g$x$data[[ii]]$text)) {
      tmpid = stringr::str_split(g$x$data[[ii]]$text,"<br />",simplify=TRUE)[,4]
      g$x$data[[ii]]$text <- newtext[match(tmpid,pp_wide$unique_id)]
    }
  }
  
  g
  
}

