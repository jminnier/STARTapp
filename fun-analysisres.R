
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
rna_volcanoplot <- function(data_results, geneids=NULL, 
                            test_sel=NULL,absFCcut=0,fdrcut=0.05) {
  print(dim(data_results))
  #group1 = group_sel[1]; group2 = group_sel[2]
  
  validate(need(mean(is.na(data_results$P.Value))<1,message = "All p-values are NA. 
                Check to make sure you have replicates or >1 groups for statistical analysis."))
  
  
  #res = data_results%>%filter(test==paste0(group1,"/",group2))
  if(test_sel%in%data_results$test) {
  res = data_results%>%filter(test==test_sel)
  }else{res = data_results}
  
  usepadj=TRUE
  pvalname = "adj-pval"
  if(is.null(res$adj.P.Val)) {
    res$adj.P.Val = res$P.Value
    usepadj = FALSE 
    pvalname = "pval"
  }
  
  
  
  res$color="None"
  res$color[which(res$adj.P.Val<fdrcut)] = paste0(pvalname,"<",fdrcut, " & abs(logfc)<",absFCcut)
  res$color[which(abs(res$logFC)>absFCcut)] = paste0(pvalname,">",fdrcut, " & abs(logfc)<",absFCcut)
  res$color[which((abs(res$logFC)>absFCcut)*(res$adj.P.Val<.05)==1)] = paste0(pvalname,"<",fdrcut," & abs(logfc)>",absFCcut)
  res$color = factor(res$color,levels = c("None",paste0("adj-pval<",fdrcut, " & abs(logfc)<",absFCcut),
                                          paste0(pvalname,">",fdrcut, " & abs(logfc)<",absFCcut),
                                          paste0(pvalname,"<",fdrcut," & abs(logfc)>",absFCcut)
                                          ))

  
  
  p <- ggplot(res,aes(x=logFC,y=-log10(P.Value),color=color,text=unique_id))+geom_point()+
    scale_color_manual(values=c("black","red","orange","green"),name="Significance")
  p <- p + theme_base() + theme(plot.margin = unit(c(2,2,2,2), "cm"))
  
  g <- plotly_build(p)

  
  #Match order of text to proper gene order
  newtext =  paste("Gene ID:",res$unique_id,"<br>",
                   "Comparison",res$test,"<br>",
                   "logFC",signif(res$logFC,3),"<br>",
                   "P.Value",signif(res$P.Value,3),"<br>",
                   "adj.P.Val",signif(res$adj.P.Val,3))

  for(ii in 1:length(g$data)) {
    tmpid = do.call(rbind,strsplit(g$data[[ii]]$text,"<br>"))[,4]
    g$data[[ii]]$text <- newtext[match(tmpid,res$unique_id)]
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
  res$color = factor(res$color,levels = c("None",paste0("adj-pval<",fdrcut),paste0("abs(logfc)>",absFCcut),paste0("adj-pval<",fdrcut," & abs(logfc)>",absFCcut)))
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

rna_scatterplot <- function(data_long, geneids=NULL, group_sel=NULL,
                            valuename="log2cpm") {
  group1 = group_sel[1]; group2 = group_sel[2]
  
  data_long$value = data_long[,valuename]
  
  pp = data_long%>%filter(group%in%group_sel)
  pp_sum = pp%>%group_by(unique_id,group)%>%summarise("Ave_value"=mean(value))
  
  pp_wide = pp_sum%>%spread(key = group,Ave_value)
  pp_wide$id = 1:nrow(pp_wide)
  
  colnames(pp_wide)[c(match(group1,colnames(pp_wide)),match(group2,colnames(pp_wide)))] = c("g1","g2")
  pp_wide = pp_wide%>%mutate(diff = g1-g2,color=1*(g1>=g2))
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
  p <- ggplot(pp_wide,aes(x=g1,y=g2,
                          color=factor(color),text=unique_id))+geom_point()
  p <- p + xlab(paste0(group1,"_Ave",valuename)) + ylab(paste0(group2,"_Ave",valuename))+
    scale_color_manual(values=c("darkred","darkorange"))
  p <- p + theme_base() + #ggtitle(paste0("Number of genes: ",nrow(pp_wide))) + 
    theme(legend.position="none",plot.margin = unit(c(2,2,2,2), "cm"))
  
  g <- plotly_build(p)
  
  #Match order of text to proper gene order
  newtext =  paste("Gene ID:",pp_wide$unique_id,"<br>",
                   paste0(group1,"_Ave",valuename,":"),round(pp_wide$g1,3),"<br>",
                   paste0(group2,"_Ave",valuename,":"),round(pp_wide$g2,3),"<br>",
                   "Difference:",round(pp_wide$diff,3))
  
  
  tmpid = do.call(rbind,strsplit(g$data[[1]]$text,"<br>"))[,4]
  g$data[[1]]$text <- newtext[match(tmpid,pp_wide$unique_id)]
  
  tmpid = do.call(rbind,strsplit(g$data[[2]]$text,"<br>"))[,4]
  g$data[[2]]$text <- newtext[match(tmpid,pp_wide$unique_id)]
  
  g
  
}

