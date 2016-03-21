
## ==================================================================================== ##
## Scatter plot of log2 fold changes
## ==================================================================================== ##  	


#change rna_volcanoplot to have an input function that depends on type of variables and then a plotting function


rna_volcanoplot <- function(data_results, geneids=NULL, 
                            test_sel=NULL,absFCcut=0,fdrcut=0.05) {
  print(dim(data_results))
  
  #group1 = group_sel[1]; group2 = group_sel[2]
  
  #res = data_results%>%filter(test==paste0(group1,"/",group2))
  
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
  pp_sum = pp%>%unite(grp,unique_id,group,sep=":")%>%
    group_by(grp)%>%summarise("Ave_value"=mean(value))%>%
    separate(grp,into=c("unique_id","group"),sep = ":")
  
  pp_wide = pp_sum%>%spread(key = group,Ave_value)
  pp_wide$id = 1:nrow(pp_wide)
  
  colnames(pp_wide)[c(match(group1,colnames(pp_wide)),match(group2,colnames(pp_wide)))] = c("g1","g2")
  pp_wide = pp_wide%>%mutate(diff = g1-g2,color=1*(g1>=g2))
  #  pp_wide = pp_wide%>%filter(value>=valuecut[1],value<=valuecut[2])
  
  all_values <- function(x){
    if(is.null(x)) return(NULL)
    row <- pp_wide[pp_wide$id==x$id,]
    show <- c("unique_id","g1","g2","diff")
    showname <- c("Gene ID",
                  paste0(group1,"_Ave",valuename),paste0(group2,"_Ave",valuename),
                  "difference")
    tmpout = paste0(showname,": ",format(row[,show],digits=3),collapse="<br />")
    tmpout
  }
  pp_wide%>%ggvis(~g1,~g2,fill=~factor(color),key := ~id)%>%
    layer_points()%>%add_axis("x",title=paste0(group1,"_Ave",valuename))%>%
    add_axis("x",orient = "top",title=paste0("Number of genes: ",nrow(pp_wide)),
             ticks=0)%>%
    add_axis("y",title=paste0(group2,"_Ave",valuename))%>%
    add_tooltip(all_values, "hover")%>%hide_legend("fill")
  
}

