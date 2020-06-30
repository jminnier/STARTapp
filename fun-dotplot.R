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
## DOTPLOT FOR RNA-SEQ
##================================================================================##

#Get data for dotplot
#data_long needs to be in long format
dotplot_dat <- function(data_long,
                        geneids,
                        sel_group=NULL,
                        sel_gene=NULL,
                        #log2y=FALSE,
                        ytype="log2expr") {
  
  ll = length(sel_gene)
  if(ll==0) {return()} #NO GENES SELECTED
  sel_group = sort(sel_group)
  sel_group_labels = sel_group
  
  #condense ids down into a unique list based on input
  tmpids = geneids[unique(na.omit(c(apply(geneids,2,function(k) match(sel_gene,k))))),]
  
  #get subset of data
  subdat_rna = filter(data_long,unique_id%in%tmpids$unique_id,group%in%sel_group)
  geneids_tmp = filter(geneids,unique_id%in%tmpids$unique_id)
  
  subdat_rna = subdat_rna %>% rename(y=ytype)
  
  subdat_rna$group = factor(subdat_rna$group)
  subdat_rna = subdat_rna%>%filter(!is.na(y))
  subdat_rna = merge(geneids_tmp,subdat_rna)
  
  return(subdat_rna)
}

#Make dotplot with ggplot2
if(FALSE){
  genelabel="MGIsymbol"
  ytype="count"
  sel_group=group_names
  sel_gene="ENSMUSG00000026177_Slc11a1"
  
  dotplot_fun(data_long = data_long,geneids = geneids,
              genelabel="MGIsymbol",
              sel_group=group_names,sel_gene="ENSMUSG00000026177_Slc11a1",
              ytype="count")
}

dotplot_fun <- function(data_long,
                        geneids,
                        genelabel="unique_id",
                        sel_group=NULL,
                        sel_gene="Gnai3",
                        #log2y=TRUE,
                        ytype="log2expr") {
  
  ll = length(sel_gene) 
  if((ll==0)||(length(sel_group)==0)) {return(NULL)}else{
    sel_group = sort(sel_group)
    subdat_all = dotplot_dat(data_long,geneids,sel_group,sel_gene,ytype)
    #print(subdat_all)
    
    subdat_all$labelgeneid = subdat_all[,match(genelabel,colnames(subdat_all))]

    p <- ggplot(subdat_all,aes(x=group,y=y,fill=group)) +geom_boxplot()
    p <- p + facet_grid(.~ labelgeneid,scales = "free_y")+
      geom_point(size=3,aes(text = paste("sampleid:", sampleid))) + 
      stat_summary(fun=mean,geom="point",shape=5,size=3,fill=1)
    p <- p + scale_fill_discrete(name="group",breaks=sel_group,
                                 labels=sel_group,
                                 guide=guide_legend(keyheight=4,keywidth=2))
    
    # p <- p + theme_base() + #base_family="mono") + 
    #   theme(	plot.title=element_text(face="bold",size = rel(2)),
    #          strip.text = element_text(face="bold",size=rel(2)),
    #          strip.background=element_rect(fill="lightgrey"),
    #          axis.title = element_text(size=rel(2),color="blue"),
    #          axis.text = element_text(size=rel(2)),
    #          legend.text=element_text(size=rel(2)),
    #          legend.title=element_text(size=rel(2)),
    #          legend.position="bottom",legend.direction="vertical"
    #   )+ guides(fill = guide_legend(nrow = 2))
    p <- p + theme_base() + ylab(" ") + xlab(" ")+theme(
      plot.margin = unit(c(1,1,1,1), "cm"),
      axis.text.x = element_text(angle = 45),
      legend.position="bottom")+theme(legend.position="none")
    
    
    #if count data, log scale y axis
    #This should work but is not scaling and removes y-axis labels, issue with scale_y_continuous?
    #    if(min(data_long[,ytype])>=0)    {
  #   p + scale_y_continuous(trans = log2_trans(),
  #                          breaks = trans_breaks("log2", function(x) 2^x),
  #                          labels = trans_format("log2", math_format(2^.x)))
  # }
  #   #
    
    #hack to check if counts are not logged first need to fix this better
    if(max(data_long[,ytype])>=500) {p <- p+scale_y_continuous(trans = log2_trans(),breaks=2^(0:100))}
      
    #print(p)
    g <- ggplotly(p)
    g %>% layout(yaxis = list(title=ytype))
  }
  
  #ggplot(subdat,aes(x=tissue,y=rpkm,color=tissue)) + geom_point()
}

