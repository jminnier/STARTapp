## ==================================================================================== ##
# START Shiny App for analysis and visualization of transcriptome data.
# Copyright (C) 2014-2016  Jessica Minnier
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
## 
## 
#update list of groups
observe({
  data_analyzed = analyzeCountDataReactive()
  tmpgroups = data_analyzed$group_names
  updateSelectizeInput(session,'sampleres_groups',
                       choices=tmpgroups, selected=tmpgroups)
  
  
  
})




observe({
  
  output$gene_pheatmap <- renderPlot({
    data_analyzed = analyzeCountDataReactive()
    data_results = data_analyzed$results
    geneids = data_analyzed$geneids
    sampledata = data_analyzed$sampledata
    
    tmpgroups = input$sampleres_groups
    tmpkeep = which(sampledata$group%in%tmpgroups)
    
    if(length(tmpkeep)>0) {
      tmpdat = data_analyzed$expr_data[,tmpkeep]
      gene_pheatmap(as.matrix(tmpdat),sampleid=sampledata$sampleid[tmpkeep],annotation_row = sampledata[tmpkeep,"group",drop=FALSE])
    }
  })
  
  output$pca_plot <- renderPlot({
    
    data_analyzed = analyzeCountDataReactive()
    data_results = data_analyzed$results
    geneids = data_analyzed$geneids
    sampledata = data_analyzed$sampledata
    
    
    tmpgroups = input$sampleres_groups
    tmpkeep = which(sampledata$group%in%tmpgroups)
    
    if(length(tmpkeep)>0) {
      tmpdat = data_analyzed$expr_data[,tmpkeep]
      gene_pcaplot(tmpdat,sampleid=sampledata$sampleid[tmpkeep], groupdat= sampledata[tmpkeep,"group",drop=FALSE],colorfactor="group")
    }
    
    
  })
  
  
  
  
  # rna_volcanoplot(data_results = data_results,
  #                 group_sel = input$analysisres_groups,
  #                 absFCcut = input$analysisres_fold_change_cut,
  #                 fdrcut = input$analysisres_fdrcut)%>%
  #   bind_shiny("volcanoplot_2groups_ggvis","volcanoplot_2groups_ggvisUI")
  
  
})