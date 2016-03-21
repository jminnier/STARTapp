
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