#update list of groups
observe({
  data_analyzed = analyzeCountDataReactive()
  tmpdat = data_analyzed$results
  tmpgroups = data_analyzed$group_names
  tmptests = unique(as.character(tmpdat$test))
  tmpdatlong = data_analyzed$data_long
  tmpynames = tmpdatlong%>%select(-unique_id,-sampleid,-group)%>%colnames()
  
  updateSelectizeInput(session,'analysisres_test',
                       choices=tmptests)
  
  updateSelectizeInput(session,'analysisres_groups',
                       choices=tmpgroups)
  
  updateRadioButtons(session,'scattervaluename',
                     choices=tmpynames)
  
  
})






observe({
  if(input$analysisres_test!="") {
    
    data_analyzed = analyzeCountDataReactive()
    data_results = data_analyzed$results
    geneids = data_analyzed$geneids
    
    withProgress(message = "Drawing volcano plot, please wait",
                 {
                   rna_volcanoplot(data_results = data_results,
                                   test_sel = input$analysisres_test,
                                   absFCcut = input$analysisres_fold_change_cut,
                                   fdrcut = input$analysisres_fdrcut)%>%
                     bind_shiny("volcanoplot_2groups_ggvis","volcanoplot_2groups_ggvisUI")
                   
                 })#end withProgress
  }
})



observe({
  if(length(input$analysisres_groups)==2) {
    data_analyzed = analyzeCountDataReactive()
    data_long = data_analyzed$data_long
    geneids = data_analyzed$geneids
    
    withProgress(message = "Drawing scatterplot, please wait",
                 {
                   rna_scatterplot(data_long = data_long,
                                   group_sel = input$analysisres_groups,
                                   valuename=input$scattervaluename)%>%
                     bind_shiny("scatterplot_fc_2groups_ggvis","scatterplot_fc_2groups_ggvisUI")
                   
                 })#end withProgress
  }
})