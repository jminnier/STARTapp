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
## 
## 
## #update list of groups
observe({
  
  print("server-analysisres-update")
  
  data_analyzed = analyzeCountDataReactive()
  tmpdat = data_analyzed$results
  tmpgroups = data_analyzed$group_names
  tmptests = unique(as.character(tmpdat$test))
  tmpdatlong = data_analyzed$data_long
  tmpynames = tmpdatlong%>%select(-unique_id,-sampleid,-group)%>%colnames()
  
  updateSelectizeInput(session,'analysisres_test',
                       choices=tmptests, selected=tmptests[1])
  
  updateSelectizeInput(session,'analysisres_groups',
                       choices=tmpgroups)
  
  updateRadioButtons(session,'scattervaluename',
                     choices=sort(tmpynames,decreasing = TRUE))
  
  
})






observe({
  
  print("drawing volcano plot")
  
  data_analyzed = analyzeCountDataReactive()
  data_results = data_analyzed$results
  geneids = data_analyzed$geneids
  
  output$volcanoplot <- renderPlotly({
    validate(need(input$analysisres_test!="","Select a test."))
    validate(need(data_results%>%filter(test==input$analysisres_test)%>%nrow()>0,"Test not found."))
    
    withProgress(message = "Drawing volcano plot, please wait",
                 {
                   # rna_volcanoplot(data_results = data_results,
                   #                 test_sel = input$analysisres_test,
                   #                 absFCcut = input$analysisres_fold_change_cut,
                   #                 fdrcut = input$analysisres_fdrcut)%>%
                   #   bind_shiny("volcanoplot_2groups_ggvis","volcanoplot_2groups_ggvisUI")
                   
                   rna_volcanoplot(data_results = data_results,
                                   test_sel = input$analysisres_test,
                                   absFCcut = input$analysisres_fold_change_cut,
                                   fdrcut = input$analysisres_fdrcut)
                   
                 })#end withProgress
    
  }) 
})



observe({
  
  print("drawing scatterplot")
  
  #if(length(input$analysisres_groups)==2) {
  data_analyzed = analyzeCountDataReactive()
  data_long = data_analyzed$data_long
  geneids = data_analyzed$geneids
  
  
  
  # rna_scatterplot(data_long = data_long,
  #                 group_sel = input$analysisres_groups,
  #                 valuename=input$scattervaluename)%>%
  #   bind_shiny("scatterplot_fc_2groups_ggvis","scatterplot_fc_2groups_ggvisUI")
  output$scatterplot <- renderPlotly({ 
    validate(need(length(input$analysisres_groups)==2,"Please select two groups."))
    withProgress(message = "Drawing scatterplot, please wait",{
      rna_scatterplot(data_long = data_long,
                      group_sel = input$analysisres_groups,
                      valuename=input$scattervaluename)
    })#end withProgress
  })
  
  
  
  #}
})






