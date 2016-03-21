## ==================================================================================== ##
## Gene Data/ DOT PLOT variables
## ==================================================================================== ##

# loads list of genes on the server side
#observeEvent(input$upload_data,ignoreNULL = TRUE, {
observe({
  data_analyzed = analyzeCountDataReactive()
  tmpgeneids = data_analyzed$geneids
  data_analyzedgenes = as.character(unlist(tmpgeneids))
  tmpgroups = data_analyzed$group_names
  #data_analyzedgenes = c("a","b","c")
  #
  tmpynames = data_analyzed$data_long%>%select(-unique_id,-sampleid,-group)%>%colnames()
  updateSelectizeInput(session,'sel_gene',
                       choices= data_analyzedgenes,
                       server=TRUE)
  updateCheckboxGroupInput(session,'sel_group',
                           choices=tmpgroups, selected=tmpgroups)
  updateRadioButtons(session,'sel_gene_header',
                     choices=colnames(tmpgeneids))
  updateRadioButtons(session,"ytype",
                     choices=tmpynames)
#   
# 
#   output$geneurl <- renderText({
#     genename = input$sel_gene_rna # this needs to be fixed to have a gene_name to look up
#     # need to have annotation to get gene name always
#     out = ""
#     if(length(genename)>0) {out <- paste0(out,"<p>Click for Info:</p>")}
#     for(gene in genename) {
#       tmprow = (unique(as.numeric(unlist(apply(tmpgeneids,2,function(k) match(gene,k,nomatch =0))))))
#       tmpout = tmpgeneids[tmprow,]
#       
#       #if(gene%in%tmpgeneids$gene.id) gene = as.character(mousedata$gene.name[match(gene,mousedata$gene.id)])			
#       
#       out <- paste(out,"<p>",gene,":",
#                    paste0("<a href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=", 
#                           gene,"\"  target=\"_blank\">Genecards;</a>"),
#                    paste0("<a href=\"http://www.ncbi.nlm.nih.gov/gene/?term=(", 
#                           gene,"%5BGene+Name%5D)+AND+(Mus+musculus%5BOrganism%5D+)\"  target=\"_blank\">NCBI:gene</a></p>")
#       )
#     }
#     out
#   })
#   
#   
})

#Show dotplot	
output$dotplot <- renderPlotly({
  #browser()
  validate(need(length(input$sel_gene)>0,"Please select a gene."))
  validate(need(length(input$sel_group)>0,"Please select group(s)."))
  
  data_analyzed = analyzeCountDataReactive()
  data_long = data_analyzed$data_long
  geneids = data_analyzed$geneids
  dotplot_fun(data_long = data_long,geneids = geneids,
              genelabel=input$sel_gene_header,
              sel_group=input$sel_group,sel_gene=input$sel_gene,
              #log2y=input$log2cpm_checked,
              ytype=input$ytype)
}) #renderPlot



#Based on dotplot filters create data

DataDotplotReactive <- reactive({
  data_analyzed = analyzeCountDataReactive()
  
  subdat = dotplot_dat(data_long = data_analyzed$data_long,geneids = data_analyzed$geneids,
                       sel_group=input$sel_group,sel_gene=input$sel_gene,
                       ytype=input$ytype)
  return(subdat)
})

output$dat_dotplot <- DT::renderDataTable({
  tmpdat = DataDotplotReactive()
  DT::datatable(tmpdat)
}
)			

output$downloadSubsetData <- downloadHandler(
  filename = c('dotplot_data.csv'),
  content = function(file) {write.csv(DataDotplotReactive(), file, row.names=FALSE)}
)
