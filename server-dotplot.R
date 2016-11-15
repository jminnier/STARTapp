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
## Gene Data/ DOT PLOT variables
## ==================================================================================== ##

# loads list of genes on the server side
#observeEvent(input$upload_data,ignoreNULL = TRUE, {
observe({
  
  print("server-dotplot-update")
  
  data_analyzed = analyzeDataReactive()
  tmpgeneids = data_analyzed$geneids
  #data_analyzedgenes = as.character(unlist(tmpgeneids))
  tmpgroups = data_analyzed$group_names
  #data_analyzedgenes = c("a","b","c")
  #
  tmpynames = data_analyzed$data_long%>%select(-unique_id,-sampleid,-group)%>%colnames()
  updateSelectizeInput(session,'sel_gene',
                       choices= tmpgeneids,server=TRUE)
  updateCheckboxGroupInput(session,'sel_group',
                           choices=tmpgroups, selected=tmpgroups)
  updateRadioButtons(session,'sel_gene_header',
                     choices=colnames(tmpgeneids))
  updateRadioButtons(session,"ytype",
                     choices=sort(tmpynames,decreasing = TRUE))
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
  
  print("drawing dotplot")
  
  validate(need(length(input$sel_gene)>0,"Please select a gene."))
  validate(need(length(input$sel_group)>0,"Please select group(s)."))
  
  data_analyzed = analyzeDataReactive()
  data_long = data_analyzed$data_long
  geneids = data_analyzed$geneids
  if (names(dev.cur()) != "null device") dev.off()
  pdf(NULL)
  p=dotplot_fun(data_long = data_long,geneids = geneids,
              genelabel=input$sel_gene_header,
              sel_group=input$sel_group,sel_gene=input$sel_gene,
              #log2y=input$log2cpm_checked,
              ytype=input$ytype)
}) #renderPlot


#Based on dotplot filters create data

DataDotplotReactive <- reactive({
  print("DataDotplotReactive")
  data_analyzed = analyzeDataReactive()
  
  subdat = dotplot_dat(data_long = data_analyzed$data_long,geneids = data_analyzed$geneids,
                       sel_group=input$sel_group,sel_gene=input$sel_gene,
                       ytype=input$ytype)
  return(subdat)
})

output$dat_dotplot <- renderDataTable({
  tmpdat = DataDotplotReactive()
  tmpdat[,sapply(tmpdat,is.numeric)] <- signif(tmpdat[,sapply(tmpdat,is.numeric)],3)
  datatable(tmpdat)
}
)			

output$downloadSubsetData <- downloadHandler(
  filename = c('dotplot_data.csv'),
  content = function(file) {write.csv(DataDotplotReactive(), file, row.names=FALSE)}
)
