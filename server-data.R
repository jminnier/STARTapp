## ================================================================================== ##
## ALL DATA variables
## ================================================================================== ##		

#why is this reactive?
AllRNAdatReactive <- reactive({
  outdat=mousedata[,c("gene.name","gene.id","sample","count","log2cpm","cpm")]
  colnames(outdat) = c("gene.name","gene.id","sample.id","count","log2cpm.edgeR.adjusted","cpm.bowtie.raw")
  return(outdat)		
})

output$downloadData <- downloadHandler(filename = c('all_data.csv'),
                                       content = function(file) {write.csv(AllRNAdatReactive(), file, row.names=FALSE)})



output$outdat <- DT::renderDataTable({
  tmpdat = AllRNAdatReactive()
  DT::datatable(tmpdat)
})
