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


observe({
  # Observe any changes in the input data, then update choices for column names
  # in the input data tab.
  # 
  # Check if example selected, or if not then ask to upload a file.
  validate(
    need((input$data_file_type=="examplecounts")|((!is.null(input$rdatafile))|(!is.null(input$datafile))), 
         message = "Please select a file")
  )
  inFile <- input$datafile
  if(!is.null(inFile)) {
    # update options for various analyzed data columns
    if(input$inputdat_type=="analyzed") {
      print("updating analyzed data choices")
      seqdata <- inputDataReactive()$data
      tmpcols = colnames(seqdata)
      updateSelectInput(session,"c_geneid1",choices =tmpcols)
      updateSelectInput(session,"c_geneid2",choices =tmpcols)
      updateSelectInput(session,"c_expr1",choices =tmpcols)
      updateSelectInput(session,"c_expr2",choices =tmpcols)
      updateSelectInput(session,"c_fc1",choices =tmpcols)
      updateSelectInput(session,"c_fc2",choices =tmpcols)
      updateSelectInput(session,"c_pval1",choices =tmpcols)
      updateSelectInput(session,"c_pval2",choices =tmpcols)
      updateSelectInput(session,"c_qval1",choices =tmpcols)
      updateSelectInput(session,"c_qval2",choices =tmpcols)
    }
  }
  
})


inputDataReactive <- reactive({
  
  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.
  print("inputting data")
  # Check if example selected, or if not then ask to upload a file.
  validate(
    need((input$data_file_type=="examplecounts")|((!is.null(input$rdatafile))|(!is.null(input$datafile))), 
         message = "Please select a file")
  )
  inFile <- input$datafile
  inRFile <- input$rdatafile
 # browser()
  
  if(input$data_file_type=="examplecounts") {
    # upload example data
    seqdata <- read_csv("data/mousecounts_example.csv")
    print("uploaded mousecounts data")
    return(list('data'=seqdata))
  }else if(input$data_file_type=="previousrdata"){
    if (!is.null(inRFile)) {
      load(inRFile$datapath,envir=environment())
      return(list("data"=data_results_table)) # this is so something shows in data upload window
    }else{return(NULL)}
  }else { # if uploading data
    if (!is.null(inFile)) {
      seqdata <- read_csv(inFile$datapath)
      print('uploaded seqdata')
      if(ncol(seqdata)==1) { # if file appears not to work as csv try tsv
        seqdata <- read_tsv(inFile$datapath)
        print('changed to tsv, uploaded seqdata')
      }
      validate(need(ncol(seqdata)>1,
                    message="File appears to be one column. Check that it is a comma-separated (.csv) file."))
      
      # Check for numeric columns
      not_numeric <- function(input) {
        if(sum(unlist(lapply(input,function(k) class(k)%in%c("numeric","integer"))))==0) {
          "Your data does not appear to be formatted correctly (no numeric columns). 
                        Please check your input file."
        } else if (input == "") { 
          FALSE
        } else {
          NULL
        }
      }
      
      validate(not_numeric(seqdata))
      
      
      
      return(list('data'=seqdata))}else{return(NULL)}
  }
})

# check if a file has been uploaded and create output variable to report this
output$fileUploaded <- reactive({
  return(!is.null(inputDataReactive()))
})
outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)

# after the data is uploaded or example data is selected, analyze the data
analyzeDataReactive <- 
  eventReactive(input$upload_data,
                ignoreNULL = FALSE, {
                  withProgress(message = "Analyzing RNA-seq data, please wait",{
                    
                    print("analysisCountDataReactive")
                    ## ==================================================================================== ##
                    ## Example data
                    ## ==================================================================================== ##
                    if(input$data_file_type=="examplecounts") {
                      # load('data/mousecounts_example_analysis_results.RData')
                      # load('data/mousecounts_example_analyzed.RData') #example_data_results for data_results_table
                      load_existing_rdata('data/mousecounts_example.RData')
                    
                    ## ==================================================================================== ##
                    ## Upload previously downloaded RData
                    ## ==================================================================================== ##
                    
                    }else if(input$data_file_type=="previousrdata"){
                      inRfile <- input$rdatafile
                      load_existing_rdata(inRfile$datapath)
                    }else{
                    
                    ## ==================================================================================== ##
                    ## Else, continue on with uploading csv data
                    ## ==================================================================================== ##
                      alldata <- inputDataReactive()$data
                      # remove empty columns
                      alldata = alldata %>% remove_empty(which=c("rows","cols"))
                      
                      ## ==================================================================================== ##
                      ## Count/expression data
                      ## ==================================================================================== ##
                      if(input$inputdat_type=="analyzed") {
                        tmpgenecols = seq(match(input$c_geneid1,colnames(alldata)),match(input$c_geneid2,colnames(alldata)))
                        tmpexprcols = seq(match(input$c_expr1,colnames(alldata)),match(input$c_expr2,colnames(alldata)))
                        tmpfccols = seq(match(input$c_fc1,colnames(alldata)),match(input$c_fc2,colnames(alldata)))
                        tmppvalcols = seq(match(input$c_pval1,colnames(alldata)),match(input$c_pval2,colnames(alldata)))
                        tmpqvalcols = seq(match(input$c_qval1,colnames(alldata)),match(input$c_qval2,colnames(alldata)))
                        
                        validate(need((length(tmpfccols)==length(tmppvalcols))&(length(tmpfccols)==length(tmpqvalcols)),message =
                                        "Number of fold change columns needs to be same number as 
                                      p-value and q-value columns (and in the same order)."))
                        
                        tmpres <- load_analyzed_data(alldata, 
                                                     tmpgenecols, tmpexprcols, tmpfccols, tmppvalcols, tmpqvalcols,
                                                     isfclogged = input$isfclogged)
                        return(list("countdata"=tmpres$countdata,
                                    "group_names"=tmpres$group_names,
                                    "sampledata"=tmpres$sampledata,
                                    "results"=tmpres$results,
                                    "data_long"=tmpres$data_long, 
                                    "geneids"=tmpres$geneids, 
                                    "data_results_table"=tmpres$data_results_table))
                        
                      }else if(input$inputdat_type=="expression_only") {
                        tmpres <- analyze_expression_data(alldata, analysis_method = input$analysis_method)
                        return(list("countdata"=tmpres$countdata,
                                    "group_names"=tmpres$group_names,
                                    "sampledata"=tmpres$sampledata,
                                    "results"=tmpres$results,
                                    "data_long"=tmpres$data_long, 
                                    "geneids"=tmpres$geneids, 
                                    "data_results_table"=tmpres$data_results_table))
                      }     
                      
                    }

                    return(list("countdata"=countdata,
                                "group_names"=group_names,
                                "sampledata"=sampledata,
                                "results"=results,
                                "data_long"=data_long, 
                                "geneids"=geneids, 
                                "data_results_table"=data_results_table))
                  })
                })


output$countdataDT <- renderDataTable({
  tmp <- inputDataReactive()
  if(!is.null(tmp)) tmp$data
})

observeEvent(input$upload_data, ({
  updateCollapse(session,id =  "input_collapse_panel", open="analysis_panel",
                 style = list("analysis_panel" = "success",
                              "data_panel"="primary"))
}))

observeEvent(inputDataReactive(),({
  updateCollapse(session,id =  "input_collapse_panel", open="data_panel",
                 style = list("analysis_panel" = "default",
                              "data_panel"="success"))
})
)

output$analysisoutput <- renderDataTable({
  print("output$analysisoutput")
  getresults <- analyzeDataReactive()
  res = getresults$results
  res[,sapply(res,is.numeric)] <- signif(res[,sapply(res,is.numeric)],3)
  datatable(res)
})


output$downloadResults_CSV <- downloadHandler(filename = paste0("START_results_",Sys.Date(),".csv"),
                                              content = function(file) {
                                                write.csv(analyzeDataReactive()$data_results_table, file, row.names=FALSE)})

output$downloadResults_RData <- downloadHandler(filename= paste0("START_results_",Sys.Date(),".RData"),
                                                content=function(file){
                                                  tmp = analyzeDataReactive()
                                                  
                                                  group_names = tmp$group_names
                                                  sampledata = tmp$sampledata
                                                  results = tmp$results
                                                  data_long = tmp$data_long
                                                  geneids = tmp$geneids
                                                  data_results_table = tmp$data_results_table
                                                  
                                                  save(group_names,sampledata,results,
                                                       data_long,geneids,
                                                       data_results_table,file=file)
                                                })


output$example_counts_file <- downloadHandler(filename="examplecounts_short.csv",
                                              content=function(file){
                                                file.copy("data/examplecounts_short.csv",file)
                                              })

output$example_analysis_file <- downloadHandler(filename="exampleanalysisres_short.csv",
                                                content=function(file){
                                                  file.copy("data/exampleanalysisres_short.csv",file)
                                                })



output$instructionspdf <- downloadHandler(filename="Instructions.pdf",
                                          content=function(file){
                                            file.copy("instructions/Instructions.pdf",file)
                                          })




