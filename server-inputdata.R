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
  
  validate(
    need((input$use_example_file!="upload")|((!is.null(input$rdatafile))|(!is.null(input$datafile))), "Please select a file")
  )
  inFile <- input$datafile
  if(!is.null(inFile)) {
    
    seqdata <- read.csv(inFile$datapath, header=TRUE, sep=input$sep, stringsAsFactors = FALSE)
    print('uploaded seqdata')
    
    if(ncol(seqdata)==1) {
      updateRadioButtons(session,'sep',selected="\t")
      seqdata <- read.csv(inFile$datapath, header=TRUE, sep=input$sep, stringsAsFactors = FALSE)
      print('changed to tsv, uploaded seqdata')
    }
    
    if(input$inputdat_type=="analyzed") {
      tmpcols = colnames(seqdata)
      updateSelectInput(session,"c_geneid1",choices =tmpcols)
      updateSelectInput(session,"c_geneid2",choices =tmpcols)
      updateSelectInput(session,"c_expr1",choices =tmpcols)
      updateSelectInput(session,"c_expr2",choices =tmpcols)
      updateSelectInput(session,"c_fc1",choices =tmpcols)
      updateSelectInput(session,"c_fc2",choices =tmpcols)
      updateSelectInput(session,"c_pval1",choices =tmpcols)
      updateSelectInput(session,"c_pval2",choices =tmpcols)
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
  
  validate(
    need((input$use_example_file=="examplecounts")|(!is.null(input$rdatafile))|(!is.null(input$datafile)), "Please select a file")
  )
  inFile <- input$datafile
  
  
  if (is.null(inFile)) {
    if(input$use_example_file=="examplecounts") {
      seqdata <- read.csv("data/mousecounts_example.csv",stringsAsFactors = FALSE)
      
      print("uploaded mousecounts data")
      
    }else if(input$use_example_file=="previousrdata"){
      inRfile <- input$rdatafile
      load(inRfile$datapath,envir=environment())
      return(list("data"=data_results_table)) # this is so something shows in data upload window
    }else{return(NULL)}
  }else { # if uploading data
    seqdata <- read.csv(inFile$datapath, header=TRUE, sep=input$sep, stringsAsFactors = FALSE)
  }
  
  return(list('data'=seqdata))
})


output$fileUploaded <- reactive({
  return(!is.null(inputDataReactive()))
})
outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)

analyzeCountDataReactive <- 
  eventReactive(input$upload_data,
                ignoreNULL = FALSE, {
                  withProgress(message = "Analyzing RNA-seq data, please wait",{
                    
                    print("analysisCountDataReactive")

                    
                    #if an example just return previously analyzed results
                    if(input$use_example_file=="examplecounts") {
                      
                      load('data/mousecounts_example_analysis_results.RData')
                      load('data/mousecounts_example_analyzed.RData') #example_data_results
                      return(list('group_names'=group_names,'sampledata'=sampledata,
                                  "results"=results,"data_long"=data_long, "geneids"=geneids,
                                  "expr_data"=expr_data,"data_results_table"=example_data_results))
                    }
                    
                    #if uploading own data:
                    
                    if(input$use_example_file=="previousrdata"){
                      inRfile <- input$rdatafile
                      load(inRfile$datapath,envir=environment())
                      
                      return(list('group_names'=group_names,'sampledata'=sampledata,
                                  "results"=results,"data_long"=data_long, 
                                  "geneids"=geneids, "expr_data"=expr_data,
                                  "data_results_table"=data_results_table))
                    }
                    
                    alldata <- inputDataReactive()$data
                    
                    not_numeric <- function(input) {
                      if(sum(unlist(lapply(input,function(k) class(k)%in%c("numeric","integer"))))==0) {
                        "Your data does not appear to be formatted correctly (no numeric columns). Please check your input file."
                      } else if (input == "") { 
                        FALSE
                      } else {
                        NULL
                      }
                    }
                    
                    
                    validate(
                      not_numeric(alldata)
                    )
                    
                    if(input$inputdat_type=="counts") {
                      numgeneids <- input$numgeneids
                      #catch incorrect gene id error, only works if geneids are 1:numgeneids and no other columns are characters
                      numgeneids = max(numgeneids,max(which(sapply(alldata,class)=="character")))
                      
                      tmpgenecols = 1:numgeneids
                      tmpexprcols = setdiff(1:ncol(alldata),tmpgenecols)
                      tmpfccols = NULL
                      tmppvalcols = NULL
                    }
                    if(input$inputdat_type=="analyzed") {
                      tmpgenecols = seq(match(input$c_geneid1,colnames(alldata)),match(input$c_geneid2,colnames(alldata)))
                      tmpexprcols = seq(match(input$c_expr1,colnames(alldata)),match(input$c_expr2,colnames(alldata)))
                      tmpfccols = seq(match(input$c_fc1,colnames(alldata)),match(input$c_fc2,colnames(alldata)))
                      tmppvalcols = seq(match(input$c_pval1,colnames(alldata)),match(input$c_pval2,colnames(alldata)))
                      
                      validate(need(length(tmpfccols)==length(tmppvalcols),message =
                                      "Number of fold change columns needs to be same number as p-value columns (and in the same order)."))
                    }
                    
                    
                    #split expression names into groups
                    sampleid <- colnames(alldata[,tmpexprcols])
                    tmpnames <- do.call(rbind,strsplit(sampleid,"_",fixed=TRUE))
                    group_names <- unique(tmpnames[,1])
                    group <- tmpnames[,1]
                    rep_id <- tmpnames[,2]
                    sampledata = data.frame(sampleid,group,rep_id)
                    
                    countdata <- alldata[,tmpexprcols,drop=FALSE]
                    geneids <- alldata[,tmpgenecols,drop=FALSE]
                    tmpkeep = which(apply(is.na(geneids),1,mean)<1) #remove rows with no gene identifiers
                    print(paste0("Num genes kept after removing empty geneids: ",length(tmpkeep)," of ", nrow(geneids)))
                    
                    
                    geneids = geneids[tmpkeep,,drop=FALSE]
                    countdata = countdata[tmpkeep,,drop=FALSE]
                    
                    geneids = geneids%>%unite_("unique_id",colnames(geneids),remove = FALSE)
                    
                    #if geneids not unique
                    if(length(unique(geneids$unique_id))<nrow(geneids)) {
                      geneids = geneids%>%group_by(unique_id)%>%
                        mutate(rn=row_number(unique_id),new=
                                 ifelse(rn==1,unique_id,paste(unique_id,rn,sep="_")))%>%
                        ungroup()%>%mutate(unique_id=new)%>%select(-rn,-new)
                      
                    }
                    
                    #add in catch for if length(tmpkeep) = 0
                    
                    #add filter for max # counts
                    
                    #handle NAs 
                    countdata[which(is.na(countdata),arr.ind=T)] <- 0 #allow choice of this or removal
                    rownames(countdata) = geneids$unique_id
                    
                    if(input$inputdat_type=="analyzed") {
                      
                      
                      expr_data <- alldata[tmpkeep,tmpexprcols,drop=FALSE]
                      rownames(expr_data) = geneids$unique_id
                      
                      tmpfc = alldata[tmpkeep,tmpfccols,drop=F]
                      if(input$isfclogged=="No (Log my data please)") {log2(tmpfc)}
                      
                      fcdata = cbind("unique_id"=geneids$unique_id,tmpfc)
                      pvaldata = cbind("unique_id"=geneids$unique_id,alldata[tmpkeep,tmppvalcols,drop=F])
                      
                      tmpnames = paste(colnames(fcdata),colnames(pvaldata),sep=":")[-1]
                      colnames(fcdata)[-1] = tmpnames
                      colnames(pvaldata)[-1] = tmpnames
                      
                      fcdatalong = fcdata%>%gather(key = "test",value = "logFC",-1)
                      pvaldatalong = pvaldata%>%gather(key = "test",value = "P.Value",-1)
                      tmpres = full_join(fcdatalong,pvaldatalong)
                      
                      tmpdat = cbind("unique_id"=geneids$unique_id,expr_data)
                      tmpdatlong = tmpdat%>%gather(key="sampleid",value="expr",-1)
                      data_long = left_join(tmpdatlong,sampledata%>%select(sampleid,group))
                      
                      tmpres$test = as.character(tmpres$test)
                      
                      return(list('group_names'=group_names,'sampledata'=sampledata,
                                  "results"=tmpres,"data_long"=data_long, 
                                  "geneids"=geneids,"expr_data"=expr_data,
                                  "data_results_table"=alldata))
                    }else if(input$inputdat_type=="counts") {
                      
                      #analyze data
                      
                      not_counts <- function(input) {
                        
                        remainder = sum(apply(input,2,function(k) sum(k%%1,na.rm=T)),na.rm=T)
                        if (remainder !=0) {
                          "Your data appears to not be counts, please double check your data"
                        } else if (input == "") { 
                          FALSE
                        } else {
                          NULL
                        }
                      }
                      
                      
                      validate(
                        not_counts(countdata)
                      )
                      
                      print("analyze data: counts")
                      
                      design <- model.matrix(~0+sampledata$group) # allow selection of reference group
                      colnames(design) = levels(as.factor(sampledata$group))
                      
                      #voom+limma
                      dge <- DGEList(counts=countdata) #TMM normalization first
                      dge <- calcNormFactors(dge)
                      log2cpm <- cpm(dge, prior.count=0.5, log=TRUE)
                      # v <- voom(dge,design,plot=FALSE)
                      v <- voom(dge,design,plot=FALSE,normalize.method = "cyclicloess")
                      # v <- voom(countdata,design,plot=TRUE,normalize="quantile") #use this to allow different normalization
                      
                      
                      
                      #fit <- lmFit(v,design)
                      #fit <- eBayes(fit)
                      
                      expr_data = v$E
                      
                      tmpgroup = sampledata$group
                      #contrasts(tmpgroup)
                      
                      lmobj_res = list()
                      for(ii in 1:length(group_names)) {
                        grp = relevel(tmpgroup,ref= group_names[ii] )
                        lm.obj = lm(t(expr_data) ~ grp)
                        beta.lm = t(lm.obj$coefficients)
                        pval.lm = t(lm.pval(lm.obj)$pval)
                        pval.adj.lm = apply(pval.lm,2,p.adjust,method="BH")
                        
                        colnames(beta.lm) = colnames(pval.lm) = colnames(pval.adj.lm) = 
                          gsub("grp","",colnames(beta.lm))
                        
                        tmpout = cbind(melt(beta.lm[,-1,drop=FALSE]),
                                       melt(pval.lm[,-1,drop=FALSE])$value,
                                       melt(pval.adj.lm[,-1,drop=FALSE])$value)
                        colnames(tmpout) = c("unique_id","numer_group","logFC","P.Value","adj.P.Val")
                        tmpout$denom_group = group_names[ii]
                        tmpout$test = with(tmpout, paste(numer_group,denom_group,sep="/"))
                        tmpout = tmpout[,c("unique_id","test","denom_group","numer_group",
                                           "logFC","P.Value","adj.P.Val")]
                        
                        
                        lmobj_res[[ii]] = tmpout
                      }
                      lmobj_res = do.call(rbind,lmobj_res)
                      
                      pvals = lmobj_res%>%select(unique_id,test,adj.P.Val)%>%spread(test,adj.P.Val)
                      logfcs = lmobj_res%>%select(unique_id,test,logFC)%>%spread(test,logFC)
                      
                      colnames(pvals)[-1] = paste0("padj_",colnames(pvals)[-1])
                      colnames(logfcs)[-1] = paste0("logFC_",colnames(logfcs)[-1])
                      tmpdat = cbind(geneids,log2cpm)
                      tmpdat = left_join(tmpdat,logfcs)
                      tmpdat = left_join(tmpdat,pvals)
                      
                      data_results_table = tmpdat%>%select(-unique_id) #save this into csv
                      
                      tmpexprdata = data.frame("unique_id" =geneids$unique_id,expr_data)
                      tmpcountdata = data.frame("unique_id"=geneids$unique_id,countdata)
                      
                      tmplog2cpm = data.frame("unique_id"=geneids$unique_id,log2cpm)
                      log2cpm_long = melt(tmplog2cpm,variable.name = "sampleid",value.name="log2cpm")
                      
                      countdata_long = melt(tmpcountdata,variable.name = "sampleid",value.name="count")
                      #countdata_long$log2count = log2(countdata_long$count+.25)
                      
                      exprdata_long = melt(tmpexprdata,variable.name = "sampleid",value.name="log2cpm_voom")
                      data_long = left_join(countdata_long,log2cpm_long)
                      data_long = left_join(data_long,exprdata_long)
                      data_long$group = do.call(rbind,strsplit(as.character(data_long$sampleid),"_",fixed=TRUE))[,1]
                      tmpgeneidnames = colnames(geneids%>%select(-unique_id))
                      if(length(tmpgeneidnames)>0) {
                        data_long = data_long%>%select(-one_of(tmpgeneidnames))
                        }
                      
                      #expr_data = tmplog2cpm[,-1]
                      
                      print('analyze data: done')
                      
                      
                      return(list('group_names'=group_names,'sampledata'=sampledata,
                                  "results"=lmobj_res,"data_long"=data_long, "geneids"=geneids, 
                                  "expr_data"=expr_data,
                                  "data_results_table"=data_results_table))
                      
                    }     
                    
                    
                  })
                })


output$countdataDT <- renderDataTable({
  tmp <- inputDataReactive()
  if(!is.null(tmp)) tmp$data
})

observeEvent(input$upload_data, ({
  updateCollapse(session,id =  "input_collapse_panel", open="analysis_panel",style = list("analysis_panel" = "success",
                                                                                          "data_panel"="primary"))
}))

observeEvent(inputDataReactive(),({
  updateCollapse(session,id =  "input_collapse_panel", open="data_panel",style = list("analysis_panel" = "default",
                                                                                      "data_panel"="success"))
})
)

output$analysisoutput <- DT::renderDataTable({
  print("output$analysisoutput")
  getresults <- analyzeCountDataReactive()
  res = getresults$results
  DT::datatable(res)
})


output$downloadResults_CSV <- downloadHandler(filename = paste0("START_results_",Sys.Date(),".csv"),
                                              content = function(file) {
                                                write.csv(analyzeCountDataReactive()$data_results_table, file, row.names=FALSE)})

output$downloadResults_RData <- downloadHandler(filename= paste0("START_results_",Sys.Date(),".RData"),
                                                content=function(file){
                                                  tmp = analyzeCountDataReactive()
                                                  
                                                  group_names = tmp$group_names
                                                  sampledata = tmp$sampledata
                                                  results = tmp$results
                                                  data_long = tmp$data_long
                                                  geneids = tmp$geneids
                                                  expr_data = tmp$expr_data
                                                  data_results_table = tmp$data_results_table
                                                  
                                                  save(group_names,sampledata,results,
                                                       data_long,geneids,expr_data,
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




