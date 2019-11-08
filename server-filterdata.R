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

observe({
  print("server-datafilter-update-filters")
  data_analyzed = analyzeDataReactive()
  tmpdatlong = data_analyzed$data_long
  tmpynames = tmpdatlong%>%select(-unique_id,-sampleid,-group,-one_of("rep"))%>%colnames()
  tmpgroups = data_analyzed$group_names
  tmpsamples = as.character(data_analyzed$sampledata$sampleid)
  tmpgeneids = data_analyzed$geneids
  data_analyzedgenes = as.character(unlist(tmpgeneids))
  tmpdat = data_analyzed$results
  tmptests = unique(as.character(tmpdat$test))
  
  updateSelectizeInput(session,"datafilter_groups", 
                       choices=tmpgroups,selected=tmpgroups)
  updateSelectizeInput(session,"datafilter_samples", 
                       choices=tmpsamples,selected=tmpsamples)
  
  updateSelectizeInput(session,"datafilter_gene_select",
                       choices=data_analyzedgenes,server=TRUE)
  
  updateSelectizeInput(session,"datafilter_selecttest",choices=tmptests)
  
  updateRadioButtons(session,'datafilter_selectexpr',
                     choices=sort(tmpynames,decreasing = TRUE))
  
})

# after selecting test

observe({
  print("server-datafilter-update-tests")
  data_analyzed = analyzeDataReactive()
  if(!(input$datafilter_selecttest=="")) {
    tmptest = input$datafilter_selecttest
    # get max abs fold change for this test
    tmpdat = data_analyzed$results
    tmpdat = tmpdat%>%filter(test==tmptest)
    tmpfc = abs(tmpdat$logFC)
    tmpfc = tmpfc[tmpfc<Inf]
    tmpmax = max(tmpfc,na.rm=T)
    if(tmpmax==Inf)
      
      updateNumericInput(session,"datafilter_fccut",
                         min=0,max= ceiling(tmpmax),value=0)
  }
})

# after selecting expression value
observe({
  print("server-datafilter-update-expr")
  data_analyzed = analyzeDataReactive()
  if(!(input$datafilter_selectexpr=="")) {
    exprname = input$datafilter_selectexpr
    #calculate miin and max
    tmpdat = data_analyzed$data_long # add filter by group and sample id
    tmpmin = min(tmpdat[,colnames(tmpdat)==exprname],na.rm=T)
    tmpmax = max(tmpdat[,colnames(tmpdat)==exprname],na.rm=T)
    
    updateNumericInput(session,"datafilter_exprmin",
                       min=floor(tmpmin),max= ceiling(tmpmax),value=floor(tmpmin))
    updateNumericInput(session,"datafilter_exprmax",
                       min=floor(tmpmin),max= ceiling(tmpmax),value=ceiling(tmpmax))
  }
})



filterDataReactive <- reactive({
  print("filterDataReactive")
  data_analyzed = analyzeDataReactive()
  tmpsampledata = data_analyzed$sampledata
  tmpgeneids = data_analyzed$geneids
  tmpres = data_analyzed$results
  tmpgroups = data_analyzed$group_names
  tmpdatlong = data_analyzed$data_long
  
  # tmpdatlong = data_analyzed$data_long
  # tmpynames = tmpdatlong%>%select(-unique_id,-sampleid,-group)%>%colnames()
  # 
  # tmptests = unique(as.character(tmpdat$test))
  
  mydata <- data_analyzed$data_results_table
  mydata_genes = left_join(mydata,tmpgeneids) # need also to have unique id
  
  groupids = lapply(tmpgroups,function(k) grep(k,colnames(mydata)))
  
  # filter by group
  if(!(input$datafilter_groups[1]=="")) {
    tmpselected = input$datafilter_groups
    tmprem = match(as.character(tmpsampledata$sampleid[which(!(tmpsampledata$group%in%tmpselected))]),colnames(mydata))
    tmpkeep = setdiff(1:ncol(mydata),tmprem)
    mydata = mydata[,tmpkeep]
  }
  
  # filter by sampleid
  if(!(input$datafilter_samples[1]=="")) {
    tmpselected = input$datafilter_samples
    tmpsamplesrem = setdiff(as.character(tmpsampledata$sampleid),tmpselected) # leftover samples
    tmprem = match(tmpsamplesrem,colnames(mydata))
    tmpkeep = setdiff(1:ncol(mydata),tmprem)
    mydata = mydata[,tmpkeep]
  }
  
  # filter by geneid or name
  if((input$datafilter_genelist)&(length(input$datafilter_gene_select)>0)) {
    tmpselected = input$datafilter_gene_select
    
    # find the columns with gene identifiers
    tmpmydata_genes = mydata_genes[,match(colnames(tmpgeneids),colnames(mydata_genes),nomatch=0)]
    # try to match gene names to each column, then take the union of all the indx
    tmpind = unique(na.omit(c(apply(tmpmydata_genes,2,function(k) match(tmpselected,k)))))
    
    mydata = mydata[tmpind,]
    mydata_genes = mydata_genes[tmpind,]
  }
  
  #add filter by gene name file like in heatmap
  
  if(input$datafilter_signif) {
    tmpres_filter = tmpres%>%filter(test==input$datafilter_selecttest)
    tmpres_filter = tmpres_filter%>%filter(P.Value<=input$datafilter_pvaluecut,
                                           adj.P.Val<=input$datafilter_qvaluecut)
    tmpres_up = tmpres_filter%>%filter(logFC>=input$datafilter_fccut)
    tmpres_down = tmpres_filter%>%filter(logFC<=input$datafilter_fccut)
    if(input$datafilter_logfc_dir=="up") {
      tmpgenes=as.character(tmpres_up$unique_id)
    }else if(input$datafilter_logfc_dir=="down"){
      tmpgenes=as.character(tmpres_down$unique_id)
    }else{
      tmpgenes=c(as.character(tmpres_up$unique_id),as.character(tmpres_down$unique_id))
    }
    tmpind = match(tmpgenes,mydata_genes$unique_id,nomatch=0)
    
    mydata = mydata[tmpind,]
    mydata_genes = mydata_genes[tmpind,]
  }
  
  
  if(input$datafilter_expr) {
    tmpdatlong_filter = data.table::data.table(tmpdatlong)[unique_id%in%mydata_genes$unique_id,]
    data.table::setnames(tmpdatlong_filter,input$datafilter_selectexpr, "mycol")
    tmpdatlong_filter = tmpdatlong_filter[,.(min=min(mycol,na.rm=T),max=max(mycol,na.rm=T)),by=unique_id]
    tmpdatlong_filter = tmpdatlong_filter[(min>=input$datafilter_exprmin)&(max<=input$datafilter_exprmax),]

    tmpgenes = as.character(tmpdatlong_filter$unique_id)
    tmpind = match(tmpgenes,mydata_genes$unique_id,nomatch=0)
    
    mydata = mydata[tmpind,]
    mydata_genes = mydata_genes[tmpind,]
  }
  
  
  shiny::validate(need(nrow(mydata)>0,message="No features pass this filter."))
  mydata  
  # need to add: 
  # filter only within some groups, should filter be based on above group selections? no because test is not
  # save data as file with filter settings concatinated?
  # show number of genes that pass filter like in heatmap
  # data frame display too wide, truncate columns?
  # output data should have counts too? which is shown, log2cpm? log2cpm_voom?
  #download record of filters buttons
  #DF display, make prettier?
  #data summary?
})


output$filterdataoutput <- renderDataTable({
  print("output$filterdataoutput")
  res <- filterDataReactive()
  res[,sapply(res,is.numeric)] <- signif(res[,sapply(res,is.numeric)],3)
  datatable(res, class = 'cell-border stripe', rownames = FALSE)
})

output$nrow_filterdata <- renderText({
  res <- filterDataReactive()
  tmpnum = ifelse(is.null(res),0,nrow(res))
  paste("Chosen filters result in ",tmpnum, " genes.")
})


# if datafilter_fold_change_groups selected

# observe({
#   group1 = input$datafilter_fold_change_groups[1]
#   group2 = input$datafilter_fold_change_groups[2]
#   tmpdatlong%>%filter(group==group1)
#   tmpdatlong%>%filter(group==group2)
#   updateNumericInput(session,"datafilter_log2fc_cut",min=0,max=max(abs))
# })


output$download_filtered_data_csv <- downloadHandler(filename = paste0("START_results_filtered_",Sys.Date(),".csv"),
                                              content = function(file) {
                                                write_csv(filterDataReactive(), file)})



