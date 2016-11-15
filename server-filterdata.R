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
  tmpynames = tmpdatlong%>%select(-unique_id,-sampleid,-group)%>%colnames()
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
  if(!(input$datafilter_selectexpr=="")) {
    exprname = input$datafilter_selectexpr
    #calculate miin and max
    tmpdat = data_analyzed$data_long # add filter by group and sample id
    tmpmin = min(tmpdat[,colnames(tmpdat)==exprname],na.rm=T)
    tmpmax = max(tmpdat[,colnames(tmpdat)==exprname],na.rm=T)
    
    updateNumericInput(session,"datafilter_exprmin",
                 min=tmpmin,max= tmpmax,value=tmpmin)
    updateNumericInput(session,"datafilter_exprmax",
                 min=tmpmin,max= tmpmax,value=tmpmax)
  }
})



#download buttons
#DF display, make prettier?
#data summary?






# if datafilter_fold_change_groups selected

# observe({
#   group1 = input$datafilter_fold_change_groups[1]
#   group2 = input$datafilter_fold_change_groups[2]
#   tmpdatlong%>%filter(group==group1)
#   tmpdatlong%>%filter(group==group2)
#   updateNumericInput(session,"datafilter_log2fc_cut",min=0,max=max(abs))
# })







