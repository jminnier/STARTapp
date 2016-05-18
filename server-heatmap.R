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
## HEATMAP variables
## ==================================================================================== ##		

#update list of groups
observe({
  data_analyzed = analyzeCountDataReactive()
  tmpgroups = data_analyzed$group_names
  tmpdat = data_analyzed$results
  tmptests = unique(as.character(tmpdat$test))
  tmpdatlong = data_analyzed$data_long
  tmpynames = tmpdatlong%>%select(-unique_id,-sampleid,-group,-count)%>%colnames()
  
  updateRadioButtons(session,'heatmapvaluename', choices=sort(tmpynames,decreasing = TRUE))
  updateCheckboxGroupInput(session,'view_group_heatmap',
                           choices=tmpgroups, selected=tmpgroups)
  updateSelectizeInput(session,'sel_test_heatmap',
                       choices=tmptests, selected=NULL)
  updateSelectizeInput(session,"fold_change_groups",
                       choices=tmpgroups)
  
})


inputHeatmapSubsetReactive <- reactive({
  
  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.
  
  validate(
    need((input$heatmap_subset=="all")|(!is.null(input$heatmap_file)), "Please select a file")
  )
  inFile <- input$heatmap_file
  if(!is.null(inFile)) {
    heatmap_geneids <- unlist(read.csv(inFile$datapath, header=FALSE, stringsAsFactors = FALSE))
    print('uploaded heatmap gene ids')
    
    
    data_analyzed = analyzeCountDataReactive()
    tmpgeneids = data_analyzed$geneids
    tmpmatch = apply(tmpgeneids,2,function(k) match(heatmap_geneids,k,nomatch = 0))
    subsetids = tmpgeneids$unique_id[unique(c(tmpmatch))]
    print(subsetids)
    
    validate(need(length(subsetids)>0,"No match found."))
    
    return(subsetids)
  }else{return(NULL)}
  
  
  
})

# output$filter_fc_ui <- renderUI({
#   data_analyzed = analyzeCountDataReactive()
#   tmpgroups = data_analyzed$group_names
# 
#   
#   list=list(
#     selectizeInput("fold_change_group", label="Select 2 Groups",
#                    choices=tmpgroups,
#                    selected=tmpselect,multiple=TRUE,options = list(maxItems = 2)),
#     sliderInput("fold_change_range",
#                 label="Choose Log2Fold Change Filter",min= -20, max=20,value=c(-20,20))
#   )
# })

# #withProgress(message = "Drawing heat map, please wait", #detail = "This may take a few moments...",{})

output$heatmap_rna <- renderPlot({
  if(input$action_heatmaps==0) return()
  #input$action_heatmaps
  
  data_analyzed = analyzeCountDataReactive()
  subsetids = inputHeatmapSubsetReactive()
  
  isolate({ #avoid dependency on everything else except action_heatmaps
    print("drawing heatmap rna")
    withProgress(message = "Drawing heatmap, please wait",{
      heatmap_render(
        data_analyzed=data_analyzed,
        yname = input$heatmapvaluename,
        usesubset = input$heatmap_subset=="subset",
        subsetids = subsetids,
        orderby = input$heatmap_order,
        FDRcut=input$FDRcut,
        maxgenes=input$maxgenes,
        view_group=input$view_group_heatmap,
        sel_test=input$sel_test_heatmap,
        filter_fdr=input$filter_fdr,
        filter_maxgene=input$filter_maxgene,
        # filter_cpm=input$filter_cpm,
        filter_fc=input$filter_fc,
        fold_change_range=input$fold_change_range,
        fold_change_groups=input$fold_change_groups)
    })
  })#isolate
})



observe({
  if(input$action_heatmaps==0) return()
  
  data_analyzed = analyzeCountDataReactive()
  subsetids = inputHeatmapSubsetReactive()
  
  isolate({
    #input$action_heatmaps
    #browser()
    mydat = heatmap_ggvis_data(
      data_analyzed = data_analyzed,
      yname = input$heatmapvaluename,
      usesubset = input$heatmap_subset=="subset",
      subsetids = subsetids,
      orderby = input$heatmap_order,
      FDRcut=input$FDRcut,
      maxgenes=input$maxgenes,
      view_group=input$view_group_heatmap,
      sel_test=input$sel_test_heatmap,
      filter_fdr=input$filter_fdr,
      filter_maxgene=input$filter_maxgene,
      # filter_cpm=input$filter_cpm,
      filter_fc=input$filter_fc,
      fold_change_range=input$fold_change_range,
      fold_change_groups=input$fold_change_groups)
    #print(mydat)
    #print(unique(mydat$unique_id))
    
    if(!is.null(mydat)) {
      
      yname = input$heatmapvaluename
      all_values2 <- function(x) {
        if(is.null(x)) return(NULL)
        thisrow <- mydat[mydat$id == x$id,]
        thisrow <- thisrow[,c("unique_id","sample",yname,paste0(yname,"_scaled"))]
        thisrow[,yname] <- round(thisrow[,yname],3)
        paste0(names(thisrow), ": ", format(thisrow), collapse = "<br />")
      }
      
      #browser()
      mydat%>%
        ggvis(~sample, ~unique_id, fill:=~defcolor,stroke:=~defcolor, key:=~id)%>%
        layer_rects(width = band(), height = band()) %>%
        scale_nominal("x", padding = 0, points = FALSE) %>%
        scale_nominal("y", padding = 0, points = FALSE) %>%
        add_tooltip(html=all_values2,on = "hover") %>%
        #lb$input()%>%
        set_options(width=600,height=1000)%>%bind_shiny("heatmapggvis_rna","heatmapggvisUI_rna")
      
      if(input$heatmap_subset=="subset") {
        tmptext = renderText("The heatmap is filtered based on the gene ids you uploaded.")
      }else{
        if(input$heatmap_order=="significance") {
          tmptext = renderText(paste(
            "The heatmap is filtered based on the smallest p-values
        from the following test:",input$sel_test_heatmap))
        }else{ 
          tmptext = renderText(paste(
            "The heatmap is filtered based on the largest standard deviation of the expression value selected."))
        }
      }
      
      output$heatmap_rna_title = tmptext
      output$heatmap_rna_title_int = tmptext
    }else{
      tmptext = renderText("These inputs do not produce data.")
      output$heatmap_rna_title = tmptext
      output$heatmap_rna_title_int = tmptext
    }
    
  })
})



HeatdatReactive_rna <- reactive({
  if(input$action_heatmaps==0) return()
  #input$action_heatmaps
  
  data_analyzed = analyzeCountDataReactive()
  subsetids = inputHeatmapSubsetReactive()
  
  tmp<- heatmap_data(
    data_analyzed=data_analyzed,
    yname = input$heatmapvaluename,
    usesubset = input$heatmap_subset=="subset",
    subsetids = subsetids,
    orderby = input$heatmap_order,
    FDRcut=input$FDRcut,
    maxgenes=input$maxgenes,
    view_group=input$view_group_heatmap,
    sel_test=input$sel_test_heatmap,
    filter_by_go=FALSE,
    filter_fdr=input$filter_fdr,
    filter_maxgene=input$filter_maxgene,
    filter_cpm=input$filter_cpm,
    filter_fc=input$filter_fc,
    fold_change_range=input$fold_change_range,
    fold_change_groups=input$fold_change_groups,
    group_filter_range=list(
      input$cpm_range_BE,
      input$cpm_range_HE,
      input$cpm_range_HP,
      input$cpm_range_SM),
    fixed_genes_list=NULL)
  return(tmp)
})

output$numheat <- renderText({
  if(input$action_heatmaps==0) return()
  #input$action_heatmaps
  tmp = HeatdatReactive_rna()
  #return(as.character(nrow(tmp)))
  tmpnum = ifelse(is.null(tmp),0,nrow(tmp))
  paste("Chosen filters result in ",tmpnum, " genes.")
})

output$heatdat_rna <- DT::renderDataTable({
  if(input$action_heatmaps==0) return()
  tmpdat = HeatdatReactive_rna()
  DT::datatable(tmpdat)
},options=list(sDom="ilftpr"))

output$downloadHeatmapData_rna <- downloadHandler(
  filename = c('heatmap_data_rna.csv'),
  content = function(file) {
    write.csv(HeatdatReactive_rna(), file, row.names=FALSE)}
)

#   fixedGenesReactive <- reactive({
#     if(input$fix_genes) {tmp = isolate(HeatdatReactive_rna()); return(tmp$gene.id)}else{return(NULL)}
#   })