## ==================================================================================== ##
## HEATMAP TAB 
## ==================================================================================== ## 	
tabPanel("Heatmaps",  
         sidebarLayout(sidebarPanel(
           actionButton("action_heatmaps","Generate Heatmaps"),  
           h6(textOutput("numheat")),
           radioButtons("heatmapvaluename",label="Select Value to Plot in Heatmap",choices=""),
           checkboxGroupInput("view_group_heatmap",
                              label=h5("Select Groups to View"),
                              choices="",
                              selected=""
           ),
           selectizeInput("sel_test_heatmap",
                              label=h5("Select Test to Use for Filtering"), 
                              choices="",
                              selected=""),	
           h3("Filters"),
           
           checkboxInput("filter_fdr","FDR cutoff",value = FALSE),
           conditionalPanel(condition="input.filter_fdr==true",
                            numericInput("FDRcut",label="Choose P-value 
                                         (FDR if analyzed by START) cutoff",
                                         min=0,max= 1,value=0.05)),
           
           checkboxInput("filter_maxgene",
                         "Show a maximum number of genes (recommended)",value=TRUE),
           conditionalPanel(condition="input.filter_maxgene==true",    		
                            numericInput("maxgenes",label="Choose Max # of Genes",
                                         min=1,max= 5000,value=100,step=1)),               
           
           checkboxInput("filter_fc","Filter by fold change for a pair of groups"),
           conditionalPanel(condition="input.filter_fc==true",
                            selectizeInput("fold_change_groups", label="Select 2 Groups",
                                           choices=NULL,
                                           selected=NULL,
                                           multiple=TRUE,options = list(maxItems = 2)),
                            sliderInput("fold_change_range",
                                        label="Choose Log2Fold Change Filter",min= -20, max=20,value=c(-20,20))),
br(),br(),br(),br(),br(),br(),br(),br(), 
img(src="KCardio_CMYK_4C_pos_small.jpg",height=150,width= 275,align="right")	
         ),#sidebarPanel
         mainPanel(
           tabsetPanel(
             tabPanel(title="HeatMap",
                      #textOutput("which_genes"),
                      h4(textOutput("heatmap_rna_title")),
                      plotOutput("heatmap_rna",height="800px")                        
             ),
             tabPanel(title="Interactive HeatMap",
                      h4(textOutput("heatmap_rna_title_int")),
                      uiOutput("heatmapggvisUI_rna"),
                      ggvisOutput("heatmapggvis_rna")
             ),
             tabPanel(title="Data Output",
                      downloadButton('downloadHeatmapData_rna', 
                                     'Download Heatmap Data as CSV File'),
                      DT::dataTableOutput("heatdat_rna"))
           )#tabsetPanel
         )#mainPanel
         )#sidebarLayout
) #tabPanel Heatmap