
## ==================================================================================== ##
## Analysis Plots
## ==================================================================================== ##   
tabPanel("Analysis Plots",  
         sidebarLayout(sidebarPanel( 
           conditionalPanel("input.analysisres_tabset=='Volcano Plot'",
                            selectizeInput("analysisres_test", label="Select Test for Volcano Plot",
                                           choices=NULL),
                            numericInput("analysisres_fold_change_cut",
                                         label="Choose log2(Fold Change) Threshold\n(based on your input FCs, 
                       or fitted FCs if your data has been analyzed by START)",min= 0, max=20,value=2),
                            numericInput("analysisres_fdrcut",
                                         label="Choose P-value Threshold\n(based on your input p-value, 
                        or the FDR adjusted p-value if 
                                         data has been analyzed by START)",min=0,max=1,value=0.05)
           ),#conditionalpanel
           
           conditionalPanel("input.analysisres_tabset=='Scatterplot of Fold Changes'",
                            selectizeInput("analysisres_groups",label="Select Groups for Scatterplot",
                                           choices=NULL,
                                           multiple=TRUE,options = list(maxItems = 2)),
                            radioButtons("scattervaluename",label="Select Scatterplot Value",choices="")
           ),#conditionalpanel
           br(),br(),br(),br(),br(),br(),br(),br(), br(),br(),br(), 
           img(src="KCardio_CMYK_4C_pos_small.jpg",height=150,width= 275,align="right")	
         ),#sidebarPanel
         mainPanel(
           tabsetPanel(id="analysisres_tabset",
                       tabPanel(title="Volcano Plot",
                                #h5(textOutput("corPR")),
                                uiOutput("volcanoplot_2groups_ggvisUI"),
                                ggvisOutput("volcanoplot_2groups_ggvis")  
                       ),#tabPanel
                       tabPanel(title="Scatterplot of Fold Changes",
                                #h5(textOutput("corPR")),
                                uiOutput("scatterplot_fc_2groups_ggvisUI"),
                                ggvisOutput("scatterplot_fc_2groups_ggvis")  
                       )#tabPanel
           )#tabsetPanel
         )#mainPanel
         )#sidebarLayout
) #END tabPanel