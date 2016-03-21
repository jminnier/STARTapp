
## ==================================================================================== ##
## Group Plots
## ==================================================================================== ##   
tabPanel("Group Plots",  
         sidebarLayout(sidebarPanel( 
           selectizeInput("sampleres_groups", label="Select Groups",
                          choices=NULL,
                          multiple=TRUE),
           br(),br(),br(),br(),br(),br(),br(),br(), br(),br(),br(), 
           img(src="KCardio_CMYK_4C_pos_small.jpg",height=150,width= 275,align="right")	 
           
         ),#sidebarPanel
         mainPanel(
           tabsetPanel(
             tabPanel(title="PCA Plot",
                      plotOutput("pca_plot")
             ),#tabPanel
             tabPanel(title="Sample Correlation Heatmap",
                      plotOutput("gene_pheatmap") 
             )#tabPanel
           )#tabsetPanel
         )#mainPanel
         )#sidebarLayout
) #END tabPanel