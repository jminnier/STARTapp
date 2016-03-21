# This tab is used to input the count or normalized data files

tabPanel("Input Data", 
         sidebarLayout(
           sidebarPanel(
             radioButtons('use_example_file','Use example file or upload your own data',
                          c('Example: RNA-seq gene counts'="examplecounts",
                            'RData from previous START upload'="previousrdata",
                            'Upload Data'="upload")),
             conditionalPanel(condition="input.use_example_file=='previousrdata'",
                              fileInput('rdatafile','Upload Previously Downloaded RData File')
                              ),
             conditionalPanel(condition="input.use_example_file=='upload'",
                              radioButtons("inputdat_type","Input Data Type:",
                                           c("Raw data: Gene Counts"="counts",
                                             "Analyzed data: Expression Values, p-values, fold changes"="analyzed")),
                                               fileInput('datafile', 'Choose CSV or TSV File Containing Data',
                                                         accept=c('text/csv', 
                                                                  'text/comma-separated-values,text/plain', 
                                                                  '.csv')),
                                               radioButtons('sep', 'Separator',
                                                            c(Comma=',',
                                                              Tab='\t'),
                                                            ',',inline = TRUE),
                                               tags$hr(),
                              
             conditionalPanel(condition="input.inputdat_type=='counts'",
                              downloadLink("example_counts_file",label="Example Count File"),
                              p("File must have a header row. First column(s) must be gene identifiers. Column names of 
                                                 count data must include group names and replicate number in this format:
                                Group1_1, Group1_2, Group2_1, Group2_2..."),
                              numericInput("numgeneids",label="Number of columns with gene identifiers",
                                           min=1,max= 50,value=1,step=1)
                              ),
             conditionalPanel(condition="input.inputdat_type=='analyzed'",
                              downloadLink("example_analysis_file",label="Example Analysis Results File"),
                              p("File must have a header row. Column names of expression values must be in the format
                                Group1_1, Group1_2, Group2_1, Group2_2... Number/order of fold changes must match p-values."),
                              #checkboxInput('header', 'Header', TRUE),
                        selectInput("c_geneid1",label="First column # with gene IDs",choices=NULL),
                        selectInput("c_geneid2",label="Last column # with gene IDs",choices=NULL),
                        selectInput("c_expr1",label="First column # with expression values",choices=NULL),
                        selectInput("c_expr2",label="Last column # with expression values",choices=NULL),
                        selectInput("c_fc1",label="First column # with fold changes",choices=NULL),
                        selectInput("c_fc2",label="Last column # with fold changes",choices=NULL),
                        radioButtons("isfclogged",label="Is FC logged? (if false, expression values will be log2-transformed for visualization)",choices=c("Yes (Leave it alone)","No (Log my data please)"),selected="No (Log my data please)"),
                        selectInput("c_pval1",label="First column # with p-values",choices=NULL),
                        selectInput("c_pval2",label="Last column # with p-values",choices=NULL)
                              )
             ),
             conditionalPanel("output.fileUploaded",actionButton("upload_data","Submit Data")),
             # add reference group selection
             # add instructions
             # missing value character?
             # allow to input counts and analyze via standard limma, OR
             # input counts, fitted logcpm, and p-values? if multiple comparisons
             # p-values would be complicated
             br(),br(),br(),br(),br(),br(),br(),br(), br(),br(),br(), 
             img(src="KCardio_CMYK_4C_pos_small.jpg",height=150,width= 275,align="right")	
           ),
           mainPanel(
             bsCollapse(id="input_collapse_panel",open="data_panel",multiple = FALSE,
                        bsCollapsePanel(title="Data Contents: Check Before `Submit`",value="data_panel",
                                        dataTableOutput('countdataDT')                       
                        ),
                        bsCollapsePanel(title="Analysis Results: Ready to View Other Tabs",value="analysis_panel",
                                        downloadButton('downloadResults_CSV','Save Results as CSV File'),
                                        downloadButton('downloadResults_RData','Save Results as RData File for Future Upload'),
                                        dataTableOutput('analysisoutput')
                        )
             )#tabsetPanel
           )#mainPanel
         )
)
