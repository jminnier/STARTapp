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
## 
## # This tab is used to input the count or normalized data files

tabPanel("Input Data", 
         fluidRow(column(4,wellPanel(
           downloadLink("instructionspdf",label="Download Instructions (pdf)"),
           radioButtons('use_example_file','Use example file or upload your own data',
                        c('Example: RNA-seq gene counts'="examplecounts",
                          'RData from previous START upload'="previousrdata",
                          'Upload Data'="upload")),
           conditionalPanel(condition="input.use_example_file=='previousrdata'",
                            fileInput('rdatafile','Upload Previously Downloaded RData File'),
                            conditionalPanel("output.fileUploaded",h4(strong("Check data contents then click:")))
           ),
           conditionalPanel(condition="input.use_example_file=='upload'",
                            radioButtons("inputdat_type","Input Data Type:",
                                         c("Count data: Gene Counts"="counts",
                                           #"Microarray expression data"="microarray",
                                           "Analyzed data: Expression Values, p-values, fold changes"="analyzed")),
                            # radioButtons('sep', 'Separator',
                            #              c(Comma=',',
                            #                Tab='\t'),
                            #              ',',inline = TRUE),
                            tags$hr(),
                            
                            conditionalPanel(condition="input.inputdat_type=='counts'",
                                             downloadLink("example_counts_file",label="Download Example Count File"),
                                             p(""),
                                             img(src="examplecounts.png",width="100%"),
                                             p("File must have a header row. First column(s) must be gene identifiers. 
Column names of 
                                                 count data must include group names and replicate number in this format:
                                Group1_1, Group1_2, Group2_1, Group2_2...")#,
                                             # numericInput("numgeneids",label="Number of columns with gene identifiers",
                                             #              min=1,max= 50,value=1,step=1)
                            ),
                            conditionalPanel(condition="input.inputdat_type=='analyzed'",
                                             downloadLink("example_analysis_file",label="Download Example Analysis Results File"),
                                             p(""),
                                             img(src="exampleanalysisdata.png",width="100%"),
                                             #a("<img src=`exampleanalysisdata.tiff` width=100%>"),href="www.google.com",target="_blank"),
                                             p("File must have a header row. Column names of expression values must be in the format
                                Group1_1, Group1_2, Group2_1, Group2_2... Number/order of fold changes must match p-values.")),
                            fileInput('datafile', 'Choose File Containing Data (.CSV)',
                                      accept=c('text/csv', 
                                               'text/comma-separated-values,text/plain', 
                                               '.csv')),
                            conditionalPanel(condition="input.inputdat_type=='analyzed'",
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
           conditionalPanel("output.fileUploaded",actionButton("upload_data","Submit Data",
                                                               style="color: #fff; background-color: #CD0000; border-color: #9E0000"))
         )#,
         # add reference group selection
         # add instructions
         # missing value character?
         # allow to input counts and analyze via standard limma, OR
         # input counts, fitted logcpm, and p-values? if multiple comparisons
         # p-values would be complicated
         #img(src="KCardio_CMYK_4C_pos_small.jpg",height=150,width= 275,align="right")	
         ),#column
         column(8,
                bsCollapse(id="input_collapse_panel",open="data_panel",multiple = FALSE,
                           bsCollapsePanel(title="Data Contents: Check Before `Submit`",value="data_panel",
                                           dataTableOutput('countdataDT')                       
                           ),
                           bsCollapsePanel(title="Analysis Results: Ready to View Other Tabs",value="analysis_panel",
                                           downloadButton('downloadResults_CSV','Save Results as CSV File'),
                                           downloadButton('downloadResults_RData','Save Results as RData File for Future Upload'),
                                           dataTableOutput('analysisoutput')
                           )
                )#bscollapse
         )#column
         )#fluidrow
)#tabpanel
