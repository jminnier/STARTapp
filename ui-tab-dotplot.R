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
## GENE DATA TAB 
## ==================================================================================== ## 

tabPanel("Gene Expression Boxplots",
         fluidRow(column(4,wellPanel(
           tags$head( tags$style("body {background-color: white; }")),
           selectizeInput("sel_gene",
                          label="RNA-Seq Gene Name (Select 1 or more)",#or Ensembl ID",
                          choices = NULL,
                          multiple=TRUE,
                          options = list(
                            placeholder = 
                              'Start typing to search for a gene name'# or ID',
                          ) #,
           ),#end selectInput
           #h5(htmlOutput("geneurl")),
           checkboxGroupInput("sel_group",
                              label="Select Group",
                              choices="", selected=""
           ),	
           radioButtons("sel_gene_header",label="Select Gene Identifier Label",
                        choices=" "),
           radioButtons("ytype","Y axis:",choices="")
                        #c("fitted cpm"="cpm","count"="count")),
#            radioButtons("log2cpm_checked","Y axis transformation:",
#                         c("log2"=TRUE,"raw value"=FALSE))
         ),#br(),br(),br(),br(),br(),br(),br(),br(), 
         img(src="KCardio_CMYK_4C_pos_small.jpg",height=150,width= 275,align="right")	
         ), 
         ## ==================================================================================== ##
         ## GENE DATA: DOT PLOT
         ## Search gene name to view dotplot(s) of expression    
         ## ==================================================================================== ## 
         column(8,
                tabsetPanel(
                  tabPanel(title="DotPlot",
                           plotlyOutput("dotplot",height=600)
                  ),#end tabPanel
                  tabPanel(title="Info",
                           h5("This panel constructs box and whisker plots of log2(CPM) or CPM values with dot plots
                              superimposed to show the raw data. When there are three data points the median and 
                              interquartile ranges are precisely the data values. Medians are denoted by horizontal
                              lines and averages are denoted by open diamonds."),
                           br(),br(),br(),br(),br(),br(),br(),br(), br(),br(),br(), 
                           img(src="KCardio_CMYK_4C_pos_small.jpg",height=150,width= 275,align="right")	
                  ),#end tabPanel
                  ## ==================================================================================== ##
                  ## GENE DATA: DOT PLOT DATA
                  ## ==================================================================================== ## 
                  tabPanel(title="Data Output",
                           downloadButton('downloadSubsetData', 'Download Data Subset as CSV File'),
                           DT::dataTableOutput("dat_dotplot")
                  )#tabsetPanel
                ))#column
         
         )#sidebarLayout	
) #end tabPanel Gene Data