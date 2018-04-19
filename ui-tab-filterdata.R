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

## ==================================================================================== ##
## This tab is used to filter the data, can be saved in various formats after filtering
## ==================================================================================== ##

tabPanel("Filter Data", 
         ## ==================================================================================== ##
         ## Left hand column has the filtera settings and options
         ## ==================================================================================== ##
         fluidRow(column(4,wellPanel(
           h5("Note: this does not redo any analysis, it merely filters the data already analyzed based on gene identifiers or sample/group names."),
           selectizeInput("datafilter_groups", label="Select Groups",
                          choices=NULL,
                          multiple=TRUE),
           selectizeInput("datafilter_samples", label="Select Samples",
                          choices=NULL,
                          multiple=TRUE),
           checkboxInput("datafilter_genelist","Filter by gene names",value=FALSE), #upload a file of gene names
           conditionalPanel("input.datafilter_genelist==true",
                            selectizeInput("datafilter_gene_select","Search gene names.",choices=NULL,multiple=TRUE),
                            fileInput('datafilter_gene_file', 'Or, upload file containing gene IDs\n(one row per gene)',
                                      accept=c('text/csv', 
                                               'text/comma-separated-values,text/plain', 
                                               '.csv'))
           ),
           # this is complicated since need to calculate fold change based on expression, but we may not know if that
           # expression value is logged or not, so how do we get fold change vs log fold change from an unknown value? needs more work
           # checkboxInput("datafilter_fc","Filter by fold change for a pair of groups",
           #               value=FALSE),
           # conditionalPanel(condition="input.datafilter_fc==true",
           #                  selectizeInput("datafilter_fold_change_groups", label="Select 2 Groups (numerator first, denominator second)",
           #                                 choices=NULL,
           #                                 selected=NULL,
           #                                 multiple=TRUE,options = list(maxItems = 2)),
           #                  numericInput("datafilter_log2fc_cut",
           #                               label="Choose abs(Log2FC) Minimum",
           #                               min= 0, max=Inf,value=0),
           #                  radioButtons("datafilter_log2fc_direction",label = "Restrict Direction of FC",
           #                               choices=c("No","Up in first group","Down in first group"))
           # ),#conditionalpanel
           checkboxInput("datafilter_signif","Filter by significance",
                         value=FALSE),
           conditionalPanel(condition="input.datafilter_signif==true",
                            selectizeInput("datafilter_selecttest","Select Test",choices=NULL),
                            numericInput("datafilter_fccut",label="Minimum absolute value of Log2 Fold Change",
                                         min=0,max= 20,value=1),
                            radioButtons("datafilter_logfc_dir",label="Direction of Fold Change",
                                         choices=c("Both"="both","Up in first group"="up","Down in first group"="down")),
                            numericInput("datafilter_qvaluecut",label="Maximim Q-value (FDR)",
                                         min=0,max= 1,value=1),
                            numericInput("datafilter_pvaluecut",label="Maximim P-value",
                                         min=0,max= 1,value=1)
           ),
           checkboxInput("datafilter_expr","Filter by expression",
                         value=FALSE),
           conditionalPanel(condition="input.datafilter_expr==true",
                            radioButtons("datafilter_selectexpr",label="Select Expression Value",choices=""),
                            numericInput("datafilter_exprmin",label="Filter by minimum expression (remove genes with minimum expression lower than this value)",
                                         min=0,max= 1000,value=1),
                            numericInput("datafilter_exprmax",label="Filter by maximum expression (remove genes with maximum expression lower than this value)",
                                         min=0,max= 1000,value=1)
           )
           
           
         )#wellpanel
         ),#column
         ## ==================================================================================== ##
         ## Right hand column shows data input DT and data analysis result DT
         ## ==================================================================================== ##
         column(8,wellPanel(
           downloadButton('download_filtered_data_csv','Save Filtered Results as CSV File'),
           dataTableOutput('filterdataoutput')
         )#wellpanel
         )#column
         )#fluidRow
         
         #download results, download START data with current results but filtered gene lists; download count data with filtered gene lists
)