
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
## Analysis Plots
## ==================================================================================== ##   
tabPanel("Analysis Plots",  
         fluidRow(column(4,wellPanel(
           conditionalPanel("input.analysisres_tabset=='Volcano Plot'",
                            selectizeInput("analysisres_test", label="Select Test for Volcano Plot",
                                           choices=NULL),
                            numericInput("analysisres_fold_change_cut",
                                         label="Choose log2(Fold Change) Threshold\n(based on your input FCs, 
                       or fitted FCs if your data has been analyzed by START)",min= 0, max=20,value=1),
                            numericInput("analysisres_pvalcut",
                                         label="Choose P-value Threshold",min=0,max=1,value=0.05),
                            numericInput("analysisres_fdrcut",
                                         label="Choose adjusted P-value (or FDR) Threshold",min=0,max=1,value=0.05)
           ),#conditionalpanel
           
           conditionalPanel("input.analysisres_tabset=='Scatterplot of Fold Changes'",
                            selectizeInput("analysisres_groups",label="Select Groups for Scatterplot",
                                           choices=NULL,
                                           multiple=TRUE,options = list(maxItems = 2)),
                            radioButtons("scattervaluename",label="Select Scatterplot Value",choices=""),
                            radioButtons("scattercolor",label="Select Color Factor Value",
                                         choices=c("Sign of FC","logFC","P.Value","adj.P.Val"),
                                         selected = "P.Value"),
                            radioButtons("scatterresultsname",label="Select Test for Color Factor",
                                         choices="")
                            
           ),#conditionalpanel
           selectizeInput("analysisres_genes",label="Select Genes to Highlight",
                          choices=NULL,
                          multiple=TRUE)
         )#,
         #img(src="KCardio_CMYK_4C_pos_small.jpg",height=150,width= 275,align="right")	
         ),#column
         column(8,
                tabsetPanel(id="analysisres_tabset",
                            tabPanel(title="Volcano Plot",
                                     #h5(textOutput("corPR")),
                                     #uiOutput("volcanoplot_2groups_ggvisUI"),
                                     #ggvisOutput("volcanoplot_2groups_ggvis")  
                                     plotlyOutput("volcanoplot",height=600)
                            ),#tabPanel
                            tabPanel(title="Scatterplot of Fold Changes",
                                     #h5(textOutput("corPR")),
                                     # uiOutput("scatterplot_fc_2groups_ggvisUI"),
                                     # ggvisOutput("scatterplot_fc_2groups_ggvis")  
                                     plotlyOutput("scatterplot",height=600)
                            )#tabPanel
                )#tabsetPanel
         )#column
         )#fluidrow
) #END tabPanel