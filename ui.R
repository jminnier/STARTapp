## ==================================================================================== ##
# START Shiny App for analysis and visualization of transcriptome data.
# Copyright (C) 2014-2016  Jessica Minnier
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
# ui.R
source("helpers.R")
customHeaderPanel <- function(title,windowTitle=title){
  tagList(
    tags$head(
      tags$title(windowTitle),
      tags$link(rel="stylesheet", type="text/css",
                href="app.css"),
      tags$h1(a(href="http://www.ohsu.edu/xd/health/services/heart-vascular/"))
    )
  )
}
# collects all of the tab UIs
#shinyUI(
navbarPage(
  fluid = TRUE,theme = "bootstrap.min.united.css",
  #United theme from http://bootswatch.com/
  #customHeaderPanel(title="START: RNAseq Analysis and Visualization Resource"),#img(src="KCardio_CMYK_4C_pos_small.jpg",height=50,width= 92,align="right")	,
  title="START: Shiny Transcriptome Analysis Resource Tool", 
  ## =========================================================================== ##
  ## DOWNLOAD DATA TABS
  ## =========================================================================== ##
  source("ui-tab-inputdata.R",local=TRUE)$value,
  ## =========================================================================== ##
  ## Visualization TABS
  ## =========================================================================== ##
  source("ui-tab-samplegroupplots.R",local=TRUE)$value,
  source("ui-tab-analysisres.R",local=TRUE)$value,
  source("ui-tab-dotplot.R",local=TRUE)$value,
  source("ui-tab-heatmap.R",local=TRUE)$value,
  ## ============================================================================ ##
  ## INFO TAB
  ## ============================================================================ ##   
  tabPanel("Info",
           includeMarkdown("instructions/Instructions.md")),#end definitions of tabs, now footer
  
  ## ==================================================================================== ##
  ## FOOTER
  ## ==================================================================================== ##              
  footer=p(hr(),p("ShinyApp created by ", strong("{Jessica Minnier + Jiri Sklenar + Jonathan Nelson}")," of ",align="center",width=4),
           p(("Knight Cardiovascular Institute, Oregon Health & Science University"),align="center",width=4),
           p(("Copyright (C) 2014-2016, code licensed under GPLv3"),align="center",width=4),
           p(("Code available on Github:"),a("https://github.com/jminnier/STARTapp",href="https://github.com/jminnier/STARTapp"),align="center",width=4)
  )
  
  ## ==================================================================================== ##
  ## end
  ## ==================================================================================== ## 
) #end Navbar

