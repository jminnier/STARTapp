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
    title="START: RNAseq Analysis and Visualization Resource", 
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
    tabPanel("Info"),#end definitions of tabs, now footer
    
    ## ==================================================================================== ##
    ## FOOTER
    ## ==================================================================================== ##              
    footer=p(hr(),p("ShinyApp created by ", 
                    strong("{Jessica Minnier + Jiri Sklenar + Jonathan Nelson}")," of ",align="center",width=4),
                    p(strong("Knight Cardiovascular Institute, Oregon Health & Science University"),align="center",width=4),
                    p(strong("Copyright (C) 2014-2016, code licenesed under GPLv3"),align="center",width=4
    )
    )
    
    ## ==================================================================================== ##
    ## end
    ## ==================================================================================== ## 
  ) #end Navbar
#)#end shiny UI
