options(shiny.maxRequestSize = 100*1024^2)

source("helpers.R")
print(sessionInfo())

shinyServer(function(input, output,session) {
  ## Server functions are divided by tab
  ## 
  
  
  source("server-inputdata.R",local = TRUE)
  source("server-dotplot.R",local = TRUE)
  source("server-heatmap.R",local = TRUE)
  source("server-samplegroupplots.R",local=TRUE)
  source("server-analysisres.R",local = TRUE)
  source("server-data.R",local = TRUE)
  
})