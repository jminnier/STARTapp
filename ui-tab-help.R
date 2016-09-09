tabPanel("Instructions",
         fluidRow(
           column(4,wellPanel(
             h4("Instructions"),
             a("Input Data", href = "#inputdata"), br(),
             a("Data Formats", href = "#dataformat"), br(),
             a("Save Data For Future START Sessions", href="#rdata"), br(),
             a("Visualizations", href="#vis"), br(),
             a("PCA Plots", href="#pcaplots"), br(),
             a("Analysis Plots", href="#analysisplots"), br(),
             a("Volcano Plots", href="#volcano"), br(),
             a("Scatterplots", href="#scatterplots"), br(),
             a("Gene Expression Boxplots", href="#boxplots"), br(),
             a("Heatmaps", href="#heatmaps"), br()
           )
           ),#column
           column(8,
                  includeMarkdown("instructions/Instructions.md"))
         ))