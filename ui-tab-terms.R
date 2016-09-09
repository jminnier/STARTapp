tabPanel("Terms & Conditions",
         fluidRow(
           column(4,wellPanel(
             h4("Shinyapps.io Terms & Conditions")
           )
           ),#column
           column(8,
                  includeMarkdown("instructions/terms.md"))
         ))