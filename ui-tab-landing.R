tabPanel("Getting Started",
         fluidRow(
           column(4,wellPanel(
             h4("Getting Started with START"),
             a("Features", href="#features"),br(),
             a("Data Formats", href = "#dataformats"), br(),
             a("Save Data for Future Upload", href="#savedata"), br(),
             a("More Help", href = "#help"), br()
           )
           ),#column
           column(8,
                  includeMarkdown("instructions/landing.md"))
         ))