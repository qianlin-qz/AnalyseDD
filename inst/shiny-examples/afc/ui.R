
library(FactoMineR)
library(DT)

navbarPage("Analyses factorielles des correspondances",
  tabPanel("Original data",
   DT::dataTableOutput("table")
    ),

  tabPanel("Plot",
   sidebarLayout(
     sidebarPanel(
       radioButtons("plot_type", "Plot type",
         c("Types of Accommodations", "Socio-Professional Categories", "Comparison")
       )
   ),
   mainPanel(
    plotOutput("plot")
   )
  )
 ),

 tabPanel("Details",
  verbatimTextOutput("summary")
 )
)


