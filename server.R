#server
server <- function(input, output) {  
  output$distPlot1 <- renderPlot({plot(input$qtype, type=input$type)+
      ggtitle("How are feeling with following statements today?")}, height = 1000)
  scale <- reactive({
    get(input$qtype)
  })
  
  output$dat <- renderPrint({
    scale()
  })
  
}

# server()
# 
# shiny::runApp()

# And then do plot(scale() will show selected plot.
