library(scales)
library(HH)
library(shiny)
input <- 0
ouput <- 10
a <- sample(rep((1:5),5))
b <- sample(rep((1:5),5))
c <- data.frame(sapply(data.frame(a), factor))
d <- data.frame(sapply(data.frame(b), factor))

# scaledc <- likert(a, as.percent = T)
# scaledd <- likert(b)

# ?likert

Pop <- rbind(a=c(3,2,4,9), b=c(6,10,12,10))
dimnames(Pop)[[2]] <- c("Very Low", "Low", "High", "Very High")

scaledc <- scaledd <- likert(as.listOfNamedMatrices(Pop),
       as.percent=FALSE,
       resize.height="rowSums",
       strip=FALSE,
       strip.left=FALSE,
       main=paste("Area and Height are proportional to 'Row Count Totals'.",
                  "Width is exactly 100%.", sep="\n"))

# data(ProfChal)
# ## Percent plot calculated automatically from Count data
# a <- likert(Question ~ . , ProfChal[ProfChal$Subtable=="Employment sector",],
#        as.percent=TRUE,
#        main='Is your job professionally challenging?',
#        ylab=NULL,
#        sub="This plot looks better in a 9in x 4in window.")
# 
# str(a)
# The shiny codes are:
  
ui <- fluidPage(
  titlePanel("Survey"),
  sidebarLayout(
    sidebarPanel(
      selectInput("type",
                  "Plot Type",
                  choices = c("Likert"="bar",
                              "Density"="density",
                              "Heatmap"="heat"), selected="Likert"),
      radioButtons("qtype", 
                   "Question type:",
                   c("Agreement"="scaledc", "Helpfulness"="scaledd"),
                   selected="scaledc")
    ),
    
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Yearly Data", plotOutput("distPlot1"))
      )
    )
  )
)

# ?tabPanel()
