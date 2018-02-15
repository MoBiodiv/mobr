#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(mobr)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("MOBR"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput("comMat", "Choose Community Matrix"),
      
      fileInput("siteTab", "Choose Site Attribute Table")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
       plotOutput("rarePlot")
    )
  )
))
