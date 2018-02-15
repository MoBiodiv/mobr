#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(mobr)
library(dplyr)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  output$rarePlot <- renderPlot({
    
    outfile <- tempfile(fileext = ".png")
    
    file1 = write.csv(input$comMat, "file1.csv", row.names=FALSE)
    comMat = read.csv(file = "file1.csv", header=TRUE, sep=",")

    file2 = write.csv(input$siteTab, "file2.csv", row.names=FALSE)
    siteTab = read.csv(file = "file2.csv", header=TRUE, sep=",")
    #siteTab$group = as.numeric(factor(siteTab$group , levels=c("uninvaded", "invaded")))
    
    inv_mob_in = make_mob_in(comMat, siteTab)
    
    plot_rarefaction(inv_mob_in, 'group', 'spat', lwd=4)
    
  })
  
})

