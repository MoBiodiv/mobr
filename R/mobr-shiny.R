# Module UI function
csvFileInput <- function(id, label = "CSV file") {
  # Create a namespace function using the provided id
  ns <- NS(id)

  tagList(
    fileInput(ns("file"), label)
  )
}

# Module server function
csvFile <- function(input, output, session, stringsAsFactors = TRUE) {
  # The selected file, if any
  userFile <- reactive({
    # If no file is selected, don't do anything
    validate(need(input$file, message = FALSE))
    input$file
  })

  # The user's data, parsed into a data frame
  dataframe <- reactive({
    read.csv(userFile()$datapath,
      header = TRUE,
      stringsAsFactors = stringsAsFactors)
  })

  # We can run observers in here if we want to
  observe({
    msg <- sprintf("File %s was uploaded", userFile()$name)
    cat(msg, "\n")
  })

  # Return the reactive that yields the data frame
  return(dataframe)
}

# Define UI for the mobr application
ui <- fluidPage(

    # App title ----
    titlePanel("MoB-R"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(

            # Input: Select a file ----
            csvFileInput("comm", "Upload community data (.csv format)"),
        
            # Horizontal line ----
            tags$hr(),

            # Input: Select a file ----
            csvFileInput("plot_attr", "Upload plot attribute data (.csv format)"),
            
            selectInput("graphType", "Select Graph Type",
                        c("Spacial Rarefaction", "Individual Rarefaction - Unpooled", "Individual Rarefaction - Pooled", "Unpooled Abundance", "Pooled Abundance", "All MoB Metrics", "MoB Metrics - S", "MoB Metrics - N", "MoB Metrics - S_n", "MoB Metrics - S_PIE", "MoB Delta Stats"))
            
            #textInput("filename", "Enter filename:"),
            
            #downloadLink("downloadData", "Download Graph"),
        
            #actionButton("do", "Submit")
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
          tabsetPanel(
            tabPanel("Plot", plotOutput('plot')),
            tabPanel("Stats", verbatimTextOutput('stats'))
          )
    )
  )
)

# define the server for the mobr application
server <- function(input, output) {
        
    comm <- callModule(csvFile, "comm")
    
    plot_attr <- callModule(csvFile, "plot_attr")
    
    mob_in <- reactive(make_mob_in(comm(), plot_attr()))
         
    output$mob_in <- renderPrint(mob_in())
            
    output$plot <- renderPlot({
      if(input$graphType == "Spacial Rarefaction"){
        plot_rarefaction(mob_in(), 'group', 'spat', lwd = 4, leg_loc = 'topright')
      }
      else if(input$graphType == "Individual Rarefaction - Unpooled"){
        plot_rarefaction(mob_in(), 'group', 'indiv', pooled = F, lwd = 2)
      }
      else if(input$graphType == "Individual Rarefaction - Pooled"){
        plot_rarefaction(mob_in(), 'group', 'indiv', pooled = T, lwd = 2)
      }
      else if(input$graphType == "Unpooled Abundance"){
        plot_abu(mob_in(), 'group', type = 'rad', pooled = F, log='x')
      }
      else if(input$graphType == "Pooled Abundance"){
        plot_abu(mob_in(), 'group', type = 'rad', pooled = T, log='x')
      }
      else if(input$graphType == "MoB Delta Stats"){
        delta_stats <- get_delta_stats(mob_in(), 'group',
                                       ref_group = 'uninvaded',
                                       type='discrete', log_scale=TRUE,
                                       n_perm=20)
        plot(delta_stats, 'invaded', 'uninvaded')
      }
      else if(input$graphType == "MoB Metrics - S"){
        mob_stats <- get_mob_stats(mob_in(), 'group')
        plot(mob_stats, 'S')
      }
      else if(input$graphType == "MoB Metrics - N"){
        mob_stats <- get_mob_stats(mob_in(), 'group')
        plot(mob_stats, 'N')
      }
      else if(input$graphType == "MoB Metrics - S_n"){
        mob_stats <- get_mob_stats(mob_in(), 'group')
        plot(mob_stats, 'S_n')
      }
      else if(input$graphType == "MoB Metrics - S_PIE"){
        mob_stats <- get_mob_stats(mob_in(), 'group')
        plot(mob_stats, 'S_PIE')
      }
      else if(input$graphType == "All MoB Metrics"){
        mob_stats <- get_mob_stats(mob_in(), 'group')
        plot(mob_stats, multi_panel = TRUE)
      }
})
    
    output$stats = renderPrint({
      if(input$graphType == "MoB Delta Stats"){
        invisible(capture.output(stats = get_delta_stats(mob_in(), 'group',
                                                         ref_group = 'uninvaded',
                                                         type='discrete', log_scale=TRUE,
                                                         n_perm=20)))
        #Large outputs - showing null
        
        list("Indiv_rare", stats$indiv_rare,
        "Sample_rare", stats$sample_rare,
        "SAD", stats$SAD,
        "N", stats$N,
        "Agg", stats$agg)
      }
      else if(input$graphType == "MoB Metrics - S" 
              || input$graphType == "MoB Metrics - N"
              || input$graphType == "MoB Metrics - S_n"
              || input$graphType == "MoB Metrics - S_PIE"
              || input$graphType == "All MoB Metrics"){
        invisible(capture.output(stats = get_mob_stats(mob_in(), 'group')))
        list("Group Stats", stats$groups_stats,
             "Samples Tests", stats$samples_tests,
             "Group Tests", stats$groups_tests)
      }
      
      else(
        "There are no associated stats for this model"
      )
      
      # else if(input$graphType == "MoB Metrics - N"){
      #   stats = get_mob_stats(mob_in(), 'group')
      #   stats$groups_stats
      #   stats$samples_tests
      #   stats$groups_tests
      # }
      # else if(input$graphType == "MoB Metrics - S_n"){
      #   stats = get_mob_stats(mob_in(), 'group')
      #   stats$groups_stats
      #   stats$samples_tests
      #   stats$groups_tests
      # }
      # else if(input$graphType == "MoB Metrics - S_PIE"){
      #   stats = get_mob_stats(mob_in(), 'group')
      #   stats$groups_stats
      #   stats$samples_tests
      #   stats$groups_tests
      # }
      # else if(input$graphType == "All MoB Metrics"){
      #   stats = get_mob_stats(mob_in(), 'group')
      #   stats$groups_stats
      #   stats$samples_tests
      #   stats$groups_tests
      # }
      
    })
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(input$filename, '.png', sep='')
      },
      content=function(file){
        png(file)
        print(plot_rarefaction(mob_in(), 'group', 'spat', lwd = 4, leg_loc = 'topright'))
        dev.off()
      },
      contentType='image/png')
}


# Create a Shiny app object
#shinyApp(ui = ui, server = server)

#' mobr package Graphic User Interface
#'
#' User interface of the mobr package.
#'
#' @param port char. The TCP port that the application should listen on (see
#'   \code{\link[shiny]{runApp}} for more details).
#' @param host char. The IPv4 address that the application should listen on (see
#'   \code{\link[shiny]{runApp}} for more details).
#' @param working.directory char. Directory in which the application will run.
#'
#' @return Open a window with a shiny app to use the soar package with an
#'   user-friendly interface.
#'
#' @examples
#' \dontrun{
#' gui()
#' }
#'
#' @export
gui <- function(port = getOption("shiny.port"),
                host = getOption("shiny.host", "127.0.0.1")) {

  shiny::runApp(shinyApp(ui = ui, server = server),
                display.mode = "normal", port = port, host = host)
  rm(ui, server, envir = .GlobalEnv)
}
