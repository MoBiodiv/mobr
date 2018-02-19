
# Module UI function
csvFileInput <- function(id, label = "CSV file") {
  # Create a namespace function using the provided id
  ns <- NS(id)

  tagList(
    fileInput(ns("file"), label),
    checkboxInput(ns("heading"), "Has heading"),
    selectInput(ns("quote"), "Quote", c(
      "None" = "",
      "Double quote" = "\"",
      "Single quote" = "'"
    ))
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
      header = input$heading,
      quote = input$quote,
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
    titlePanel("Uploading Files"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(

            # Input: Select a file ----
            csvFileInput("comm", "Upload community data (.csv format)"),
        
            # Horizontal line ----
            tags$hr(),

            # Input: Select a file ----
            csvFileInput("plot_attr", "Upload plot attribute data (.csv format)")
        ),
        # Main panel for displaying outputs ----
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Data",
                                 htmlOutput("mob_in")),
                        tabPanel("MoB Metrics", 
                                 plotOutput("mob_stats")),
                        tabPanel("MoB Delta Stats",
                                plotOutput("delta_stats"))
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
            
    output$mob_stats <- renderPlot({
        mob_stats <- get_mob_stats(mob_in(), 'group')
        plot(mob_stats)
    })
    
    output$delta_stats <- renderPlot({
        delta_stats <- get_delta_stats(mob_in(), 'group',
                                       ref_group = 'uninvaded',
                                       type='discrete', log_scale=TRUE,
                                       n_perm=20)
        plot(delta_stats, 'invaded', 'uninvaded')
    })
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