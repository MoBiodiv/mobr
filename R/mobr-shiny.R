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

ui <- dashboardPage(
  dashboardHeader(title = "MoB-R"),

  
  dashboardSidebar(
    sidebarMenu(style = "position: fixed; overflow: visible;",
      menuItem("Home", tabName = "Home", icon = icon("home")),
      menuItem("Data Tab", tabName = "DataTab", icon = icon("table")),
      menuItem("Plot Rarefaction", tabName = "plot_rarefaction", icon = icon("angle-right")),
      menuItem("MoB Metrics", tabName = "mob_metrics", icon = icon("angle-right"),
        menuSubItem("All MoB Metrics", tabName = "all_mob_tab", icon = icon("angle-right")),
        menuSubItem("Individual MoB Metrics", tabName = "ind_mob_tab", icon = icon("angle-right"))),
      menuItem("Delta Stats", tabName = "delta_stats", icon = icon("angle-right")),
      menuItem("MoB-R GitHub Page", icon = icon("github"), 
               href = "https://github.com/MoBiodiv/mobr"),
      menuItem("Got an Issue? Submit it!", icon = icon("exclamation-triangle"), 
               href = "https://github.com/MoBiodiv/mobr/issues")
    )
  ),
  dashboardBody(
    
    
    tabItems(
      # First tab content
      tabItem(tabName = "Home",
              h2("Welcome to the MoB-R App!"),
              h4("NOTE: Graphics work best if viewed in full window OR when side panel is closed"),
              fluidRow(
              box(
                hieght = "300px",
                width = 6,
                h3(icon("globe"), "Navigation Tabs"),
                #tags$b(h3("Navigation")),
                br(),
                h4(icon("table"),"DataTab"),
                "- Where you will enter you Community Matrix and Plot Data as CSV files", 
                br(), 
                br(), 
                h4(icon("angle-right"),"Plot Rarefaction"),
                "- Graphs and Code for Plot Rarefaction data",
                br(),
                br(), 
                h4(icon("angle-right"),"MoB Metrics"),
                "- A pull down tab containing 'All MoB Metrics' and 'Individual MoB Metrics'", 
                br(), 
                "- Stats for both can be found in the 'Individual MoB Metrics' tab",
                br(),
                br(), 
                h4(icon("angle-right"),"Delta Stats"),
                "- Graph and Code for Delta Stats data"
              ),
              box(
                width = 5,
                h3(icon("link"), "Link Tabs"),
                #tags$b(h3("Link Tabs")),
                br(),
                "The last two tabs represent links to the MoB-R Github page",
                br(),
                h3(icon("github")),
                "If you want to learn more about MoB-R fromthe source check out the first link. It will take you to MoB-R's main GitHub page.",
                br(),
                br(),
                h4(icon("exclamation-triangle")),
                "If you run into a problem while using this app, please let the developers know! Use the second link to navigate to the issues page on the MoB-R GitHub.",
                br(),
                "There you can submit an issue describing your problem. We appreciate your feedback!"
              )
              ),
              fluidRow(
                box(
                  width = 4,
                  h3(icon("users"), "The people behind MoB-R:"),
                  tags$b("Author and Creator: "),
                  br(),
                  "Xiao Xiao - email: xiao@weecology.org",
                  br(),
                  tags$b("Authors: "),
                  br(),
                  "Daniel McGlinn - email: danmcglinn@gmail.com",
                  br(),
                  "Felix May - email: felix.may@idiv.de",
                  br(),
                  "Caroline Oliver - email: olivercs@g.cofc.edu"
                )
               #box(
               #   imageOutput("image")
               # )
                
              )
      ),
      
      tabItem(tabName = "DataTab",
              h2("Enter your data below in CSV file format:"),
              # Input: Select a file ----
              csvFileInput("comm", "Upload community data"),
              
              # Horizontal line ----
              tags$hr(),
              
              # Input: Select a file ----
              csvFileInput("plot_attr", "Upload plot attribute data")
              
      ),
      
      # Second tab content
      tabItem(tabName = "plot_rarefaction",
              h2("Plot Rarefaction"),
              h4("SR = Spacial Rarefaction, IR = Indiviudal Rarefaction, Abu = Abundance"),
              fluidRow(
                tabBox(
                  title = "Plot Rarefaction",
                  side = "left", width = "10",
                  selected = "SR",
                  tabPanel("SR", withSpinner(plotOutput('s_rare'))),
                  tabPanel("IR Unpooled", withSpinner(plotOutput('i_rare_up'))),
                  tabPanel("IR Pooled", withSpinner(plotOutput('i_rare_p'))),
                  tabPanel("Unpooled Abu", withSpinner(plotOutput('up_abu'))),
                  tabPanel("Pooled Abu", withSpinner(plotOutput('p_abu')))
                )
              ),
              fluidRow(
                tabBox(
                  title = "Plot Rarefaction Code",
                  side = "left", width = "10",
                  selected = "SR",
                  tabPanel("SR", withSpinner(htmlOutput('s_rare_code'))),
                  tabPanel("IR Unpooled", withSpinner(htmlOutput('ir_up_code'))),
                  tabPanel("IR Pooled", withSpinner(htmlOutput('ir_p_code'))),
                  tabPanel("Unpooled Abu", withSpinner(htmlOutput('up_abu_code'))),
                  tabPanel("Pooled Abu", withSpinner(htmlOutput('p_abu_code')))
                )
                
              )
      ),
      
      tabItem(tabName = "all_mob_tab",
              h2("All MoB Metrics"),
              fluidRow(
                box(
                  width = 7,
                  height = "750px",
                  title = "All MoB Metrics",
                  withSpinner(plotOutput('mob_all', height = "695px"))
                ),
                fluidRow(
                  box(
                    width = 7,
                    title = "All MoB Metrics Code",
                    withSpinner(htmlOutput('all_mob_code'))
                  )
                )
              )
        ),
              
      
      tabItem(tabName = "ind_mob_tab",
              h2("MoB Metrics"),
              fluidRow(
                tabBox(
                  title = "Individual MoB Metrics",
                  side = "left", width = "10",
                  selected = "S",
                  tabPanel("S", withSpinner(plotOutput('mob_s'))),
                  tabPanel("N", withSpinner(plotOutput('mob_n'))),
                  tabPanel("S_n", withSpinner(plotOutput('mob_Sn'))),
                  tabPanel("S_PIE", withSpinner(plotOutput('mob_Spie')))
                )
              ),
              fluidRow(
                box(
                  title = "Groups Stats",
                  width = 4,
                  withSpinner(verbatimTextOutput('mob_groups_stats'))
                ),
                box(
                  title = "Samples Tests",
                  width = 4,
                  withSpinner(verbatimTextOutput('mob_samples_tests'))
                ),
                box(
                  title = "Groups Tests",
                  width = 4,
                  withSpinner(verbatimTextOutput('mob_groups_tests'))
                )
              ),
              fluidRow(
                tabBox(
                  title = "MoB Metrics Code",
                  side = "left", width = "10",
                  selected = "S",
                  tabPanel("S", withSpinner(htmlOutput('s_code'))),
                  tabPanel("N", withSpinner(htmlOutput('n_code'))),
                  tabPanel("S_n", withSpinner(htmlOutput('sn_code'))),
                  tabPanel("S_PIE", withSpinner(htmlOutput('spie_code')))
                )
                
              )
      ),
      
      # Forth tab content
      tabItem(tabName = "delta_stats",
              h2("Delta Stats"),
              fluidRow(
                box(
                  title = "Delta Stats",
                  width = 8,
                  withSpinner(plotOutput('delta_plot'))
                  
                )
              ),
              fluidRow(
                box(
                  title = "Delta Stats Code",
                  width = 8,
                  withSpinner(htmlOutput('delta_code'))
                )
              )
      )
    )
    )
)
  



server <- function(input, output) {

  comm <- callModule(csvFile, "comm")
  
  plot_attr <- callModule(csvFile, "plot_attr")
  
  mob_in <- reactive(make_mob_in(comm(), plot_attr()))

  mob_stats <- reactive(get_mob_stats(mob_in(), 'group'))
  
  output$mob_s <- renderPlot({
      #mob_stats <- get_mob_stats(mob_in(), 'group')
      plot(mob_stats(), 'S')
  })
  output$mob_n <- renderPlot({
    #mob_stats <- get_mob_stats(mob_in(), 'group')
    plot(mob_stats(), 'N')
  })
  output$mob_Sn <- renderPlot({
    #mob_stats <- get_mob_stats(mob_in(), 'group')
    plot(mob_stats(), 'S_n')
  })
  output$mob_Spie <- renderPlot({
    #mob_stats <- get_mob_stats(mob_in(), 'group')
    plot(mob_stats(), 'S_PIE')
  })
  output$mob_all <- renderPlot({
    #mob_stats <- get_mob_stats(mob_in(), 'group')
    plot(mob_stats(), multi_panel = TRUE)
  })
  
  
  output$mob_groups_stats = renderPrint({
    #invisible(capture.output(stats = get_mob_stats(mob_in(), 'group')))
    mob_stats()$groups_stats
  })
  
  output$mob_samples_tests = renderPrint({
    #invisible(capture.output(stats = get_mob_stats(mob_in(), 'group')))
    mob_stats()$samples_tests
  })
  
  output$mob_groups_tests = renderPrint({
    #invisible(capture.output(stats = get_mob_stats(mob_in(), 'group')))
    mob_stats()$groups_tests
  })
  
  output$s_code <- renderText({
    paste("mob_in <- make_mob_in(community_matrix_FILENAMEHERE, plot_attribute_FILENAMEHERE)",
          "<br>", 
          "mob_stats <- get_mob_stats(mob_in(), 'group')",
          "<br>",
          "plot(mob_stats, 'S')")
  })
  output$n_code <- renderText({
    paste("mob_in <- make_mob_in(community_matrix_FILENAMEHERE, plot_attribute_FILENAMEHERE)",
          "<br>", 
          "mob_stats <- get_mob_stats(mob_in(), 'group')",
          "<br>",
          "plot(mob_stats, 'N')")
  })
  output$sn_code <- renderText({
    paste("mob_in <- make_mob_in(community_matrix_FILENAMEHERE, plot_attribute_FILENAMEHERE)",
          "<br>", 
          "mob_stats <- get_mob_stats(mob_in(), 'group')",
          "<br>",
          "plot(mob_stats, 'S_n')")
  })
  output$spie_code <- renderText({
    paste("mob_in <- make_mob_in(community_matrix_FILENAMEHERE, plot_attribute_FILENAMEHERE)",
          "<br>", 
          "mob_stats <- get_mob_stats(mob_in(), 'group')",
          "<br>",
          "plot(mob_stats, 'S_PIE')")
  
  })
  
  output$all_mob_code <- renderText({
    paste("mob_in <- make_mob_in(community_matrix_FILENAMEHERE, plot_attribute_FILENAMEHERE)",
          "<br>", 
          "mob_stats <- get_mob_stats(mob_in(), 'group')",
          "<br>",
          "plot(mob_stats, multi_panel = TRUE)")
    
  })
  
  output$s_rare <- renderPlot({
    plot_rarefaction(mob_in(), 'group', 'spat', lwd = 4, leg_loc = 'topright')
  })
  
  output$i_rare_up <- renderPlot({
    plot_rarefaction(mob_in(), 'group', 'indiv', pooled = F, lwd = 2)
  })
  
  output$i_rare_p <- renderPlot({
    plot_rarefaction(mob_in(), 'group', 'indiv', pooled = T, lwd = 2)
  })
  
  output$up_abu <- renderPlot({
    plot_abu(mob_in(), 'group', type = 'rad', pooled = F, log='x')
  })
  
  output$p_abu <- renderPlot({
    plot_abu(mob_in(), 'group', type = 'rad', pooled = T, log='x')
  })
  
  
  output$s_rare_code <- renderText({
    paste("mob_in <- make_mob_in(community_matrix_FILENAMEHERE, plot_attribute_FILENAMEHERE)",
          "<br>", 
          "plot_rarefaction(mob_in(), 'group', 'spat', lwd = 4, leg_loc = 'topright')")
  })
  output$ir_up_code <- renderText({
    paste("mob_in <- make_mob_in(community_matrix_FILENAMEHERE, plot_attribute_FILENAMEHERE)",
          "<br>", 
          "plot_rarefaction(mob_in(), 'group', 'indiv', pooled = F, lwd = 2)")
  })
  output$ir_p_code <- renderText({
    paste("mob_in <- make_mob_in(community_matrix_FILENAMEHERE, plot_attribute_FILENAMEHERE)",
          "<br>", 
          "plot_rarefaction(mob_in(), 'group', 'indiv', pooled = T, lwd = 2)")
  })
  output$up_abu_code <- renderText({
    paste("mob_in <- make_mob_in(community_matrix_FILENAMEHERE, plot_attribute_FILENAMEHERE)",
          "<br>", 
          "plot_abu(mob_in(), 'group', type = 'rad', pooled = F, log='x')")
    
  })
  output$p_abu_code <- renderText({
    paste("mob_in <- make_mob_in(community_matrix_FILENAMEHERE, plot_attribute_FILENAMEHERE)",
          "<br>", 
          "plot_abu(mob_in(), 'group', type = 'rad', pooled = T, log='x')")
    
  })
  
  
  output$delta_plot <- renderPlot({
    delta_stats <- get_delta_stats(mob_in(), 'group',
                                   ref_group = 'uninvaded',
                                   type='discrete', log_scale=TRUE,
                                   n_perm=20)
    plot(delta_stats, 'invaded', 'uninvaded')
    
  })
  
  output$delta_code <- renderText({
    paste("mob_in <- make_mob_in(community_matrix_FILENAMEHERE, plot_attribute_FILENAMEHERE)",
          "<br>", 
          "delta_stats <- get_delta_stats(mob_in, 'group', ref_group = 'uninvaded', type='discrete', log_scale=TRUE, n_perm=20)",
          "<br>",
          "plot(delta_stats, 'invaded', 'uninvaded')")
    
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
