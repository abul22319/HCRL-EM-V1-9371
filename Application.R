# Monoexponeial Modeling Application
source("monoexp_model_HCRL.R")
# Load packages
packages <- c(
  "shiny",
  "readxl",
  "ggplot2",
  "dplyr",
  "minpack.lm",
  "bslib",
  "signal",
  "zoo"
)

invisible(lapply(packages, function(x){
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
}))

# UI
ui <- page_sidebar(
  
  title = div(
    style = "font-weight:700; font-size:20px;",
    "MonoExpFitLab"
  ),
  
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#1f3c88"
  ),
  
  sidebar = sidebar(
    
    width = 320,
    
    fileInput("file", "Upload Excel File", accept = c(".xlsx", ".xls")),
    
    radioButtons(
      "direction",
      "Model Type:",
      choices = c("Rise" = 1, "Decay" = 2),
      inline = TRUE
    ),
    
    tags$hr(),
    
    div(style = "font-weight:600; font-size:13px; margin-bottom:4px;",
        "Column Mapping"),
    helpText(style = "font-size:11px; margin-top:-6px;",
             "Choose which spreadsheet column belongs to each tab."),
    
    selectInput("bf_col", "Blood Flow column:", choices = NULL),
    selectInput("vc_col", "Vascular Conductance column:", choices = NULL),
    selectInput("vo2_col", "VO2 column:", choices = NULL),
    
    tags$hr(),
    
    div(style = "font-weight:600; font-size:13px; margin-bottom:4px;",
        "Filter Settings"),
    
    checkboxInput("apply_filter", "Apply Butterworth filter", value = TRUE),
    
    sliderInput(
      "cutoff",
      "Cutoff frequency (normalized, 0-1):",
      min = 0.05, max = 0.95, value = 0.3, step = 0.05
    ),
    
    radioButtons(
      "order",
      "Filter order:",
      choices = c("2" = 2, "4" = 4, "6" = 6, "8" = 8, "10" = 10),
      selected = 2,
      inline = TRUE
    ),
    
    tags$hr(),
    
    div(
      style = "font-size:11px; color:gray; text-align:center;",
      "MonoExpFitLab — v1.0"
    )
  ),
  
  navset_tab(
    nav_panel("Blood Flow",
              br(),
              plotOutput("bf_plot", height = "400px"),
              plotOutput("bf_resid", height = "300px"),
              tableOutput("bf_param"),
              tableOutput("bf_cor")
    ),
    
    nav_panel("Vascular Conductance",
              br(),
              plotOutput("vc_plot", height = "400px"),
              plotOutput("vc_resid", height = "300px"),
              tableOutput("vc_param"),
              tableOutput("vc_cor")
    ),
    
    nav_panel("VO2",
              br(),
              plotOutput("vo2_plot", height = "400px"),
              plotOutput("vo2_resid", height = "300px"),
              tableOutput("vo2_param"),
              tableOutput("vo2_cor")
    )
  )
)

# Server
server <- function(input, output, session){
  
data_reactive <- reactive({
    req(input$file)
    df <- readxl::read_excel(input$file$datapath)
    
    # Auto-generate Time if missing (assumes 2-second sampling interval)
    if(!"Time" %in% names(df)){
      df$Time <- seq(0, by = 2, length.out = nrow(df))
    }
    
    df
  })
  

  observeEvent(data_reactive(), {
    df <- data_reactive()
    col_choices <- setdiff(names(df), "Time")
    
    guess_col <- function(patterns, fallback_index = 1){
      hit <- col_choices[grepl(patterns, col_choices, ignore.case = TRUE)]
      if(length(hit) > 0) hit[1] else col_choices[min(fallback_index, length(col_choices))]
    }
    
    bf_guess  <- guess_col("blood ?flow|\\bbf\\b|\\bflow\\b", fallback_index = 1)
    vc_guess  <- guess_col("conductance|\\bvc\\b", fallback_index = 2)
    vo2_guess <- guess_col("vo2|oxygen", fallback_index = 3)
    
    updateSelectInput(session, "bf_col",  choices = col_choices, selected = bf_guess)
    updateSelectInput(session, "vc_col",  choices = col_choices, selected = vc_guess)
    updateSelectInput(session, "vo2_col", choices = col_choices, selected = vo2_guess)
  })
  
# Models: one per tab
  model_bf <- reactive({
    req(data_reactive(), input$bf_col)
    MonoExpModel(
      data_reactive(), input$bf_col, as.numeric(input$direction),
      filter = input$apply_filter,
      cutoff = input$cutoff,
      order = as.numeric(input$order)
    )
  })
  
  model_vc <- reactive({
    req(data_reactive(), input$vc_col)
    MonoExpModel(
      data_reactive(), input$vc_col, as.numeric(input$direction),
      filter = input$apply_filter,
      cutoff = input$cutoff,
      order = as.numeric(input$order)
    )
  })
  
  model_vo2 <- reactive({
    req(data_reactive(), input$vo2_col)
    MonoExpModel(
      data_reactive(), input$vo2_col, as.numeric(input$direction),
      filter = input$apply_filter,
      cutoff = input$cutoff,
      order = as.numeric(input$order)
    )
  })
  
# Output for Blood Flow
  output$bf_plot  <- renderPlot({ model_bf()$Exp.Model })
  output$bf_resid <- renderPlot({ model_bf()$RefLine.Model })
  output$bf_param <- renderTable({ model_bf()$Parameters })
  output$bf_cor   <- renderTable({ model_bf()$Cor.Result })
# Output for Vascular Conductance
  output$vc_plot  <- renderPlot({ model_vc()$Exp.Model })
  output$vc_resid <- renderPlot({ model_vc()$RefLine.Model })
  output$vc_param <- renderTable({ model_vc()$Parameters })
  output$vc_cor   <- renderTable({ model_vc()$Cor.Result })
  
# Output for VO2
  output$vo2_plot  <- renderPlot({ model_vo2()$Exp.Model })
  output$vo2_resid <- renderPlot({ model_vo2()$RefLine.Model })
  output$vo2_param <- renderTable({ model_vo2()$Parameters })
  output$vo2_cor   <- renderTable({ model_vo2()$Cor.Result })
}

# Run app
shinyApp(ui, server)