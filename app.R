
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

    radioButtons(
      "model_type",
      "Number of Components:",
      choices = c("1" = 1, "2" = 2, "3" = 3),
      selected = 1,
      inline = TRUE
    ),

    tags$hr(),

    div(style = "font-weight:600; font-size:13px; margin-bottom:4px;",
        "Column to Model"),
    helpText(style = "font-size:11px; margin-top:-6px;",
             "Choose which spreadsheet column to fit."),

    selectInput("column", "Column:", choices = NULL),

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

  # Single output area (no tabs)
  br(),
  plotOutput("fit_plot",   height = "400px"),
  plotOutput("fit_resid",  height = "300px"),
  tableOutput("fit_param"),
  tableOutput("fit_cor")
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

  # Populate the column selector whenever a file is loaded
  observeEvent(data_reactive(), {
    df <- data_reactive()
    col_choices <- setdiff(names(df), "Time")

    updateSelectInput(
      session, "column",
      choices  = col_choices,
      selected = col_choices[1]
    )
  })

  # Single model for the selected column
  model_fit <- reactive({
    req(data_reactive(), input$column)
    MonoExpModel(
      data_reactive(), input$column, as.numeric(input$direction),
      n_comp = as.numeric(input$model_type),
      filter = input$apply_filter,
      cutoff = input$cutoff,
      order  = as.numeric(input$order)
    )
  })

  # Outputs
  output$fit_plot  <- renderPlot({ model_fit()$Exp.Model })
  output$fit_resid <- renderPlot({ model_fit()$RefLine.Model })
  output$fit_param <- renderTable({ model_fit()$Parameters })
  output$fit_cor   <- renderTable({ model_fit()$Cor.Result })
}

# Run app
shinyApp(ui, server)
