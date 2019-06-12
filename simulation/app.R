#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny);library(shinythemes);library(dplyr);library(purrr);library(MASS);library(ggplot2)
source("sim_functions.R")
results = readRDS("../results/03012019_calibration_sim_results.rds")
results = results %>% filter(
  !(y_true == "x" & y_spec %in% c("x2","z-x","z-x2","z-x2-full")),
  !(y_true == "x2" & y_spec %in% c("z-x2","z-x2-full")),
  !(y_true == "z-x" & y_spec %in% c("z-x2","z-x2-full")),
  !(y_true == "z-x2" & y_spec %in% c("z-x2-full"))
)

ui <- fluidPage(
  
  # Application title
  #titlePanel("Calibrating Validation Samples: Simulation"),
  withMathJax(),
  navbarPage("Calibrating Validation Samples: Simulation Results",
             id = "tabs",
             theme = shinytheme("flatly"),
             tabPanel("Results",
                      sidebarPanel(
                        checkboxGroupInput("y_spec", "Specified Y model",
                                     choices = list("\\(Y \\sim Z\\)" = "simple",
                                                    "\\(Y \\sim Z + X\\)" = "x",
                                                    "\\(Y \\sim Z + X + X^2\\)" = "x2",
                                                    "\\(Y \\sim Z + X + ZX\\)" = "z-x",
                                                    "\\(Y \\sim Z + X + ZX + X^2\\)" = "z-x2",
                                                    "\\(Y \\sim Z + X + ZX + X^2 + ZX^2\\)" = "z-x2-full"
                                     ), selected = "simple"),
                        radioButtons("b0","\\(\\beta_0\\)", 
                                     choices = list("0" = 0,
                                                    "1" = 1,
                                                    "5" = 5,
                                                    "10" = 10,
                                                    "20" = 20
                                     )),
                        radioButtons("g0","\\(\\gamma_0\\)", 
                                     choices = list("0" = 0,
                                                    "1" = 1,
                                                    "5" = 5,
                                                    "10" = 10,
                                                    "20" = 20
                                     ))
                      ),
                      # Show a plot of the generated distribution
                      mainPanel(
                        plotOutput("bias_plot"),
                        plotOutput("dom_plot")
                      ))
  ))

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$bias_plot = renderPlot({plot_bias(results, input$b0, input$g0, input$y_spec)})
  
  output$dom_plot = renderPlot({plot_dom(results, input$b0, input$g0, input$y_spec)})
}

# Run the application 
shinyApp(ui = ui, server = server)

