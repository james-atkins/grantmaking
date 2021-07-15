library(shiny)
library(ggplot2)
library(glue)
library(purrr)
library(assertthat)

source("model.R")
source("ge.R")


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Grantmaking"),

    tabsetPanel(
        tabPanel(
            "Constant Payline",

            sidebarLayout(
                sidebarPanel(
                    selectInput(
                        "type_dist",
                        "Type Distribution",
                        choices = structure(names(DISTRIBUTIONS), names = stringr::str_to_title(names(DISTRIBUTIONS))),
                        selected = "UNIFORM"
                    ),

                    selectInput(
                        "noise_dist",
                        "Noise Distribution",
                        choices = structure(names(DISTRIBUTIONS), names = stringr::str_to_title(names(DISTRIBUTIONS))),
                        selected = "NORMAL"
                    ),

                    sliderInput(
                        "gamma",
                        withMathJax("\\(\\gamma\\)"),
                        min = 0,
                        max = 1,
                        value = 0.5,
                        step = 0.01
                    ),

                    sliderInput(
                        "sigma",
                        withMathJax("\\(\\sigma\\)"),
                        min = 0,
                        max = 4,
                        value = 1,
                        step = 0.01
                    ),

                    sliderInput(
                        "p",
                        withMathJax("\\(p\\)"),
                        min = 0,
                        max = 1,
                        value = 0.6,
                        step = 0.01
                    )
                ),

                mainPanel(
                    plotOutput("demandSupplyPlot")
                )
            )
        ),

        tabPanel("General Equilibrium", ge_ui())

    )
)

server <- function(input, output) {

    ge_server(input, output)

    output$demandSupplyPlot <- renderPlot({
        field <- new_field(
            type_dist = DISTRIBUTIONS[[input$type_dist]],
            noise_dist = DISTRIBUTIONS[[input$noise_dist]],
            gamma = input$gamma,
            sigma = input$sigma
        )

        payline <- constant_payline(input$p)

        eq <- compute_equilibrium(list(field), payline)

        ggplot() +
            xlim(0, 1) +
            geom_function(fun = x_D, args = list(field = field), colour = "blue") +
            geom_segment(aes(x = 0, xend = 0, y = x_D(0, field), yend = x_D(0, field) + 0.2), colour = "blue") +
            geom_function(fun = x_P, args = list(p = input$p, field = field), colour = "red") +
            geom_point(data = as.data.frame(eq), aes(x = as, y = xs), size = 2) +
            labs(x = "Applications", y = "Standard") +
            theme_minimal()
    })
}

shinyApp(ui = ui, server = server)
