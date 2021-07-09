library(shiny)
library(ggplot2)

source("model.R")


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
                        "gamma",
                        min = 0,
                        max = 1,
                        value = 0.5,
                        
                    ),
                    
                    sliderInput(
                        "sigma",
                        "sigma",
                        min = 0,
                        max = 4,
                        value = 1,
                        step = 0.01
                    ),
                    
                    sliderInput(
                        "p",
                        "p",
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
        )
        
    )
)

server <- function(input, output) {
    
    output$demandSupplyPlot <- renderPlot({
        
        field <- new_field(
            type_dist = DISTRIBUTIONS[[input$type_dist]],
            noise_dist = DISTRIBUTIONS[[input$noise_dist]],
            gamma = input$gamma,
            sigma = input$sigma
        )
        
        payline <- constant_payline(input$p)
        
        ggplot() +
            xlim(0.001, 1) +
            geom_function(fun = x_D, args = list(field = field), colour = "blue") +
            geom_function(fun = x_P, args = list(p = input$p, field = field), colour = "red") +
            labs(x = "Applications", y = "Standard") +
            theme_minimal()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

