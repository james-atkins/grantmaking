

ge_ui <- function() {
  sidebarLayout(
    sidebarPanel(
      numericInput("ge_num_fields", "Number of fields", value = 2, min = 2, max = 4),
      sliderInput("ge_total_budget", "Total Budget", value = 1, min = 0.1, max = 2, step = 0.1),
      uiOutput("ge_field_controls")
    ),
    mainPanel(
      uiOutput("ge_plots")
    )
  )
}

# General equilibrium

make_id <- function(id, prefix = NA, suffix = NA) {
  stopifnot(length(id) == 1)
  if (!is.na(prefix)) {
    prefix <- paste0(prefix, "_")
  }

  if (!is.na(suffix)) {
    suffix <- paste0("_", suffix)
  }

  paste0(prefix, id, suffix)
}



field_controls <- function(prefix = NA, suffix = NA) {
  tagList(
    selectInput(
      make_id("type_dist", prefix, suffix),
      "Type Distribution",
      choices = structure(names(DISTRIBUTIONS), names = stringr::str_to_title(names(DISTRIBUTIONS))),
      selected = "UNIFORM"
    ),

    selectInput(
      make_id("noise_dist", prefix, suffix),
      "Noise Distribution",
      choices = structure(names(DISTRIBUTIONS), names = stringr::str_to_title(names(DISTRIBUTIONS))),
      selected = "UNIFORM"
    ),

    sliderInput(
      make_id("gamma", prefix, suffix),
      withMathJax("\\(\\gamma\\)"),
      min = 0,
      max = 1,
      value = 0.5,
      step = 0.01
    ),

    sliderInput(
      make_id("sigma", prefix, suffix),
      withMathJax("\\(\\sigma\\)"),
      min = 0,
      max = 4,
      value = 1,
      step = 0.01
    )
  )
}


ge_server <- function(input, output) {

  output$ge_field_controls <- renderUI({
    num_fields <- as.integer(input$ge_num_fields)

    map(1:num_fields, function(i) {
      tagList(
        h3(glue("Field {i}")),
        field_controls("ge", i)
      )
    })
  })

  fields <- reactive({
    num_fields <- as.integer(input$ge_num_fields)

    map(1:num_fields, function(i) {
      new_field(
        type_dist = DISTRIBUTIONS[[input[[make_id("type_dist", "ge", i)]] %||% "UNIFORM"]],
        noise_dist = DISTRIBUTIONS[[input[[make_id("noise_dist", "ge", i)]] %||% "UNIFORM"]],
        gamma = input[[make_id("gamma", "ge", i)]] %||% 0.5,
        sigma = input[[make_id("sigma", "ge", i)]] %||% 1
      )
    })
  })

  eq <- reactive({
    f <- fields()
    payline <- qpa(rep.int(1, length(f)), input$ge_total_budget)
    compute_equilibrium(f, payline)
  })

  output$ge_plots <- renderUI({
    eq <- eq()
    plots <- imap(fields(), function(field, i) {
      tagList(
        h3(glue("Field {i}")),
        renderPlot(plot_ge(eq, field, i))
      )
    })

    do.call(tagList, plots)
  })
}

plot_ge <- function(eq, field, i) {
  a <- eq$as[[i]]
  x <- eq$xs[[i]]
  p <- eq$ps[[i]]

  x_D_1 <- x_D(1, field)

  ggplot() +
    xlim(0, 1) +
    geom_function(fun = x_D, args = list(field = field), colour = "blue") +
    geom_segment(aes(x = 0, xend = 0, y = x_D(0, field), yend = x_D(0, field) + 0.2), colour = "blue") +
    geom_segment(aes(x = 1, xend = 1, y = x_D_1, yend = x_D_1 - 0.2), colour = "blue") +
    geom_function(fun = x_P, args = list(p = p, field = field), colour = "red") +
    geom_point(aes(x = a, y = x), size = 2) +
    labs(x = "Applications", y = "Standard") +
    theme_minimal()
}
