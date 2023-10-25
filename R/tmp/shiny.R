# Load the necessary libraries
library(shiny)
library(shinyalert)
library(jsonlite)
library(dplyr)
library(shinyjs)
library(ggraph)
library(igraph)

species_layout_prompt = function()
{
  sidebarLayout(
    sidebarPanel(
      textInput("name", "Name", value = 'A'),
      numericInput("methylation", "Methylation Rate", value = 0.01, min = 0),
      numericInput("demethylation", "Demethylation Rate", value = 0.01, min = 0),
      numericInput("death_rate_minus", "Death Rate (-)", value = 0.1, min = 1e-9),
      numericInput("duplication_rate_minus", "Duplication Rate (-)", value = 0.2, min = 1e-9),
      numericInput("death_rate_plus", "Death Rate (+)", value = 0.05, min = 1e-9),
      numericInput("duplication_rate_plus", "Duplication Rate (+)", value = 0.2, min = 1e-9),
      selectInput("ancestor", "Ancestor", choices = c("None")),
      selectInput("driver_event", "Driver Event", choices = c("Event A", "Event B", "Event C")),
      numericInput("time_of_origin", "Time of Origin", value = 0, min = 0),  # Add time_of_origin field
      actionButton("add_species", "Add Species"),
      actionButton("submit", "Submit")
    ),
    mainPanel(
      div(id = "message_bar"),
      tableOutput("species_table"),
      tableOutput("clones_table"),
      plotOutput("tree_plot"),
      verbatimTextOutput("json_output")
    )
  )
}

# Define the UI for the Shiny app
ui <- fluidPage(
  useShinyalert(),
  useShinyjs(),

  titlePanel("Clonal Evolution Model"),
  species_layout_prompt(),
)

# Define the server logic
server <- function(input, output, session) {
  species_data <- reactiveVal(data.frame(
    name = character(0),
    methylation = numeric(0),
    demethylation = numeric(0),
    status = character(0),
    death_rate = numeric(0),
    duplication_rate = numeric(0),
    ancestor = character(0),
    driver_event = character(0),
    time_of_origin = numeric(0)
  ))

  observeEvent(input$add_species, {
    runjs("document.getElementById('message_bar').innerHTML = '';")

    if (nchar(input$name) == 0) {
      shinyalert(
        title = "Invalid Input",
        text = "Please enter a non-empty name for the species.",
        type = "error"
      )
      return()
    }

    if (input$methylation <= 0 || input$demethylation <= 0 ||
        input$death_rate_minus <= 0 || input$duplication_rate_minus <= 0 ||
        input$death_rate_plus <= 0 || input$duplication_rate_plus <= 0) {
      shinyalert(
        title = "Invalid Input",
        text = "Please enter positive values for numeric inputs.",
        type = "error"
      )
      return()
    }

    existing_species <- species_data()
    if (input$name %in% existing_species$name) {
      shinyalert(
        title = "Species Already Exists",
        text = paste0("Species ", input$name, " already exists and will be substituted."),
        type = "warning"
      )
      existing_species <- existing_species %>% filter(name != input$name)
    }

    time_of_origin_val <- input$time_of_origin
    new_species_minus <- data.frame(
      name = input$name,
      methylation = input$methylation,
      demethylation = input$demethylation,
      status = "-",
      death_rate = input$death_rate_minus,
      duplication_rate = input$duplication_rate_minus,
      ancestor = input$ancestor,
      driver_event = input$driver_event,
      time_of_origin = time_of_origin_val
    )

    new_species_plus <- data.frame(
      name = input$name,
      methylation = input$methylation,
      demethylation = input$demethylation,
      status = "+",
      death_rate = input$death_rate_plus,
      duplication_rate = input$duplication_rate_plus,
      ancestor = input$ancestor,
      driver_event = input$driver_event,
      time_of_origin = time_of_origin_val
    )

    current_species_data <- existing_species
    updated_species_data <- rbind(current_species_data, new_species_minus, new_species_plus)
    species_data(updated_species_data)

    updateSelectInput(session, "ancestor", choices = unique(updated_species_data$name))
    # updateSelectInput(session, "time_of_origin", disabled = FALSE)


    # Increment time_of_origin by 1 for subsequent species
    if (time_of_origin_val == 0) {
      updateNumericInput(session, "time_of_origin", value = 1)
    } else {
      updateNumericInput(session, "time_of_origin", value = time_of_origin_val + 1)
    }

    runjs("document.getElementById('message_bar').innerHTML = '<div style=\"background-color: green; color: white; text-align: center;\">Species added successfully.</div>'")
  })

  output$species_table <- renderTable({
    species_data()
  })

  output$clones_table <- renderTable({
    species_data() %>%
      select(name, ancestor, driver_event, time_of_origin) %>%
      distinct()
  })

  output$tree_plot <- renderPlot({
    nsp = species_data() %>% nrow()
    if (nsp == 0) {
      print(ggplot())
    }
    else {

      tree_data <- species_data() %>%
        distinct(name, ancestor, driver_event, time_of_origin) %>%
        select(ancestor, name, driver_event, time_of_origin)

      library(ggraph)

      pch = tree_data %>%
        select(ancestor,name) %>%
        distinct()

      # print(tree_data)

      too = tree_data %>%
        select(name, time_of_origin) %>%
        distinct()

      too_v = too$time_of_origin
      names(too_v) = too$name

      # Plot call
      layout <- create_layout(pch, layout = 'tree')

      # layout$y = too_v[layout$name]
      # layout$y[1] = -3
      # layout$y = max(layout$y) - layout$y
      #
      # layout$name[1] = ' '

      # print(layout)

      col_names = layout$name %>% unique()
      col = rep("gray", length(col_names))
      names(col) = col_names
      col[' '] = 'white'

      p = ggraph::ggraph(layout) +
        theme_void() +
        ggraph::geom_edge_diagonal(
          # aes(label = name),
          # linetype = 'dashed',
          arrow = arrow(length = unit(2.5, 'mm')),
          end_cap = circle(2.5, 'mm'),
          start_cap  = circle(2.5, 'mm'),
          color = 'steelblue'
        ) +
        ggraph::geom_node_point(
          aes(color = name),
          size = 2.5) +
        ggraph::geom_node_text(
          aes(label = name),
          color = 'black',
          repel = TRUE,
          point.padding = unit(1, "lines"),
          size = 3) +
        coord_cartesian(clip = 'off') +
        scale_color_manual(values = col) +
        guides(color = 'none')

      print(p)
    }
  })

  output$json_output <- renderText({
    toJSON(species_data())
  })
}

shinyApp(ui, server)
