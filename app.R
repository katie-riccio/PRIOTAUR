# ///////////////////////////////////////// ----

# initial setup ----

## libraries ----
library(bslib)
library(echo)
library(ggplot2)
library(ggetho)
library(ggstats)
library(shiny)
library(shinycssloaders)
library(rain)
library(readxl)
library(tibble)
library(tidyverse)
library(viridis)
library(writexl)

## loading files ----
nuc_proteome <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM4_ESM.xlsx", sheet = 2) # nuclear-encoded proteome data
org_proteome <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM4_ESM.xlsx", sheet = 3) # organellar-encoded proteome data
rhythm_pvalues <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM5_ESM.xlsx") # pre-calculated eJTK/rain/echo LD and LL p-values

head(nuc_proteome) # check nuclear proteome
head(org_proteome) # check organellar proteome
head(rhythm_pvalues) # check rhythmicity values

nuc_abundance_LD <- select(nuc_proteome, "Identifier", "description", "0":"24.5") # extract protein abundance values at each LD time point
nuc_abundance_LL <- select(nuc_proteome, "Identifier", "description", "48":"121.5") # extract protein abundance values at each LL time point

# ion abundance values
ions_table <- read_excel("ions.xlsx", range = "A3:M41", col_names = c(
  "Time",
  "Mg_1", "Mg_2", "Mg_3", "Mg_4",
  "Ca_1", "Ca_2", "Ca_3", "Ca_4",
  "K_1", "K_2", "K_3", "K_4"
)) 


# ///////////////////////////////////////// ----
# UI ----
# ///////////////////////////////////////// ----

ui <- fluidPage(
  #theme = bs_theme(bootswatch = "flatly"),
  mainPanel(
    tabsetPanel(
      id = "tabset",

      ## Ionomics ----
      tabPanel(
        title = "Ionome",
        sidebarLayout(
          
          ### Sidebar ----
          sidebarPanel(
            style = "position: fixed; overflow-y: auto",
            checkboxGroupInput(
              input = "ions", label = "Select ions", selected = NULL,
              choices = c("Ca2+" = "Ca", "K+" = "K", "Mg2+" = "Mg")
            ),
            radioButtons("abs_rel", label = "Plot abundance as:",
                         choices = c("Absolute" = "abs", "Relative" = "rel")),
            checkboxInput(input = "ion_raw", label = "Show all values?")
            
          ),
          
          ### Main page ----
          mainPanel(
            
            conditionalPanel(
              condition = "output.chosen_ions == ''",
              p("Select one or more ions to get started..")
            ),
            
            # p("Show all values?", textOutput("ion_show_raw")),
            
            conditionalPanel(
              condition = "input.ions != ''",
              # p("test message - will display if graph can display"),
              # p("You have selected:", textOutput("chosen_ions")),
              withSpinner(plotOutput("ion_graph", width = 600, height = 300)),
              downloadButton("ion_graph_download", "Download graph view"),
              
              #p("Source: Gil and Crosby, unpublished")
            
            conditionalPanel(
              condition = "output.ion_show_raw == 'TRUE'",
              downloadButton("ion_table_download", "Download table of values"),
              withSpinner(tableOutput("ion_table")),
            )
            ),
          )
        )
      ),

      ## Proteomics ----
      tabPanel(
        title = "Proteome",
        sidebarLayout(
          ### Sidebar ----
          sidebarPanel(
            selectizeInput(input = "protein", label = "Select a protein", choices = NULL),
            conditionalPanel(
              condition = "input.protein != ''",
              checkboxGroupInput(
                input = "rhythm_method", label = "Show rhythmicity p-value?", selected = NULL,
                choices = c("eJTK" = "ejtk", "RAIN" = "rain", "ECHO" = "echo")
              )
            )
          ),
          ### Main page ----
          mainPanel(
            conditionalPanel(
              condition = "input.protein == ''",
              p("Select a protein to get started..")
            ),
            conditionalPanel(
              condition = "input.protein != ''",
              withSpinner(plotOutput("abundance_graph", width = 600, height = 300)),
              p("Source: Kay et al. Comms Bio 2021")
            ),
            conditionalPanel(
              condition = "input.rhythm_method != ''",
              p("Rhythmicity p-values:", textOutput("protein_with_desc")),
              tableOutput("rhythm_pvalues"),
              p("Source: Kay et al. Comms Bio 2021")
            )
          )
        )
      ),

      ## Table display page ----
      tabPanel(
        title = "Comparison",
        mainPanel(
          conditionalPanel(
            condition = "input.protein == '' OR input.chosen_ions == ''",
            p("Select ions and proteins to get started..")
          ),
          conditionalPanel(
            condition = "input.protein != ''",
            tableOutput("protein_table"),
            p("Source: Kay et al. Comms Bio 2021")
          )
        )
      )
    )
  )
)

# ///////////////////////////////////////// ----
# SERVER ----
# ///////////////////////////////////////// ----

server <- function(input, output, session) {

  # Ionomics ----
  
  
  
  ## Get values for selected ions to go into graph ----
  ion_data <- reactive({
    req(input$ions)
    ion_data <- select(ions_table, "Time", contains(input$ions)) # create table with all necessary columns
    ion_data <- pivot_longer(ion_data, cols = !Time, names_to = "Ion", values_to = "Abundance") %>% drop_na()
    ion_data$Ion <- str_remove(ion_data$Ion, "\\_.*") # clean up ion names
    ion_data <- mutate(ion_data,
                       Condition = case_when(ion_data$Time <= 24 ~ "LD", ion_data$Time >= 48 ~ "LL"),
                       .after = Time
    ) # add condition column
    ion_data %>%
      group_by(Time, Condition, Ion) %>%
      summarise(across(), mean = mean(Abundance), sd = sd(Abundance)) %>% 
      group_by(Ion) %>% mutate(rel_abun = mean/max(mean)) # add other columns
  })
  
  output$ion_show_raw <- renderText({
    paste(input$ion_raw)
    })
  outputOptions(output, "ion_show_raw", suspendWhenHidden = FALSE)
  
  # get values for selected ions as an output? (may not need)
  output$ion_data <- reactive({
    req(input$ions)
    ion_data()
  })
  
  # show values for selected ions
  output$ion_table <- renderTable({
    req(input$ions)
    ion_data()
    
    
    
  })
  
  # show which ions are selected
  output$chosen_ions <- renderText({
    paste(input$ions)
  })
  
  ## Graph of selected ion abundance changes over time ----
  output$ion_graph <- renderPlot({
    ggplot(data = ion_data(), aes(x = Time, y = rel_abun, group = Ion)) +
      ggtitle("Mean ion abundance") +
      geom_stripped_cols(width = 12, nudge_x = 6, odd = "#e5e5e5", even = "#ffffff", colour = "grey") +
      #geom_stripped_cols(width = 12, nudge_x = 6, odd = "darkgray", even = "#ffffff", colour = "black") +
      geom_line(aes(y = rel_abun, group = interaction(Ion, Condition))) +
      #geom_errorbar(aes(y = mean, ymin = mean - sd, ymax = mean + sd)) +
      geom_point(size = 2, aes(fill = Ion), shape = 21, colour = "black") +
      
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 12)) +
      labs(x = "Time (h)", y = "Ion abundance (ug/L)") +
      scale_x_continuous(breaks = scales::breaks_width(12), minor_breaks = NULL) +
      # scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96), minor_breaks = NULL)
      scale_colour_brewer(palette = "Set1")
  })
  
  ## Download table values for currently-selected ions ----
  ion_table_download <- downloadHandler(
    filename = function() {paste0(input$chosen_ions, ".csv")},
    content = function(file) {write.csv(ion_data(), file)}
  )
  
  ## Download current graph view ----
  ion_graph_download <- downloadHandler(
    filename = "graph.png",
    content = function(file) {
      ggsave(file, plot = ion_graph())
    }
  )
  
  # Proteomics ----
  
  ## select protein ----
  updateSelectizeInput(session, "protein", choices = c("Select a protein" = "", nuc_proteome["Identifier"]), server = TRUE)

  # get selected protein from list
  protein_with_desc <- reactive({
    req(input$protein)
    paste(select(nuc_proteome, Identifier, description) %>% filter(Identifier == input$protein), collapse = " - ")
  })

  # list selected protein
  output$protein_with_desc <- reactive({
    req(input$protein)
    protein_with_desc()
  })

  ## get relevant values from proteome spreadsheet based on selected protein ----
  protein_data <- reactive({
    req(input$protein)

    protein_data_LD <- (subset(nuc_abundance_LD, Identifier == input$protein)
    ) %>%
      pivot_longer(
        cols = !Identifier, names_to = "Time", names_transform = list("Time" = as.numeric),
        values_to = "Abundance", values_transform = list("Abundance" = as.numeric)
      ) %>%
      add_column(Condition = "LD")

    protein_data_LL <- (subset(nuc_abundance_LL, Identifier == input$protein)
    ) %>%
      pivot_longer(
        cols = !Identifier, names_to = "Time", names_transform = list("Time" = as.numeric),
        values_to = "Abundance", values_transform = list("Abundance" = as.numeric)
      ) %>%
      add_column(Condition = "LL")

    protein_data <- rbind(protein_data_LD, protein_data_LL) %>% drop_na() # merge two tables together
    
    
  })


  # show table of values for selected protein
  output$protein_table <- renderTable({
    protein_data()
  })


  ## generate a graph for selected protein ----
  output$abundance_graph <- renderPlot({
    ggplot(data = protein_data(), aes(x = Time, y = Abundance, group = 1)) +
      ggtitle(protein_with_desc()) +
      # geom_stripped_cols(aes(group = "LD"), width = 12, nudge_x = 6, odd = "#ffffff", even = "#a5a5a5") +
      geom_stripped_cols(aes(group = "LL"), width = 12, nudge_x = 6, odd = "#e5e5e5", even = "#ffffff") +
      geom_line(aes(group = Condition)) +
      geom_point(size = 2, color = "red") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 12)) +
      labs(x = "Time (h)", y = "Protein abundance") +
      scale_x_continuous(breaks = scales::breaks_width(12), minor_breaks = NULL)
  })

  # rhythmicity calculations ----
  ## rain analysis ----

  # rain_data_LD <- #to add
  # rain_analysis_LD <- rain(nuc_abundance_LD, period = 24, period.delta = 0, method = "independent", deltat = 3.5, nr.series = 1)
  #
  #
  # rain_data_LL <- (subset(nuc_abundance_LL, Identifier == input$protein)
  # ) %>% pivot_longer(
  #   cols = everything(), names_to = "Time", names_transform = list("Time" = as.numeric),
  #   values_to = "Abundance", values_transform = list("Abundance" = as.numeric)
  # ) %>% drop_na() %>% as.matrix()
  #
  # rain_analysis_LL <- rain(rain_data_LL, period = 24, period.delta = 4, method = "longitudinal",
  #                          deltat = 3.5, nr.series = 1)

  ## echo analysis ----

  # echo_data_LD <- (subset(nuc_abundance_LD, select = -c(description)))
  # echo_analysis_LD <- echo_find(nuc_abundance_LD, low = 24, high = 24, is_de_linear_trend = FALSE, begin = 0, end = 24.5, resol = 3.5, num_reps = 1)
  #
  # echo_data_LL <- (subset(nuc_abundance_LL, select = -c(description)))
  # echo_analysis_LL <- echo_find(nuc_abundance_LL, low = 20, high = 28, is_de_linear_trend = TRUE, begin = 0, end = 24.5, resol = 3.5, num_reps = 1)

  ## show rhythmicity ----
  output$rhythmicity <- renderTable(
    rain_analysis_LD
  )

  # get rhythmicity p-values
  output$rhythm_pvalues <- renderTable(
    {
      req(input$rhythm_method)
      protein_rhythm_pvalues <- (filter(rhythm_pvalues, Identifier == input$protein)
      ) %>% select(contains(input$rhythm_method))
    },
    width = "600px",
    digits = 5
  )

  # Misc. ----

  output$page <- renderText({
    paste("Currently on", input$tabset)
  })
}

# ///////////////////////////////////////// ----

# Run Shiny app ----
shinyApp(ui, server)
#runApp(display.mode = "showcase")

# ///////////////////////////////////////// ----