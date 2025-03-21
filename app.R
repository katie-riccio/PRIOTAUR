# ///////////////////////////////////////// ----

# initial setup ----

## libraries ----
library(bslib)
library(echo)
library(ggplot2)
# library(ggetho)
library(ggbreak)
library(ggstats) # for stripped columns
library(shiny)
library(shinycssloaders)
library(plotly)
library(rain)
library(readxl)
# library(tibble)
library(tidyverse)
library(viridis)
library(writexl)

## raw files ----
nuc_proteome <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM4_ESM.xlsx", sheet = 2, 
                           na = "ND", .name_repair = function(x) gsub("\\s+", "_", x)) # nuclear-encoded proteome data
org_proteome <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM4_ESM.xlsx", sheet = 3,
                           na = "ND", .name_repair = function(x) gsub("\\s+", "_", x)) # organellar-encoded proteome data
protein_rhythm_pvalues <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM5_ESM.xlsx",
                                     na = "ND") %>% mutate(across(!Identifier, as.numeric))  # pre-calculated eJTK/rain/echo LD and LL p-values (names don't need fixing here)

nuc_proteome # check nuclear proteome
org_proteome # check organellar proteome
protein_rhythm_pvalues # check rhythmicity values

ion_abundance <- read_excel("ions.xlsx", range = "A3:M41", col_names = c(
  "Time",
  "Mg_1", "Mg_2", "Mg_3", "Mg_4",
  "Ca_1", "Ca_2", "Ca_3", "Ca_4",
  "K_1", "K_2", "K_3", "K_4"
)) # ion abundance values - columns aren't named in the actual raw file so that's done here

ion_abundance # check ion abundance values

ion_rhythmicity <- "ion_rhythmicity.xlsx"
ion_rhythmicity <- ion_rhythmicity %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = ion_rhythmicity) %>% 
  list2env(.GlobalEnv) # load in pre-calculated RAIN and ECHO outputs from excel sheet


## one-time-only processes ----

protein_abundance <- select(nuc_proteome, "Identifier", "0":"24.5", "48":"121.5") # extract necessary columns for abundance graphs

double_rhythmics <- filter(protein_rhythm_pvalues, LL_rain < 0.05 & LL_echo < 0.05) # find proteins rhythmic in both RAIN and ECHO in LL
double_rhythmics # preview outcome

double_rhythmics_full <- filter(nuc_proteome, Identifier %in% double_rhythmics$Identifier) # get full information for those proteins
double_rhythmics_full # preview outcome

mg_co_rhythmics <- read_excel("mg_co_rhythmics.xlsx") # load in list of rhythmic proteins that use Mg as a co-factor
mg_co_rhythmics

# rhythmic_ion_binders_full <- filter(double_rhythmics_full, 
#                                Identifier %in% mg_co_rhythmics$Identifier) %>%
#   select(Identifier, Amino_acid_sequence, description, Absolute_Phase_LL, 
#          Circadian_Phase_LL, protein_domains, everything()) %>% 
#   arrange(Absolute_Phase_LL) # get and reorder full table for the 16

rhythmic_ion_binders_full <- filter(double_rhythmics_full, 
                                    Identifier %in% mg_co_rhythmics$Identifier) # just get 16

rhythmic_ion_binders_full # preview outcome


# ///////////////////////////////////////// ----
# UI ----
# ///////////////////////////////////////// ----

ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly"),
  mainPanel(
    textOutput("page"),
    tabsetPanel(
      id = "tabset",
    
      ## Ion visualisation ----
      tabPanel(
        title = "Ion visualisation",
        sidebarLayout(
          ### Sidebar ----
          sidebarPanel(
            # style = "position: fixed; overflow-y: auto",
            checkboxGroupInput(
              input = "ions", label = "Select ions:", selected = NULL,
              choices = c("Calcium" = "Ca", "Potassium" = "K", "Magnesium" = "Mg")
            ),
            actionButton("ion_go", label = "Go!"),
            radioButtons("ion_abs_rel",
              label = "Plot abundance as:",
              choices = c("Absolute" = "abs", "Relative" = "rel")
            ),
            checkboxInput(input = "ion_raw", label = "View raw values?"),
            radioButtons("ion_rain_echo",
              label = "Show LL values from:",
              choices = c("RAIN" = "rain", "ECHO" = "echo")
            )
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
              conditionalPanel(
                condition = "input.ion_abs_rel == 'abs'",
                withSpinner(plotOutput("ion_graph_abs", width = 700, height = 350))
              ),
              conditionalPanel(
                condition = "input.ion_abs_rel == 'rel'",
                withSpinner(plotOutput("ion_graph_rel", width = 700, height = 350))
              ),
              # withSpinner(plotOutput("ion_graph", width = 600, height = 300)),
              downloadButton("ion_graph_download", "Download graph view"),
    
              # p("You have selected:", textOutput("ion_rain_echo")),
              tableOutput("ion_rhythm_pvalues"),
              conditionalPanel(
                condition = "output.ion_show_raw == 'TRUE'",
                downloadButton("ion_table_download", "Download table of values"),
                withSpinner(tableOutput("ion_table"))
              )
            )
          )
        )
      ),
    
      ## Protein visualisation ----
      tabPanel(
        title = "Protein visualisation",
        sidebarLayout(
          ### Sidebar ----
          sidebarPanel(
            selectizeInput(
              input = "proteins", label = "Select up to 5 proteins:", choices = NULL,
              multiple = TRUE, options = list(maxItems = 5)
            ),
            radioButtons("protein_abs_rel",
              label = "Plot abundance as:",
              choices = c("Absolute" = "abs", "Relative" = "rel")
            ),
            conditionalPanel(
              condition = "input.proteins != ''",
              checkboxGroupInput(
                input = "rhythm_method", label = "Show rhythmicity p-value?", selected = NULL,
                choices = c("eJTK" = "ejtk", "RAIN" = "rain", "ECHO" = "echo")
              )
            )
          ),
          ### Main page ----
          mainPanel(
            conditionalPanel(
              condition = "input.proteins == ''",
              p("Select a protein to get started..")
            ),
            conditionalPanel(
              condition = "input.proteins != ''",
              conditionalPanel(
                condition = "input.protein_abs_rel == 'abs'",
                withSpinner(plotOutput("protein_graph_abs", width = 700, height = 350))
              ),
              conditionalPanel(
                condition = "input.protein_abs_rel == 'rel'",
                withSpinner(plotOutput("protein_graph_rel", width = 700, height = 350))
              ),
              conditionalPanel(
                condition = "input.rhythm_method != ''",
                #p("Rhythmicity p-values for: ", textOutput("protein_with_desc", inline = TRUE)),
                tableOutput("protein_rhythm_pvalues")
              ),
              conditionalPanel(
                condition = "input.proteins != ''",
                tableOutput("protein_table")
              ),
              # withSpinner(plotOutput("abundance_graph", width = 700, height = 350)),
              p("Source: Kay et al. Comms Bio 2021")
            )
          )
        )
      ),
    
      ## Rhythmic Mg proteins ----
      tabPanel(
        title = "Rhythmic Mg2+ binding proteins",
        sidebarLayout(
          ### Sidebar ----
          sidebarPanel(
            selectizeInput(
              input = "rhythmics", label = "Select proteins:", choices = NULL,
              multiple = TRUE
            ),
            actionButton("rhythmics_go", label = "Go!"),
            radioButtons("rhythmics_abs_rel",
              label = "Plot abundance as:",
              choices = c("Absolute" = "abs", "Relative" = "rel")
            )
          ),
          ### Main page ----
          mainPanel(
            conditionalPanel(
              condition = "input.rhythmics == ''",
              p("Select a protein to get started..")
            ),
            conditionalPanel(
              condition = "input.rhythmics != ''",
              conditionalPanel(
                condition = "input.rhythmics_abs_rel == 'abs'",
                withSpinner(plotOutput("rhythmics_graph_abs", width = 700, height = 350))
              ),
              conditionalPanel(
                condition = "input.rhythmics_abs_rel == 'rel'",
                withSpinner(plotOutput("rhythmics_graph_rel", width = 700, height = 350))
              ),
              conditionalPanel(
                condition = "input.rhythmics != ''",
                withSpinner(tableOutput("rhythmics_table"))
              ),
              # withSpinner(plotOutput("abundance_graph", width = 700, height = 350)),
              p("Source: Kay et al. Comms Bio 2021")
            )
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
  # Ions ----

  ## Get relevant values for selected ions, to go into graph/table ----
  ion_data <- eventReactive(input$ion_go, {
    req(input$ions)
    ion_data <- select(ion_abundance, "Time", contains(input$ions)) # create table with all necessary columns
    ion_data <- pivot_longer(ion_data, cols = !Time, names_to = "Ion", values_to = "Abundance") %>% drop_na()
    ion_data$Ion <- str_remove(ion_data$Ion, "\\_.*") # clean up ion names
    ion_data <- mutate(ion_data,
      Condition = case_when(ion_data$Time <= 24 ~ "LD", ion_data$Time >= 48 ~ "LL"),
      .after = Time
    ) # add condition column
    ion_data %>%
      group_by(Time, Condition, Ion) %>%
      summarise(across(), mean = mean(Abundance), sd = sd(Abundance)) %>%
      group_by(Ion) %>%
      mutate(rel_abun = mean / max(mean)) # add other/stats columns
    
  })

  # determine if raw values want to be shown (needed for raw values table to be visible)
  output$ion_show_raw <- renderText({
    paste(input$ion_raw)
  })
  outputOptions(output, "ion_show_raw", suspendWhenHidden = FALSE)

  # get values for selected ions as an output? (may not need)
  # output$ion_data <- reactive({
  #   req(input$ions)
  #   ion_data()
  # })

  # show table of values for selected ions
  output$ion_table <- renderTable({
    req(input$ions)
    ion_data() # get full table
    distinct(select(ion_data(), -c(Abundance, rel_abun))) # condense down for better visibility
  })
  
  # show which ions have been selected
  output$chosen_ions <- renderText({
    paste(input$ions)
  })

  # get whether abundance or relative are selected
  output$ion_abs_rel <- reactive({
    input$ion_abs_rel
  })

  ## Graph selected ion abundance changes over time ----

  ### Absolute abundance (best for single ions) ----
  output$ion_graph_abs <- renderPlot({
    ggplot(data = ion_data(), aes(x = Time, y = Abundance, group = Ion)) +
      ggtitle("Ion abundance") +
      geom_stripped_cols(width = 12, nudge_x = 6, odd = "#e5e5e5", even = "#ffffff", colour = "grey") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 12)) +
      scale_x_continuous(breaks = scales::breaks_width(12), minor_breaks = NULL) +
      scale_x_break(c(26, 48)) +
      geom_smooth(aes(colour = Ion, fill = Ion), span = 0.2, se = FALSE) +
      ggborderline::geom_borderline(
        linewidth = 0.5, aes(
          y = mean, colour = Ion,
          group = interaction(Ion, Condition)
        ),
        bordercolour = "black"
      ) +
      geom_errorbar(aes(y = mean, ymin = mean - sd, ymax = mean + sd)) +
      geom_point(size = 2, aes(x = Time, y = Abundance, group = Ion, fill = Ion), shape = 21, colour = "black") +
      
      labs(x = "Time (h)", y = "Ion abundance (ug/L)") 
  })


  ### Relative abundance (best for multiple ions) ----
  output$ion_graph_rel <- renderPlot({
    ggplot(data = ion_data(), aes(x = Time, y = rel_abun, group = Ion)) +
      ggtitle("Ion abundance") +
      geom_stripped_cols(width = 12, nudge_x = 6, odd = "#e5e5e5", even = "#ffffff", colour = "grey") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 12)) +
      scale_x_continuous(breaks = scales::breaks_width(12), minor_breaks = NULL) +
      scale_x_break(c(26, 48)) +
      geom_smooth(aes(colour = Ion, fill = Ion), span = 0.2, se = FALSE) +
      ggborderline::geom_borderline(
        linewidth = 0.5, aes(colour = Ion, group = interaction(Ion, Condition)),
        bordercolour = "black"
      ) +
      geom_point(size = 2, aes(fill = Ion), shape = 21, colour = "black") +
      labs(x = "Time (h)", y = "Relative ion abundance")
  })

  ## Consistent colour palette ----
  ion_graph_palette <- reactive({
    
  })
  
  ## Rhythmicity display ----

  output$ion_rain_echo <- reactive({
    input$ion_rain_echo
  })
  #outputOptions(output, "ion_rain_echo", suspendWhenHidden = FALSE)

  output$ion_rhythm_pvalues <- renderTable({
    req(input$ion_rain_echo)
    
    if (input$ion_rain_echo == "rain") {
      ion_rhythm_pvalues <- filter(ions_rain_LL, Ion %in% input$ions) 
    }
    if (input$ion_rain_echo == "echo") {
      ion_rhythm_pvalues <- filter(ions_echo_LL, Ion %in% input$ions) %>% select("Ion", "P-Value":"BY Adj P-Value")
    } 
    ion_rhythm_pvalues

  })



  ## Download table values for currently-selected ions ----
  ion_table_download <- downloadHandler(
    filename = function() {
      paste(input$chosen_ions, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(ion_data(), file)
    }
  )

  ## Download current ion graph view ----
  ion_graph_download <- downloadHandler(
    filename = "graph.png",
    content = function(file) {
      ggsave(file, plot = ion_graph())
    }
  )


  # Proteins ----

  ## update input options  ----
  updateSelectizeInput(session, "proteins", choices = c("Select a protein" = "", protein_abundance["Identifier"]), server = TRUE)
  

  # get selected protein from list
  protein_with_desc <- reactive({
    req(input$proteins)
    # paste(select(nuc_proteome, Identifier, description) %>% filter(Identifier == input$proteins), collapse = " - ")
  })

  # list selected protein
  output$protein_with_desc <- reactive({
    req(input$proteins)
    # protein_with_desc()
  })

  ## Get relevant values from proteome spreadsheet based on selected protein(s) ----
  protein_data <- reactive({
    req(input$proteins) # wait for proteins to be selected
    
    protein_data <- filter(protein_abundance, Identifier %in% (input$proteins)) # extract only chosen proteins
    
    protein_data <- pivot_longer(protein_data,
      cols = -Identifier,
      names_to = "Time", names_transform = as.numeric,
      values_to = "Abundance"
    ) # format like ions
    
    protein_data <- mutate(protein_data,
      Condition = case_when(
        protein_data$Time <= 24.5 ~ "LD",
        protein_data$Time >= 48 ~ "LL"
      ),
      .after = Time
    ) %>%
      group_by(Identifier) %>%
      mutate(rel_abun = Abundance / max(Abundance)) %>%
      ungroup() # add extra columns
  })

  # show more-nicely-formatted table of values for selected protein
  output$protein_table <- renderTable({
    protein_data() %>%
      arrange(Time) %>%
      relocate(Time, Condition) %>%
      rename("Relative abundance" = rel_abun) # %>% pivot_wider(names_from = Identifier, values_from = Abundance)
  })


  ## Graph for selected protein ----

  ### Absolute abundance ----
  output$protein_graph_abs <- renderPlot({
    ggplot(data = protein_data(), aes(x = Time, y = Abundance, group = Identifier)) +
      # ggtitle(protein_with_desc()) +
      ggtitle("Protein abundance") +
      # geom_stripped_cols(aes(group = "LD"), width = 12, nudge_x = 6, odd = "#ffffff", even = "#a5a5a5") +
      geom_stripped_cols(aes(group = "LL"), width = 12, nudge_x = 6, odd = "#e5e5e5", even = "#ffffff") +
      # geom_line(aes(group = Condition)) +
      geom_smooth(aes(colour = Identifier, fill = Identifier), span = 0.2, se = FALSE) +
      ggborderline::geom_borderline(
        linewidth = 0.5, aes(
          y = Abundance, colour = Identifier,
          group = interaction(Identifier, Condition)
        ),
        bordercolour = "black"
      ) +
      geom_point(size = 2, aes(x = Time, y = Abundance, group = Identifier, fill = Identifier), shape = 21, colour = "black") +
      
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 12)) +
      labs(x = "Time (h)", y = "Protein abundance") +
      scale_x_continuous(breaks = scales::breaks_width(12)) +
      scale_x_break(c(26, 48))
  })

  ### Relative abundance ----
  output$protein_graph_rel <- renderPlot({
    ggplot(data = protein_data(), aes(x = Time, y = rel_abun, group = Identifier)) +
      # ggtitle(protein_with_desc()) +
      ggtitle("Relative protein abundance") +
      # geom_stripped_cols(aes(group = "LD"), width = 12, nudge_x = 6, odd = "#ffffff", even = "#a5a5a5") +
      geom_stripped_cols(aes(group = "LL"), width = 12, nudge_x = 6, odd = "#e5e5e5", even = "#ffffff") +
      # geom_line(aes(group = Condition)) +
      geom_smooth(aes(colour = Identifier, fill = Identifier), span = 0.2, se = FALSE) +
      ggborderline::geom_borderline(
        linewidth = 0.5, aes(
          y = rel_abun, colour = Identifier,
          group = interaction(Identifier, Condition)
        ),
        bordercolour = "black"
      ) +
      geom_point(size = 2, aes(x = Time, y = rel_abun, group = Identifier, fill = Identifier), shape = 21, colour = "black") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 12)) +
      labs(x = "Time (h)", y = "Relative abundance") +
      scale_x_continuous(breaks = scales::breaks_width(12)) +
      scale_x_break(c(26, 48)) +
      ylim(NA,1)
  })


  ## show rhythmicity ----
  output$protein_rhythmicity <- renderTable(
    NULL
  )

  # get rhythmicity p-values for selected proteins
  output$protein_rhythm_pvalues <- renderTable({
    req(input$rhythm_method)
    filter(protein_rhythm_pvalues, Identifier %in% input$proteins) %>%
      select(Identifier, contains(input$rhythm_method))
  })

# Rhythmic proteins ----
  
## update input options ----
  updateSelectizeInput(session, "rhythmics", choices = c("Select a protein" = "", rhythmic_ion_binders_full["Identifier"]), server = TRUE)
  
## show table ----  
output$rhythmics_table <- renderTable({
  req(input$rhythmics)
  rhythmics_table <- filter(rhythmic_ion_binders_full, Identifier %in% (input$rhythmics)) %>% 
    select("Identifier", "description", "Absolute_Phase_LL", "Circadian_Phase_LL")
})
  
## get relevant values, to plot ----
rhythmics_data <- eventReactive(input$rhythmics_go, {
  req(input$rhythmics)
  rhythmics_data <- filter(rhythmic_ion_binders_full, Identifier %in% (input$rhythmics)) %>%
    select("Identifier", "0":"24.5", "48":"121.5") %>%
    pivot_longer(
      cols = -Identifier,
      names_to = "Time", names_transform = as.numeric,
      values_to = "Abundance"
    )
    rhythmics_data <- mutate(rhythmics_data,
      Condition = case_when(
        rhythmics_data$Time <= 24.5 ~ "LD",
        rhythmics_data$Time >= 48 ~ "LL"
      ),
      .after = Time
    ) 
    rhythmics_data %>%
    group_by(Identifier) %>%
    mutate(rel_abun = Abundance / max(Abundance)) %>%
    ungroup()
})





## Absolute abundance graph ----
  output$rhythmics_graph_abs <- renderPlot({
    ggplot(data = rhythmics_data(), aes(x = Time, y = Abundance, group = Identifier)) +
      # ggtitle(protein_with_desc()) +
      ggtitle("Protein abundance") +
      # geom_stripped_cols(aes(group = "LD"), width = 12, nudge_x = 6, odd = "#ffffff", even = "#a5a5a5") +
      geom_stripped_cols(aes(group = "LL"), width = 12, nudge_x = 6, odd = "#e5e5e5", even = "#ffffff") +
      # geom_line(aes(group = Condition)) +
      geom_smooth(aes(colour = Identifier, fill = Identifier), span = 0.2, se = FALSE) +
      ggborderline::geom_borderline(
        linewidth = 0.5, aes(
          y = Abundance, colour = Identifier,
          group = interaction(Identifier, Condition)
        ),
        bordercolour = "black"
      ) +
      geom_point(size = 2, aes(x = Time, y = Abundance, group = Identifier, fill = Identifier), shape = 21, colour = "black") +
      
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 12)) +
      labs(x = "Time (h)", y = "Protein abundance") +
      scale_x_continuous(breaks = scales::breaks_width(12)) +
      scale_x_break(c(26, 48))
  })

## Relative abundance ----  
output$rhythmics_graph_rel <- renderPlot({
  ggplot(data = rhythmics_data(), aes(x = Time, y = rel_abun, group = Identifier)) +
    # ggtitle(protein_with_desc()) +
    ggtitle("Relative protein abundance") +
    # geom_stripped_cols(aes(group = "LD"), width = 12, nudge_x = 6, odd = "#ffffff", even = "#a5a5a5") +
    geom_stripped_cols(aes(group = "LL"), width = 12, nudge_x = 6, odd = "#e5e5e5", even = "#ffffff") +
    # geom_line(aes(group = Condition)) +
    #geom_smooth(aes(colour = Identifier, fill = Identifier), span = 0.2, se = FALSE) +
    ggborderline::geom_borderline(
      linewidth = 0.5, aes(
        y = rel_abun, colour = Identifier,
        group = interaction(Identifier, Condition)
      ),
      bordercolour = "black"
    ) +
    geom_point(size = 2, aes(x = Time, y = rel_abun, group = Identifier, fill = Identifier), shape = 21, colour = "black") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 12)) +
    labs(x = "Time (h)", y = "Relative abundance") +
    scale_x_continuous(breaks = scales::breaks_width(12)) +
    scale_x_break(c(26, 48)) +
    ylim(NA,1)
})
  
  
# Misc. ----

  output$page <- renderText({
    paste("Currently on", input$tabset)
  })
}

# ///////////////////////////////////////// ----

# Run Shiny app ----
shinyApp(ui, server)
# runApp(display.mode = "showcase")

# ///////////////////////////////////////// ----
