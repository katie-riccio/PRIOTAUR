# ///////////////////////////////////////// ----

# initial setup ----

## libraries ----
library(bslib)
library(DT) # for nicer tables
library(ggplot2)
library(ggbreak)
library(ggstats) # for stripped columns
library(ggnewscale) # for setting colour palette
library(shiny)
library(shinycssloaders) # for spinners
library(readxl)
library(tidyverse)
library(viridis)
library(writexl)

## raw files ----
nuc_proteome <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM4_ESM.xlsx", sheet = 2, 
                           na = "ND", .name_repair = function(x) gsub("\\s+", "_", x)) # nuclear-encoded proteome data
org_proteome <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM4_ESM.xlsx", sheet = 3,
                           na = "ND", .name_repair = function(x) gsub("\\s+", "_", x)) # organellar-encoded proteome data
protein_pvalue_list <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM5_ESM.xlsx",
                                     na = "ND") %>% mutate(across(!Identifier, as.numeric))  # pre-calculated eJTK/rain/echo LD and LL p-values (names don't need fixing here)

nuc_proteome # check nuclear proteome
org_proteome # check organellar proteome
protein_pvalue_list # check rhythmicity values

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

double_rhythmics <- filter(protein_pvalue_list, LL_rain < 0.05 & LL_echo < 0.05) # find proteins rhythmic in both RAIN and ECHO in LL
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

time_shading <- tibble(Time = seq(0, 122, 1)) %>%
  mutate(
    Condition = case_when(Time <= 24 ~ "LD", Time >= 48 ~ "LL"),
    ToD = case_when(
      ((Time %/% 12) %% 2) == 1 & Condition == "LD" ~ "Night",
      ((Time %/% 12) %% 2) == 1 & Condition == "LL" ~ "Subjective night",
      ((Time %/% 12) %% 2) == 0 ~ "Day"
    ),
    colour = case_when(
      ToD == "Day" ~ "#ffffff",
      ToD == "Night" ~ "#a5a5a5",
      ToD == "Subjective night" ~ "#e5e5e5"
    )
  ) %>% drop_na() # set up graph background 

#ion_colours <- c("Mg" = "tomato", "Ca" = "seagreen", "K" = "royalblue") # hard code palette test

protein_history <- c() #initialise blank history
protein_colour_table <- tibble("Identifier" = NA, "Colour" = NA) # blank

# ///////////////////////////////////////// ----
# UI ----
# ///////////////////////////////////////// ----

ui <- page_fillable(
  theme = bs_theme(bootswatch = "flatly"),
  mainPanel(
    textOutput("page"),
  #  input_dark_mode(id = "mode"), #dark mode toggle
    card(
    tabsetPanel(
      id = "tabset",
      
      ## Ion visualisation ----
      tabPanel(
        title = "Ion visualisation",
        sidebarLayout(
          ### Sidebar ----
          sidebarPanel(
            # style = "position: fixed; overflow-y: auto",
            card(
              checkboxGroupInput(
                input = "ions", label = "Select ions:", selected = NULL,
                choices = c("Calcium" = "Ca", "Potassium" = "K", "Magnesium" = "Mg")
              ),
              radioButtons("ion_abs_rel",
                label = "Plot abundance as:",
                choices = c("Absolute" = "abs", "Relative" = "rel")
              ),
              radioButtons("ion_rain_echo",
                label = "Show LL values from:",
                choices = c("RAIN" = "rain", "ECHO" = "echo")
              ),
              actionButton("ion_go", label = "Draw graph")
            ),
            checkboxInput(input = "ion_raw", label = "View raw values?"),
            selectInput("ion_graph_palette", label = "Select palette",
                        choices = c("viridis", "magma", "inferno", "plasma",
                                    "cividis", "rocket", "mako", "turbo"))
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
              card(
              conditionalPanel(
                condition = "input.ion_abs_rel == 'abs'",
                withSpinner(plotOutput("ion_graph_abs")),
                downloadButton("ion_graph_abs_dl", "Download graph view")
              
              ),
              
              conditionalPanel(
                condition = "input.ion_abs_rel == 'rel'",
                withSpinner(plotOutput("ion_graph_rel")),
                downloadButton("ion_graph_rel_dl", "Download graph view")
              )
              ),
              # withSpinner(plotOutput("ion_graph", width = 600, height = 300))
            ),
            # p("You have selected:", textOutput("ion_rain_echo")),
            dataTableOutput("ion_rhythm_table"),
            conditionalPanel(
              condition = "output.ion_show_raw == 'TRUE'",
              downloadButton("ion_table_dl", "Download table of values"),
              withSpinner(dataTableOutput("ion_table"))
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
          selectInput("protein_graph_palette", label = "Select palette",
                      choices = c("viridis", "magma", "inferno", "plasma",
                                  "cividis", "rocket", "mako", "turbo")),
          actionButton("protein_go", label = "Draw graph"),
          conditionalPanel(
            condition = "input.proteins != ''",
            checkboxGroupInput(
              input = "rhythm_method", label = "Show rhythmicity p-value?", selected = NULL,
              choices = c("eJTK" = "ejtk", "RAIN" = "rain", "ECHO" = "echo")
            )
          ),
          textOutput("protein_history")
        ),
        
        ### Main page ----
        mainPanel(
          textOutput("protein_colours"),
          conditionalPanel(
            condition = "input.proteins == ''",
            p("Select a protein to get started..")
          ),
          conditionalPanel(
            condition = "input.proteins != ''",
            card(
            conditionalPanel(
              condition = "input.protein_abs_rel == 'abs'",
              withSpinner(plotOutput("protein_graph_abs")),
              downloadButton("protein_graph_abs_dl", "Download graph view")
            ),
            conditionalPanel(
              condition = "input.protein_abs_rel == 'rel'",
              withSpinner(plotOutput("protein_graph_rel")),
              downloadButton("protein_graph_rel_dl", "Download graph view")
            )
            ),
            conditionalPanel(
              condition = "input.rhythm_method != ''",
              # p("Rhythmicity p-values for: ", textOutput("protein_with_desc", inline = TRUE)),
              dataTableOutput("protein_rhythm_table")
            ),
            conditionalPanel(
              condition = "input.proteins != ''",
              downloadButton("protein_table_dl", "Download table of values"),
              dataTableOutput("protein_table")
            ),
            # withSpinner(plotOutput("abundance_graph")),
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
          actionButton("rhythmics_go", label = "Draw graph"),
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
              withSpinner(plotOutput("rhythmics_graph_abs")),
              downloadButton("rhythmics_graph_abs_dl", "Download graph view")
            ),
            conditionalPanel(
              condition = "input.rhythmics_abs_rel == 'rel'",
              withSpinner(plotOutput("rhythmics_graph_rel")),
              downloadButton("rhythmics_graph_rel_dl", "Download graph view")
            ),
            conditionalPanel(
              condition = "input.rhythmics != ''",
              downloadButton("rhythmics_table_dl", "Download table of values"),
              withSpinner(dataTableOutput("rhythmics_table"))
            ),
            # withSpinner(plotOutput("abundance_graph")),
            p("Source: Kay et al. Comms Bio 2021")
          )
        ) # mg-binding mainpanel
      ) # mg-binding sidebarlayout
    ) # mg-binding tabpanel
    
  ) # tabsetpanel
) # card
) # mainpanel
) # ui


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
  output$ion_table <- renderDataTable({
    req(input$ion_go)
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
  ion_abs_data <- reactive({
    ggplot(data = ion_data(), aes(x = Time, y = Abundance, group = Ion)) +
      
      geom_rect(data = time_shading,
                aes(y = NULL, group = NULL,
                    xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour), 
                show.legend = FALSE) + # set up background
      scale_fill_identity(time_shading$colour) + # shade appropriately
      new_scale_fill() + # allow rest of graph to use different shading method
      
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
      scale_colour_manual(values = ion_colours()) +
      scale_fill_manual(values = ion_colours()) +
      
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), 
            axis.title = element_text(size = 12)) +
      scale_x_continuous(breaks = scales::breaks_width(12), minor_breaks = NULL,
                         limits = c(0, 96)) +
      scale_x_break(c(26, 48)) +
      
      labs(x = "Time (h)", y = "Ion abundance (ug/L)",
           title = "Ion abundance") 
  })
  
  output$ion_graph_abs <- renderPlot({
    ion_abs_data()
  })


  ### Relative abundance (best for multiple ions) ----
  ion_rel_data <- reactive({
    ggplot(data = ion_data(), aes(x = Time, y = rel_abun, group = Ion)) +
      ggtitle("Ion abundance") +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 12)) +
      scale_x_continuous(breaks = scales::breaks_width(12), minor_breaks = NULL,
                         limits = c(0, 96)) +
      scale_x_break(c(26, 48)) +
      
      geom_rect(data = time_shading,
                aes(y = NULL, group = NULL,
                    xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour), 
                show.legend = FALSE) + # set up background
      scale_fill_identity(time_shading$colour) + # shade appropriately
      new_scale_fill() + # allow rest of graph to use different shading method
      
      geom_smooth(aes(colour = Ion, fill = Ion), span = 0.2, se = FALSE) +
      ggborderline::geom_borderline(
        linewidth = 0.5, aes(colour = Ion, group = interaction(Ion, Condition)),
        bordercolour = "black"
      ) +
      geom_point(size = 2, aes(fill = Ion), shape = 21, colour = "black") +
      scale_colour_manual(values = ion_colours()) +
      scale_fill_manual(values = ion_colours()) +
      
      labs(x = "Time (h)", y = "Relative ion abundance")
  })
  
  output$ion_graph_rel <- renderPlot({
    ion_rel_data()
  })

  ## Consistent colour palette ----
  ion_colours <- reactive({
    ion_palette <- viridis(n = 3, option = input$ion_graph_palette) # set up based on user selection
    ion_colours <- c("Mg" = ion_palette[1], "Ca" = ion_palette[2],"K" = ion_palette[3]) # assign colours specifically to each ion based on current palette, so colours are consistent across plots
  })
  
  #output$ion_colours <- renderText({ion_colours()})
  
  ## Rhythmicity display ----

  output$ion_rain_echo <- reactive({
    input$ion_rain_echo
  })
  #outputOptions(output, "ion_rain_echo", suspendWhenHidden = FALSE)

  ion_rhythm_pvalues <- eventReactive(input$ion_go, {
    if (input$ion_rain_echo == "rain") {
      ion_rhythm_pvalues <- filter(ions_rain_LL, Ion %in% input$ions) %>% 
        select("Ion", "P-Value" = "pVal", "Phase" = "phase", "Period" = "period")
    }
    if (input$ion_rain_echo == "echo") {
      ion_rhythm_pvalues <- filter(ions_echo_LL, Ion %in% input$ions) %>% 
        select("Ion", "P-Value":"BY Adj P-Value")
    } 
    
    ion_rhythm_pvalues
    
  })

  output$ion_rhythm_table <- renderDataTable({
    datatable(ion_rhythm_pvalues(), options = list(dom = "t"), rownames = FALSE) %>% 
      formatStyle(-c(1), 
                  fontStyle = styleInterval(c(0.05, 1), c("normal", "italic", "normal")),
                  colour = styleInterval(0.05, c('black', 'darkgrey'))
                  ) %>%
      formatSignif(-c(1), 5)
  })
  
  ## Downloads ----
  
  ### absolute abundance ----  
  output$ion_graph_abs_dl <- downloadHandler(
    filename = function() {paste((format(Sys.time(), "%F_%H%M%S")),
                                 "-ions_absolute", 
                                 ".png", 
                                 sep = "")},
    content = function(file) {
      ggsave(file, plot = ion_abs_data(), device = "png",
             width = 8, height = 5) #this becomes 2400x1500px
    }
  )
  ### relative abundance ----
  output$ion_graph_rel_dl <- downloadHandler(
    filename = function() {paste((format(Sys.time(), "%F_%H%M%S")),
                                 "-ions_relative", ".png", 
                                 sep = "")},
    content = function(file) {
      ggsave(file, plot = ion_rel_data(), device = "png",
             width = 8, height = 5) #this becomes 2400x1500px
    }
  )
  ### table of values ----
  output$ion_table_dl <- downloadHandler(
    filename = function() {paste((format(Sys.time(), "%F_%H%M%S")), 
                                 "-ions", ".tsv", 
                                 sep = "")},
    content = function(file) {write_tsv(ion_data(), file)}
  )


  # Proteins ----

  ## update input options  ----
  updateSelectizeInput(session, "proteins", 
                       choices = c("Select a protein" = "", protein_abundance["Identifier"]), 
                       server = TRUE) # fill options
  

  # get selected proteins from list (old)
  protein_with_desc <- reactive({
    req(input$proteins)
    # paste(select(nuc_proteome, Identifier, description) %>% filter(Identifier == input$proteins), collapse = " - ")
  })

  # list selected proteins
  output$protein_with_desc <- reactive({
    req(input$proteins)
    # protein_with_desc()
  })

  ## Get relevant values from proteome spreadsheet based on selected protein(s) ----
  protein_data <- eventReactive(input$protein_go, {
    req(input$proteins) # wait for proteins to be selected
    
    protein_data <- filter(protein_abundance, Identifier %in% (input$proteins)) # extract only chosen proteins
    
    protein_data <- pivot_longer(protein_data,
      cols = -Identifier,
      names_to = "Time", names_transform = as.numeric,
      values_to = "Abundance"
    ) # format to plot
    protein_data <- mutate(protein_data,
      Condition = case_when(
        protein_data$Time <= 24.5 ~ "LD",
        protein_data$Time >= 48 ~ "LL"
      ), .after = Time ) %>%
      group_by(Identifier) %>%
      mutate(rel_abun = Abundance / max(Abundance)) %>%
      ungroup() # add extra columns
  })

  ### store history of selected proteins ----
  
  history <- reactiveValues(protein = NULL) # initialise empty list on startup
  
  observeEvent(input$protein_go, {
    history$protein <- unique(c(history$protein, protein_colours())) # update list when graph is drawn
  })
  
  # output list as a test
  output$protein_history <- renderText({
    paste(history$protein, collapse = ", ")
    }) 
  
  
  # show more-nicely-formatted table of values for selected protein
  output$protein_table <- renderDataTable({
    protein_data() %>%
      arrange(Time) %>%
      relocate(Time, Condition) %>%
      rename("Relative abundance" = rel_abun) # %>% pivot_wider(names_from = Identifier, values_from = Abundance)
  })


  ## Graph for selected proteins ----
  
  ### Consistent colours across plots ----
  protein_colours <- eventReactive(input$protein_go, {
    protein_palette <- viridis(n = length(input$proteins), option = input$protein_graph_palette) # set up colours to assign based on user selection
    
    for (id_num in 1:length(input$proteins)) {
      protein_colour_table <- add_row(protein_colour_table,
                                 Identifier = (input$proteins)[id_num], 
                                 Colour = protein_palette[id_num]
                                 ) # create palette one row at a time
    }
    protein_colours <- unlist(split(protein_colour_table$Colour, 
                                    protein_colour_table$Identifier))
  })

  output$protein_colours <- renderText({
    protein_colours()
    })
  
  ### Absolute abundance ----
  
  protein_abs_data <- eventReactive(input$protein_go, {
    ggplot(data = protein_data(), aes(x = Time, y = Abundance, group = Identifier)) +
      # ggtitle(protein_with_desc()) +
      ggtitle("Protein abundance") +
      
      geom_rect(data = time_shading,
                aes(y = NULL, group = NULL,
                    xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour), 
                show.legend = FALSE) + # set up background
      scale_fill_identity(time_shading$colour) + # shade appropriately
      new_scale_fill() + # allow rest of graph to use a different shading method
      
      geom_smooth(aes(colour = Identifier, fill = Identifier), span = 0.2, se = FALSE) +
      ggborderline::geom_borderline(
        linewidth = 0.5, aes(
          y = Abundance, colour = Identifier,
          group = interaction(Identifier, Condition)
        ),
        bordercolour = "black"
      ) +
      geom_point(size = 2, aes(x = Time, y = Abundance, group = Identifier, fill = Identifier), shape = 21, colour = "black") +
      
      scale_colour_manual(values = protein_colours()) +
      scale_fill_manual(values = protein_colours()) +
      
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 12)) +
      labs(x = "Time (h)", y = "Protein abundance") +
      scale_x_continuous(breaks = scales::breaks_width(12)) +
      scale_x_break(c(26, 48))
  })
  
  output$protein_graph_abs <- renderPlot({
    protein_abs_data()
  })

  ### Relative abundance ----
  protein_rel_data <- eventReactive(input$protein_go, {
    
    ggplot(data = protein_data(), aes(x = Time, y = rel_abun, group = Identifier)) +
      # ggtitle(protein_with_desc()) +
      ggtitle("Relative protein abundance") +
      
      geom_rect(data = time_shading,
                aes(y = NULL, group = NULL,
                    xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour), 
                show.legend = FALSE) + # set up background
      scale_fill_identity(time_shading$colour) + # shade appropriately
      new_scale_fill() + # allow rest of graph to use different shading method
      
      geom_smooth(aes(colour = Identifier, fill = Identifier), span = 0.2, se = FALSE) +
      ggborderline::geom_borderline(
        linewidth = 0.5, aes(
          y = rel_abun, colour = Identifier,
          group = interaction(Identifier, Condition)
        ),
        bordercolour = "black"
      ) +
      geom_point(size = 2, aes(x = Time, y = rel_abun, group = Identifier, fill = Identifier), shape = 21, colour = "black") +
      
      scale_colour_manual(values = protein_colours()) +
      scale_fill_manual(values = protein_colours()) +
      
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 12)) +
      labs(x = "Time (h)", y = "Relative abundance") +
      scale_x_continuous(breaks = scales::breaks_width(12)) +
      scale_x_break(c(26, 48)) +
      ylim(NA,1)
  })
  
  output$protein_graph_rel <- renderPlot({
    protein_rel_data()
  })


## show rhythmicity ----

  # get rhythmicity p-values for selected proteins
  protein_rhythm_pvalues <- reactive({
    req(input$rhythm_method)
    protein_rhythm_pvalues <- filter(protein_pvalue_list, Identifier %in% input$proteins) %>%
      rename(c("eJTK (LD)" = "LD_ejtk",
               "RAIN (LD)" = "LD_rain",
               "ECHO (LD)" = "LD_echo",
               "eJTK (LL)" = "LL_ejtk",
               "RAIN (LL)" = "LL_rain",
               "ECHO (LL)" = "LL_echo")
             ) %>% 
      select(Identifier, contains(input$rhythm_method)) %>%
    relocate(ends_with("(LD)"), .after = Identifier) #reorder
  })

  # show values for selected proteins
  output$protein_rhythm_table <- renderDataTable({
    datatable(protein_rhythm_pvalues(), options = list(dom = "t"), rownames = FALSE) %>% 
      formatStyle(-c(1), 
                  #fontStyle = styleInterval(c(0.05, 1), c("normal", "italic", "normal")),
                  fontWeight = styleInterval(0.05, c("bold", "normal")),
                  colour = styleInterval(0.05, c('black', 'darkgrey'))
      ) %>%
      formatSignif(-c(1), 5)
  })
  
  ## Downloads ----
  
  ### absolute abundance ----  
  output$protein_graph_abs_dl <- downloadHandler(
    filename = function() {paste((format(Sys.time(), "%F_%H%M%S")),
                                 "-proteins_absolute", 
                                 ".png", 
                                 sep = "")},
    content = function(file) {
      ggsave(file, plot = protein_abs_data(), device = "png",
             width = 8, height = 5) #this becomes 2400x1500px
    }
  )
  ### relative abundance ----
  output$protein_graph_rel_dl <- downloadHandler(
    filename = function() {paste((format(Sys.time(), "%F_%H%M%S")),
                                 "-proteins_relative", ".png", 
                                 sep = "")},
    content = function(file) {
      ggsave(file, plot = protein_rel_data(), device = "png",
             width = 8, height = 5) #this becomes 2400x1500px
    }
  )
  ### table of values ----
  output$protein_table_dl <- downloadHandler(
    filename = function() {paste((format(Sys.time(), "%F_%H%M%S")), 
                                 "-proteins", ".tsv", 
                                 sep = "")},
    content = function(file) {
      write_tsv(protein_data(), file)
    }
  )

  
# Rhythmic Mg2+ binding proteins ----
  
## update input options ----
  updateSelectizeInput(session, "rhythmics", 
                       choices = c("Select a protein" = "", rhythmic_ion_binders_full["Identifier"]), 
                       server = TRUE)
  
## show table ----  
rhythmics_table <- eventReactive(input$rhythmics_go, {
  filter(rhythmic_ion_binders_full, Identifier %in% (input$rhythmics)) %>% 
  select("Identifier", "description", "Absolute_Phase_LL", "Circadian_Phase_LL")
})
  
output$rhythmics_table <- renderDataTable({
  req(input$rhythmics)
  rhythmics_table()
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

## plot relevant values ----

### Absolute abundance ----
  rhythmics_abs_data <- eventReactive(input$rhythmics_go, {
    ggplot(data = rhythmics_data(), aes(x = Time, y = Abundance, group = Identifier)) +

      ggtitle("Protein abundance") +
      
      geom_rect(data = time_shading,
                aes(y = NULL, group = NULL,
                    xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour), 
                show.legend = FALSE) + # set up background
      scale_fill_identity(time_shading$colour) + # shade appropriately
      new_scale_fill() + # allow rest of graph to use different shading method
      
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
  
  output$rhythmics_graph_abs <- renderPlot({
    rhythmics_abs_data()
  })
  

### Relative abundance ----  
  rhythmics_rel_data <- eventReactive(input$rhythmics_go, {
    ggplot(data = rhythmics_data(), aes(x = Time, y = rel_abun, group = Identifier)) +
      # ggtitle(protein_with_desc()) +
      ggtitle("Relative protein abundance") +
      
      geom_rect(data = time_shading,
                aes(y = NULL, group = NULL,
                    xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour), 
                show.legend = FALSE) + # set up background
      scale_fill_identity(time_shading$colour) + # shade appropriately
      new_scale_fill() + # allow rest of graph to use different shading method
      
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
  
output$rhythmics_graph_rel <- renderPlot({
  rhythmics_rel_data()
})

## Downloads ----

### absolute abundance ----  
output$rhythmics_graph_abs_dl <- downloadHandler(
  filename = function() {paste((format(Sys.time(), "%F_%H%M%S")),
                               "-rhythmics_absolute", 
                               ".png", 
                               sep = "")},
    content = function(file) {
      ggsave(file, plot = rhythmics_abs_data(), device = "png",
             width = 8, height = 5) #this becomes 2400x1500px
    }
)
### relative abundance ----
output$rhythmics_graph_rel_dl <- downloadHandler(
  filename = function() {paste((format(Sys.time(), "%F_%H%M%S")),
                               "-rhythmics_relative", ".png", 
                               sep = "")},
  content = function(file) {
    ggsave(file, plot = rhythmics_rel_data(), device = "png",
           width = 8, height = 5) #this becomes 2400x1500px
  }
)
### table of values ----
output$rhythmics_table_dl <- downloadHandler(
  filename = function() {paste((format(Sys.time(), "%F_%H%M%S")), 
                               "-rhythmics", ".tsv", 
                               sep = "")},
  content = function(file) {
    write_tsv(rhythmics_table(), file)
  }
)
  
# Misc. ----

  output$page <- renderText({
    paste("Current panel:", input$tabset)
  })
}

# ///////////////////////////////////////// ----

# Run Shiny app ----
shinyApp(ui, server)
# runApp(display.mode = "showcase")

# ///////////////////////////////////////// ----
