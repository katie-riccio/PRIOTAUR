# ///////////////////////////////////////// ----

# PRIOTAUR ----

# ///////////////////////////////////////// ----



# initial setup ----

## libraries ----
library(bslib) # for layout options
library(bsicons) # for tooltip icons
library(DT) # for nicer tables
library(ggplot2) # for nicer graphs
library(ggbreak) # for axis breaks
library(ggborderline) # for easier-to-see line display
library(ggnewscale) # for background time shading
library(scales) # for setting graph scales
library(shiny)
library(shinycssloaders) # for spinners
library(readxl)
library(tidyverse)
library(viridis) # for colour-blind-safer palette options
library(writexl)

## loading files ----
nuc_proteome <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM4_ESM.xlsx",
  sheet = 2,
  na = "ND", .name_repair = function(x) gsub("\\s+", "_", x)
) # nuclear-encoded proteome data
org_proteome <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM4_ESM.xlsx",
  sheet = 3,
  na = "ND", .name_repair = function(x) gsub("\\s+", "_", x)
) # organellar-encoded proteome data

nuc_pvals <- read_excel("Kay et al Comms Bio 2021/42003_2021_2680_MOESM5_ESM.xlsx",
  na = "ND"
) %>% mutate(across(!Identifier, as.numeric)) # pre-calculated eJTK/rain/echo LD and LL p-values (names don't need fixing here)
org_pvals <- select(org_proteome, "Identifier", "LD_ejtk" = "p_val_LD", "LL_ejtk" = "p_val_LL") # pre-calculated eJTK values for organellar proteome
my_org_pvals <- read_excel("organelle_rhythmicity_pvals.xlsx") # RAIN and ECHO p-values for organellar proteome, that I calculated

ion_abundance <- read_excel("ions.xlsx", range = "A3:M41", col_names = c(
  "Time",
  "Mg_1", "Mg_2", "Mg_3", "Mg_4",
  "Ca_1", "Ca_2", "Ca_3", "Ca_4",
  "K_1", "K_2", "K_3", "K_4"
)) # ion abundance values - columns aren't named in the actual raw file so that's done here

ion_rhythmicity <- "ion_rhythmicity.xlsx"
ion_rhythmicity <- ion_rhythmicity %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = ion_rhythmicity) %>%
  list2env(.GlobalEnv) # load in pre-calculated RAIN and ECHO outputs from excel sheet


## one-time-only processes ----

all_org_pvals <- left_join(org_pvals, my_org_pvals, by = "Identifier") %>%
  relocate(LL_ejtk, .after = LD_echo) # combine and set correct column order
protein_pvals <- bind_rows(nuc_pvals, all_org_pvals) # combine into one table

protein_abundance <- rbind(
  select(nuc_proteome, "Identifier", "0":"24.5", "48":"121.5"),
  select(org_proteome, "Identifier", "0":"24.5", "48":"121.5")
) # combine abundance values into one table

double_rhythmics <- filter(protein_pvals, LL_rain < 0.05 & LL_echo < 0.05) # list of proteins rhythmic in both RAIN and ECHO in LL

double_rhythmics_nuc <- filter(nuc_proteome, Identifier %in% double_rhythmics$Identifier) # get full information for rhythmic proteins from main table
double_rhythmics_org <- filter(org_proteome, Identifier %in% double_rhythmics$Identifier) %>%
  rename(GeneType = Gene, description = Name) # get full information (renaming columns so they bind properly)
double_rhythmics_full <- bind_rows(double_rhythmics_nuc, double_rhythmics_org) # bind to one table

mg_co_rhythmics <- read_excel("mg_co_rhythmics.xlsx") # load in info of rhythmic proteins that use Mg as a co-factor

rhythmic_ion_binders_full <- filter(
  double_rhythmics_full,
  Identifier %in% mg_co_rhythmics$Identifier
) # get full proteome table for 17


# ion_colours <- c("Mg" = "tomato", "Ca" = "seagreen", "K" = "royalblue") # hard code palette test

### table of identifiers and colours to use for proteins
protein_colour_table <- tibble("Identifier" = NA, "Colour" = NA) # blank list
rhythmics_colour_table <- tibble("Identifier" = NA, "Colour" = NA) # blank list


### time shading ----
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
  ) %>%
  drop_na() # set up graph background

#### new time test (WIP) ----
time_spacing <- tibble(Time = seq(0, 120, 12))


### plot theme setup ----
plot_theme <- theme_set(
  theme_linedraw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
)


# ///////////////////////////////////////// ----
# UI ----
# ///////////////////////////////////////// ----

ui <- page_navbar(
  # Settings ----
  title = "PRIOTAUR",
  id = "tab",
  # position = "fixed-top",
  theme = bs_theme(bootswatch = "flatly", secondary = "steelblue"),

  # Sidebars ----
  sidebar = sidebar(
    mobile = "always-above",
    conditionalPanel(
      condition = "input.tab == 'Home'",
      p("Select a tab to get started!")
    ),



    ## Ions ----
    conditionalPanel(
      condition = "input.tab == 'Ions'",
      card( # put everything that depends on the action button in its own card
        checkboxGroupInput(
          input = "ions", label = "Ions to plot:", selected = NULL,
          choices = c("Calcium" = "Ca", "Potassium" = "K", "Magnesium" = "Mg")
        ),
        radioButtons("ion_abs_rel",
          label = list(
            "Plot abundance as:",
            tooltip(
              trigger = bs_icon("info-circle", size = "0.8em"),
              "Relative is recommended for multiple ions"
            )
          ),
          choices = c("Absolute" = "abs", "Relative" = "rel")
        ),
        radioButtons("ion_rain_echo",
          label = "Show LL rhythmicity values from:",
          choices = c("RAIN" = "rain", "ECHO" = "echo")
        ),
        selectizeInput("ion_graph_palette",
          label = "Palette:",
          choices = c(
            "viridis", "magma", "inferno", "plasma",
            "cividis", "rocket", "mako", "turbo"
          ),
          options = list(dropdownParent = "body")
        ), # viridis palette names

        actionButton("ion_go", label = "Plot graph")
      ), # card end

      conditionalPanel(
        condition = "input.tab == 'Proteins'", # hide for the moment
        checkboxInput(input = "ion_raw", label = "View raw values?", value = TRUE)
      )
    ),

    ## Proteins ----
    conditionalPanel(
      condition = "input.tab == 'Proteins'",
      card( # put everything that relies on the action button in its own card
        selectizeInput(
          input = "proteins",
          label = list(
            "Proteins to plot:",
            tooltip(
              trigger = bs_icon("info-circle", size = "0.8em"),
              "Up to 5 can be selected at once"
            )
          ),
          choices = NULL,
          multiple = TRUE, options = list(maxItems = 5)
        ),
        radioButtons("protein_abs_rel",
          label = list(
            "Plot abundance as:",
            tooltip(
              trigger = bs_icon("info-circle", size = "0.8em"),
              "Relative is recommended for multiple proteins"
            )
          ),
          choices = c("Absolute" = "abs", "Relative" = "rel")
        ),
        selectizeInput("protein_graph_palette",
          label = "Palette:",
          choices = c(
            "viridis", "magma", "inferno", "plasma",
            "cividis", "rocket", "mako", "turbo"
          ), # viridis palette names
          selected = "plasma",
          options = list(dropdownParent = "body")
        ),
        actionButton("protein_go", label = "Plot graph")
      ) # card end
      ,
      card(
        checkboxGroupInput(
          input = "rhythm_method", label = "View rhythmicity p-value from:", selected = c("ejtk", "rain", "echo"),
          choices = c("eJTK" = "ejtk", "RAIN" = "rain", "ECHO" = "echo")
        ),
        conditionalPanel(
          condition = "input.protein_go != ''",
          card(
            p("Previously-selected proteins:"), textOutput("protein_history")
          )
        )
      )
    ),
    ## Rhythmic ion-binding proteins ----
    conditionalPanel(
      condition = "input.tab == 'Rhythmic ion-binding proteins'",
      card( # put everything that relies on the action button in a card with it
        selectizeInput(
          input = "rhythmics", label = list(
            "Proteins to plot:",
            tooltip(
              trigger = bs_icon("info-circle", size = "0.8em"),
              "Up to 5 can be selected at once"
            )
          ),
          choices = NULL,
          multiple = TRUE, options = list(maxItems = 5)
        ),
        radioButtons("rhythmics_abs_rel",
          label = list(
            "Plot abundance as:",
            tooltip(
              trigger = bs_icon("info-circle", size = "0.8em"),
              "Relative is recommended for multiple proteins"
            )
          ),
          choices = c("Absolute" = "abs", "Relative" = "rel")
        ),
        conditionalPanel(
          condition = "input.rhythmics_abs_rel == 'rel'",
          

          checkboxInput("add_mg_trace",
                        label = HTML("Add Mg2+ trace? <br> (NB: Limits view to 2 days)")
          )
        ),
        selectizeInput("rhythmics_graph_palette",
          label = "Palette:",
          choices = c(
            "viridis", "magma", "inferno", "plasma",
            "cividis", "rocket", "mako", "turbo"
          ), # viridis palette names
          selected = "plasma",
          options = list(dropdownParent = "body")
        ),
        actionButton("rhythmics_go", label = "Plot graph")
      ) # card end
    )
  ),

  # Main pages ----

  ## Homepage ----
  nav_panel(
    title = "Home",
    HTML("<h1> PRIOTAUR: Proteome and Rhythmic Ions of <i>Ostreococcus TAURi</i> </h1>
         <p> Welcome to PRIOTAUR! This is an online interface I developed for my undergraduate honours project.
         <br> It allows for the exploration and critical integration of a deep-coverage proteome with multi-ion abundance data for <i>Ostreococcus tauri</i>, a key model microalga species.
         <br> faciliating comparative analyses at multiple omic scales through an easy-to-use interface.
         <br>
         <br> Use the tabs above to navigate between pages.
         <br> If you have any questions about the tool, please feel free to get in contact!
         </p>")
  ),

  ## Ions ----
  nav_panel(
    title = "Ions",
    conditionalPanel(
      condition = "input.ions == ''",
      card(p("Select one or more ions to get started!"))
    ),
    # p("Show all values?", textOutput("ion_show_raw")),
    conditionalPanel(
      condition = "input.ion_go != 0",
      conditionalPanel(
        condition = "input.ions != ''",
        # p("test message - will display if graph can display"),
        # p("You have selected:", textOutput("chosen_ions")),
        card(
          min_height = 650,
          card_header("Here is your graph!"),
          conditionalPanel(
            condition = "input.ion_abs_rel == 'abs'",
            withSpinner(plotOutput(width = 800, height = 500, outputId = "ion_graph_abs"))
          ),
          conditionalPanel(
            condition = "input.ion_abs_rel == 'rel'",
            withSpinner(plotOutput(width = 800, height = 500, outputId = "ion_graph_rel"))
          ),

          # p("Ion abundance values over 1 day LD and 2 days LL."),

          card_footer(list(
            conditionalPanel(
              "input.ion_abs_rel == 'abs'",
              downloadButton("ion_graph_abs_dl", "Download graph view")
            ),
            conditionalPanel(
              "input.ion_abs_rel == 'rel'",
              downloadButton("ion_graph_rel_dl", "Download graph view")
            )
          ))
          # withSpinner(plotOutput(width = 800, height = 500, outputId = "ion_graph", width = 600, height = 300))
        ),
        # p("You have selected:", textOutput("ion_rain_echo")),

        card(
          min_height = 300,
          card_header("Rhythmicity of selected ion(s):"),
          withSpinner(dataTableOutput("ion_rhythm_table")),
          card_footer(HTML("<i>Italics</i> indicate values not significantly rhythmic at the <i>p</i> < 0.05 threshold."))
        )
      ),
      conditionalPanel(
        condition = "output.ion_show_raw == 'TRUE'",
        card(
          min_height = 500, fill = FALSE,
          card_header("Table of values:"),
          withSpinner(dataTableOutput("ion_table")),
          card_footer(downloadButton("ion_table_dl", "Download table of values"))
        )
      ),
      p("Source for ion abundance data: Gil Rodríguez et al. bioRXiv 2024")
    )
  ),

  ## Proteins ----
  nav_panel(
    title = "Proteins",
    conditionalPanel(
      condition = "input.proteins == ''",
      card(p("Select one or more proteins to get started!"))
    ),
    conditionalPanel(
      condition = "input.protein_go != 0",
      conditionalPanel(
        condition = "input.proteins != ''",
        card(
          min_height = 650,
          card_header("Here is your graph!"),
          conditionalPanel(
            condition = "input.protein_abs_rel == 'abs'",
            withSpinner(plotOutput(width = 800, height = 500, outputId = "protein_graph_abs"))
          ),
          conditionalPanel(
            condition = "input.protein_abs_rel == 'rel'",
            withSpinner(plotOutput(width = 800, height = 500, outputId = "protein_graph_rel"))
          ),
          card_footer(
            conditionalPanel(
              condition = "input.protein_abs_rel == 'abs'",
              downloadButton("protein_graph_abs_dl", "Download graph view")
            ),
            conditionalPanel(
              condition = "input.protein_abs_rel == 'rel'",
              downloadButton("protein_graph_rel_dl", "Download graph view")
            )
          )
        ),
        conditionalPanel(
          condition = "input.rhythm_method != ''",
          # p("Rhythmicity p-values for: ", textOutput("protein_with_desc", inline = TRUE)),
          card(
            min_height = 400, fill = FALSE,
            card_header("Rhythmicity values for the selected protein(s):"),
            dataTableOutput("protein_rhythm_table"),
            card_footer(HTML("<b>Bolded</b> values indicate significant rhythmicity at the <i>p</i> < 0.05 threshold."))
          )
        ),
        conditionalPanel(
          condition = "input.proteins != ''",
          card(
            min_height = 500, fill = FALSE,
            card_header("Abundance values for the selected protein(s):"),
            dataTableOutput("protein_table"),
            card_footer(downloadButton("protein_table_dl", "Download table of values"))
          )
        ),
        p("Proteome data source: Kay et al. Comms Bio 2021")
      )
    )
  ),

  ## Rhythmic ion-binding proteins ----
  nav_panel(
    title = "Rhythmic ion-binding proteins",
    conditionalPanel(
      condition = "input.rhythmics == ''",
      card(p("Select one or more proteins to get started!"))
    ),
    conditionalPanel(
      condition = "input.rhythmics_go != 0",
      card(
        min_height = 650,
        card_header("Here is your graph!"),
        conditionalPanel(
          condition = "input.rhythmics != ''",
          conditionalPanel(
            condition = "input.rhythmics_abs_rel == 'abs'",
            withSpinner(plotOutput(width = 800, height = 500, outputId = "rhythmics_graph_abs"))
          ),
          conditionalPanel(
            condition = "input.rhythmics_abs_rel == 'rel'",
            withSpinner(plotOutput(width = 800, height = 500, outputId = "rhythmics_graph_rel"))
          )
        ),
        card_footer(
          conditionalPanel(
            condition = "input.rhythmics_abs_rel == 'abs'",
            downloadButton("rhythmics_graph_abs_dl", "Download graph view")
          ),
          conditionalPanel(
            condition = "input.rhythmics_abs_rel == 'rel'",
            downloadButton("rhythmics_graph_rel_dl", "Download graph view")
          )
        )
      ),
      conditionalPanel(
        condition = "input.rhythmics_go != 0",
        card(
          min_height = 400,
          card_header("Information on selected protein(s):"),
          conditionalPanel(
            condition = "input.rhythmics != ''",
            withSpinner(dataTableOutput("rhythmics_info_table"))
          )
        ),
        card(
          min_height = 400,
          card_header("Abundance values for selected protein(s):"),
          conditionalPanel(
            condition = "input.rhythmics != ''",
            withSpinner(dataTableOutput("rhythmics_table"))
          ),
          card_footer(downloadButton("rhythmics_table_dl", "Download table of values"))
        )
      ),
      # withSpinner(plotOutput(width = 800, height = 500, outputId = "abundance_graph")),
      p("Source for the proteome data: Kay et al. Comms Bio 2021")
    )
  ) # last navpanel
  # Other menu items ----

  # nav_item(input_dark_mode(id = "mode")) #dark mode toggle
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
    distinct(select(ion_data(), -c(Abundance))) %>% # condense down for better visibility
      arrange(Time) %>% rename("Relative abundance" = rel_abun)
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
  ion_abs_plot <- reactive({
    ggplot(data = ion_data(), aes(x = Time, y = Abundance, group = Ion)) +
      geom_rect(
        data = time_shading,
        aes(
          y = NULL, group = NULL,
          xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour
        ),
        show.legend = FALSE
      ) + # set up background
      scale_fill_identity(time_shading$colour) + # shade appropriately
      new_scale_fill() + # allow rest of graph to use different shading method
      
      # geom_smooth(aes(colour = Ion, fill = Ion), span = 0.2, se = FALSE) +
      geom_borderline(
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
      scale_x_continuous(
        breaks = breaks_width(12), minor_breaks = NULL,
        limits = c(0, 96)
      ) +
      scale_x_break(c(26, 48)) +
      labs(
        x = "Time (h)", y = "Ion abundance (ug/L)",
        title = "Ion abundance"
      )
  })
  
  output$ion_graph_abs <- renderPlot({
    ion_abs_plot()
  })
  
  
  ### Relative abundance (best for multiple ions) ----
  ion_rel_plot <- reactive({
    ggplot(data = ion_data(), aes(x = Time, y = rel_abun, group = Ion)) +
      geom_rect(
        data = time_shading,
        aes(
          y = NULL, group = NULL,
          xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour
        ),
        show.legend = FALSE
      ) + # set up background
      scale_fill_identity(time_shading$colour) + # shade appropriately
      new_scale_fill() + # allow rest of graph to use different shading method
      
      # geom_smooth(aes(colour = Ion, fill = Ion), span = 0.2, se = FALSE) +
      geom_borderline(
        linewidth = 0.5, aes(colour = Ion, group = interaction(Ion, Condition)),
        bordercolour = "black"
      ) +
      geom_point(size = 2, aes(fill = Ion), shape = 21, colour = "black") +
      scale_colour_manual(values = ion_colours()) +
      scale_fill_manual(values = ion_colours()) +
      scale_x_continuous(
        breaks = breaks_width(12), minor_breaks = NULL,
        limits = c(0, 96)
      ) +
      scale_x_break(c(26, 48)) +
      labs(
        x = "Time (h)", y = "Relative ion abundance",
        title = "Ion abundance"
      )
  })
  
  output$ion_graph_rel <- renderPlot({
    ion_rel_plot()
  })
  
  ## Consistent colour palette ----
  ion_colours <- eventReactive(input$ion_go, {
    ion_palette <- viridis(n = 3, option = input$ion_graph_palette) # set up based on user selection
    ion_colours <- c("Mg" = ion_palette[1], "Ca" = ion_palette[2], "K" = ion_palette[3]) # assign colours specifically to each ion based on current palette, so colours are consistent across plots
  })
  
  # output$ion_colours <- renderText({ion_colours()})
  
  
  ## Rhythmicity display ----
  
  output$ion_rain_echo <- reactive({
    input$ion_rain_echo
  })
  # outputOptions(output, "ion_rain_echo", suspendWhenHidden = FALSE)
  
  ion_rhythm_pvalues <- eventReactive(input$ion_go, {
    if (input$ion_rain_echo == "rain") {
      ion_rhythm_pvalues <- filter(ions_rain_LL, Ion %in% input$ions) %>%
        select("Ion", "P-Value" = "pVal", "Phase" = "phase", "Period" = "period")
    }
    if (input$ion_rain_echo == "echo") {
      ion_rhythm_pvalues <- filter(ions_echo_LL, Ion %in% input$ions) %>%
        select("Ion", "Oscillation Type", "Period", "P-Value":"BY Adj P-Value")
    }
    
    ion_rhythm_pvalues
  })
  
  output$ion_rhythm_table <- renderDataTable({
    datatable(ion_rhythm_pvalues(), options = list(dom = "t"), rownames = FALSE) %>%
      formatStyle(-c(1),
                  fontStyle = styleInterval(c(0.05, 1), c("normal", "italic", "normal")),
                  colour = styleInterval(0.05, c("black", "darkgrey")) # format based on significance
      ) %>%
      formatSignif(-c(1), 5)
  })
  
  ## Downloads ----
  
  ### absolute abundance ----
  output$ion_graph_abs_dl <- downloadHandler(
    filename = function() {
      paste((format(Sys.time(), "%F_%H%M%S")),
            "-ions_absolute",
            ".png",
            sep = ""
      )
    },
    content = function(file) {
      ggsave(file,
             plot = ion_abs_plot(), device = "png",
             width = 8, height = 5
      ) # this becomes 2400x1500px
    }
  )
  ### relative abundance ----
  output$ion_graph_rel_dl <- downloadHandler(
    filename = function() {
      paste((format(Sys.time(), "%F_%H%M%S")),
            "-ions_relative", ".png",
            sep = ""
      )
    },
    content = function(file) {
      ggsave(file,
             plot = ion_rel_plot(), device = "png",
             width = 8, height = 5
      ) # this becomes 2400x1500px
    }
  )
  ### table of values ----
  output$ion_table_dl <- downloadHandler(
    filename = function() {
      paste((format(Sys.time(), "%F_%H%M%S")),
            "-ions", ".tsv",
            sep = ""
      )
    },
    content = function(file) {
      write_tsv(ion_data(), file)
    }
  )
  
  ## Mg2+ abundance trace for rhythmic ion-binder overlay ----
  mg_trace <- eventReactive(input$rhythmics_go, {
    mg_trace <- select(ion_abundance, "Time", contains("Mg")) %>% # create table with all necessary columns
      pivot_longer(cols = !Time, names_to = "Ion", values_to = "Abundance") %>% drop_na()
    mg_trace$Ion <- str_remove(mg_trace$Ion, "\\_.*") # clean up ion names
    mg_trace <- mutate(mg_trace,
                       Condition = case_when(Time <= 24 ~ "LD", Time >= 48 ~ "LL"),
                       .after = Time
    ) %>% # add condition column
      group_by(Time, Condition, Ion) %>%
      summarise(across(), mean = mean(Abundance), sd = sd(Abundance)) %>%
      group_by(Ion) %>%
      mutate(rel_abun = mean / max(mean)) # add other/stats columns
  })
  
  
  
  # Proteins ----
  
  ## update input options  ----
  updateSelectizeInput(session, "proteins",
                       choices = c("Select a protein" = "", protein_abundance["Identifier"]),
                       server = TRUE
  ) # fill options
  
  
  # get selected proteins from list (old)
  protein_with_desc <- reactive({
    req(input$proteins)
    paste(select(nuc_proteome, Identifier, description) %>% filter(Identifier == input$proteins), collapse = " - ")
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
                           ), .after = Time
    ) %>%
      group_by(Identifier) %>%
      mutate(rel_abun = Abundance / max(Abundance, na.rm = TRUE)) %>%
      ungroup() # add extra columns
  })
  
  ### store history of selected proteins ----
  
  history <- reactiveValues(protein = NULL, rhythmics = NULL) # initialise empty list on startup
  
  observeEvent(input$protein_go, {
    history$protein <- unique(c(history$protein, names(protein_colours()))) # update list when graph is drawn
  })
  
  # output list as a test
  output$protein_history <- renderText({
    paste("Previously-selected proteins: \n")
    paste(history$protein, collapse = ", \n") # display list of previous proteins
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
    protein_colours <- unlist(split(
      protein_colour_table$Colour,
      protein_colour_table$Identifier
    ))
  })
  
  output$protein_colours <- renderText({
    protein_colours()
  })
  
  ### Absolute abundance ----
  
  protein_abs_plot <- eventReactive(input$protein_go, {
    ggplot(data = protein_data(), aes(x = Time, y = Abundance, group = Identifier)) +
      geom_rect(
        data = time_shading,
        aes(
          y = NULL, group = NULL,
          xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour
        ),
        show.legend = FALSE
      ) + # set up background
      scale_fill_identity(time_shading$colour) + # shade appropriately
      new_scale_fill() + # allow rest of graph to use a different shading method
      
      # geom_smooth(aes(colour = Identifier, fill = Identifier), span = 0.2, se = FALSE) +
      geom_borderline(
        aes(
          y = Abundance, colour = Identifier,
          group = interaction(Identifier, Condition)
        ),
        linewidth = 0.5, bordercolour = "black"
      ) +
      geom_point(
        size = 2, aes(x = Time, y = Abundance, group = Identifier, fill = Identifier),
        shape = 21, colour = "black"
      ) +
      scale_colour_manual(values = protein_colours()) +
      scale_fill_manual(values = protein_colours()) +
      scale_x_continuous(breaks = breaks_width(12)) +
      scale_x_break(c(26, 48)) +
      labs(
        x = "Time (h)", y = "Protein abundance",
        title = "Protein abundance"
      )
  })
  
  output$protein_graph_abs <- renderPlot({
    protein_abs_plot()
  })
  
  ### Relative abundance ----
  protein_rel_plot <- eventReactive(input$protein_go, {
    ggplot(data = protein_data(), aes(x = Time, y = rel_abun, group = Identifier)) +
      geom_rect(
        data = time_shading,
        aes(
          y = NULL, group = NULL,
          xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour
        ),
        show.legend = FALSE
      ) + # set up background
      scale_fill_identity(time_shading$colour) + # shade appropriately
      new_scale_fill() + # allow rest of graph to use different shading method
      
      # geom_smooth(aes(colour = Identifier, fill = Identifier), span = 0.2, se = FALSE) +
      geom_borderline(
        aes(
          y = rel_abun, colour = Identifier,
          group = interaction(Identifier, Condition)
        ),
        linewidth = 0.5, bordercolour = "black"
      ) +
      geom_point(
        size = 2, aes(x = Time, y = rel_abun, group = Identifier, fill = Identifier),
        shape = 21, colour = "black"
      ) +
      scale_colour_manual(values = protein_colours()) +
      scale_fill_manual(values = protein_colours()) +
      scale_x_continuous(breaks = breaks_width(12)) +
      scale_x_break(c(26, 48)) +
      labs(
        x = "Time (h)", y = "Relative abundance",
        title = "Protein abundance"
      ) +
      ylim(NA, 1)
  })
  
  output$protein_graph_rel <- renderPlot({
    protein_rel_plot()
  })
  
  
  ## show rhythmicity ----
  
  # get rhythmicity p-values for selected proteins
  protein_chosen_pvals <- reactive({
    req(input$rhythm_method)
    protein_chosen_pvals <- filter(protein_pvals, Identifier %in% input$proteins) %>%
      rename(c(
        "eJTK (LD)" = "LD_ejtk",
        "RAIN (LD)" = "LD_rain",
        "ECHO (LD)" = "LD_echo",
        "eJTK (LL)" = "LL_ejtk",
        "RAIN (LL)" = "LL_rain",
        "ECHO (LL)" = "LL_echo"
      )) %>%
      select(Identifier, contains(input$rhythm_method)) %>%
      relocate(ends_with("(LD)"), .after = Identifier) # reorder
  })
  
  # show values for selected proteins
  output$protein_rhythm_table <- renderDataTable({
    datatable(protein_chosen_pvals(), options = list(dom = "t"), rownames = FALSE) %>%
      formatStyle(-c(1),
                  # fontStyle = styleInterval(c(0.05, 1), c("normal", "italic", "normal")),
                  fontWeight = styleInterval(0.05, c("bold", "normal")),
                  colour = styleInterval(0.05, c("black", "darkgrey"))
      ) %>%
      formatSignif(-c(1), 5)
  })
  
  ## Downloads ----
  
  ### absolute abundance ----
  output$protein_graph_abs_dl <- downloadHandler(
    filename = function() {
      paste((format(Sys.time(), "%F_%H%M%S")),
            "-proteins_absolute",
            ".png",
            sep = ""
      )
    },
    content = function(file) {
      ggsave(file,
             plot = protein_abs_plot(), device = "png",
             width = 8, height = 5
      ) # this becomes 2400x1500px
    }
  )
  ### relative abundance ----
  output$protein_graph_rel_dl <- downloadHandler(
    filename = function() {
      paste((format(Sys.time(), "%F_%H%M%S")),
            "-proteins_relative", ".png",
            sep = ""
      )
    },
    content = function(file) {
      ggsave(file,
             plot = protein_rel_plot(), device = "png",
             width = 8, height = 5
      ) # this becomes 2400x1500px
    }
  )
  ### table of values ----
  output$protein_table_dl <- downloadHandler(
    filename = function() {
      paste((format(Sys.time(), "%F_%H%M%S")),
            "-proteins", ".tsv",
            sep = ""
      )
    },
    content = function(file) {
      write_tsv(protein_data(), file)
    }
  )
  
  
  
  # Rhythmic Mg2+ binding proteins ----
  
  ## update input options ----
  updateSelectizeInput(session, "rhythmics",
                       choices = c("Select a protein" = "", rhythmic_ion_binders_full["Identifier"]),
                       server = TRUE
  )
  
  ## show info table ----
  rhythmics_info_table <- eventReactive(input$rhythmics_go, {
    filter(mg_co_rhythmics, Identifier %in% (input$rhythmics)) %>%
      select("Identifier", "Protein_names", "Peak_time_of_day", "Gene_Ontology_IDs")
  })
  
  output$rhythmics_info_table <- renderDataTable({
    req(input$rhythmics)
    rhythmics_info_table()
  })
  
  ## show abundance values ----
  rhythmics_table <- eventReactive(input$rhythmics_go, {
    rhythmics_data()
  })
  
  output$rhythmics_table <- renderDataTable({
    req(input$rhythmics)
    rhythmics_table() %>% 
      arrange(Time) %>% 
      rename("Relative abundance" = rel_abun)
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
      mutate(rel_abun = Abundance / max(Abundance, na.rm = TRUE)) %>%
      ungroup()
  })
  
  ### consistent graph colour palette ----
  rhythmics_colours <- eventReactive(input$rhythmics_go, {
    rhythmics_palette <- viridis(
      n = length(input$rhythmics),
      option = input$rhythmics_graph_palette
    ) # set up colours to assign based on user selection
    
    for (id_num in 1:length(input$rhythmics)) {
      rhythmics_colour_table <- add_row(rhythmics_colour_table,
                                        Identifier = (input$rhythmics)[id_num],
                                        Colour = rhythmics_palette[id_num]
      ) # create palette one row at a time
    }
    rhythmics_colours <- unlist(split(
      rhythmics_colour_table$Colour,
      rhythmics_colour_table$Identifier
    ))
  })
  
  ### store history of previous selections ----
  # observeEvent(input$rhythmics_go, {
  #   history$rhythmics <- unique(c(history$rhythmics, names(rhythmics_colours()) )) # update list when graph is drawn
  # })
  
  # output list as a test
  output$rhythmics_history <- renderText({
    paste("Previously-selected proteins: \n")
    paste(history$rhythmics, collapse = ", \n") # display list of previous proteins
  })
  
  ## plot relevant values ----
  
  ### Absolute abundance ----
  rhythmics_abs_plot <- eventReactive(input$rhythmics_go, {
    ggplot(data = rhythmics_data(), aes(x = Time, y = Abundance, group = Identifier)) +
      geom_rect(
        data = time_shading,
        aes(
          y = NULL, group = NULL,
          xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour
        ),
        show.legend = FALSE
      ) + # set up background
      scale_fill_identity(time_shading$colour) + # shade appropriately
      new_scale_fill() + # allow rest of graph to use different shading method
      
      # geom_smooth(aes(colour = Identifier, fill = Identifier), span = 0.2, se = FALSE) +
      geom_borderline(
        aes(
          y = Abundance, colour = Identifier,
          group = interaction(Identifier, Condition)
        ),
        linewidth = 0.5, bordercolour = "black"
      ) +
      geom_point(
        size = 2, aes(x = Time, y = Abundance, group = Identifier, fill = Identifier),
        shape = 21, colour = "black"
      ) +
      scale_colour_manual(values = rhythmics_colours()) +
      scale_fill_manual(values = rhythmics_colours()) +
      scale_x_continuous(breaks = breaks_width(12)) +
      scale_x_break(c(26, 48)) +
      labs(
        x = "Time (h)", y = "Abundance",
        title = "Rhythmic Mg2+ binding proteins"
      )
  })
  
  output$rhythmics_graph_abs <- renderPlot({
    rhythmics_abs_plot()
  })
  
  
  
  ### Relative abundance ----
  rhythmics_rel_plot <- eventReactive(input$rhythmics_go, {
    if (input$add_mg_trace == FALSE) {
      rhythmics_rel_plot <- ggplot(data = rhythmics_data(), aes(x = Time, y = rel_abun, group = Identifier)) +
        geom_rect(
          data = time_shading, aes(
            y = NULL, group = NULL,
            xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour
          ),
          show.legend = FALSE
        ) + # set up background
        scale_fill_identity(time_shading$colour) + # shade appropriately
        new_scale_fill() + # allow rest of graph to use different shading method
        
        ## geom_smooth(aes(colour = Identifier, fill = Identifier), span = 0.2, se = FALSE) +
        geom_borderline(
          aes(
            y = rel_abun, colour = Identifier,
            group = interaction(Identifier, Condition)
          ),
          linewidth = 0.5, bordercolour = "black"
        ) +
        geom_point(
          size = 2, aes(x = Time, y = rel_abun, group = Identifier, fill = Identifier),
          shape = 21, colour = "black"
        ) +
        labs(x = "Time (h)", y = "Relative abundance",
             title = "Rhythmic Mg2+ binding proteins") +
        scale_colour_manual(values = rhythmics_colours()) +
        scale_fill_manual(values = rhythmics_colours()) +
        scale_x_continuous(breaks = breaks_width(12)) +
        scale_x_break(c(26, 48)) +
        ylim(NA, 1)
      return(rhythmics_rel_plot)
    }
    
    #### with trace ----
    if (input$add_mg_trace == TRUE) {
      rhythmics_rel_plot <- ggplot(data = rhythmics_data(), aes(x = Time, y = rel_abun, group = Identifier)) +
        geom_rect(
          data = time_shading, aes(
            y = NULL, group = NULL,
            xmin = Time, xmax = Time + 1, ymin = -Inf, ymax = Inf, fill = colour
          ),
          show.legend = FALSE
        ) + # set up background
        scale_fill_identity(time_shading$colour) + # shade appropriately
        new_scale_fill() + # allow rest of graph to use different shading method
        
        geom_line(
          data = mg_trace(), aes(x = Time, y = rel_abun, group = NULL, colour = Ion),
          linetype = "dashed", alpha = 0.8
        ) +
        geom_point(
          data = mg_trace(), aes(x = Time, y = rel_abun, group = NULL, colour = Ion),
          shape = 4
        ) +
        scale_colour_manual(values = "black") +
        new_scale_colour() +
        
        ## geom_smooth(aes(colour = Identifier, fill = Identifier), span = 0.2, se = FALSE) +
        geom_borderline(
          aes(
            y = rel_abun, colour = Identifier,
            group = interaction(Identifier, Condition)
          ),
          linewidth = 0.5, bordercolour = "black"
        ) +
        geom_point(
          size = 2, aes(x = Time, y = rel_abun, group = Identifier, fill = Identifier),
          shape = 21, colour = "black"
        ) +
        labs(x = "Time (h)", y = "Relative abundance",
             title = "Rhythmic Mg2+ binding proteins") +
        scale_colour_manual(values = rhythmics_colours()) +
        scale_fill_manual(values = rhythmics_colours()) +
        scale_x_continuous(limits = c(0, 96), breaks = breaks_width(12)) +
        scale_x_break(c(26, 48)) +
        ylim(NA, 1)
      return(rhythmics_rel_plot)
    }
  })
  
  output$rhythmics_graph_rel <- renderPlot({
    rhythmics_rel_plot()
  })
  
 
  
  
  
  ## Downloads ----
  
  ### absolute abundance ----
  output$rhythmics_graph_abs_dl <- downloadHandler(
    filename = function() {
      paste((format(Sys.time(), "%F_%H%M%S")),
            "-rhythmics_absolute",
            ".png",
            sep = ""
      )
    },
    content = function(file) {
      ggsave(file,
             plot = rhythmics_abs_plot(), device = "png",
             width = 8, height = 5
      ) # this becomes 2400x1500px
    }
  )
  ### relative abundance ----
  output$rhythmics_graph_rel_dl <- downloadHandler(
    filename = function() {
      paste((format(Sys.time(), "%F_%H%M%S")),
            "-rhythmics_relative", ".png",
            sep = ""
      )
    },
    content = function(file) {
      ggsave(file,
             plot = rhythmics_rel_plot(), device = "png",
             width = 8, height = 5
      ) # this becomes 2400x1500px
    }
  )
  ### table of values ----
  output$rhythmics_table_dl <- downloadHandler(
    filename = function() {
      paste((format(Sys.time(), "%F_%H%M%S")),
            "-rhythmics", ".tsv",
            sep = ""
      )
    },
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
