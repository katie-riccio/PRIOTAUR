library(bslib)
library(ggplot2)
library(shiny)
library(readxl)
library(tibble)
library(tidyverse)

nuc_proteome <- read_excel(path = "Kay et al Comms Bio 2021/42003_2021_2680_MOESM4_ESM.xlsx", sheet = 2) # nuclear-encoded proteome
org_proteome <- read_excel(path = "Kay et al Comms Bio 2021/42003_2021_2680_MOESM4_ESM.xlsx", sheet = 3) # organellar-encoded proteome

head(nuc_proteome) # check nuclear proteome
head(org_proteome) # check organellar proteome

nuc_abundance <- nuc_proteome[23:44] # extract abundance values at each time point
nuc_identifiers <- nuc_proteome[1] # extract gene identifiers
nuc_abundance <- add_column(nuc_abundance, nuc_identifiers, .before = 1) # add gene identifiers 

nuc_abundance[1] # check that this is gene identifiers like it should be

test_gene <- "ostta15g00570" # TEST

test_frame <- ((subset(nuc_abundance, Identifier == test_gene)) # extract desired row of data
               %>% pivot_longer(cols = !Identifier, names_to = "Time", values_to = "Abundance")) # make long
test_frame

g <- ggplot(test_frame) + geom_point(mapping = aes(x = Time, y = Abundance))

g


# SHINY

# Define UI
ui <- fluidPage(

  titlePanel("Kay et al. Comms Bio 2021"),
  
  sidebarLayout(
    
    sidebarPanel(
      selectizeInput(input = "gene", label = "Select a gene", c(Choose='', nuc_identifiers) )
      ),
    
    mainPanel(
      plotOutput("graph"),
      verbatimTextOutput("gene")
    )
    
  )
  
)

# Define server logic
server <- function(input, output) {
  
  # Generate a graph for selected gene (WIP)
  output$graph <- renderPlot({ # below: extract desired row of plot and make longer
    data <- (subset(nuc_abundance, Identifier == input$gene)
      ) %>% pivot_longer(cols = !Identifier, names_to = "Time", values_to = "Abundance")
    
    data %>% ggplot(aes(x = Time, y = Abundance))  + geom_point() # plot
    
  })
  
  output$gene <- renderPrint(
    input$gene
    )
  
}

# Create Shiny app
shinyApp(ui = ui, server = server)

