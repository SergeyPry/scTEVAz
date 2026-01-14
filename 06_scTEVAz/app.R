library(shiny)
library(bslib)
library(shinyjs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(paletteer)

# define relevant data frames

#################################### Daniocell files #######################
Daniocell_genes <- read.csv("daniocell_map_df.csv", header = TRUE)

# Daniocell tissues
Daniocell_tissues <- read.csv("Daniocell_tissues.txt", header = TRUE )

# Daniocell stages
Daniocell_stages <- read.csv("Daniocell_stages.txt", header = TRUE )


#################################### ZCL files ############################
zcl_genes <- read.csv("ZCL_map_df.csv", header = TRUE)

zcl_tissues <- read.csv("zcl_cell-lineages.txt", header = TRUE )

zcl_cell_types <- read.csv("zcl_cell-types.txt", header = TRUE )

zcl_stages <- read.csv("zcl_stages.txt", header = TRUE )


#################################### zebrahub files ############################
zhub_genes <- read.csv("zhub_map_df.csv", header = TRUE)

zhub_tissues <- read.csv("zhub_cell-lineages.txt", header = TRUE )

zhub_stages <- read.csv("zhub_stages.txt", header = TRUE )



ui <- fluidPage(
  useShinyjs(), # Must include this line in the UI
  
  HTML('<div class="row"><div style="background: #152437; margin-left: 5px; margin-right: 5px; height: 100px; text-indent: 10px; line-height: 80px; font-size: 35px; text-align: center; color: white; text-transform: none; "><img src="logo.png" style="float:left; height="100px"; width="148px"> scTEVAz - single-cell Tissue Expression Value Aggregator for zebrafish</div></div>'),
  
  navbarPage(
    
    id = "nav", 
    
    title = "scTEVAz App",
    
    
    # The Landing Page (Home)
    tabPanel(
      
      title = "Home",
      icon = icon("home"),
      
      div(
        style = "padding: 20px; text-align: center;",
        h1("Welcome to the scTEVAz (single-cell Tissue Expression Value Aggregator for zebrafish)"),
        h3("This tool contains data from three scRNA-seq datasets:"),
        
        br(),
        actionButton("goToTool1", h4("Daniocell"), class = "btn-lg btn-primary"),
        actionButton("goToTool2", h4("Zebrafish Cell Landscape"), class = "btn-lg btn-success"),
        actionButton("goToTool3", h4("Zebrahub"), class = "btn-lg btn-info")
      ),
      
      
      div(
        style = "text-align: center;",
        img(src = "logo_video.gif", height = 500, width = 500)
      )

    ),
    
    # Daniocell panel
    tabPanel(
      
      title = "Daniocell", 
      
      tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"),
                tags$style(HTML(".shiny-output-error-validation {color: green;font-size: 20px;}")),
                tags$style(HTML(".checkbox { font-size: 18px; }")),
                tags$style(HTML(".radio { font-size: 18px; }")),
                tags$style(HTML(".radio-inline { font-size: 18px; }"))
                ),
      
      tags$style("body {
                          -moz-transform: scale(0.95 0.95); /* Moz-browsers */
                          zoom: 0.95; /* Other non-webkit browsers */
                          zoom: 95%; /* Webkit browsers */
                          }"),
      
      
      # Sidebar layout for inputs and output
      sidebarLayout(
        # Sidebar Panel - Inputs
        sidebarPanel(
          width = 4,
          class = "sidebar",
          h3(strong("Input Parameters"), style = "color: #1e3a8a;"),
          
          tags$label(h4(strong("Please input gene symbol:"))),
          # Input 1: Gene Selection
          selectizeInput(
            inputId = "daniocell_gene_id", 
            label = NULL,
            choices = NULL, 
            selected = NULL,
            multiple = FALSE,
            width = '200px'
          ),
      
          
          # Input 2: Cell Type Filtering
          tags$label(h4(strong("Filter by tissue:"))),
          
          fluidRow(column(width = 5, 
                          actionLink("deselectall_daniocell","Deselect All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);"),
                          actionLink("selectall_daniocell","Select All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);")
                          )),
          tags$br(),
          
          fluidRow(
            column(width = 3,
                   
                   checkboxGroupInput(
                     inputId = "cell_types_1",
                     label = NULL,
                     choices = Daniocell_tissues$tissue[1:11],
                     selected = Daniocell_tissues$tissue[1:11],
                   ) 
                ),
            column(width = 3,
                   
                   checkboxGroupInput(
                     inputId = "cell_types_2",
                     label = NULL,
                     choices = Daniocell_tissues$tissue[12:22],
                     selected = Daniocell_tissues$tissue[12:22],
                   ) 
            ),
            
            column(width = 3,
                   
                   checkboxGroupInput(
                     inputId = "cell_types_3",
                     label = NULL,
                     choices = Daniocell_tissues$tissue[23:33],
                     selected = Daniocell_tissues$tissue[23:33],
                   )
            ),
            
            column(width = 3,
                   
                   checkboxGroupInput(
                     inputId = "cell_types_4",
                     label = NULL,
                     choices = Daniocell_tissues$tissue[34:43],
                     selected = Daniocell_tissues$tissue[34:43],
                   )
            )
          ),
          
          
          
  
          # Input 3: Stage Filtering
          tags$label(h4(strong("Filter by stage:"))),
  
          fluidRow(column(width = 5, 
                          actionLink("deselect_dcell_stages","Deselect All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);"),
                          actionLink("select_dcell_stages","Select All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);")
          )),
          tags$br(),
          
          
                  
          fluidRow(
            column(width = 4,
                   
                   checkboxGroupInput(
                     inputId = "stages_1",
                     label = NULL,
                     choices = Daniocell_stages$stage[1:5],
                     selected = Daniocell_stages$stage[1:5],
                   )
            ),
            column(width = 4,
                   
                   checkboxGroupInput(
                     inputId = "stages_2",
                     label = NULL,
                     choices = Daniocell_stages$stage[6:10],
                     selected = Daniocell_stages$stage[6:10],
                   )
            ),
            
            column(width = 4,
                   
                   checkboxGroupInput(
                     inputId = "stages_3",
                     label = NULL,
                     choices = Daniocell_stages$stage[11:14],
                     selected = Daniocell_stages$stage[11:14],
                   )
            )
          ),
      
      radioButtons("data_type", h4(strong("Choose the type of data:")),
                       c("Sum of normalized read counts" = "count_sum",
                         "Normalized cell counts" = "cell_count"),
                   selected = "count_sum"),
      
      radioButtons("y_scale", h4(strong("Choose the y-axis scale:")),
                   c("Linear" = "linear",
                     "Square root" = "sqrt"),
                   selected = "linear",
                   inline = TRUE),
      
            
          
      # Action Button
      actionButton("Daniocell_update_plot", "Generate Plot", class = "bg-blue-600 text-white hover:bg-blue-700")
          
      ), # End of sidebarPanel
        
        # Main Panel - Output
        mainPanel(
          style = "margin-top: -20px;",
          class = "main-panel",
          tags$strong(h3("Daniocell Aggregated expression values plot", style = "width: 100%; color: #0c4a6e; background: #f8f5f0; text-align: center; border: 1px solid black; padding: 10px; display: inline-block;")),
   
          # plotting the resulting graph
          plotOutput("daniocell_expr_plot")
          
          
        ) # End of mainPanel
      ) # End of sidebarLayout
      
    ),
    
    # ZCL Panel
    tabPanel(
      
      title = "Zebrafish Cell Landscape",
      
      tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"),
                tags$style(HTML(".shiny-output-error-validation {color: green;font-size: 20px;}"))),
      
      tags$style(HTML("
      .multicol {
        -webkit-column-count: 5; /* Chrome, Safari, Opera */
        -moz-column-count: 5;    /* Firefox */
        column-count: 5;
        -moz-column-fill: balanced;
        column-fill: balanced;
      }
      /* Optional: Adjust vertical spacing */
      .checkbox {
        margin-top: 0px !important;
        -webkit-margin-after: 0px !important;
      }
    ")),
      
      tags$style("body {
                          -moz-transform: scale(0.95 0.95); /* Moz-browsers */
                          zoom: 0.95; /* Other non-webkit browsers */
                          zoom: 95%; /* Webkit browsers */
                          }"),
      
      
      # Sidebar layout for inputs and output
      sidebarLayout(
        # Sidebar Panel - Inputs
        sidebarPanel(
          class = "sidebar",
          h3(strong("Input Parameters"), style = "color: #1e3a8a;"),
          
          tags$label(h4(strong("Please input gene symbol:"))),
          # Input 1: Gene Selection
          selectizeInput(
            inputId = "zcl_gene_id", 
            label = NULL,
            choices = NULL, 
            selected = NULL,
            multiple = FALSE,
            width = '200px'
          ),
          
          
          # Action buttons to control visibility
          actionButton("show_lineage", "Set up lineages plot"),
          actionButton("show_celltypes", "Set up cell types plot"),
          
          # Wrap checkbox groups in div tags with unique IDs
          # Initially hide the second group
          hidden(
            div(id = "lineage_wrapper",
              
              # Input 2: Cell Type Filtering
              tags$label(h4(strong("Filter by lineage:"))),
              
              fluidRow(column(width = 5, 
                              actionLink("deselectall_zcl","Deselect All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);"),
                              actionLink("selectall_zcl","Select All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);")
              )),
              tags$br(),
              
                            
              fluidRow(
                column(width = 3,
                       
                       checkboxGroupInput(
                         inputId = "zcl_lineages_1",
                         label = NULL,
                         choices = zcl_tissues$tissue[1:5],
                         selected = zcl_tissues$tissue[1:5],
                       ) 
                ),
                column(width = 3,
                       
                       checkboxGroupInput(
                         inputId = "zcl_lineages_2",
                         label = NULL,
                         choices = zcl_tissues$tissue[6:10],
                         selected = zcl_tissues$tissue[6:10],
                       ) 
                )
              ),
              
              tags$br(),
              
              # Input 3: Stage Filtering
              tags$label(h4(strong("Filter by stage:"))),
              
              fluidRow(
                column(width = 4,
                       
                       checkboxGroupInput( 
                         inputId = "lineage_zcl_stages",
                         label = NULL,
                         choices = zcl_stages$stage[1:5],
                         selected = zcl_stages$stage[1:5],
                       )
                )
              ),
              
              radioButtons("lineage_data_type", h4(strong("Choose the type of data:")),
                           c("Sum of normalized read counts" = "count_sum",
                             "Normalized cell counts" = "cell_count"),
                           selected = "count_sum"),
              
              radioButtons("lineage_y_scale", h4(strong("Choose the y-axis scale:")),
                           c("Linear" = "linear",
                             "Square root" = "sqrt"),
                           selected = "linear",
                           inline = TRUE),
              
              
              # Action Button
              actionButton("lineages_plot_zcl", "Generate Lineages Plot", class = "bg-blue-600 text-white hover:bg-blue-700")
              
              
            ) #end of lineage div
          ),
          
          hidden(
            div(id = "celltypes_wrapper",
                
                # Input 2: Cell Type Filtering
                tags$label(h4(strong("Filter by cell types:"))),
                
                fluidRow(column(width = 5, 
                                actionLink("deselectall_zcl_celltypes","Deselect All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);"),
                                actionLink("selectall_zcl_celltypes","Select All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);")
                )),
                tags$br(),
                
                
                fluidRow(
                  column(width = 3,
                         
                         checkboxGroupInput(
                           inputId = "zcl_cell_types_1",
                           label = NULL,
                           choices = zcl_cell_types$cell_type[1:11],
                           selected = zcl_cell_types$cell_type[1:11],
                         ) 
                  ),
                  column(width = 3,
                         
                         checkboxGroupInput(
                           inputId = "zcl_cell_types_2",
                           label = NULL,
                           choices = zcl_cell_types$cell_type[12:22],
                           selected = zcl_cell_types$cell_type[12:22],
                         ) 
                  ),
                  
                  column(width = 3,
                         
                         checkboxGroupInput(
                           inputId = "zcl_cell_types_3",
                           label = NULL,
                           choices = zcl_cell_types$cell_type[23:33],
                           selected = zcl_cell_types$cell_type[23:33],
                         ) 
                  ),
                  
                  column(width = 3,
                         
                         checkboxGroupInput(
                           inputId = "zcl_cell_types_4",
                           label = NULL,
                           choices = zcl_cell_types$cell_type[34:41],
                           selected = zcl_cell_types$cell_type[34:41],
                         ) 
                  )
                ),
                
                tags$br(),
                
                # Input 3: Stage Filtering
                tags$label(h4(strong("Filter by stage:"))),
                
                fluidRow(
                  column(width = 4,
                         
                         checkboxGroupInput(
                           inputId = "celltypes_zcl_stages",
                           label = NULL,
                           choices = zcl_stages$stage[1:5],
                           selected = zcl_stages$stage[1:5],
                         )
                  )
                ),
                
                radioButtons("celltypes_data_type", h4(strong("Choose the type of data:")),
                             c("Sum of normalized read counts" = "count_sum",
                               "Normalized cell counts" = "cell_count"),
                             selected = "count_sum"),
                
                radioButtons("celltypes_y_scale", h4(strong("Choose the y-axis scale:")),
                             c("Linear" = "linear",
                               "Square root" = "sqrt"),
                             selected = "linear",
                             inline = TRUE),
                
                
                # Action Button
                actionButton("celltypes_plot_zcl", "Generate Cell types Plot", class = "bg-blue-600 text-white hover:bg-blue-700")
                
                
            ) # end of celltypes_wrapper div
          ), 
          
        ), # End of sidebarPanel
        
        # Main Panel - Output
        mainPanel(
          style = "margin-top: -20px;",
          class = "main-panel",
          tags$strong(h3("ZCL Aggregated expression values plot", style = "width: 100%; color: #0c4a6e; background: #f8f5f0; text-align: center; border: 1px solid black; padding: 10px; display: inline-block;")),
          
          # plotting the resulting graph
          plotOutput("zcl_expr_plot")
          
          # 
          # # Panel visible 
          # conditionalPanel(
          #   condition = "input.plotType == 'Histogram'",
          #   h3("Histogram Panel"),
          #   plotOutput("histPlot")
          # ),
          # # Panel visible only if "Scatter Plot" is selected
          # conditionalPanel(
          #   condition = "input.plotType == 'Scatter Plot'",
          #   h3("Scatter Plot Panel"),
          #   plotOutput("scatterPlot")
          # )
          
          
        ) # End of mainPanel
      ) # End of sidebarLayout
      
    ),
    
    # Zebrahub Panel
    tabPanel(
      
      title = "Zebrahub",
      
      tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"),
                tags$style(HTML(".shiny-output-error-validation {color: green;font-size: 20px;}"))),
      
      tags$style(HTML("
      .multicol {
        -webkit-column-count: 5; /* Chrome, Safari, Opera */
        -moz-column-count: 5;    /* Firefox */
        column-count: 5;
        -moz-column-fill: balanced;
        column-fill: balanced;
      }
      /* Optional: Adjust vertical spacing */
      .checkbox {
        margin-top: 0px !important;
        -webkit-margin-after: 0px !important;
      }
    ")),
      
      tags$style("body {
                          -moz-transform: scale(0.95 0.95); /* Moz-browsers */
                          zoom: 0.95; /* Other non-webkit browsers */
                          zoom: 95%; /* Webkit browsers */
                          }"),
      
      
      # Sidebar layout for inputs and output
      sidebarLayout(
        # Sidebar Panel - Inputs
        sidebarPanel(
          class = "sidebar",
          h3(strong("Input Parameters"), style = "color: #1e3a8a;"),
          
          tags$label(h4(strong("Please input gene symbol:"))),
          # Input 1: Gene Selection
          selectizeInput(
            inputId = "zhub_gene_id", 
            label = NULL,
            choices = NULL, 
            selected = NULL,
            multiple = FALSE,
            width = '200px'
          ),
          
          
          # Input 2: Cell Type Filtering
          tags$label(h4(strong("Filter by lineage:"))),
          
          fluidRow(column(width = 5, 
                          actionLink("deselectall_zhub","Deselect All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);"),
                          actionLink("selectall_zhub","Select All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);")
          )),
          tags$br(),
          
          fluidRow(
            column(width = 5,
                   
                   checkboxGroupInput(
                     inputId = "zhub_cell_types_1",
                     label = NULL,
                     choices = zhub_tissues$tissue[1:5],
                     selected = zhub_tissues$tissue[1:5],
                   ) 
            ),
            column(width = 5,
                   
                   checkboxGroupInput(
                     inputId = "zhub_cell_types_2",
                     label = NULL,
                     choices = zhub_tissues$tissue[6:10],
                     selected = zhub_tissues$tissue[6:10],
                   ) 
            )
          ),
          
          # Input 3: Stage Filtering
          tags$label(h4(strong("Filter by stage:"))),
          
          # select_zhub_stages
          fluidRow(column(width = 5, 
                          actionLink("deselect_zhub_stages","Deselect All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);"),
                          actionLink("select_zhub_stages","Select All", style = "font-size: 17px; background-color: #337ab7; color: white; padding: 10px 15px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 2px rgba(0,0,0,0.1);")
          )),
          tags$br(),
          
          
          fluidRow(
            column(width = 4,
                   
                   checkboxGroupInput(
                     inputId = "zhub_stages_1",
                     label = NULL,
                     choices = zhub_stages$stage[1:5],
                     selected = zhub_stages$stage[1:5],
                   )
            ),
            
            column(width = 4,
                   
                   checkboxGroupInput(
                     inputId = "zhub_stages_2",
                     label = NULL,
                     choices = zhub_stages$stage[6:10],
                     selected = zhub_stages$stage[6:10],
                   )
            )
          ),
          
          
          radioButtons("zhub_data_type", h4(strong("Choose the type of data:")),
                       c("Sum of normalized read counts" = "count_sum",
                         "Normalized cell counts" = "cell_count"),
                       selected = "count_sum"),
          
          radioButtons("zhub_y_scale", h4(strong("Choose the y-axis scale:")),
                       c("Linear" = "linear",
                         "Square root" = "sqrt"),
                       selected = "linear",
                       inline = TRUE),
          
          
          # Action Button
          actionButton("update_plot_zhub", "Generate Plot", class = "bg-blue-600 text-white hover:bg-blue-700")
          
        ), # End of sidebarPanel
        
        # Main Panel - Output
        mainPanel(
          style = "margin-top: -20px;",
          class = "main-panel",
          tags$strong(h3("Zebrahub Aggregated Expression Values plot", style = "width: 100%; color: #0c4a6e; background: #f8f5f0; text-align: center; border: 1px solid black; padding: 10px; display: inline-block;")),
          
          # plotting the resulting graph
          plotOutput("zebrahub_expr_plot")
        
          ) # End of mainPanel
      ) # End of sidebarLayout
      
    ),
  )
) 


server <- function(input, output, session) {
  
  # Update the selectizeInput with the data frame and enable server-side processing
  updateSelectizeInput(session, "daniocell_gene_id", 
                       choices = unique(Daniocell_genes$orig_gene), 
                       server = TRUE)
  
  # Update the selectizeInput with the data frame and enable server-side processing
  updateSelectizeInput(session, "zcl_gene_id", 
                       choices = unique(zcl_genes$orig_gene), 
                       server = TRUE)

  # Update the selectizeInput with the data frame and enable server-side processing
  updateSelectizeInput(session, "zhub_gene_id", 
                       choices = unique(zhub_genes$orig_gene), 
                       server = TRUE)
  

  # --- Navigation Logic ---
  
  # When 'Start with Tool 1' button is clicked
  observeEvent(input$goToTool1, {
    # inputId MUST match the id given to navbarPage in the UI ("nav")
    # selected MUST match the title of the target tabPanel ("Tool 1: Data Input")
    updateNavbarPage(session, inputId = "nav", selected = "Daniocell")
  })
  
  # When 'Explore Tool 2' button is clicked
  observeEvent(input$goToTool2, {
    updateNavbarPage(session, inputId = "nav", selected = "Zebrafish Cell Landscape")
  })
  
  # When 'Analyze Tool 3' button is clicked
  observeEvent(input$goToTool3, {
    updateNavbarPage(session, inputId = "nav", selected = "Zebrahub")
  })
  
  ########## UI modification ##########################################
  
  # Observe the button for Group 1
  observeEvent(input$show_lineage, {
    shinyjs::show("lineage_wrapper") # Show Group 1 wrapper
    shinyjs::hide("celltypes_wrapper") # Hide Group 2 wrapper
  })
  
  # Observe the button for Group 2
  observeEvent(input$show_celltypes, {
    shinyjs::hide("lineage_wrapper") # Hide Group 1 wrapper
    shinyjs::show("celltypes_wrapper") # Show Group 2 wrapper
  })

  #------------------------------- Lineages/cell types - select all / deselect all --------------
  
  observe({
    if(is.null(input$deselectall_daniocell)) return(NULL)
    
    if (input$deselectall_daniocell > 0) {
  
      updateCheckboxGroupInput(session, "cell_types_1", label = NULL,
        choices = Daniocell_tissues$tissue[1:11],selected = NULL)
      updateCheckboxGroupInput(session, "cell_types_2", label = NULL,
                               choices = Daniocell_tissues$tissue[12:22],selected = NULL)
      updateCheckboxGroupInput(session, "cell_types_3", label = NULL,
                               choices = Daniocell_tissues$tissue[23:33],selected = NULL)
      updateCheckboxGroupInput(session, "cell_types_4", label = NULL,
                               choices = Daniocell_tissues$tissue[34:43],selected = NULL)
      
    }
  })
  
  
  observe({
    if(is.null(input$selectall_daniocell)) return(NULL)
    
    if (input$selectall_daniocell > 0) {
      
      updateCheckboxGroupInput(session, "cell_types_1", label = NULL,
                               choices = Daniocell_tissues$tissue[1:11],selected = Daniocell_tissues$tissue[1:11] )
      updateCheckboxGroupInput(session, "cell_types_2", label = NULL,
                               choices = Daniocell_tissues$tissue[12:22],selected = Daniocell_tissues$tissue[12:22] )
      updateCheckboxGroupInput(session, "cell_types_3", label = NULL,
                               choices = Daniocell_tissues$tissue[23:33],selected = Daniocell_tissues$tissue[23:33] )
      updateCheckboxGroupInput(session, "cell_types_4", label = NULL,
                               choices = Daniocell_tissues$tissue[34:43],selected = Daniocell_tissues$tissue[34:43] )
      
    }
  })
  
  # zcl_lineages
  
  observe({
    if(is.null(input$deselectall_zcl)) return(NULL)
    
    if (input$deselectall_zcl > 0) {
      
      updateCheckboxGroupInput(session, "zcl_lineages_1", label = NULL,
                               choices = zcl_tissues$tissue[1:5],selected = NULL)
      updateCheckboxGroupInput(session, "zcl_lineages_2", label = NULL,
                               choices = zcl_tissues$tissue[6:10],selected = NULL)

    }
  })
  
  
  observe({
    if(is.null(input$selectall_zcl)) return(NULL)
    
    if (input$selectall_zcl > 0) {
      
      updateCheckboxGroupInput(session, "zcl_lineages_1", label = NULL,
                               choices = zcl_tissues$tissue[1:5],selected = zcl_tissues$tissue[1:5])
      
      updateCheckboxGroupInput(session, "zcl_lineages_2", label = NULL,
                               choices = zcl_tissues$tissue[6:10],selected = zcl_tissues$tissue[6:10])
      
    }
  })  
  
  

  # zcl cell types
  observe({
    if(is.null(input$deselectall_zcl_celltypes)) return(NULL)
    
    if (input$deselectall_zcl_celltypes > 0) {
      
      updateCheckboxGroupInput(session, "zcl_cell_types_1", label = NULL,
                               choices = zcl_cell_types$cell_type[1:11],selected = NULL)
      
      updateCheckboxGroupInput(session, "zcl_cell_types_2", label = NULL,
                               choices = zcl_cell_types$cell_type[12:22],selected = NULL)
      
      updateCheckboxGroupInput(session, "zcl_cell_types_3", label = NULL,
                               choices = zcl_cell_types$cell_type[23:33],selected = NULL)
      
      updateCheckboxGroupInput(session, "zcl_cell_types_4", label = NULL,
                               choices = zcl_cell_types$cell_type[34:41],selected = NULL)
      
    }
  })
  
  
  observe({
    if(is.null(input$selectall_zcl_celltypes)) return(NULL)
    
    if (input$selectall_zcl_celltypes > 0) {
      
      updateCheckboxGroupInput(session, "zcl_cell_types_1", label = NULL,
                               choices = zcl_cell_types$cell_type[1:11],selected = zcl_cell_types$cell_type[1:11])
      
      updateCheckboxGroupInput(session, "zcl_cell_types_2", label = NULL,
                               choices = zcl_cell_types$cell_type[12:22],selected = zcl_cell_types$cell_type[12:22])
      
      updateCheckboxGroupInput(session, "zcl_cell_types_3", label = NULL,
                               choices = zcl_cell_types$cell_type[23:33],selected = zcl_cell_types$cell_type[23:33])
      
      updateCheckboxGroupInput(session, "zcl_cell_types_4", label = NULL,
                               choices = zcl_cell_types$cell_type[34:41],selected = zcl_cell_types$cell_type[34:41])
      
    }
  })  
    

  # Zebrahub lineages
  
  observe({
    if(is.null(input$deselectall_zhub)) return(NULL)
    
    if (input$deselectall_zhub > 0) {
      
      updateCheckboxGroupInput(session, "zhub_cell_types_1", label = NULL,
                               choices = zhub_tissues$tissue[1:5],selected = NULL)
      updateCheckboxGroupInput(session, "zhub_cell_types_2", label = NULL,
                               choices = zhub_tissues$tissue[6:10],selected = NULL)      
    }
  })
  
  
  observe({
    if(is.null(input$selectall_zhub)) return(NULL)
    
    if (input$selectall_zhub > 0) {
      
      updateCheckboxGroupInput(session, "zhub_cell_types_1", label = NULL,
                               choices = zhub_tissues$tissue[1:5],selected = zhub_tissues$tissue[1:5])
      updateCheckboxGroupInput(session, "zhub_cell_types_2", label = NULL,
                               choices = zhub_tissues$tissue[6:10],selected = zhub_tissues$tissue[6:10])      
    }
  })
  
  
  #------------------------------- Stages - select all / deselect all --------------
  
  # Daniocell stages
  
  observe({
    if(is.null(input$deselect_dcell_stages)) return(NULL)
    
    if (input$deselect_dcell_stages > 0) {
      
      updateCheckboxGroupInput(session, "stages_1", label = NULL,
                               choices = Daniocell_stages$stage[1:5],selected = NULL)
      
      updateCheckboxGroupInput(session, "stages_2", label = NULL,
                               choices = Daniocell_stages$stage[6:10],selected = NULL)
      
      updateCheckboxGroupInput(session, "stages_3", label = NULL,
                               choices = Daniocell_stages$stage[11:14],selected = NULL)
      

    }
  })
  
  
  observe({
    if(is.null(input$select_dcell_stages)) return(NULL)
    
    if (input$select_dcell_stages > 0) {
      
      updateCheckboxGroupInput(session, "stages_1", label = NULL,
                               choices = Daniocell_stages$stage[1:5],selected = Daniocell_stages$stage[1:5])
      
      updateCheckboxGroupInput(session, "stages_2", label = NULL,
                               choices = Daniocell_stages$stage[6:10],selected = Daniocell_stages$stage[6:10])
      
      updateCheckboxGroupInput(session, "stages_3", label = NULL,
                               choices = Daniocell_stages$stage[11:14],selected = Daniocell_stages$stage[11:14])
      
      
    }
  })
  
  # Zebrafish Cell Landscape - no need to filter because there are only 5 of them
  
  
  # Zebrahub stages
  
  observe({
    if(is.null(input$deselect_zhub_stages)) return(NULL)
    
    if (input$deselect_zhub_stages > 0) {
      
      updateCheckboxGroupInput(session, "zhub_stages_1", label = NULL,
                               choices = zhub_stages$stage[1:5],selected = NULL)
      
      updateCheckboxGroupInput(session, "zhub_stages_2", label = NULL,
                               choices = zhub_stages$stage[6:10],selected = NULL)
      
      
    }
  })
  
  
  observe({
    if(is.null(input$select_zhub_stages)) return(NULL)
    
    if (input$select_zhub_stages > 0) {
      
      updateCheckboxGroupInput(session, "zhub_stages_1", label = NULL,
                               choices = zhub_stages$stage[1:5],selected = zhub_stages$stage[1:5])
      
      updateCheckboxGroupInput(session, "zhub_stages_2", label = NULL,
                               choices = zhub_stages$stage[6:10],selected = zhub_stages$stage[6:10])
      
      
    }
  })
  
  
  # ------------- Daniocell data ----------------------------------------------
  
  # 1. Reactive filtering of data based on user inputs
  daniocell_data <- reactive({
    
    # validate the inputs
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
    req(input$Daniocell_update_plot)
    
    # 2. Combined validation checks with user-friendly messages
    validate(
      need(input$daniocell_gene_id != "", "Please select or type a gene symbol to begin."),
      #Daniocell_genes
      need(input$daniocell_gene_id %in% Daniocell_genes$orig_gene, "Your input gene needs to be present in the dataset."),
      need(c(input$cell_types_1, input$cell_types_2, input$cell_types_3, input$cell_types_4), 
        "Please select at least one tissue." ),
      need(c(input$stages_1, input$stages_2, input$stages_3), 
        "Please select at least one developmental stage." )
    )
    
    # 3. Store the input items in variables
    
    #collect the gene ID
    gene_id = input$daniocell_gene_id
    
    # collect all tissue terms
    selected_tissues <- c(input$cell_types_1, input$cell_types_2, input$cell_types_3, input$cell_types_4)
    
    # collect stages    
    selected_stages <- c(input$stages_1, input$stages_2, input$stages_3)
    
    # 4. Load the relevant data frame based on the other options
    
    # obtain the processed gene name
    proc_gene <- unique(Daniocell_genes[Daniocell_genes$orig_gene == gene_id,]$proc_gene)
    
    # read counts option
    if(input$data_type == "count_sum"){
      
      df_filename <- paste0("./Daniocell_dataset/single_genes_dfs/daniocell_", proc_gene, "_count_sums.csv")
      
      # Validate existence
      validate(
        need(file.exists(df_filename), paste("Error: No data file found for gene", input$daniocell_gene_id))
      )
      
      # 5. Read the data frame    
      data <- read.csv(df_filename)
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "counts_sum")
      data_long <- data_long[data_long$counts_sum > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
      
            
    }else{
      
      df_filename <- paste0("./Daniocell_dataset/single_genes_cell_counts/daniocell_", proc_gene, "_cell_counts.csv")
      
      # Validate existence
      validate(
        need(file.exists(df_filename), paste("Error: No data file found for gene", input$daniocell_gene_id))
      )
      
      # 5. Read the data frame    
      data <- read.csv(df_filename)
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "cell_counts")
      data_long <- data_long[data_long$cell_counts > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
      
    }
    
    ############################# Plotting part of the function ###############
    
    # define relevant factors
    data_long$stage <- factor(data_long$stage, levels = Daniocell_stages$stage)
    data_long$tissue <- factor(data_long$tissue, levels = Daniocell_tissues$tissue)
    
    
    # 7. filter the data frame by the selected tissues and stages
    data <- filter(data_long, tissue %in% selected_tissues & stage %in% selected_stages)
    gene <- gene_id
    
    
    # Check if the filtered dataset is empty
    if (nrow(data) == 0) {
      
      # Display a friendly message instead of an error
      return(ggplot() + 
               geom_text(aes(x=0.5, y=0.5, label="No data points match the current filters."), 
                         size=6, color="#dc2626") +
               theme_void())
    } else{
      
      # specify the y-axis label name
      if(input$data_type == "count_sum"){
        y_value = "counts_sum"
        y_label = "Normalized cell expression value"
      } else{
        y_value = "cell_counts"
        y_label = "Normalized cell count"
      }
      
      # define the plot depending on the y scale that was specified in the input
      if(input$y_scale == "linear"){
        
        return(ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
                 geom_col() +
                 scale_y_continuous(n.breaks = 6) +
                 facet_wrap(~tissue, scales = "free_x") +
                 ylab(y_label)+
                 ggtitle(paste("Summarised expression plot for", gene , "in zebrafish"))+
                 theme_bw() + paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
                 theme(plot.title = element_text(size = 14),
                       strip.text.x = element_text(size = 14, margin = margin(0.12,0,0.12,0, "cm")),
                       axis.text.y = element_text(size = 12),
                       axis.title.y = element_text(size = 12),
                       axis.title.x = element_text(size = 12),
                       axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 12, face="plain"),
                       legend.title = element_text( size = 12, face = "bold"),
                       legend.text = element_text( size = 12, face = "plain"),
                       legend.key.size = unit(0.5, "cm"),
                       panel.spacing = unit(0.1, "lines")) 
        )
        
      } else{
        
        return(ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
                 geom_col() +
                 scale_y_sqrt(n.breaks = 6) +
                 facet_wrap(~tissue, scales = "free_x") +
                 ylab( y_label )+
                 ggtitle(paste("Summarised expression plot for", gene , "in zebrafish"))+
                 theme_bw() + paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
                 theme(plot.title = element_text(size = 14),
                       strip.text.x = element_text(size = 14, margin = margin(0.12,0,0.12,0, "cm")),
                       axis.text.y = element_text(size = 12),
                       axis.title.y = element_text(size = 12),
                       axis.title.x = element_text(size = 12),
                       axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 12, face="plain"),
                       legend.title = element_text( size = 12, face = "bold"),
                       legend.text = element_text( size = 12, face = "plain"),
                       legend.key.size = unit(0.5, "cm"),
                       panel.spacing = unit(0.1, "lines"))
        )
        
      } # else for specific scale plot
    } # else for the plot code
    
  })
  
  
  # 2. Render the plot
  output$daniocell_expr_plot <- renderPlot({
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
#    req(input$Daniocell_update_plot)
    
    daniocell_data()
    
  }, width = 1250, height = 1250, res = 96)

  
  ################### detecting which of the two ZCL plot buttons was clicked ##
  
  rv <- reactiveValues(lastBtn = character())
  
  observeEvent(input$lineages_plot_zcl, {
    if (input$lineages_plot_zcl > 0 ) {
      rv$lastBtn = "lineages_plot_zcl"
    }
  })
  
  
  observeEvent(input$celltypes_plot_zcl, {
    if (input$celltypes_plot_zcl > 0 ) {
      rv$lastBtn = "celltypes_plot_zcl"
    }
  })

  # output$lastButtonClicked <- renderText({
  #   paste(rv$lastBtn)
  # })
  
  
  #############################################################################
  
  
  
  # ------------- ZCL lineage data ---------------------------------------------
  
  # 1. Reactive filtering of data based on user inputs
  ZCL_lineage_data <- reactive({
    
    # validate the inputs
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
    req(input$lineages_plot_zcl)
    
    # 2. Combined validation checks with user-friendly messages
    validate(
      need(input$zcl_gene_id != "", "Please select or type a gene symbol to begin."),
      
      #zcl_genes
      need(input$zcl_gene_id %in% zcl_genes$orig_gene, "Your input gene needs to be present in the dataset."),
      
      need(c(input$zcl_lineages_1, input$zcl_lineages_2), "Please select at least one lineage." ),
      
      need(c(input$lineage_zcl_stages), "Please select at least one developmental stage." )
    )
    
    # 3. Store the input items in variables 
    #collect the gene ID
    gene_id = input$zcl_gene_id
    
    # if lineages
    selected_tissues <- c(input$zcl_lineages_1, input$zcl_lineages_2)
    
    # collect stages    
    selected_stages <- c(input$lineage_zcl_stages)
    
    # 4. Load the relevant data frame based on the other options
    
    # obtain the processed gene name
    proc_gene <- unique(zcl_genes[zcl_genes$orig_gene == gene_id,]$proc_gene)
    
    
    # read counts option
    if(input$lineage_data_type == "count_sum"){
      
      df_filename <- paste0("./ZCL_dataset_lineages/single_genes_dfs/ZCL_", proc_gene, "_count_sums.csv")
      
      # Validate existence
      validate(
        need(file.exists(df_filename), paste("Error: No data file found for gene", input$zcl_gene_id))
      )
      
      # 5. Read the data frame    
      data <- read.csv(df_filename)
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "counts_sum")
      data_long <- data_long[data_long$counts_sum > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
                                       
      
    }else{
      
      df_filename <- paste0("./ZCL_dataset_lineages/single_genes_cell_counts/ZCL_", proc_gene, "_cell_counts.csv")

      # Validate existence
      validate(
        need(file.exists(df_filename), paste("Error: No data file found for gene", input$zcl_gene_id))
      )
      
      # 5. Read the data frame    
      data <- read.csv(df_filename)
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "cell_counts")
      data_long <- data_long[data_long$cell_counts > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
      
    }
    
    ############################# Plotting part of the function ###############
    
    # define relevant factors
    data_long$stage <- factor(data_long$stage, levels = zcl_stages$stage)
    data_long$tissue <- factor(data_long$tissue, levels = zcl_tissues$tissue)
    
    
        
    # 7. filter the data frame by the selected tissues and stages
    data <- filter(data_long, tissue %in% selected_tissues & stage %in% selected_stages)
    gene <- gene_id
    
    
    # Check if the filtered dataset is empty
    if (nrow(data) == 0) {
      
      # Display a friendly message instead of an error
      return(ggplot() + 
               geom_text(aes(x=0.5, y=0.5, label="No data points match the current filters."), 
                         size=6, color="#dc2626") +
               theme_void())
    } else{
      
      # specify the y-axis label name
      if(input$lineage_data_type == "count_sum"){
        y_value = "counts_sum"
        y_label = "Normalized cell expression value"
      } else{
        y_value = "cell_counts"
        y_label = "Normalized cell count"
      }
      
      # define the plot depending on the y scale that was specified in the input
      if(input$lineage_y_scale == "linear"){
        
        return(ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
                 geom_col() +
                 scale_y_continuous(n.breaks = 6) +
                 facet_wrap(~tissue, scales = "free_x") +
                 ylab(y_label)+
                 ggtitle(paste("Summarised expression plot for", gene , "in zebrafish"))+
                 theme_bw() + paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
                 theme(plot.title = element_text(size = 14),
                       strip.text.x = element_text(size = 14, margin = margin(0.12,0,0.12,0, "cm")),
                       axis.text.y = element_text(size = 12),
                       axis.title.y = element_text(size = 12),
                       axis.title.x = element_text(size = 12),
                       axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 12, face="plain"),
                       legend.title = element_text( size = 12, face = "bold"),
                       legend.text = element_text( size = 12, face = "plain"),
                       legend.key.size = unit(0.5, "cm"),
                       panel.spacing = unit(0.1, "lines")) 
        )
        
      } else{
        
        return(ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
                 geom_col() +
                 scale_y_sqrt(n.breaks = 6) +
                 facet_wrap(~tissue, scales = "free_x") +
                 ylab( y_label )+
                 ggtitle(paste("Summarised expression plot for", gene , "in zebrafish"))+
                 theme_bw() + paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
                 theme(plot.title = element_text(size = 14),
                       strip.text.x = element_text(size = 14, margin = margin(0.12,0,0.12,0, "cm")),
                       axis.text.y = element_text(size = 12),
                       axis.title.y = element_text(size = 12),
                       axis.title.x = element_text(size = 12),
                       axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 12, face="plain"),
                       legend.title = element_text( size = 12, face = "bold"),
                       legend.text = element_text( size = 12, face = "plain"),
                       legend.key.size = unit(0.5, "cm"),
                       panel.spacing = unit(0.1, "lines"))
        )
        
      } # else for specific scale plot
    } # else for the plot code
    
  })
  
  
  # ------------- ZCL cell types data ---------------------------------------------
  
  # 1. Reactive filtering of data based on user inputs
  ZCL_celltypes_data <- reactive({
    
    # validate the inputs
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
    req(input$celltypes_plot_zcl)
    
    # 2. Combined validation checks with user-friendly messages
    validate(
      need(input$zcl_gene_id != "", "Please select or type a gene symbol to begin."),
      
      #zcl_genes
      need(input$zcl_gene_id %in% zcl_genes$orig_gene, "Your input gene needs to be present in the dataset."),
      
      need(c(input$zcl_cell_types_1, input$zcl_cell_types_2, input$zcl_cell_types_3, input$zcl_cell_types_4), "Please select at least one cell type." ),
      
      need(c(input$celltypes_zcl_stages), "Please select at least one developmental stage." )
    )
    
    # 3. Store the input items in variables 
    #collect the gene ID
    gene_id = input$zcl_gene_id
    
    # if lineages
    selected_tissues <- c(input$zcl_cell_types_1, input$zcl_cell_types_2, input$zcl_cell_types_3, input$zcl_cell_types_4)
    
    # collect stages    
    selected_stages <- c(input$celltypes_zcl_stages)
    
    # 4. Load the relevant data frame based on the other options
    
    # obtain the processed gene name
    proc_gene <- unique(zcl_genes[zcl_genes$orig_gene == gene_id,]$proc_gene)
    
    
    # read counts option
    if(input$celltypes_data_type == "count_sum"){
      
      df_filename <- paste0("./ZCL_dataset_celltypes/single_genes_dfs/ZCL_", proc_gene, "_count_sums.csv")
      
      # Validate existence
      validate(
        need(file.exists(df_filename), paste("Error: No data file found for gene", input$zcl_gene_id))
      )
      
      # 5. Read the data frame    
      data <- read.csv(df_filename)
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "counts_sum")
      data_long <- data_long[data_long$counts_sum > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
      
      
    }else{
      
      df_filename <- paste0("./ZCL_dataset_celltypes/single_genes_cell_counts/ZCL_", proc_gene, "_cell_counts.csv")
      
      # Validate existence
      validate(
        need(file.exists(df_filename), paste("Error: No data file found for gene", input$zcl_gene_id))
      )
      
      # 5. Read the data frame    
      data <- read.csv(df_filename)
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "cell_counts")
      data_long <- data_long[data_long$cell_counts > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
      
    }
    
    ############################# Plotting part of the function ###############
    # define relevant factors
    data_long$stage <- factor(data_long$stage, levels = zcl_stages$stage)
    data_long$tissue <- factor(data_long$tissue, levels = zcl_cell_types$cell_type)
    
    # 7. filter the data frame by the selected tissues and stages
    data <- filter(data_long, tissue %in% selected_tissues & stage %in% selected_stages)
    gene <- gene_id
    
    
    # Check if the filtered dataset is empty
    if (nrow(data) == 0) {
      
      # Display a friendly message instead of an error
      return(ggplot() + 
               geom_text(aes(x=0.5, y=0.5, label="No data points match the current filters."), 
                         size=6, color="#dc2626") +
               theme_void())
    } else{
      
      # specify the y-axis label name
      if(input$celltypes_data_type == "count_sum"){
        y_value = "counts_sum"
        y_label = "Normalized cell expression value"
      } else{
        y_value = "cell_counts"
        y_label = "Normalized cell count"
      }
      
      # define the plot depending on the y scale that was specified in the input
      if(input$celltypes_y_scale == "linear"){
        
        return(ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
                 geom_col() +
                 scale_y_continuous(n.breaks = 6) +
                 facet_wrap(~tissue, scales = "free_x") +
                 ylab(y_label)+
                 ggtitle(paste("Summarised expression plot for", gene , "in zebrafish"))+
                 theme_bw() + paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
                 theme(plot.title = element_text(size = 14),
                       strip.text.x = element_text(size = 14, margin = margin(0.12,0,0.12,0, "cm")),
                       axis.text.y = element_text(size = 12),
                       axis.title.y = element_text(size = 12),
                       axis.title.x = element_text(size = 12),
                       axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 12, face="plain"),
                       legend.title = element_text( size = 12, face = "bold"),
                       legend.text = element_text( size = 12, face = "plain"),
                       legend.key.size = unit(0.5, "cm"),
                       panel.spacing = unit(0.1, "lines")) 
        )
        
      } else{
        
        return(ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
                 geom_col() +
                 scale_y_sqrt(n.breaks = 6) +
                 facet_wrap(~tissue, scales = "free_x") +
                 ylab( y_label )+
                 ggtitle(paste("Summarised expression plot for", gene , "in zebrafish"))+
                 theme_bw() + paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
                 theme(plot.title = element_text(size = 14),
                       strip.text.x = element_text(size = 14, margin = margin(0.12,0,0.12,0, "cm")),
                       axis.text.y = element_text(size = 12),
                       axis.title.y = element_text(size = 12),
                       axis.title.x = element_text(size = 12),
                       axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 12, face="plain"),
                       legend.title = element_text( size = 12, face = "bold"),
                       legend.text = element_text( size = 12, face = "plain"),
                       legend.key.size = unit(0.5, "cm"),
                       panel.spacing = unit(0.1, "lines"))
        )
        
      } # else for specific scale plot
    } # else for the plot code
    
  })
  
  
  # # lineages plot  
  # output$zcl_expr_plot_lineage <- renderPlot({
  #   
  #   ZCL_lineage_data()
  #   
  # }, width = 1250, height = 1250, res = 96)
  # 
  # 
  # # Cell types plot
  # output$zcl_expr_plot_celltypes <- renderPlot({
  #   
  #   ZCL_celltypes_data()
  #   
  # }, width = 1250, height = 1250, res = 96)
  
  
  # 
  # 2. Render the plot
  output$zcl_expr_plot <- renderPlot({
    req(input$lineages_plot_zcl | input$celltypes_plot_zcl )

    # test
    if( rv$lastBtn == "lineages_plot_zcl" ){
      return(ZCL_lineage_data())
    }

    if( rv$lastBtn == "celltypes_plot_zcl" ){
       return(ZCL_celltypes_data())
    }

  }, width = 1000, height = 1000, res = 96)
    
  
  # ------------- Zebrahub data ----------------------------------------------
  
  # 1. Reactive filtering of data based on user inputs
  zhub_data <- reactive({
    
    # validate the inputs
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
    req(input$update_plot_zhub)
    
    # 2. Combined validation checks with user-friendly messages
    validate(
      need(input$zhub_gene_id != "", "Please select or type a gene symbol to begin."),
      #zhub_genes
      need(input$zhub_gene_id %in% zhub_genes$orig_gene, "Your input gene needs to be present in the dataset."),
      
      need(c(input$zhub_cell_types_1, input$zhub_cell_types_2), "Please select at least one cell type." ),
      need(c(input$zhub_stages_1, input$zhub_stages_2), "Please select at least one developmental stage." )
    )

    # 3. Store the input items in variables
    
    #collect the gene ID
    gene_id = trimws(input$zhub_gene_id)
    
    # collect all tissue terms
    selected_tissues <- c(input$zhub_cell_types_1, input$zhub_cell_types_2)
    
    # collect stages    
    selected_stages <- c(input$zhub_stages_1, input$zhub_stages_2)
    
    # 4. Load the relevant data frame based on the other options
    
    # obtain the processed gene name
    proc_gene <- unique(zhub_genes[zhub_genes$orig_gene == gene_id,]$proc_gene)
    
    # read counts option
    if(input$zhub_data_type == "count_sum"){
      
      df_filename <- paste0("./Zebrahub_dataset/single_genes_dfs/zebrahub_", proc_gene, "_count_sums.csv")
      
      # Validate existence
      validate(
        need(file.exists(df_filename), paste("Error: No data file found for gene", input$zhub_gene_id))
      )
      
      # 5. Read the data frame    
      data <- read.csv(df_filename)
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "counts_sum")
      data_long <- data_long[data_long$counts_sum > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
      
      
    }else{
      
      df_filename <- paste0("./Zebrahub_dataset/single_genes_cell_counts/zebrahub_", proc_gene, "_cell_counts.csv")
      
      # Validate existence
      validate(
        need(file.exists(df_filename), paste("Error: No data file found for gene", input$zhub_gene_id))
      )
      
      # 5. Read the data frame    
      data <- read.csv(df_filename)
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "cell_counts")
      data_long <- data_long[data_long$cell_counts > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
      
    }
    
    ############################# Plotting part of the function ###############

    # define relevant factors
    data_long$stage <- factor(data_long$stage, levels = zhub_stages$stage)
    data_long$tissue <- factor(data_long$tissue, levels = zhub_tissues$tissue)
    
    # 7. filter the data frame by the selected tissues and stages
    data <- filter(data_long, tissue %in% selected_tissues & stage %in% selected_stages)
    gene <- gene_id
    
    
    # Check if the filtered dataset is empty
    if (nrow(data) == 0) {
      
      # Display a friendly message instead of an error
      return(ggplot() + 
               geom_text(aes(x=0.5, y=0.5, label="No data points match the current filters."), 
                         size=6, color="#dc2626") +
               theme_void())
    } else{
      
      # specify the y-axis label name
      if(input$zhub_data_type == "count_sum"){
        y_value = "counts_sum"
        y_label = "Normalized cell expression value"
      } else{
        y_value = "cell_counts"
        y_label = "Normalized cell count"
      }
      
      # define the plot depending on the y scale that was specified in the input
      if(input$zhub_y_scale == "linear"){
        
        return(ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
                 geom_col() +
                 scale_y_continuous(n.breaks = 6) +
                 facet_wrap(~tissue, scales = "free_x") +
                 ylab(y_label)+
                 ggtitle(paste("Summarised expression plot for", gene , "in zebrafish"))+
                 theme_bw() + paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
                 theme(plot.title = element_text(size = 14),
                       strip.text.x = element_text(size = 14, margin = margin(0.12,0,0.12,0, "cm")),
                       axis.text.y = element_text(size = 12),
                       axis.title.y = element_text(size = 12),
                       axis.title.x = element_text(size = 12),
                       axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 12, face="plain"),
                       legend.title = element_text( size = 12, face = "bold"),
                       legend.text = element_text( size = 12, face = "plain"),
                       legend.key.size = unit(0.5, "cm"),
                       panel.spacing = unit(0.1, "lines")) 
        )
        
      } else{
        
        return(ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
                 geom_col() +
                 scale_y_sqrt(n.breaks = 6) +
                 facet_wrap(~tissue, scales = "free_x") +
                 ylab( y_label )+
                 ggtitle(paste("Summarised expression plot for", gene , "in zebrafish"))+
                 theme_bw() + paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
                 theme(plot.title = element_text(size = 14),
                       strip.text.x = element_text(size = 14, margin = margin(0.12,0,0.12,0, "cm")),
                       axis.text.y = element_text(size = 12),
                       axis.title.y = element_text(size = 12),
                       axis.title.x = element_text(size = 12),
                       axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 12, face="plain"),
                       legend.title = element_text( size = 12, face = "bold"),
                       legend.text = element_text( size = 12, face = "plain"),
                       legend.key.size = unit(0.5, "cm"),
                       panel.spacing = unit(0.1, "lines"))
        )
        
      } # else for specific scale plot
    } # else for the plot code
    
  })
  
  
  # 2. Render the plot
  output$zebrahub_expr_plot <- renderPlot({
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
    #    req(input$update_plot_zhub)
    
    zhub_data()
    
  }, width = 1250, height = 1250, res = 96)
  
  
  
}

# Run the complete application
shinyApp(ui = ui, server = server)