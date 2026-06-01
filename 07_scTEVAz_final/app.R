library(shiny)
library(readr)
library(shinyjs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(paletteer)
library(shinycssloaders)
library(bsicons)
library(stringr)

# packages for interactive plots 
library(ggiraph)

################################## CONSTANT VARIABLES #####################

FACET_SIZE = 2


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
        h1("Welcome to the scTEVAz"),
        h3("Choose which dataset you want to use:"),
        
        br(),
        actionButton("goToTool1", h4("Daniocell"), class = "btn-lg btn-primary"),
        actionButton("goToTool2", h4("Zebrafish Cell Landscape"), class = "btn-lg btn-success"),
        actionButton("goToTool3", h4("Zebrahub"), class = "btn-lg btn-info")
      ),
      
      
      div(
        style = "text-align: center;",
        img(src = "logo_video.gif", height = 500, width = 500)
      ),
      br(),
      
      HTML("<footer style = 'position:fixed; bottom:0; width:100%; height:50px; color: white; padding: 10px; background-color: #152437; z-index: 1000;'>
                      <p style='font-size: 18px'> scTEVAz was developed by Sergey Prykhozhij</strong> while working at the CHEO Research Institute in Ottawa, Canada. If you encounter a problem, please send an email to <strong>Sergey Prykhozhij</strong> at <strong><font style= 'color: lightblue'>s.prykhozhij@gmail.com</font></strong>.</p>
          </footer>")
       ),
    
    # Daniocell panel
    tabPanel(
      
      title = "Daniocell", 
      
      tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"),
                tags$style(HTML(".shiny-output-error-validation {color: green;font-size: 20px;}")),
                tags$style(HTML(".checkbox { font-size: 18px; }")),
                tags$style(HTML(".radio { font-size: 18px; }")),
                tags$style(HTML(".radio-inline { font-size: 18px; }")),
                tags$style(HTML("/* Professional Header Styling */
                  .plot-header {
                    background-color: #ffffff;
                    border-bottom: 2px solid #0c4a6e;
                    padding: 20px;
                    margin-bottom: 25px;
                    border-radius: 8px 8px 0 0;
                    box-shadow: 0 4px 6px rgba(0,0,0,0.05);
                  }
                  .plot-title {
                    color: #0c4a6e;
                    font-weight: 700;
                    margin: 0;
                    letter-spacing: 0.5px;
                  }
                  /* Button refinement */
                  .btn-download {
                    margin-bottom: 15px;
                    border-radius: 20px;
                    font-weight: 500;
                    transition: all 0.3s ease;
                  }
                  .btn-download:hover {
                    transform: translateY(-2px);
                    box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                  }
                  /* Container for the plot */
                  .plot-container {
                    background: white;
                    padding: 20px;
                    border-radius: 12px;
                    box-shadow: 0 10px 15px -3px rgba(0,0,0,0.1);
                  }"))
                ),
      
      tags$style("body {
                          -moz-transform: scale(0.95 0.95); /* Moz-browsers */
                          zoom: 0.8; /* Other non-webkit browsers */
                          zoom: 80%; /* Webkit browsers */
                          }"),
      
      
      # Sidebar layout for inputs and output
      sidebarLayout(
        # Sidebar Panel - Inputs
        sidebarPanel(
          width = 4,
          class = "sidebar",
          h3(strong("Input Parameters"), style = "color: #1e3a8a;"),
          
          
          fluidRow(column(width = 5,           
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
                    )),
      
          
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
      
          fluidRow(
            column(width = 6,       
                   radioButtons("data_type", h4(strong("Choose the type of data:")),
                                                 c("Sum of normalized read counts" = "count_sum",
                                                   "Normalized cell counts" = "cell_count"),
                                                 selected = "count_sum")),
            
            column(width = 5,      
                   radioButtons("y_scale", h4(strong("Choose the y-axis scale:")),
                                                c("Linear" = "linear",
                                                  "Square root" = "sqrt"),
                                                selected = "linear",
                                                inline = TRUE))
            
          ),
          
          fluidRow(
            
            column(width = 6,      
                   radioButtons("term_sort", h4(strong("Choose how to sort the data:")),
                                c("Sum of values" = "sum",
                                  "Alphabetical" = "alpha"),
                                selected = "sum",
                                inline = TRUE)),
            column(width = 5,
                   # Action Button
                   actionButton("Daniocell_update_plot", "Generate Plot", class = "btn-lg btn-primary text-white hover:bg-blue-700")
                   
            )
          ),
          
          
      ), # End of sidebarPanel
        
        # Main Panel - Output
        mainPanel(
          style = "margin-top: -20px;",
          class = "main-panel",

          # panels for instructions 
          tabsetPanel(id="maintabset_daniocell",
                      
                      tabPanel(p(class = "panel-title",style="width: 100%, font-size: 14px; color: blue", "Instructions"), value = "Instructions", includeHTML("./www/instructions.html")),
                      
                      tabPanel(p(class = "panel-title",style="width: 100%, font-size: 14px; color: blue", "Results"), value = "Results",
                        
                        # Title Section
                        div(class = "plot-header",
                            h3(class = "plot-title", "Daniocell Aggregated Expression Values")
                        ),
                        
                        # Action Row (Buttons)
                        fluidRow(
                          column(12, style = "display: flex; gap: 10px; margin-bottom: 20px;",
                                 withSpinner(
                                   downloadButton("download_png", "Save Plot (PNG)", icon = icon("image"), class = "btn-download"),
                                   type = 4),
                                 withSpinner(
                                   downloadButton("download_data", "Export Data (CSV)", icon = icon("table"), class = "btn-download"),
                                   type = 4 ),
                                 withSpinner(
                                   downloadButton("download_code", "Save the plot code", icon = icon("code"), class = "btn-download"),
                                   type = 4 )
                          )
                        ),
                        
                        # Plot Section
                  
                        withSpinner(girafeOutput("daniocell_expr_plot"), type = 4)
                      
                      )
                                              
                      
          ) # end of maintabset
                    
          
        ) # End of mainPanel
      ) # End of sidebarLayout
      
    ),
    
    # ZCL Panel
    tabPanel(
      
      title = "Zebrafish Cell Landscape",
      
      tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"),
                tags$style(HTML(".shiny-output-error-validation {color: green;font-size: 20px;}")),
                tags$style(HTML(".checkbox { font-size: 18px; }")),
                tags$style(HTML(".radio { font-size: 18px; }")),
                tags$style(HTML(".radio-inline { font-size: 18px; }")),
                tags$style(HTML("/* Professional Header Styling */
                  .plot-header {
                    background-color: #ffffff;
                    border-bottom: 2px solid #0c4a6e;
                    padding: 20px;
                    margin-bottom: 25px;
                    border-radius: 8px 8px 0 0;
                    box-shadow: 0 4px 6px rgba(0,0,0,0.05);
                  }
                  .plot-title {
                    color: #0c4a6e;
                    font-weight: 700;
                    margin: 0;
                    letter-spacing: 0.5px;
                  }
                  /* Button refinement */
                  .btn-download {
                    margin-bottom: 15px;
                    border-radius: 20px;
                    font-weight: 500;
                    transition: all 0.3s ease;
                  }
                  .btn-download:hover {
                    transform: translateY(-2px);
                    box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                  }
                  /* Container for the plot */
                  .plot-container {
                    background: white;
                    padding: 20px;
                    border-radius: 12px;
                    box-shadow: 0 10px 15px -3px rgba(0,0,0,0.1);
                  }"))
      ),
      
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
                          zoom: 0.8; /* Other non-webkit browsers */
                          zoom: 80%; /* Webkit browsers */
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
              

              # Input 3: Stage Filtering
              fluidRow(
                column(width = 5,
                       
                       checkboxGroupInput( 
                         inputId = "lineage_zcl_stages",
                         label =  h4(strong("Filter by stage:")),
                         choices = zcl_stages$stage[1:5],
                         selected = zcl_stages$stage[1:5],
                       )
                ),
                
                column(width = 7,       
                       radioButtons("lineage_data_type", h4(strong("Choose the type of data:")),
                                    c("Sum of normalized read counts" = "count_sum",
                                      "Normalized cell counts" = "cell_count"),
                                    selected = "count_sum"),
                       radioButtons("lineage_y_scale", h4(strong("Choose the y-axis scale:")),
                                    c("Linear" = "linear",
                                      "Square root" = "sqrt"),
                                    selected = "linear",
                                    inline = TRUE)
                       
                       )
              ),
              
              fluidRow(
                
                column(width = 5,      
                       radioButtons("zcl_lin_term_sort", h4(strong("Choose how to sort the data:")),
                                    c("Sum of values" = "sum",
                                      "Alphabetical" = "alpha"),
                                    selected = "sum",
                                    inline = TRUE)),
                column(width = 7,
                       # Action Button
                       actionButton("lineages_plot_zcl", "Generate Lineages Plot", class = "btn-lg btn-primary text-white hover:bg-blue-700")
                       
                )
              )             
              
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
                

                fluidRow(
                  column(width = 5,
                         
                         checkboxGroupInput(
                           inputId = "celltypes_zcl_stages",
                           label = tags$label(h4(strong("Filter by stage:"))),
                           choices = zcl_stages$stage[1:5],
                           selected = zcl_stages$stage[1:5],
                         )
                  ),
                  
                  column(width = 7,       
                         radioButtons("celltypes_data_type", h4(strong("Choose the type of data:")),
                                      c("Sum of normalized read counts" = "count_sum",
                                        "Normalized cell counts" = "cell_count"),
                                      selected = "count_sum"),
                         
                         radioButtons("celltypes_y_scale", h4(strong("Choose the y-axis scale:")),
                                      c("Linear" = "linear",
                                        "Square root" = "sqrt"),
                                      selected = "linear",
                                      inline = TRUE)                         
                         
                         )
                  
                ),
                
                
                fluidRow(
                  
                  column(width = 5,      
                         radioButtons("zcl_celltype_term_sort", h4(strong("Choose how to sort the data:")),
                                      c("Sum of values" = "sum",
                                        "Alphabetical" = "alpha"),
                                      selected = "sum",
                                      inline = TRUE)),
                  column(width = 7,
                         # Action Button
                         actionButton("celltypes_plot_zcl", "Generate Cell types Plot", class = "btn-lg btn-primary text-white hover:bg-blue-700")
                         
                  )
                )                
                
            ) # end of celltypes_wrapper div
          ), 
          
        ), # End of sidebarPanel
        
        # Main Panel - Output
        mainPanel(
          style = "margin-top: -20px;",
          class = "main-panel",
    
          # panels for instructions and results
          tabsetPanel(id="maintabset_zcl",
                      
                      tabPanel(p(class = "panel-title",style="width: 100%, font-size: 14px; color: blue", "Instructions"), value = "Instructions", includeHTML("./www/instructions.html")),
                      
                      tabPanel(p(class = "panel-title",style="width: 100%, font-size: 14px; color: blue", "Results"), value = "Results", 
                               
                                 
                                 # plotting the resulting graph
                                 
                               # Title Section
                               div(class = "plot-header",
                                   h3(class = "plot-title", "ZCL Aggregated expression values plot")
                               ),
                               
                               # Action Row (Buttons)
                               fluidRow(
                                 column(12, style = "display: flex; gap: 10px; margin-bottom: 20px;",
                                        withSpinner(
                                          downloadButton("zcl_download_png", "Download as PNG"),
                                          type = 4),
                                        withSpinner(
                                          downloadButton("zcl_download_data", "Export Data (CSV)", icon = icon("table"), class = "btn-download"),
                                          type = 4 ),
                                        withSpinner(
                                          downloadButton("zcl_download_code", "Save the plot code", icon = icon("code"), class = "btn-download"),
                                          type = 4 )
                                 )
                               ),
                               
                               # Plot Section
                               
                               withSpinner(girafeOutput("zcl_expr_plot"), type = 4)
                               
                               
                               
                                          
                      )
                      
          ) # end of maintabset
          
          
        ) # End of mainPanel
      ) # End of sidebarLayout
      
    ),
    
    # Zebrahub Panel
    tabPanel(
      
      title = "Zebrahub",
      
      tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"),
                tags$style(HTML(".shiny-output-error-validation {color: green;font-size: 20px;}")),
                tags$style(HTML(".checkbox { font-size: 18px; }")),
                tags$style(HTML(".radio { font-size: 18px; }")),
                tags$style(HTML(".radio-inline { font-size: 18px; }")),
                tags$style(HTML("/* Professional Header Styling */
                  .plot-header {
                    background-color: #ffffff;
                    border-bottom: 2px solid #0c4a6e;
                    padding: 20px;
                    margin-bottom: 25px;
                    border-radius: 8px 8px 0 0;
                    box-shadow: 0 4px 6px rgba(0,0,0,0.05);
                  }
                  .plot-title {
                    color: #0c4a6e;
                    font-weight: 700;
                    margin: 0;
                    letter-spacing: 0.5px;
                  }
                  /* Button refinement */
                  .btn-download {
                    margin-bottom: 15px;
                    border-radius: 20px;
                    font-weight: 500;
                    transition: all 0.3s ease;
                  }
                  .btn-download:hover {
                    transform: translateY(-2px);
                    box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                  }
                  /* Container for the plot */
                  .plot-container {
                    background: white;
                    padding: 20px;
                    border-radius: 12px;
                    box-shadow: 0 10px 15px -3px rgba(0,0,0,0.1);
                  }"))
      ),
      
      
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
                          zoom: 0.8; /* Other non-webkit browsers */
                          zoom: 80%; /* Webkit browsers */
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
            column(width = 3,
                   
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
          
          
          fluidRow(
            column(width = 5,       
                   radioButtons("zhub_data_type", h4(strong("Choose the type of data:")),
                                c("Sum of normalized read counts" = "count_sum",
                                  "Normalized cell counts" = "cell_count"),
                                selected = "count_sum")),
            
            column(width = 5,      
                   radioButtons("zhub_y_scale", h4(strong("Choose the y-axis scale:")),
                                c("Linear" = "linear",
                                  "Square root" = "sqrt"),
                                selected = "linear",
                                inline = TRUE))
            
          ),
          
          
          fluidRow(
            
            column(width = 5,      
                   radioButtons("zhub_term_sort", h4(strong("Choose how to sort the data:")),
                                c("Sum of values" = "sum",
                                  "Alphabetical" = "alpha"),
                                selected = "sum",
                                inline = TRUE)),
            column(width = 6,
                   # Action Button
                   actionButton("update_plot_zhub", "Generate Plot", class = "btn-lg btn-primary text-white hover:bg-blue-700")
                   
            )
          )          
        ), # End of sidebarPanel
        
        # Main Panel - Output
        mainPanel(
          style = "margin-top: -20px;",
          class = "main-panel",
          
          # panels for instructions and results
          tabsetPanel(id="maintabset_zhub",
                      
                      tabPanel(p(class = "panel-title",style="width: 100%, font-size: 14px; color: blue", "Instructions"), value = "Instructions", includeHTML("./www/instructions.html")),
                      
                      tabPanel(p(class = "panel-title",style="width: 100%, font-size: 14px; color: blue", "Results"), value = "Results", 
                               

                               # Title Section
                               div(class = "plot-header",
                                   h3(class = "plot-title", "Zebrahub Aggregated Expression Values plot")
                               ),
                               
                               # Action Row (Buttons)
                               fluidRow(
                                 column(12, style = "display: flex; gap: 10px; margin-bottom: 20px;",
                                        withSpinner(
                                          downloadButton("zhub_download_png", "Save Plot (PNG)", icon = icon("image"), class = "btn-download"),
                                          type = 4),
                                        withSpinner(
                                          downloadButton("zhub_download_data", "Export Data (CSV)", icon = icon("table"), class = "btn-download"),
                                          type = 4 ),
                                        withSpinner(
                                          downloadButton("zhub_download_code", "Save the plot code", icon = icon("code"), class = "btn-download"),
                                          type = 4 )
                                 )
                               ),
                               
                               # Plot Section
                               
                               withSpinner(girafeOutput("zebrahub_expr_plot"), type = 4)
                               
                                  
                      )
                      
          ) # end of maintabset
          
          
        
          ) # End of mainPanel
      ) # End of sidebarLayout
      
    ),
  )
) 


server <- function(input, output, session) {
  
################################## UI inputs ######################################  
  
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

####################### switching between tabs ################################
  
  # this is the code for updating the relevant tab panel
  observeEvent(input$Daniocell_update_plot, {
    if( (input$Daniocell_update_plot >= 1) ){
      
      updateTabsetPanel(session, "maintabset_daniocell", selected = "Results")    
    }
    
  })
  
  
  # this is the code for updating the relevant tab panel
  observeEvent(input$lineages_plot_zcl, {
    if( (input$lineages_plot_zcl >= 1) ){
      
      updateTabsetPanel(session, "maintabset_zcl", selected = "Results")    
    }
    
  })
  
  # this is the code for updating the relevant tab panel
  observeEvent(input$celltypes_plot_zcl, {
    if( (input$celltypes_plot_zcl >= 1) ){
      
      updateTabsetPanel(session, "maintabset_zcl", selected = "Results")    
    }
    
  })
  
  
  # this is the code for updating the relevant tab panel
  observeEvent(input$update_plot_zhub, {
    if( (input$update_plot_zhub >= 1) ){
      
      updateTabsetPanel(session, "maintabset_zhub", selected = "Results")    
    }
    
  })

########################## Navigation Logic ###################################
  
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
  daniocell_data <- eventReactive(input$Daniocell_update_plot, {
    
    # validate the inputs
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
    #req(input$Daniocell_update_plot)
    
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
    
    # convert to the version with "_"
    selected_tissues <- str_replace_all(selected_tissues, " ","_")
    
    # collect stages    
    selected_stages <- c(input$stages_1, input$stages_2, input$stages_3)
    
    # 4. Load the relevant data frame based on the other options
    
    # obtain the processed gene name
    proc_gene <- unique(Daniocell_genes[Daniocell_genes$orig_gene == gene_id,]$proc_gene)
    
    # read counts option
    if(input$data_type == "count_sum"){
      
      ############# data reading ##############################################
      
      # name of the .rds file that contains the full list of data frames
      df_filename <- "./datasets/Daniocell_count_sums.rds"
      
      # read the rds file 
      data <- readRDS(df_filename)
      
      # generate a name key to obtain a data frame
      df_key <- paste0("daniocell_", proc_gene, "_count_sums")
      
      # Validate existence of the key inside the list
      validate(
        need(df_key %in% names(data) , paste("Error: No data file found for gene", input$daniocell_gene_id))
      )
      
      # Read the data frame    
      data <- data[[df_key]]
      ##########################################################################
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "counts_sum")
      data_long <- data_long[data_long$counts_sum > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
      

    }else{
      
      ############# data reading ##############################################
      
      # name of the .rds file that contains the full list of data frames
      df_filename <- "./datasets/Daniocell_cell_counts.rds"
      
      # read the rds file 
      data <- readRDS(df_filename)
      
      # generate a name key to obtain a data frame
      df_key <- paste0("daniocell_", proc_gene, "_cell_counts")
      
      # Validate existence of the key inside the list
      validate(
        need(df_key %in% names(data) , paste("Error: No data file found for gene", input$daniocell_gene_id))
      )
      
      # Read the data frame    
      data <- data[[df_key]]
      ##########################################################################
      
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
    
    # convert tissue and selected_tissues to the version with spaces
    data_long$tissue <- str_replace_all(data_long$tissue, "_", " ")
    selected_tissues <- str_replace_all(selected_tissues, "_", " ")
    
    # redefine tissue factor
    data_long$tissue <- factor(data_long$tissue, levels = Daniocell_tissues$tissue)
    
    # 7. filter the data frame by the selected tissues and stages
    data <- filter(data_long, tissue %in% selected_tissues & stage %in% selected_stages)
    gene <- gene_id
    
    # 8. Sorting the data 
    if(input$data_type == "count_sum"){
      
      if(input$term_sort == "sum"){

        # sort the tissues by their total value
        data <- data %>%
          group_by(tissue) %>%
          mutate(Total_sum = sum(counts_sum)) %>%
          arrange(desc(Total_sum))
        
      }
    
      data$tissue <- factor(data$tissue, levels = unique(data$tissue))
      
    } else {
      
      if(input$term_sort == "sum"){
        # sort the tissues by their total value
        data <- data %>%
          group_by(tissue) %>%
          mutate(Total_sum = sum(cell_counts)) %>%
          arrange(desc(Total_sum))
      }
      
      data$tissue <- factor(data$tissue, levels = unique(data$tissue))
      
    }
    
  
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
     
      ########################### save the data ################################
      write.csv(data, "./Daniocell_dataset/temp.csv", row.names = FALSE)
      
       
      # define the plot depending on the y scale that was specified in the input
      if(input$y_scale == "linear"){
        
        p_ggiraph <- ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
          geom_col_interactive(  aes(tooltip = paste0("Stage: ", stage, 
                                 "<br>", y_label, ": ", round(.data[[y_value]], 2),
                                 "<br>Tissue: ", tissue),
                data_id = paste(stage, tissue, sep = "_") )
          ) +
          scale_y_continuous(n.breaks = 5) +
          facet_wrap(~tissue, scales = "fixed", axes = "all_x", ncol = 5) +
          ylab(y_label) +
          ggtitle(paste("Summarised expression plot for", gene, "in zebrafish")) +
          theme_bw() + 
          scale_fill_paletteer_d("colorBlindness::paletteMartin") +
          theme(plot.title = element_text(size = 8),
                strip.text.x = element_text(size = 8, margin = margin(0.1,0,0.1,0, "cm")),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 10),
                axis.title.x = element_text(size = 10),
                axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 7.5, face="plain"),
                legend.title = element_text(size = 8, face = "bold"),
                legend.text = element_text(size = 8, face = "plain"),
                legend.key.size = unit(0.15, "cm"),
                panel.spacing = unit(0.1, "lines"))   
        
        return(p_ggiraph) 
        
        
      } else{
        
        p_ggiraph <- ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
          geom_col_interactive(  aes(tooltip = paste0("Stage: ", stage, 
                                 "<br>", y_label, ": ", round(.data[[y_value]], 2),
                                 "<br>Tissue: ", tissue),
                data_id = paste(stage, tissue, sep = "_"))
          ) +
          scale_y_sqrt(n.breaks = 5) +
          facet_wrap(~tissue, scales = "fixed", axes = "all_x", ncol = 5) +
          ylab(y_label) +
          ggtitle(paste("Summarised expression plot for", gene, "in zebrafish")) +
          theme_bw() + 
          scale_fill_paletteer_d("colorBlindness::paletteMartin") +
          theme(plot.title = element_text(size = 8),
                strip.text.x = element_text(size = 8, margin = margin(0.1,0,0.1,0, "cm")),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 10),
                axis.title.x = element_text(size = 10),
                axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 7.5, face="plain"),
                legend.title = element_text(size = 8, face = "bold"),
                legend.text = element_text(size = 8, face = "plain"),
                legend.key.size = unit(0.15, "cm"),
                panel.spacing = unit(0.1, "lines"))       
        
        return(p_ggiraph)
      
      } # else for specific scale plot
    } # else for the plot code
    
  })
  
  
  # 2. Render the plot
  output$daniocell_expr_plot <- renderGirafe({
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
    req(input$Daniocell_update_plot)
    
    # run the graph output function to generate the graph and to save the data
    # that can be used
    output_obj <- daniocell_data()
    
    # 2. read the data and define how many facets there will be
    df <- read.csv("./Daniocell_dataset/temp.csv")
    n_facets = length(unique(df$tissue))
    
    if(n_facets < 5){
      calc_width = n_facets * FACET_SIZE + 0.5
      calc_height = FACET_SIZE + 0.5
    } else{
      # do the final calculation
      calc_width = 5 * FACET_SIZE
      calc_height = ceiling(n_facets/5)*FACET_SIZE + 1
    }
    
    # Convert the ggplot object to a girafe object
    girafe(ggobj = output_obj, width_svg = calc_width, height_svg = calc_height, 
           options = list(opts_sizing(rescale = TRUE), 
                          opts_zoom = opts_zoom(min = 0.25, max = 0.5) )  ) 
    
  })

  # Handle the PNG download using the base_plot
  output$download_png <- downloadHandler(
    filename = function() {
      paste0("Daniocell_plot_", input$daniocell_gene_id, '_', Sys.Date(), ".png")
    },
    content = function(file) {
      
      # obtain the data for the size parameters of the plot
      df <- read.csv("./Daniocell_dataset/temp.csv")
      n_facets = length(unique(df$tissue))
      
      if(n_facets < 5){
        calc_width = n_facets * FACET_SIZE + 0.5
        calc_height = FACET_SIZE + 0.5
      } else{
        # do the final calculation
        calc_width = 5 * FACET_SIZE
        calc_height = ceiling(n_facets/5)*FACET_SIZE + 1
      }
      
      # Use ggsave on the ggplot object      
      ggsave(file, plot = daniocell_data(), device = "png", dpi = 300, 
             width = calc_width, height = calc_height)
    }
  )
  
  # Handle data download
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("Daniocell_data_", input$daniocell_gene_id, '_', Sys.Date(), ".csv")
    },
    content = function(file) {
      
      data <- read.csv("./Daniocell_dataset/temp.csv")
      
      # Write the data to a temporary file
      write.csv(data, file, row.names = FALSE)
    }
  )
  
  # Handle data download
  output$download_code <- downloadHandler(
    
    filename = function() {
      paste0("Plot_code_Daniocell_", input$daniocell_gene_id, '_', Sys.Date(), ".R")  
    },
    content = function(file) {
      
      # generate a file name for the data in this session
      # use this file name to insert it into the code file
      data_filename <- paste0("Daniocell_data_", input$daniocell_gene_id, '_', Sys.Date(), ".csv")
      
      # load the templates as appropriate

      if(input$data_type == "count_sum"){
        
        # choose the scale to load the templates
        if(input$y_scale == "linear"){
          
          # count_sum - linear
          file_path <- "./templates/counts_sum_linear.R" 
          script_code <- read_file(file_path)
          
        } else{
          # count_sum - sqrt
          file_path <- "./templates/counts_sum_sqrt.R" 
          script_code <- read_file(file_path)          
          
        }
        

      } else {

        # choose the scale to load the templates
        if(input$y_scale == "linear"){
          # count_sum - linear
          file_path <- "./templates/cell_counts_linear.R" 
          script_code <- read_file(file_path)
          
        } else{
          # count_sum - sqrt
          file_path <- "./templates/cell_counts_sqrt.R" 
          script_code <- read_file(file_path)
                    
        }
        
      }
      
      # replace temp.csv with the current data file name
      script_code <- str_replace(script_code, 'temp.csv', data_filename)
      
      # insert gene name
      script_code <- str_replace(script_code, 'current_gene', input$daniocell_gene_id)
      
      # Write the data to a temporary file
      writeChar(script_code, file)
    }
  )
  
  
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
  ZCL_lineage_data <- eventReactive(input$lineages_plot_zcl, {
    
    # validate the inputs
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
    #req(input$lineages_plot_zcl)
    
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
      
      ############# data reading ##############################################
      
      # name of the .rds file that contains the full list of data frames
      df_filename <- "./datasets/ZCL_lineages_count_sums.rds"
      
      # read the rds file 
      data <- readRDS(df_filename)
      
      # generate a name key to obtain a data frame
      df_key <- paste0("ZCL_", proc_gene, "_count_sums")
      
      # Validate existence of the key inside the list
      validate(
        need(df_key %in% names(data) , paste("Error: No data file found for gene", input$zcl_gene_id))
      )
      
      # Read the data frame    
      data <- data[[df_key]]
      ##########################################################################
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "counts_sum")
      data_long <- data_long[data_long$counts_sum > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
                                       
      
    }else{
      
      ############# data reading ##############################################
      
      # name of the .rds file that contains the full list of data frames
      df_filename <- "./datasets/ZCL_lineages_cell_counts.rds"
      
      # read the rds file 
      data <- readRDS(df_filename)
      
      # generate a name key to obtain a data frame
      df_key <- paste0("ZCL_", proc_gene, "_cell_counts")
      
      # Validate existence of the key inside the list
      validate(
        need(df_key %in% names(data) , paste("Error: No data file found for gene", input$zcl_gene_id))
      )
      
      # Read the data frame    
      data <- data[[df_key]]
      ##########################################################################
      
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

    
    # 8. Sorting the data 
    if(input$lineage_data_type == "count_sum"){
      
      if(input$zcl_lin_term_sort == "sum"){
        # sort the tissues by their total value
        data <- data %>%
          group_by(tissue) %>%
          mutate(Total_sum = sum(counts_sum)) %>%
          arrange(desc(Total_sum))
      }
      
      data$tissue <- factor(data$tissue, levels = unique(data$tissue))
      
    } else {
      
      if(input$zcl_lin_term_sort == "sum"){
        # sort the tissues by their total value
        data <- data %>%
          group_by(tissue) %>%
          mutate(Total_sum = sum(cell_counts)) %>%
          arrange(desc(Total_sum))
      }
      
      data$tissue <- factor(data$tissue, levels = unique(data$tissue))
      
    }
    
    
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
      
      ########################### save the data ################################
      write.csv(data, "./ZCL_dataset_lineages/temp.csv", row.names = FALSE)
      
      
      # define the plot depending on the y scale that was specified in the input
      if(input$lineage_y_scale == "linear"){
        
        p_ggiraph <- ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
          geom_col_interactive(  aes(tooltip = paste0("Stage: ", stage, 
                                 "<br>", y_label, ": ", round(.data[[y_value]], 2),
                                 "<br>Tissue: ", tissue),
                data_id = paste(stage, tissue, sep = "_"))
          ) +
          scale_y_continuous(n.breaks = 5) +
          facet_wrap(~tissue, scales = "fixed", axes = "all_x", ncol = 5) +
          ylab(y_label) +
          ggtitle(paste("Summarised expression plot for", gene, "in zebrafish")) +
          theme_bw() + 
          scale_fill_paletteer_d("colorBlindness::paletteMartin") +
          theme(plot.title = element_text(size = 8),
                strip.text.x = element_text(size = 8, margin = margin(0.1,0,0.1,0, "cm")),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 10),
                axis.title.x = element_text(size = 10),
                axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 7.5, face="plain"),
                legend.title = element_text(size = 8, face = "bold"),
                legend.text = element_text(size = 8, face = "plain"),
                legend.key.size = unit(0.15, "cm"),
                panel.spacing = unit(0.1, "lines")) 
        
        return(p_ggiraph)
        
      } else{
        
        p_ggiraph <- ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
          geom_col_interactive( aes(tooltip = paste0("Stage: ", stage, 
                                 "<br>", y_label, ": ", round(.data[[y_value]], 2),
                                 "<br>Tissue: ", tissue),
                data_id = paste(stage, tissue, sep = "_"))
          ) +
          scale_y_sqrt(n.breaks = 5) +
          facet_wrap(~tissue, scales = "fixed", axes = "all_x", ncol = 5) +
          ylab(y_label) +
          ggtitle(paste("Summarised expression plot for", gene, "in zebrafish")) +
          theme_bw() + 
          scale_fill_paletteer_d("colorBlindness::paletteMartin") +
          theme(plot.title = element_text(size = 8),
                strip.text.x = element_text(size = 8, margin = margin(0.1,0,0.1,0, "cm")),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 10),
                axis.title.x = element_text(size = 10),
                axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 7.5, face="plain"),
                legend.title = element_text(size = 8, face = "bold"),
                legend.text = element_text(size = 8, face = "plain"),
                legend.key.size = unit(0.15, "cm"),
                panel.spacing = unit(0.1, "lines"))      
        
        return(p_ggiraph)
        
      } # else for specific scale plot
    } # else for the plot code
    
  })
  
  
  # ------------- ZCL cell types data ---------------------------------------------
  
  # 1. Reactive filtering of data based on user inputs
  ZCL_celltypes_data <- eventReactive(input$celltypes_plot_zcl, {
    
    # validate the inputs
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
    #req(input$celltypes_plot_zcl)
    
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
    
    # make a vector of cell types - these can be used directly since they are not
    # column names
    selected_tissues <- c(input$zcl_cell_types_1, input$zcl_cell_types_2, input$zcl_cell_types_3, input$zcl_cell_types_4)
    

    # collect stages    
    selected_stages <- c(input$celltypes_zcl_stages)
    
    # 4. Load the relevant data frame based on the other options
    
    # obtain the processed gene name
    proc_gene <- unique(zcl_genes[zcl_genes$orig_gene == gene_id,]$proc_gene)
    
    
    # read counts option
    if(input$celltypes_data_type == "count_sum"){
      
      ############# data reading ##############################################
      
      # name of the .rds file that contains the full list of data frames
      df_filename <- "./datasets/ZCL_celltypes_count_sums.rds"
      
      # read the rds file 
      data <- readRDS(df_filename)
      
      # generate a name key to obtain a data frame
      df_key <- paste0("ZCL_", proc_gene, "_count_sums")
      
      # Validate existence of the key inside the list
      validate(
        need(df_key %in% names(data) , paste("Error: No data file found for gene", input$zcl_gene_id))
      )
      
      # Read the data frame    
      data <- data[[df_key]]
      ##########################################################################
      
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "counts_sum")
      data_long <- data_long[data_long$counts_sum > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
      
      
    }else{
      
      ############# data reading ##############################################
      
      # name of the .rds file that contains the full list of data frames
      df_filename <- "./datasets/ZCL_celltypes_cell_counts.rds"
      
      # read the rds file 
      data <- readRDS(df_filename)
      
      # generate a name key to obtain a data frame
      df_key <- paste0("ZCL_", proc_gene, "_cell_counts")
      
      # Validate existence of the key inside the list
      validate(
        need(df_key %in% names(data) , paste("Error: No data file found for gene", input$zcl_gene_id))
      )
      
      # Read the data frame    
      data <- data[[df_key]]
      ##########################################################################
      
      
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
    
    # convert tissue and selected_tissues 
    # previous values are names and have '_' and values are new ones 
    map_tissues <- c("Cardiomyocyte" = "Cardiomyocyte", "Endothelial_cell" = "Endothelial cell", 
                     "Enterocyte" = "Enterocyte", "Epithelial_cell" = "Epithelial cell", 
                     "Epithelial_cell__Brain_" = "Epithelial cell (Brain)", "Erythrocyte" = "Erythrocyte", 
                     "Erythrocyte__Liver_" = "Erythrocyte (Liver)", "Erythroid_progenitor_cell" = "Erythroid progenitor cell", 
                     "Fibroblast" = "Fibroblast", "Goblet_cell" = "Goblet cell", "Granulocyte" = "Granulocyte",
                     "Granulosa_cell" = "Granulosa cell", "Hatching_gland" = "Hatching gland",  "Hepatocyte" = "Hepatocyte", 
                     "Immune_cell" = "Immune cell","Immune_progenitor_cell" = "Immune progenitor cell", 
                     "Intestinal_bulb_cell" = "Intestinal bulb cell", "Ionocyte" = "Ionocyte", "Keratinocyte" = "Keratinocyte",
                     "Macrophage" = "Macrophage", "Mesenchymal_cell" = "Mesenchymal cell",
                     "Mesenchymal_cell__caudal fin_" = "Mesenchymal cell (caudal fin)", "Mt-rich_cell" = "Mt-rich cell",
                     "Muscle_cell" = "Muscle cell", "Nephron_cell" = "Nephron cell",
                     "Neural_cell" = "Neural cell", "Neural_progenitor_cell" = "Neural progenitor cell",
                     "Neurosecretory_cell" = "Neurosecretory cell", "Oligodendrocyte" = "Oligodendrocyte",
                     "Oocyte" = "Oocyte", "Osteoblast" = "Osteoblast", "Pancreatic_cell" = "Pancreatic cell", 
                     "Pancreatic_macrophage" = "Pancreatic macrophage", "Primordial_germ_cell" = "Primordial germ cell", 
                     "Radial_glia" = "Radial glia", "Retinal_cell" = "Retinal cell", "Retinal_cone_cell" = "Retinal cone cell",
                     "Retinal_pigment_epithelial_cell" = "Retinal pigment epithelial cell", "Smooth_muscle_cell" = "Smooth muscle cell",
                     "Spermatocyte" = "Spermatocyte", "T_cell" = "T cell") 
      
    
    # Replace values in the column using the named vector for indexing
    data_long$tissue <- map_tissues[as.character(data_long$tissue)]
    
    # refresh the tissue factor
    data_long$tissue <- factor(data_long$tissue, levels = zcl_cell_types$cell_type)
    
    
    # 7. filter the data frame by the selected tissues and stages
    data <- filter(data_long, tissue %in% selected_tissues & stage %in% selected_stages)
    gene <- gene_id
    
    # 8. Sorting the data 
    if(input$celltypes_data_type == "count_sum"){
      
      if(input$zcl_celltype_term_sort == "sum"){
        # sort the tissues by their total value
        data <- data %>%
          group_by(tissue) %>%
          mutate(Total_sum = sum(counts_sum)) %>%
          arrange(desc(Total_sum))
      }
        
      data$tissue <- factor(data$tissue, levels = unique(data$tissue))
      
    } else {
      
      if(input$zcl_celltype_term_sort == "sum"){
        # sort the tissues by their total value
        data <- data %>%
          group_by(tissue) %>%
          mutate(Total_sum = sum(cell_counts)) %>%
          arrange(desc(Total_sum))
      }
      
      data$tissue <- factor(data$tissue, levels = unique(data$tissue))
      
    }
    
    
    
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
      
      ########################### save the data ################################
      write.csv(data, "./ZCL_dataset_celltypes/temp.csv", row.names = FALSE)
      
      # define the plot depending on the y scale that was specified in the input
      if(input$celltypes_y_scale == "linear"){
        
        p_ggiraph <- ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
          geom_col_interactive( aes(tooltip = paste0("Stage: ", stage, 
                                 "<br>", y_label, ": ", round(.data[[y_value]], 2),
                                 "<br>Tissue: ", tissue),
                data_id = paste(stage, tissue, sep = "_"))
          ) +
          scale_y_continuous(n.breaks = 5) +
          facet_wrap(~tissue, scales = "fixed", axes = "all_x", ncol = 5) +
          ylab(y_label) +
          ggtitle(paste("Summarised expression plot for", gene, "in zebrafish")) +
          theme_bw() + 
          scale_fill_paletteer_d("colorBlindness::paletteMartin") +
          theme(plot.title = element_text(size = 8),
                strip.text.x = element_text(size = 8, margin = margin(0.1,0,0.1,0, "cm")),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 10),
                axis.title.x = element_text(size = 10),
                axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 7.5, face="plain"),
                legend.title = element_text(size = 8, face = "bold"),
                legend.text = element_text(size = 8, face = "plain"),
                legend.key.size = unit(0.15, "cm"),
                panel.spacing = unit(0.1, "lines"))   
        
        return(p_ggiraph)
        
        
      } else{
        
        p_ggiraph <- ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
          geom_col_interactive( aes(tooltip = paste0("Stage: ", stage, 
                                 "<br>", y_label, ": ", round(.data[[y_value]], 2),
                                 "<br>Tissue: ", tissue),
                data_id = paste(stage, tissue, sep = "_"))
          ) +
          scale_y_sqrt(n.breaks = 5) +
          facet_wrap(~tissue, scales = "fixed", axes = "all_x", ncol = 5) +
          ylab(y_label) +
          ggtitle(paste("Summarised expression plot for", gene, "in zebrafish")) +
          theme_bw() + 
          scale_fill_paletteer_d("colorBlindness::paletteMartin") +
          theme(plot.title = element_text(size = 8),
                strip.text.x = element_text(size = 8, margin = margin(0.1,0,0.1,0, "cm")),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 10),
                axis.title.x = element_text(size = 10),
                axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 7.5, face="plain"),
                legend.title = element_text(size = 8, face = "bold"),
                legend.text = element_text(size = 8, face = "plain"),
                legend.key.size = unit(0.15, "cm"),
                panel.spacing = unit(0.1, "lines"))        
        
        return(p_ggiraph)
        
      } # else for specific scale plot
    } # else for the plot code
    
  })
  
  
  # 
  # 2. Render the plot
  output$zcl_expr_plot <- renderGirafe({
    req(input$lineages_plot_zcl | input$celltypes_plot_zcl )

    # test
    if( rv$lastBtn == "lineages_plot_zcl" ){
      
      # run the graph output function to generate the graph and to save the data
      # that can be used
      output_obj <- ZCL_lineage_data()
      
      # 2. read the data and define how many facets there will be
      df <- read.csv("./ZCL_dataset_lineages/temp.csv")
      n_facets = length(unique(df$tissue))
      
      if(n_facets < 5){
        calc_width = n_facets * FACET_SIZE + 0.5
        calc_height = FACET_SIZE + 0.5
      } else{
        # do the final calculation
        calc_width = 5 * FACET_SIZE
        calc_height = ceiling(n_facets/5)*FACET_SIZE + 1
      }
      
      
      # Convert the ggplot object to a girafe object
      return( girafe(ggobj = output_obj, width_svg = calc_width, height_svg = calc_height, 
                     options = list(opts_sizing(rescale = TRUE), opts_zoom = opts_zoom(min = 0.25, max = 0.5) ) ) )
    }

    if( rv$lastBtn == "celltypes_plot_zcl" ){

      # run the graph output function to generate the graph and to save the data
      # that can be used
      output_obj <- ZCL_celltypes_data()
      
      # 2. read the data and define how many facets there will be
      df <- read.csv("./ZCL_dataset_celltypes/temp.csv")
      n_facets = length(unique(df$tissue))
      
      if(n_facets < 5){
        calc_width = n_facets * FACET_SIZE + 0.5
        calc_height = FACET_SIZE + 0.5
      } else{
        # do the final calculation
        calc_width = 5 * FACET_SIZE
        calc_height = ceiling(n_facets/5)*FACET_SIZE + 1
      }
      
      # Convert the ggplot object to a girafe object
      return( girafe(ggobj = output_obj, width_svg = calc_width, height_svg = calc_height, 
                     options = list(opts_sizing(rescale = TRUE), opts_zoom = opts_zoom(min = 0.25, max = 0.5) ) ) )
      
    }

  })
    
  
  # Handle the PNG download using the base_plot
  output$zcl_download_png <- downloadHandler(
    filename = function() {
      paste0("ZCL_plot_", input$zcl_gene_id, '_', Sys.Date(), ".png")
    },
    content = function(file) {
      
      # test which button was last pressed
      if( rv$lastBtn == "lineages_plot_zcl" ){
        result <- ZCL_lineage_data()
        # obtain the data for the size parameters of the plot
        df <- read.csv("./ZCL_dataset_lineages/temp.csv")
      }
      
      if( rv$lastBtn == "celltypes_plot_zcl" ){
        result <-  ZCL_celltypes_data()
        # obtain the data for the size parameters of the plot
        df <- read.csv("./ZCL_dataset_celltypes/temp.csv")
      }
      
      # get a measure of the data size
      n_facets = length(unique(df$tissue))
      
      if(n_facets < 5){
        calc_width = n_facets * FACET_SIZE + 0.5
        calc_height = FACET_SIZE + 0.5
      } else{
        # do the final calculation
        calc_width = 5 * FACET_SIZE
        calc_height = ceiling(n_facets/5)*FACET_SIZE + 1
      }
      
      # Use ggsave on the ggplot object      
      ggsave(file, plot = result, device = "png", dpi = 300, 
             width = calc_width, height = calc_height)
      
    }
  )

    # Handle data download
  output$zcl_download_data <- downloadHandler(
    filename = function() {
      paste0("ZCL_data_", input$zcl_gene_id, '_', Sys.Date(), ".csv")
    },
    content = function(file) {

      
      if( rv$lastBtn == "lineages_plot_zcl" ){
        data <- read.csv("./ZCL_dataset_lineages/temp.csv")
      }
      
      if( rv$lastBtn == "celltypes_plot_zcl" ){
        data <- read.csv("./ZCL_dataset_celltypes/temp.csv")
      }
      
      
      # Write the data to a temporary file
      write.csv(data, file, row.names = FALSE)
    }
  )

  
  # Handle data download
  output$zcl_download_code <- downloadHandler(
    
    filename = function() {
      paste0("Plot_code_ZCL_", input$zcl_gene_id, '_', Sys.Date(), ".R")  
    },
    content = function(file) {
      
      # generate a file name for the data in this session
      # use this file name to insert it into the code file
      data_filename <- paste0("ZCL_data_", input$zcl_gene_id, '_', Sys.Date(), ".csv")
      
      if( rv$lastBtn == "lineages_plot_zcl" ){
        
        # load the templates as appropriate
        if(input$lineage_data_type == "count_sum"){
          
          # choose the scale to load the templates
          if(input$lineage_y_scale == "linear"){
            
            # count_sum - linear
            file_path <- "./templates/counts_sum_linear.R" 
            script_code <- read_file(file_path)
            
          } else{
            # count_sum - sqrt
            file_path <- "./templates/counts_sum_sqrt.R" 
            script_code <- read_file(file_path)          
            
          }
          
          
        } else {
          
          # choose the scale to load the templates
          if(input$lineage_y_scale == "linear"){
            # count_sum - linear
            file_path <- "./templates/cell_counts_linear.R" 
            script_code <- read_file(file_path)
            
          } else{
            # count_sum - sqrt
            file_path <- "./templates/cell_counts_sqrt.R" 
            script_code <- read_file(file_path)
            
          }
          
        }
        
      }
      
      
      if( rv$lastBtn == "celltypes_plot_zcl" ){
        
        # load the templates as appropriate
        if(input$celltypes_data_type == "count_sum"){
          
          # choose the scale to load the templates
          if(input$celltypes_y_scale == "linear"){
            
            # count_sum - linear
            file_path <- "./templates/counts_sum_linear.R" 
            script_code <- read_file(file_path)
            
          } else{
            # count_sum - sqrt
            file_path <- "./templates/counts_sum_sqrt.R" 
            script_code <- read_file(file_path)          
            
          }
          
          
        } else {
          
          # choose the scale to load the templates
          if(input$celltypes_y_scale == "linear"){
            # count_sum - linear
            file_path <- "./templates/cell_counts_linear.R" 
            script_code <- read_file(file_path)
            
          } else{
            # count_sum - sqrt
            file_path <- "./templates/cell_counts_sqrt.R" 
            script_code <- read_file(file_path)
            
          }
          
        }
        
      }
      
      # replace temp.csv with the current data file name
      script_code <- str_replace(script_code, 'temp.csv', data_filename)
      
      # insert gene name
      script_code <- str_replace(script_code, 'current_gene', input$zcl_gene_id)
      
      
      # insert code for sorting stages 
      script_code <- str_replace(script_code, '# sorting code', '#sort the data\n data$stage <- factor(data$stage, levels = c("24hpf", "72hpf","21Day", "3Month","22Month"))')
      
      # Write the data to a temporary file
      writeChar(script_code, file)
    }
  )
  
  
  # ------------- Zebrahub data ----------------------------------------------
  
  # 1. Reactive filtering of data based on user inputs
  zhub_data <- eventReactive(input$update_plot_zhub, {
    
    # validate the inputs
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
    #req(input$update_plot_zhub)
    
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
      
      ############# data reading ##############################################
      
      # name of the .rds file that contains the full list of data frames
      df_filename <- "./datasets/zhub_count_sums.rds"
      
      # read the rds file 
      data <- readRDS(df_filename)
      
      # generate a name key to obtain a data frame
      df_key <- paste0("zebrahub_", proc_gene, "_count_sums")
      
      # Validate existence of the key inside the list
      validate(
        need(df_key %in% names(data) , paste("Error: No data file found for gene", input$zhub_gene_id))
      )
      
      # Read the data frame    
      data <- data[[df_key]]
      ##########################################################################
      
      # 6. Convert the data frame to the long format
      data_long <- pivot_longer(data, cols = colnames(data)[2:ncol(data)], names_to = "tissue", values_to = "counts_sum")
      data_long <- data_long[data_long$counts_sum > 0, ]
      
      data_long <- data_long |> 
        group_by(tissue, stage)  |> 
        arrange(tissue, stage, .by_group = TRUE)
      
      
    }else{
      
      ############# data reading ##############################################
      
      # name of the .rds file that contains the full list of data frames
      df_filename <- "./datasets/zhub_cell_counts.rds"
      
      # read the rds file 
      data <- readRDS(df_filename)
      
      # generate a name key to obtain a data frame
      df_key <- paste0("zebrahub_", proc_gene, "_cell_counts")
      
      # Validate existence of the key inside the list
      validate(
        need(df_key %in% names(data) , paste("Error: No data file found for gene", input$zhub_gene_id))
      )
      
      # Read the data frame    
      data <- data[[df_key]]
      ##########################################################################
      
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
  
    
    # 8. Sorting the data 
    if(input$zhub_data_type == "count_sum"){
      
      
      if(input$zhub_term_sort == "sum"){
        # sort the tissues by their total value
        data <- data %>%
          group_by(tissue) %>%
          mutate(Total_sum = sum(counts_sum)) %>%
          arrange(desc(Total_sum))
      }
      
      data$tissue <- factor(data$tissue, levels = unique(data$tissue))
      
    } else {
      
      if(input$zhub_term_sort == "sum"){
        # sort the tissues by their total value
        data <- data %>%
          group_by(tissue) %>%
          mutate(Total_sum = sum(cell_counts)) %>%
          arrange(desc(Total_sum))
      }
              
      data$tissue <- factor(data$tissue, levels = unique(data$tissue))
      
    }
    
      
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
      
      ########################### save the data ################################
      write.csv(data, "./Zebrahub_dataset/temp.csv", row.names = FALSE)
      
      # define the plot depending on the y scale that was specified in the input
      if(input$zhub_y_scale == "linear"){
        
        p_ggiraph <- ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
          geom_col_interactive( aes(tooltip = paste0("Stage: ", stage, 
                                 "<br>", y_label, ": ", round(.data[[y_value]], 2),
                                 "<br>Tissue: ", tissue),
                data_id = paste(stage, tissue, sep = "_"))
          ) +
          scale_y_continuous(n.breaks = 5) +
          facet_wrap(~tissue, scales = "fixed", axes = "all_x", ncol = 5) +
          ylab(y_label) +
          ggtitle(paste("Summarised expression plot for", gene, "in zebrafish")) +
          theme_bw() + 
          scale_fill_paletteer_d("colorBlindness::paletteMartin") +
          theme(plot.title = element_text(size = 8),
                strip.text.x = element_text(size = 8, margin = margin(0.1,0,0.1,0, "cm")),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 10),
                axis.title.x = element_text(size = 10),
                axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 7.5, face="plain"),
                legend.title = element_text(size = 8, face = "bold"),
                legend.text = element_text(size = 8, face = "plain"),
                legend.key.size = unit(0.15, "cm"),
                panel.spacing = unit(0.1, "lines")) 
        
        return(p_ggiraph) 
        
      } else{
        
        p_ggiraph <- ggplot(data, aes(x = stage, y = .data[[y_value]], fill = stage)) +
          geom_col_interactive( aes(tooltip = paste0("Stage: ", stage, 
                                 "<br>", y_label, ": ", round(.data[[y_value]], 2),
                                 "<br>Tissue: ", tissue),
                data_id = paste(stage, tissue, sep = "_"))
          ) +
          scale_y_sqrt(n.breaks = 5) +
          facet_wrap(~tissue, scales = "fixed", axes = "all_x", ncol = 5) +
          ylab(y_label) +
          ggtitle(paste("Summarised expression plot for", gene, "in zebrafish")) +
          theme_bw() + 
          scale_fill_paletteer_d("colorBlindness::paletteMartin") +
          theme(plot.title = element_text(size = 8),
                strip.text.x = element_text(size = 8, margin = margin(0.1,0,0.1,0, "cm")),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 10),
                axis.title.x = element_text(size = 10),
                axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 7.5, face="plain"),
                legend.title = element_text(size = 8, face = "bold"),
                legend.text = element_text(size = 8, face = "plain"),
                legend.key.size = unit(0.15, "cm"),
                panel.spacing = unit(0.1, "lines"))       
        
        return(p_ggiraph)
        
      } # else for specific scale plot
    } # else for the plot code
    
  })
  
  
  # 2. Render the plot
  output$zebrahub_expr_plot <- renderGirafe({
    
    # 1. Ensure the 'Generate Plot' button has been clicked at least once
    req(input$update_plot_zhub)

    # run the graph output function to generate the graph and to save the data
    # that can be used
    output_obj <- zhub_data()
    
    # 2. read the data and define how many facets there will be
    df <- read.csv("./Zebrahub_dataset/temp.csv")
    n_facets = length(unique(df$tissue))
    
    if(n_facets < 5){
      calc_width = n_facets * FACET_SIZE + 0.5
      calc_height = FACET_SIZE + 0.5
    } else{
      # do the final calculation
      calc_width = 5 * FACET_SIZE
      calc_height = ceiling(n_facets/5)*FACET_SIZE + 1
    }
    
    # Convert the ggplot object to a girafe object
    return( girafe(ggobj = output_obj, width_svg = calc_width, height_svg = calc_height, 
                   options = list(opts_sizing(rescale = TRUE), opts_zoom = opts_zoom(min = 0.25, max = 0.5) )) )
    
  })
  
  
  
  
  # Handle data download
  output$zhub_download_data <- downloadHandler(
    filename = function() {
      paste0("Zebrahub_data_", input$zhub_gene_id, '_', Sys.Date(), ".csv")
    },
    content = function(file) {
      
      data <- read.csv("./Zebrahub_dataset/temp.csv")
      
      # Write the data to a temporary file
      write.csv(data, file, row.names = FALSE)
    }
  )
  
  
  # Handle the PNG download using the base_plot
  output$zhub_download_png <- downloadHandler(
    filename = function() {
      paste0("Zebrahub_plot_", input$zhub_gene_id, '_', Sys.Date(), ".png")
    },
    content = function(file) {
      
      # obtain the data for the size parameters of the plot
      df <- read.csv("./Zebrahub_dataset/temp.csv")
      n_facets = length(unique(df$tissue))
      
      if(n_facets < 5){
        calc_width = n_facets * FACET_SIZE + 0.5
        calc_height = FACET_SIZE + 0.5
      } else{
        # do the final calculation
        calc_width = 5 * FACET_SIZE
        calc_height = ceiling(n_facets/5)*FACET_SIZE + 1
      }
      
      # Use ggsave on the ggplot object      
      ggsave(file, plot = zhub_data(), device = "png", dpi = 300, 
             width = calc_width, height = calc_height)
      
      }
  )
  
  
  # Handle data download
  output$zhub_download_code <- downloadHandler(
    
    filename = function() {
      paste0("Plot_code_Zebrahub_", input$zhub_gene_id, '_', Sys.Date(), ".R")  
    },
    content = function(file) {
      
      # generate a file name for the data in this session
      # use this file name to insert it into the code file
      data_filename <- paste0("Zebrahub_data_", input$zhub_gene_id, '_', Sys.Date(), ".csv")
      
      # load the templates as appropriate
      
      if(input$zhub_data_type == "count_sum"){
        
        # choose the scale to load the templates
        if(input$zhub_y_scale == "linear"){
          
          # count_sum - linear
          file_path <- "./templates/counts_sum_linear.R" 
          script_code <- read_file(file_path)
          
        } else{
          # count_sum - sqrt
          file_path <- "./templates/counts_sum_sqrt.R" 
          script_code <- read_file(file_path)          
          
        }
        
        
      } else {
        
        # choose the scale to load the templates
        if(input$zhub_y_scale == "linear"){
          # count_sum - linear
          file_path <- "./templates/cell_counts_linear.R" 
          script_code <- read_file(file_path)
          
        } else{
          # count_sum - sqrt
          file_path <- "./templates/cell_counts_sqrt.R" 
          script_code <- read_file(file_path)
          
        }
        
      }
      
      # replace temp.csv with the current data file name
      script_code <- str_replace(script_code, 'temp.csv', data_filename)
      
      # insert gene name
      script_code <- str_replace(script_code, 'current_gene', input$zhub_gene_id)
      
      # insert code for sorting stages 
      script_code <- str_replace(script_code, '# sorting code', '#sort the data\n data$stage <- factor(data$stage, levels = c("10hpf","12hpf","14hpf","16hpf","19hpf","24hpf","2dpf","3dpf","5dpf","10dpf"))')
      
      # Write the data to a temporary file
      writeChar(script_code, file)
    }
  )
  
  
  
}

# Run the complete application
shinyApp(ui = ui, server = server)