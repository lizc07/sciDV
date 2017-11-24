############
# ui
#########
# header ---------------------------------------------------------------------------------------------------------------------
header  <- dashboardHeader(title = "sciDV", titleWidth = 320, disable = F)
# sidebar ---------------------------------------------------------------------------------------------------------------------
sidebar <- dashboardSidebar(width = 320, collapsed = F,
            sidebarMenu(
              # Setting id makes input$tabs give the tabName of currently-selected tab
              id = "tabs", 
              menuItem("Load Data", tabName = "fileInput", icon = icon("dashboard"),badgeLabel = "Step0"),
              menuItem("Heatmap", icon = icon("th"),
                       textInput("tabHeatmap.select_genes", label = ("Input Genes"), placeholder = "Input gene names..."),
                       menuItem("Plot the Heatmap", icon = icon("photo"), tabName = "tabHeatmap",badgeLabel = "Action", badgeColor = "green"),
                       menuItem("Graph Size Parameters", icon = icon("sliders"),
                                sliderInput("tabHeatmap.height", label = "Graph Height", min = 0, max = 2000, value = 600),
                                sliderInput("tabHeatmap.treeheight", label = "Dendrogram Depth", min = 0, max = 500, value = 100)
                       ),
                       menuItem("Clustering Parameters",icon = icon("check-square"),
                                
                                checkboxInput("tabHeatmap.rowClust", label = "Gene clustering?", value = FALSE),
                                radioButtons("tabHeatmap.rowClust.dist", label = NA, choices = list("cor.pearson","cor.spearman","euclidean"), selected = "cor.pearson",inline = T),
                                checkboxInput("tabHeatmap.colClust", label = "Sample clustering?", value = FALSE),
                                radioButtons("tabHeatmap.colClust.dist", label = NA, choices = list("cor.pearson","cor.spearman","euclidean"), selected = "cor.pearson",inline = T),
                                radioButtons("tabHeatmap.clust.linkage", label = "Clustering linkage:", choices = list("ward.D","single","complete","average"), selected = "ward.D",inline = T)
                       ),
                       menuItem("Display Parameters", icon = icon("list-ul"),
                                uiOutput("menuHeatmap.displayPara", inline = T)
                       )
                       
              ),
              menuItem("Dimension Plot", icon = icon("bar-chart-o"), 
                       menuItem("Plot the Dimensions", icon = icon("photo"),tabName = "tabDimPlot",badgeLabel = "Action", badgeColor = "green"),
                       menuItem("Display Parameters", icon = icon("list-ul"), 
                                uiOutput("menuDimPlot.displayPara", inline = T)
                                )
              ),
              menuItem("Violin Plot", icon = icon("bar-chart"),
                       menuItem("Plot the vlnplot", icon = icon("photo"), tabName = "tabGeneVlnplot",badgeLabel = "Action", badgeColor = "green"),
                       menuItem("Display Parameters", icon = icon("list-ul"), 
                                uiOutput("menuVlnPlot.displayPara", inline = T)
                       )
              ),
              menuItem("Bar Plot", icon = icon("bar-chart"),
                       menuItem("Plot the barplot", icon = icon("photo"), tabName = "tabGeneBarplot",badgeLabel = "Action", badgeColor = "green"),
                       menuItem("Display Parameters", icon = icon("list-ul"), 
                                uiOutput("menuBarPlot.displayPara", inline = T)
                       )
              ),
              menuItem("Subset Cells", icon = icon("hourglass-end"),
                       uiOutput("showSubcat")
              )
              
            ))
# body ---------------------------------------------------------------------------------------------------------------------
body <- dashboardBody(
  tabItems(
    tabItem(tabName = "fileInput", 
            box(title = "Loading RData File", width = 12, footer = "Based on shiny.", status = "primary", solidHeader = T, 
                fileInput("fileInput.rdata", "Saved RData, including Expression and Annotation Data", width = "100%",
                          accept = c(".RData")
                          ),
                actionButton(inputId = "fileInput.action", width = "100%", label = "Submit"),
                hr(),
                infoBoxOutput("infoBox.expr", width = 6),
                infoBoxOutput("infoBox.annot", width = 6),
                uiOutput("fileInput.submitted")
            ),
            
            uiOutput("fileInput.error")
            
    ),
    tabItem(tabName = "tabHeatmap",
            fluidRow(
              tabBox(title = h6("Created  by Zongcheng Li, using pheatmap"), width = 12, side = "left",
                     tabPanel(title = "Heatmap Graph", icon = icon("photo"),
                              uiOutput("heatmap.ui")),
                     tabPanel(title = "Heatmap Info", icon = icon("info-circle"),
                              tableOutput("heatmap.info.tab"),htmlOutput("heatmap.info.text"))
              )
            )
    ),
    tabItem(tabName = "tabDimPlot", 
            uiOutput("showDimPlot.gene"),
            uiOutput("showDimPlot.cat")
            
      
    ),
    tabItem(tabName = "tabGeneVlnplot",
            fluidRow(
              box(title = "Gene-centric Violin Plot", width = 12,side = "left", footer = "Created  by Zongcheng Li, using ggplot2",
                  status = "success", solidHeader = T, collapsible = T, collapsed = F,
                  uiOutput("vlnplot.ui")
              )
            )
    ),
    tabItem(tabName = "tabGeneBarplot",
            fluidRow(
              box(title = "Gene-centric Bar Plot", width = 12,side = "left", footer = "Created  by Zongcheng Li, using ggplot2",
                  status = "success", solidHeader = T, collapsible = T, collapsed = F,
                  uiOutput("barplot.ui"),downloadButton('barplotDL', 'Download')
              )
            )
    ),
    # Cluster -------------------------------------------------------------------------
    tabItem(tabName = "tSNE and Clustering"
    ),
    tabItem(tabName = "pointPlot"
      # DEG
    )
    
  )
)
# page ---------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(title = "Interactive Data Visualization for single cell RNA-Seq", skin = "purple",
              header = header, sidebar = sidebar, body = body)
