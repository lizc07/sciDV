############
# server.R
############
server <- function(input, output){
  
  # Load Data
  #incProgress(amount = 0.1, message = "Display infoBox")
  source(file = "upload_data.R", local = T)
  # subcat
  #incProgress(amount = 0.2, message = "prepare subcat")
  source(file = "subcat.R", local = T)
  # Heatmap
  #incProgress(amount = 0.2, message = "show Heatmap")
  source(file = "heatmap.R", local = T)
  # DEG
  source(file = "deg.R", local = T)
  # DimPlot
  #incProgress(amount = 0.1, message = "prepare DimPlot")
  source(file = "dimplot.R", local = T)

  # cell cycle
  source(file = "cellcycle.R", local = T)
  # VlnPlot
  #incProgress(amount = 0.1, message = "prepare VlnPlot")
  source(file = "vlnplot.R", local = T)
  # barplot
  #incProgress(amount = 0.1, message = "prepare BarPlot")
  source(file = "barplot.R", local = T)

}