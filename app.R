#' @include server.R
#' @include ui.R
#' @include global.R
NULL

#'
#' iteractive Data Visualization
#'
#' Interactive Data Visualization for UMI-based tagging scRNA-Seq data. 
#'
#' @param NULL
#'
#' @export
#'
#' @examples
#' # runApp
#' # iDA()
iDV <- function(rdata = NULL, expr = NULL, annot = NULL, colors = NULL, ai = NULL, maxDataSize = 1000){
  options(shiny.maxRequestSize=maxDataSize*1024^2)
  
  if(is.null(rdata))
  pal_heatmap1 <<- colorpanel(100,"blue","white","red")
  pal_heatmap2 <<- colorpanel(100,"darkblue","white","red")
  
  runApp(shinyApp(ui = ui, server = server))
}
