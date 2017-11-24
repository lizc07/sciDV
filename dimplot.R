output$showDimPlot.gene <- renderUI({
  uis <- names(dataSubmit$coord)
  
  tabs <- lapply(1:length(uis), 
                 function(x) {
                   tabPanel(title = uis[x], 
                            uiOutput(paste0(uis[x],".gene.plot.ui")),
                            icon = icon("send")
                            )
                   }
                 )
  tabs$title <- "Gene-centric"
  do.call(tabBox, tabs)
})
#outputOptions(output, "showDimPlot", suspendWhenHidden = T, priority = -)

output$showDimPlot.cat <- renderUI({
  uis <- names(dataSubmit$coord)
  
  tabs <- lapply(1:length(uis), 
                 function(x) {
                   tabPanel(title = uis[x], 
                            uiOutput(paste0(uis[x],".cat.plot.ui")),
                            icon = icon("send")
                   )
                 }
  )
  tabs$title <- "Category"
  do.call(tabBox, tabs)
})