output$vlnplot <- renderPlot({
  
  vlnplot_cols <- subCol()
  tmp <- cbind(dataSubmit$annot[vlnplot_cols, input$vlnplot.cat, drop = F], t(dataSubmit$expr[input$vlnplot.gene, vlnplot_cols, drop =F]))
  data <- melt(tmp, id.vars = input$vlnplot.cat, measure.vars = input$vlnplot.gene, variable.names = input$vlnplot.cat)
  p <- ggplot(data = data, mapping = aes(x = get(input$vlnplot.cat), y = value)) + 
    geom_violin(mapping = aes(fill = get(input$vlnplot.cat)), trim = T, scale = "width") +
    geom_boxplot(width = 0.05) + 
    theme_bw() + theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      #,panel.border = element_blank()
      ,axis.text.x = element_text(angle = 45, hjust = 1)
      ,plot.title = element_text(hjust = 0.5, face = "bold")
    ) + xlab(input$vlnplot.cat) + ylab(paste0("Expression level")) + ggtitle(input$vlnplot.gene)
  p <- p + scale_fill_manual(values = c(dataSubmit$color[[input$vlnplot.cat]], rep(pal_dolphin, 100)), name = input$vlnplot.cat) 
  return(p)
  #ggplotly()
})
output$vlnplot.ui <- renderUI({
  plotOutput("vlnplot",height = input$vlnplot.height, width = input$vlnplot.width)
})