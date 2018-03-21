### barplot
barplot <- reactive({
  barplot_cols <- subCol()
  
  tmp <- cbind(dataSubmit$annot[barplot_cols, input$barplot.cat, drop = F], t(dataSubmit$expr[input$barplot.gene, barplot_cols, drop =F]))
  tmp <- tmp[order(tmp[,input$barplot.cat]),]
  data <- melt(tmp, id.vars = input$barplot.cat, measure.vars = input$barplot.gene, variable.names = input$barplot.cat)
  #save(data, file = "tmp.Rdata")
  p <- ggplot()+ geom_col(data = data, 
                          mapping = aes(x = rep(1:nrow(tmp), times = length(input$barplot.gene)), y = value, fill = get(input$barplot.cat)), 
                          position = "identity") +
    facet_wrap(~variable, scales = "free_y") + 
    theme_classic() + theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,legend.position = "bottom"
      ,axis.text.x = element_text(angle = 45, hjust = 1)
      ,plot.title = element_text(hjust = 0.5, face = "bold")
    ) + guides(fill = guide_legend(nrow = 1)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    xlab(input$barplot.cat) + ylab(paste0("Expression level")) + ggtitle("Barplot of Gene Expression")
  p <- p + scale_fill_manual(values = c(dataSubmit$color[[input$barplot.cat]], rep(pal_dolphin, 100)), name = input$barplot.cat) 
  return(p)
  # ggplotly()
})
output$barplot <- renderPlot({
  barplot()
})
output$barplot.ui <- renderUI({
  plotOutput("barplot",height = input$barplot.height, width = input$barplot.width)
})
output$barplotDL <- downloadHandler(
  filename = function() {
    paste(input$barplot.gene, Sys.Date(), 'pdf', sep='.')
  },
  content = function(file){
    ggsave(filename = file, plot = barplot(), width = round(input$barplot.width/65),height = round(input$barplot.height/80))
  }
)