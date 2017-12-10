output$ccplot <- renderPlot({
  ccData <- dataSubmit$annot
  ccth <- input$ccplot.th
  ccx <- ceiling(max(ccData$idv.g1sScore))
  ccy <- ceiling(max(ccData$idv.g2mScore))
  
  specify_geom_point(data = ccData, 
                     one_gene = NULL, aes.x = "idv.g1sScore", aes.y = "idv.g2mScore", aes.color = input$ccplot.cat, 
                     size = input$ccplot.psize, coord_fixed.ratio = NULL, 
                     scale_color_manual = if(input$ccplot.cat %in% names(dataSubmit$color)) dataSubmit$color[[input$ccplot.cat]] else NULL
  ) + scale_x_continuous(expand = c(0, 0), limits = c(0, ccx)) + scale_y_continuous(expand = c(0, 0), limits = c(0,ccy)) +
    geom_linerange(mapping = aes_(x = ccth, ymin = 0, ymax = ccth)) + 
    geom_abline(mapping = aes_(slope = 0, intercept = ccth)) + 
    geom_segment(mapping = aes_(x = ccth, y = ccth, xend = ccx, yend = ccy))
})
output$ccplot.ui <- renderUI({plotOutput("ccplot", height = input$ccplot.height)})

output$ccbar <- renderPlot({
  ccData <- dataSubmit$annot
  ccth <- input$ccplot.th
  ccData$idv.ccphase[ccData$idv.g1sScore < ccth & ccData$idv.g2mScore < ccth] <- "Q"
  ccData$idv.ccphase[!(ccData$idv.g1sScore < ccth & ccData$idv.g2mScore < ccth) & (ccData$idv.g1sScore < ccData$idv.g2mScore)] <- "G2/M"
  ccData$idv.ccphase[!(ccData$idv.g1sScore < ccth & ccData$idv.g2mScore < ccth) & (ccData$idv.g1sScore > ccData$idv.g2mScore) & ccData$idv.g2mScore < ccth] <- "G1"
  ccData$idv.ccphase[!(ccData$idv.g1sScore < ccth & ccData$idv.g2mScore < ccth) & (ccData$idv.g1sScore > ccData$idv.g2mScore) & ccData$idv.g2mScore > ccth] <- "S"
  ccData$idv.ccphase <- factor(ccData$idv.ccphase, levels = c("Q","G1","S","G2/M"))
  ggplot(ccData %>% count_(list(input$ccplot.cat, "idv.ccphase")) %>% group_by_(input$ccplot.cat) %>% arrange(desc(idv.ccphase)) %>%
           mutate(pct = n/sum(n),
                  ypos = cumsum(pct) - 0.5*pct)) + 
    geom_bar(mapping = aes(x = get(input$ccplot.cat), y = pct, fill = idv.ccphase), stat = "identity", width = 0.8) + 
    geom_text(mapping = aes(x = get(input$ccplot.cat), y = ypos, label = paste0(n,"\n",sprintf("%1.1f", pct*100),"%"))) +
    ggtitle(label = "Cell cycle distribution of cells") + 
    theme_bw() + theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      #,panel.border = element_blank()
      ,axis.text.x = element_text(angle = 45, hjust = 1)
      ,plot.title = element_text(hjust = 0.5, face = "bold")
    ) + xlab(input$ccplot.cat) + ylab(paste0("Proportion"))
})
output$ccbar.ui <- renderUI({plotOutput("ccbar", height = input$ccplot.height)})