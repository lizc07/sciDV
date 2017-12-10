observeEvent(eventExpr = input$go.action,
             handlerExpr = {
               withProgress(message = "Calculation in progress ...",{
                 gogenes <- textInput2genes(input$go.genes, org)
                 engo <- list()
                 for(i in input$go.enrichdb){
                   local({
                     engo[[i]] <- enrichGO(gene         = gogenes,
                                           OrgDb         = org.Mm.eg.db,
                                           keytype       = 'SYMBOL',
                                           ont           = gsub("GO:","",i),
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.01,
                                           qvalueCutoff  = 0.05,
                                           readable = T)
                     i_name <- gsub(":","_", i)
                     engo.plotName <- paste0(i_name,".engo.plot")
                     engo.downloadName <- paste0(i_name,".engo.download")
                     engo.uiName <- paste0(i_name, ".engo.ui")
                     output[[engo.plotName]] <- renderPlot({
                       print(dotplot(engo[[i]], title = paste0("Enriched ",i," Terms"))+ scale_y_discrete(labels = function(x) str_wrap(x, width = 60)))
                     })
                     
                     # Downloadable csv of selected dataset ----
                     output[[engo.downloadName]] <- downloadHandler(
                       filename = function() {
                         suffixFileName(paste0(i_name,".tables"), ".csv")
                       },
                       content = function(file) {
                         write.csv(engo[[i]]@result, file, row.names = FALSE)
                       }
                     )
                     
                     output[[engo.uiName]] <- renderUI({
                       fluidRow(plotOutput(engo.plotName), downloadButton(engo.downloadName, paste0("Download ", i," Results")))
                     })
                     
                   })
                 }
               })
               
               output$go.ui <- renderUI({
                 engo.names <- names(engo)
                 tabs <- lapply(1:length(engo.names), 
                                function(x) {
                                  tabPanel(title = engo.names[x], 
                                           uiOutput(paste0(engo.names[x],".engo.ui")),
                                           icon = icon("send")
                                  )
                                }
                 )
                 tabs$title <- "Gene Ontology Enrichment Analysis"
                 do.call(tabBox, tabs)
               })
             })