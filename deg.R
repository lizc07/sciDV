observeEvent(eventExpr = input$seurat.deg.action,
             handlerExpr = {
               select_cols <- subCol() # using subcat to subset cells
               withProgress(message = "Creating Seurat Object ...",
                            pbmc <- CreateSeuratObject(project = "iDA", raw.data = dataSubmit$expr[,select_cols,drop=F], 
                                                       meta.data = dataSubmit$annot,
                                                       normalization.method = NULL, scale.factor = 100000, 
                                                       do.scale = T, do.center = T, save.raw = F))
               
               withProgress(message = "Calculating DEGs ...",{
                 tmp <- try({
                   ### diff
                   pbmc <- SetAllIdent(pbmc, id = input$seurat.deg.cat)
                   pbmc.markers <- FindAllMarkers(pbmc, only.pos = F, min.pct = input$seurat.deg.pct, 
                                                   thresh.use = input$seurat.deg.th,
                                                   test.use = input$seurat.deg.test, latent.vars = NULL)
                   if(input$seurat.deg.test == "roc"){
                     sig.markers <- subset(pbmc.markers, myAUC >= 0.7)
                     sig.markers$p_val_adj <- -sig.markers$myAUC
                   }else{
                     sig.markers <- subset(pbmc.markers, p_val_adj <= 0.01)
                   }
                   sig.markers %>% group_by(cluster) %>% top_n(input$seurat.deg.num, power) -> topn
                 }) 
               })
               
               output$seurat.deg.heatmap <- renderPlot({
                 DoHeatmap(pbmc, genes.use = topn$gene,#col.low = "darkblue", col.mid = "white",col.high = "red",
                           use.scaled = T, slim.col.label = TRUE, remove.key = TRUE, cex.col = 1)
               })
               output$seurat.deg.heatmap.u <- renderPlot({
                 pheatmap(pbmc@data[topn$gene, order(pbmc@meta.data[,input$seurat.deg.cat])], 
                          cluster_rows = F, cluster_cols = F, annotation_col = pbmc@meta.data[,c(input$seurat.deg.cat), drop=F],
                          annotation_colors = dataSubmit$color[input$seurat.deg.cat],
                          show_colnames = F, color = pal_heatmap2, border_color = NA)
               })
               output$showDEGPlot.heatmap <- renderUI({
                 tabBox(width = 12,
                        tabPanel(title = "Expression Heatmap",
                                 plotOutput("seurat.deg.heatmap.u", height = input$tabDEG.heatmap.height)),
                        tabPanel(title = "Scaled Expression Heatmap",
                                 plotOutput("seurat.deg.heatmap", height = input$tabDEG.heatmap.height))
                 )
               })
               
               # Downloadable csv of selected dataset ----
               output$seurat.deg.download <- downloadHandler(
                 filename = function() {
                   suffixFileName(paste0(input$seurat.deg.cat, ".deg."), ".csv")
                 },
                 content = function(file) {
                   write.csv(pbmc.markers, file, row.names = FALSE)
                 }
               )
             })