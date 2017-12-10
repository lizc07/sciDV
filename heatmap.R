selectRC <- reactive({
  if(input$tabHeatmap.select_genes == "")
    select_rows <- sample(rownames(dataSubmit$expr),size = 50,replace = F)
  else{
    select_rows <- textInput2genes(input$tabHeatmap.select_genes, org)
  }

  select_cols <- subCol()
  
  # check expression for clustering
  rows_valid <- intersect(select_rows,rownames(dataSubmit$expr))
  cols_valid <- intersect(select_cols,colnames(dataSubmit$expr))
  rows_invalid <- setdiff(select_rows,rownames(dataSubmit$expr))
  cols_invalid <- setdiff(select_cols,colnames(dataSubmit$expr))
  
  validMat <- dataSubmit$expr[rows_valid, cols_valid]
  
  validMat.rowSD <- apply(validMat, 1, sd)
  validMat.rowMax <- apply(validMat,1,max)
  rows_sd0 <- rows_valid[validMat.rowSD == 0]
  rows_zero <- rows_valid[validMat.rowMax == 0 & validMat.rowSD == 0]
  rows_sdna <- rows_valid[is.na(validMat.rowSD)]
  
  validMat.colSD <- apply(validMat, 2, sd)
  validMat.colMax <- apply(validMat,2,max)
  cols_sd0 <- cols_valid[validMat.colSD == 0]
  cols_zero <- cols_valid[validMat.colMax == 0 & validMat.colSD == 0]
  cols_sdna <- cols_valid[is.na(validMat.colSD)]
  
  return(list(select_rows, select_cols, rows_valid, cols_valid, rows_invalid, cols_invalid,
                   rows_sd0, rows_zero, rows_sdna, cols_sd0, cols_zero, cols_sdna))
})

heatMat <- reactive({
  rows_heat <- selectRC()[[3]]
  cols_heat <- selectRC()[[4]]
  
  if(input$tabHeatmap.rowClust)  rows_heat <- Reduce(setdiff,list(rows_heat, selectRC()[[7]], selectRC()[[9]]))
  if(input$tabHeatmap.colClust)  cols_heat <- Reduce(setdiff,list(cols_heat, selectRC()[[10]],selectRC()[[12]]))
  
  selectMat <- dataSubmit$expr[rows_heat, cols_heat]
  heatMat <- selectMat
  #print(dim(heatMat))
  # if(!exprScale()) heatMat <- 2^selectMat - 1
  return(heatMat)
})


output$heatmap <- renderPlot({
  #breaks <- range(heatMat) 
  heat_rows <- rownames(heatMat())
  heat_cols <- colnames(heatMat())
  sort_rows <- switch (input$tabHeatmap.sortRow,
                       "original" = heat_rows,
                       "alphabet" = sort(heat_rows))
  
  # sort_cols
  #print(heat_cols)
  heat_annot <- cbind(dataSubmit$annot[heat_cols,,drop = F], 
                      original = 1:length(heat_cols))
  sort_annot <- input$tabHeatmap.sortCol
  if(is.null(sort_annot) || sort_annot == ""){
    sort_annot <- "original"
  }
  if("Gene" %in% input$tabHeatmap.sortCol){
    #
    sort_cols.genes <- intersect(textInput2genes(input$tabHeatmap.sortColByGene, org), rownames(dataSubmit$expr))
    if(!is.null(sort_cols.genes)){
      heat_annot <- cbind(heat_annot, t(dataSubmit$expr[sort_cols.genes, rownames(heat_annot), drop = F]))
    }
    ind_Gene <- match("Gene", input$tabHeatmap.sortCol)
    if(ind_Gene > 1){
      if(length(input$tabHeatmap.sortCol) > (ind_Gene + 1) )
        sort_annot <- c(input$tabHeatmap.sortCol[1:(ind_Gene-1)], sort_cols.genes, input$tabHeatmap.sortCol[(ind_Gene+1):length(input$tabHeatmap.sortCol)])
      else
        sort_annot <- c(input$tabHeatmap.sortCol[1:(ind_Gene-1)], sort_cols.genes)
    }else{
      if(length(input$tabHeatmap.sortCol) > (ind_Gene + 1) )
        sort_annot <- c(sort_cols.genes, input$tabHeatmap.sortCol[(ind_Gene+1):length(input$tabHeatmap.sortCol)])
      else
        sort_annot <- c(sort_cols.genes)
    }
  }
  sort_cols <- rownames(heat_annot)[do.call(order,as.data.frame(heat_annot[,sort_annot, drop = F]))]
  
  phData <- heatMat()[sort_rows,sort_cols]
  clustData <- heatMat()[sort_rows,sort_cols]
  if(input$tabHeatmap.rowClust){
    cluster_distance_row <- switch(input$tabHeatmap.rowClust.dist,
                                   "cor.pearson" = as.dist((1-cor(t(clustData),use = "pairwise.complete.obs"))/2),
                                   "cor.spearman"= as.dist((1-cor(t(clustData),use = "pairwise.complete.obs",method = "spearman"))/2),
                                   dist(clustData,method = input$tabHeatmap.rowClust.dist))
  }
  
  if(input$tabHeatmap.colClust){
    cluster_distance_col <- switch(input$tabHeatmap.rowClust.dist,
                                   "cor.pearson" = as.dist((1-cor(clustData,use = "pairwise.complete.obs"))/2),
                                   "cor.spearman"= as.dist((1-cor(clustData,use = "pairwise.complete.obs",method = "spearman"))/2),
                                   dist(t(clustData),method = input$tabHeatmap.rowClust.dist))
  }
  
  if(input$tabHeatmap.rowScale != "non-scale"){
    phData <- t(scale(t(phData),center = T, scale = T))
    if(gsub("scale","", input$tabHeatmap.rowScale) != ""){
      cap <- as.numeric(gsub("scale","", input$tabHeatmap.rowScale))
      phData[phData > cap] <- cap
      phData[phData < -cap] <- -cap
    }
  }

  if(is.null(input$tabHeatmap.annotCol))
    phAnnot_col <- NA
  else
    phAnnot_col <- dataSubmit$annot[,input$tabHeatmap.annotCol,drop = F]
  
  ph_show_colnames <- ifelse(ncol(phData) > 100, F, T)
  pheatmap(phData,cluster_rows = input$tabHeatmap.rowClust ,cluster_cols = input$tabHeatmap.colClust, border_color = F, 
           clustering_distance_rows = cluster_distance_row, clustering_distance_cols= cluster_distance_col, 
           clustering_method = input$tabHeatmap.clust.linkage, annotation_col = phAnnot_col,annotation_colors = dataSubmit$color, 
           treeheight_row = input$tabHeatmap.treeheight, treeheight_col = input$tabHeatmap.treeheight, show_colnames = ph_show_colnames)
})
output$heatmap.ui <- renderUI({
  plotOutput("heatmap",height = input$tabHeatmap.height)
})

output$heatmap.info.tab <- renderTable({
  infoRC <- selectRC()
  nrowzero <- length(infoRC[[8]])
  nrowcc <- length(infoRC[[7]]) - nrowzero
  nrowsdna <- length(infoRC[[9]])
  nrowvariable <- length(infoRC[[3]]) - nrowcc - nrowzero
  nrowinvalid <- length(infoRC[[5]])
  nrowinput <- length(infoRC[[1]])
  
  ncolzero <- length(infoRC[[11]])
  ncolcc <- length(infoRC[[10]]) - ncolzero
  ncolsdna <- length(infoRC[[12]])
  ncolvariable <- length(infoRC[[4]]) - ncolcc - ncolzero
  ncolinvalid <- length(infoRC[[6]])
  ncolinput <- length(infoRC[[2]])
  
  infoTab <- matrix(c(nrowzero, nrowcc, nrowsdna, nrowvariable, nrowinvalid, nrowinput, ncolzero, ncolcc, ncolsdna, ncolvariable, ncolinvalid, ncolinput),nrow = 2,byrow = T)
  rownames(infoTab) <- c("Gene","Sample")
  colnames(infoTab) <- c("Zero constant","Non-zero constant", "SD N/A","Variable", " Invalid", "Total")
  
  return(infoTab)
}, striped = T, hover = T, bordered = T, spacing = "l", align = "c", rownames = T)
output$heatmap.info.text <- renderUI({
  #text1 <- strong(paste("Heatmap for ",heatMat.dim[1]," genes and ",heatMat.dim[2]," samples!",sep = ""))
  text1 <- paste(span("Variational expression genes:",style = "color:green"),paste(Reduce(setdiff,list(selectRC()[[3]],selectRC()[[7]],selectRC()[[9]])),collapse = ", "))
  text2 <- ""
  exgenes <- union(selectRC()[[7]],selectRC()[[9]])
  if(length(exgenes) > 0) text2 <- paste(span("Remove genes when Gene-clustering:",style = "color:red"),paste(exgenes,collapse = ", "))
  
  text3 <- ""
  if(length(selectRC()[[5]]) > 0) text3 <- paste(span("INVALID genes:",style = "color:grey"),paste(selectRC()[[5]],collapse = ", "))
  
  text4 <- ""
  exsamples <- union(selectRC()[[10]],selectRC()[[12]])
  if(length(exsamples) > 0 & length(exsamples) <= 10) text4 <- paste(span("Remove samples when Sample-clustering:",style = "color:red"),paste(exsamples,collapse = ", "))
  else if(length(exsamples) > 10) text4 <- paste(span("Remove samples when Sample-clustering:",style = "color:red"),paste(length(exsamples),"samples"))
  text5 <- ""
  if(length(selectRC()[[6]]) > 0) text5 <- paste(span("INVALID samples:",style = "color:grey"),paste(selectRC()[[6]],collapse = ", "))
  #text4 <- paste(span("Constant genes are:",style = "color:grey"),paste(heatMat.rownames[heatMat.rowValid],collapse = " "),sep = " ")
  HTML(paste(text1,text2,text3,text4,text5,sep = "<br/>"))
})





