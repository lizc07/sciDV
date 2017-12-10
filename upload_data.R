############
# Load Data
############
upload_data <- observeEvent(input$fileInput.action, {
  
  if(!is.null(input$fileInput.rdata)) { # if file selected
   tmp <- try(load(file = input$fileInput.rdata$datapath[1]))
    #print(pbmc)
    if(class(tmp) == "try-error"){ # if failed to load
      output$fileInput.error <- renderUI({box(renderText({tmp}),width = 8, status = "danger",solidHeader=T,
                                              title = "Error", footer = "Please Contact lizc07@vip.qq.com")
      })
      return()
    }else{ # if scucceed to load
      
      expr <- data_expr
      annot <- data_annot
      color <- data_color
      subcat <- data_subcat
      coord <- data_coord
      
      if(!is.list(color)){
        color <- list()
        for(i in data_subcat){
          tmp.cat <- sort(unique(annot[,i]))
          color[[i]] <- setNames(data_color[1:length(tmp.cat)], tmp.cat)
        }
      }
      ## cell cycle scores
      g1sGene <- sapply(c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"), textCapitalize, org = org)
      g2mGene <- sapply(c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa"), textCapitalize, org = org)
      annot$idv.g1sScore <- colMeans(expr[intersect(g1sGene, rownames(expr)), rownames(annot)], na.rm = T)
      annot$idv.g2mScore <- colMeans(expr[intersect(g2mGene, rownames(expr)), rownames(annot)], na.rm = T)
      
      dataSubmit <<- list(expr = expr, annot = annot, color = color, subcat = subcat, coord = coord)
    }
  }else{ # # if file is null
    checkdata.showModal(NULL, message = "Please load RData file...")
    return()
  }
  
  
  output$infoBox.expr <- renderInfoBox({
    if(checkdata.showModal(dataSubmit$expr, "Expression Data")){
      if(is.null(dataSubmit$annot)){
        text.submitted <- "Warning: Annotation Data is not uploaded!"
        status.submitted <- "warning"
      }else if(length(c(setdiff(rownames(dataSubmit$annot), colnames(dataSubmit$expr)), setdiff(colnames(dataSubmit$expr),rownames(dataSubmit$annot)))) != 0){
        text.submitted <- "Warning: Cell Names is not identical between Expression data and Annotation data!"
        status.submitted <- "warning"
      }else{
        text.submitted <- "Please continue to process"
        status.submitted <- "success"
      }
      output$fileInput.submitted <- renderUI({box(text.submitted, width = 12, status = status.submitted, solidHeader=T,
                                                  title = "Submitted", footer = "Produced by Zongcheng Li, using shiny")
      })
      infoBox(title = "Expression", value = paste0("Gene: ",nrow(dataSubmit$expr),", \n\rCell: ", ncol(dataSubmit$expr)), fill = F)
    }else{
      infoBox(title = "Expression", value = NULL,fill = T,subtitle="Please Load Data")
    }
  })
  outputOptions(output, "infoBox.expr", suspendWhenHidden = F, priority = 0)
  
  output$infoBox.annot <- renderInfoBox({
    if(!is.null(dataSubmit$annot)){
      checkdata.showModal(dataSubmit$annot, "Annotation Data")
      infoBox(title = "Annotation", value = paste0("Cell: ",nrow(dataSubmit$annot),
                                                   ", \n\rCategory: ", ncol(dataSubmit$annot)), fill = F)
    }else{
      infoBox(title = "Annotation", value = NULL,fill = T,subtitle="Please Load Data")
    }
  })
  outputOptions(output, "infoBox.annot", suspendWhenHidden = F, priority = 0)
  
  
  # menuItem Sub cells
  output$showSubcat <- renderUI({

    if(checkdata.showModal(dataSubmit$annot[,dataSubmit$subcat, drop = F], show = F)){
      return(getSubcatUI(dataSubmit$annot[,dataSubmit$subcat, drop = F]))
    }
  })
  outputOptions(output, "showSubcat", suspendWhenHidden = F, priority = -1)
  
  
  
  
  # menuItem Heatmap
  output$menuHeatmap.displayPara <- renderUI({
    fluidRow(
      radioButtons("tabHeatmap.rowScale", label = "Expression level scaled ?", choices = list("non-scale","scale","scale3", "scale5"), selected = "non-scale",inline = T),
      
      selectInput("tabHeatmap.annotCol", label = ("Samples' annotation"), choices = colnames(dataSubmit$annot),
                  selected = dataSubmit$subcat, multiple = T, width = "100%"),
      
      radioButtons("tabHeatmap.sortRow", label = ("Order genes by"),
                   choices = list("original", "alphabet"),inline = T),
      
      selectInput("tabHeatmap.sortCol", label = "Order samples by", choices = c("original",colnames(dataSubmit$annot),"Gene"), 
                  selected = dataSubmit$subcat, multiple = T, width = "100%"),
      
      selectInput("tabHeatmap.sortColByGene", label = "Order samples when by Gene",multiple = T, 
                  choices = rownames(dataSubmit$expr), selected = ifelse(org == "hsa","RUNX1","Runx1"))
    )
  })
  outputOptions(output, "menuHeatmap.displayPara", suspendWhenHidden = F, priority = -2)
  
 
  # menuItem DEG 
  output$menuDEG.displayPara <- renderUI({
    selectInput(inputId = "seurat.deg.cat", label = "Find Markers Among:", 
                choices = colnames(dataSubmit$annot), selected = dataSubmit$subcat[1])  
  })
  outputOptions(output, "menuDEG.displayPara", suspendWhenHidden = F, priority = -2)
  
  ## menuItem DimPlot
  output$menuDimPlot.displayPara <- renderUI({
    fluidRow(
      selectInput("dimplot.gene", "Input a gene name", choices = c(rownames(dataSubmit$expr), colnames(dataSubmit$annot)), selected = ifelse(org == "hsa","RUNX1","Runx1"), multiple = F),
      selectInput("dimplot.cat", "Input a category", choices = c(colnames(dataSubmit$annot)), selected = dataSubmit$subcat[1], multiple = F),
      sliderInput("dimplot.height", label = "Graph Height", min = 0, max = 1000, value = 600),
      sliderInput("dimplot.psize", label = "Point Size", min = 0, max = 10, value = 5, step = 0.1,round = -1)
      
    )
  })
  outputOptions(output, "menuDimPlot.displayPara", suspendWhenHidden = F, priority = -2)
  
  
  # DimPlot
  if(!is.list(dataSubmit$coord) || length(dataSubmit$coord) < 1) {} # must be list
  else{
    for(i in 1:length(dataSubmit$coord)){
      local({
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
      tabName <- names(dataSubmit$coord)[i]
      tabDims <- dataSubmit$coord[[i]]
      if(length(tabDims) < 1) next # if no axis info 
      
      gene.plotName <- paste0(tabName,".gene.plot")
      gene.uiName <- paste0(tabName,".gene.plot.ui")
      cat.plotName <- paste0(tabName,".cat.plot")
      cat.uiName <- paste0(tabName,".cat.plot.ui")
      
      if(length(tabDims) == 2){
        output[[gene.plotName]] <- renderPlot({
          specify_geom_point(data = cbind(dataSubmit$annot, t(dataSubmit$expr[,rownames(dataSubmit$annot),drop =F])), 
                             one_gene = input$dimplot.gene, aes.x = tabDims[1], aes.y = tabDims[2], aes.color = "SELIGENE", 
                             size = input$dimplot.psize, coord_fixed.ratio = NULL)
          
        })
        output[[gene.uiName]] <- renderUI({plotOutput(gene.plotName, height = input$dimplot.height)})
        
        output[[cat.plotName]] <- renderPlot({
          specify_geom_point(data = dataSubmit$annot, 
                             one_gene = NULL, aes.x = tabDims[1], aes.y = tabDims[2], aes.color = input$dimplot.cat, 
                             size = input$dimplot.psize, coord_fixed.ratio = NULL, 
                             scale_color_manual = if(input$dimplot.cat %in% names(dataSubmit$color)) dataSubmit$color[[input$dimplot.cat]] else NULL) 
          
        })
        output[[cat.uiName]] <- renderUI({plotOutput(cat.plotName, height = input$dimplot.height)})
      }
      })
    }
  }
  
  
  ## menuItem VlnPlot
  output$menuVlnPlot.displayPara <- renderUI({
    fluidRow(
      selectInput("vlnplot.gene", "Input a gene name", choices = c(rownames(dataSubmit$expr), colnames(dataSubmit$annot)), selected = ifelse(org == "hsa","RUNX1","Runx1"), multiple = F),
      selectInput("vlnplot.cat", "Input a category", choices = c(colnames(dataSubmit$annot)), selected = dataSubmit$subcat[1], multiple = F),
      sliderInput("vlnplot.width", label = "Graph Width", min = 0, max = 1600, value = 800,round = T),
      sliderInput("vlnplot.height", label = "Graph Height", min = 0, max = 1000, value = 500,round = T),
      sliderInput("vlnplot.psize", label = "Point Size", min = 0, max = 10, value = 5, step = 0.1,round = -1)
    )
  })
  outputOptions(output, "menuVlnPlot.displayPara", suspendWhenHidden = F, priority = -2)
 
  
  ## menuItem BarPlot
  output$menuBarPlot.displayPara <- renderUI({
    fluidRow(
      selectInput("barplot.gene", "Input a gene name", choices = c(rownames(dataSubmit$expr), colnames(dataSubmit$annot)), selected = ifelse(org == "hsa","RUNX1","Runx1"), multiple = F),
      selectInput("barplot.cat", "Select a category", choices = c(colnames(dataSubmit$annot)), selected = dataSubmit$subcat[1], multiple = F),
      sliderInput("barplot.width", label = "Graph Width", min = 0, max = 1600, value = 800),
      sliderInput("barplot.height", label = "Graph Height", min = 0, max = 1000, value = 340)
    )
  })
  outputOptions(output, "menuBarPlot.displayPara", suspendWhenHidden = F, priority = -2)
  
  ## menuItem CellCyclePlot
  output$menuCellCyclePlot.displayPara <- renderUI({
    fluidRow(
      #sliderInput("ccplot.width", label = "Graph Width", min = 0, max = 1600, value = 800),
      selectInput("ccplot.cat", "Select a category", choices = c(colnames(dataSubmit$annot)), selected = dataSubmit$subcat[1], multiple = F),
      sliderInput("ccplot.height", label = "Graph Height", min = 0, max = 1000, value = 500),
      sliderInput("ccplot.psize", label = "Point Size", min = 0, max = 10, value = 5, step = 0.1,round = -1),
      sliderInput("ccplot.th", label = "Score Threshold", min = 0, max = 10, value = 2, step = 0.1,round = -1)
    )
  })
  outputOptions(output, "menuCellCyclePlot.displayPara", suspendWhenHidden = F, priority = -2)
})


