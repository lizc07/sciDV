options(stringsAsFactors = F)
rm(list = ls())
library(shiny)
library(shinydashboard)
library(DT)

options(shiny.maxRequestSize=1000*1024^2)

library(Seurat)
library(Rtsne)
library(igraph)
library(dplyr)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(gridExtra)

org = "mmu" # "mmu" for mouse, "hsa" for human
library(clusterProfiler)
library(org.Mm.eg.db)
require(stringr)
#library(org.Hs.eg.db)

dataSubmit <- NULL

pal_dolphin <- c("#FF00AE", "#A020F1", "#000000", "#0403E5", "#FF8C01", "#8B0101", 
                 "#007502", "#FE0000", "#FFFF01", "#FF99CB", "#4A95FB", "#61FE69",
                 "#9A7A01", "#017F8B", "#05FDFF")
pal_heatmap1 <- gplots::colorpanel(100,"blue","white","red")
pal_heatmap2 <- gplots::colorpanel(100,"darkblue","white","red")

is.color <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}

myFeaturePlot <- function(object, features.plot, nrow = NULL, ncol = NULL, ...){
  require(ggplot2)
  require(gridExtra)
  ggData <- as.data.frame(cbind(object@dr$tsne@cell.embeddings,FetchData(object, features.plot)))
  colnames(ggData) <- c(colnames(object@dr$tsne@cell.embeddings),gsub("-",".",features.plot))
  # print(feature.tmp)
  # ggData[,feature.tmp] <- t(object@data[feature.tmp,])
  ggl <- lapply(features.plot, function(feature){
    ggplot(ggData) + geom_point(mapping = aes_string(x = "tSNE_1", y = "tSNE_2", color = gsub("-",".",feature)), size = 1) + 
      scale_color_gradientn(colours = c("grey","yellow","red")) + 
      theme(legend.title = element_blank(),axis.title = element_blank()) + ggtitle(feature) 
  })
  grid.arrange(grobs = ggl, nrow= nrow,ncol = ncol)
}


textCapitalize <- function(x, org = "mmu"){
  if(!is.character(x)){
    return(x)
  }else{
    if(org == "mmu"){
      tmp <- strsplit(as.character(x), split = "")[[1]]
      return(paste0(toupper(tmp[1]), paste0(tmp[-1], collapse = "")))
    }else if(org == "hsa"){
      return(toupper(x))
    }else{
      return(x)
    }
  }
}
textInput2PCs <- function(x){
  tmp <- textInput2genes(x)
  x <- NULL
  for(i in tmp){
    if(grepl("-",i)){
      range = as.numeric(strsplit(i, "-")[[1]])
      x <- c(x, seq(range[1], range[2]))
    }else{
      x <- c(x, as.numeric(i))
    }
  }
}
textInput2genes <- function(x, org = "mmu"){
  x <- gsub("'", "", x)
  x <- gsub('"', "", x)
  x <- setdiff(unique(unlist(strsplit(x, split = ' |\n|,|\t'))), "")
  sapply(x, textCapitalize, org = org)
}

suffixFileName <- function(x, extension){
  paste0(x, format(Sys.time(),'_%Y%m%d_%H%M%S'), extension)
}

return.numeric <- function(x, ret = NULL){
  ifelse(is.numeric(x), x, ret)
}

checkdata.showModal <- function(x, label = NULL, message = NULL, show = T){
  #print(dim(x))
  text.null <- paste0("Please Load ", label)
  text.invalid <- paste0("Please check Your Data! It seems like 0 row/column in the ", label)
  
  if(!is.null(message)){
    text.null <- text.invalid <- message
  }
  
  if(is.null(x)){
    if(show) showModal(modalDialog(text.null,title = "Check Data", easyClose = T, fade = T))
    return(F)
  }else if(min(dim(x)) == 0){
    if(show) showModal(modalDialog(text.invalid, title = "Check Data", easyClose = T, fade = T))
    return(F)
  }else{
    return(T)
  }
}

## for guessCategory
guessNumeric <- function(x, exclude = c(NA, "NA", NULL, "NULL", ""," ","null")){
  return(!anyNA(as.numeric(setdiff(x, exclude))))
}
fillNA <- function(x, asNA = c(NA, "NA", NULL, "NULL", ""," ","null")){
  x[x %in% asNA] <- NA
  return(x)
}
formatAnnot <- function(annot, na.fill = T){
  numericCol <- apply(annot, 2, guessNumeric)
  if(na.fill){
    annot <- fillNA(annot, asNA = c(""," "))
  }
  if(any(numericCol)){
    annot[,numericCol] <- apply(annot[,numericCol], 2, as.numeric)
  }
  return(annot)
}
getSubcatUI <- function(annot){
  subcatUI <- lapply(colnames(annot), function(x){
                    cats <- as.character(sort(unique(annot[,x]), na.last = T))
                    checkboxGroupInput(inputId = paste0("subcat.",gsub(pattern = " |-",replacement = ".",x)), label = x, 
                                       choices = cats, selected = cats, inline = T, width = "100%")
                  })
  return(subcatUI)
}

specify_geom_point <- function(data, one_gene = NULL,aes.x, aes.y, aes.color, size = 1, scale_color_manual = NULL, coord_fixed.ratio = NULL){
  p <- ggplot() + 
    theme_classic() + theme(plot.title = element_text(hjust = 0.5,face = "bold"))
  if((!is.null(coord_fixed.ratio)) && is.finite(coord_fixed.ratio)){p <- p + coord_fixed(ratio = coord_fixed.ratio)}
  
  if(!is.null(one_gene) && !is.numeric(data[,one_gene])){ 
    one_gene <- NULL
    aes.color <- one_gene
  }
  
  if((!is.null(one_gene)) && (!is.na(one_gene))){
    mid.tmp <- quantile(data[,one_gene],c(0,1), na.rm = T)
    data.tmp <- data[,c(aes.x, aes.y, one_gene)]
    colnames(data.tmp) <- c(aes.x, aes.y, "SELIGENE")
    p <- p + geom_point(data = data.tmp, mapping = aes_string(x = aes.x, y = aes.y,color = "SELIGENE"), size = size) + labs(color = "SELIGENE") +
      ggtitle(paste0("Expression of ",one_gene)) + scale_color_gradient2(low = "grey", mid = "yellow",high = "red",midpoint = mean(mid.tmp))#scale_color_gradientn(colours = c("grey","yellow","red","black"), values = c(0,0.1,0.6,1))
  }else{
    p <- p + geom_point(data = data, mapping = aes_string(x = aes.x, y = aes.y,color = aes.color), size = size) + 
      ggtitle(paste0("Cell ", aes.color)) + guides(color = guide_legend(order = 2))
  }
  if((!is.null(scale_color_manual)) && (!is.na(scale_color_manual))){
    print(scale_color_manual)
    p <- p + scale_color_manual(values = scale_color_manual)
  }
  return(p)
}

