subCol <- reactive({
  subcol <- T
  for(i in dataSubmit$subcat){
    cats <- input[[paste0("subcat.",gsub(pattern = " |-",replacement = ".", i))]]
    #print(cats)
    cats[cats == ""] <- NA
    subcol <- subcol & (dataSubmit$annot[,i] %in% cats) # NA is ""
  }
#  print(length(rownames(dataSubmit$annot)[subcol]))
  return(rownames(dataSubmit$annot)[subcol])
})
