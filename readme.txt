#################################################################
# sciDV - iteractive Data Visualization for single cell RNA-Seq
#################################################################
# Rdata file should include the following datasets
# Notice: Variable names should not be changed
## data_expr    
- data.frame()
  of which rownames are gene names 
  and colnames are sample names
## data_annot   
- data.frame(), 
  of which rownames are sample names (identical to colnames of data_expr) 
  and colnames are annotation names 
## data_color   
- list(), 
  of which names should be among annotation names of data_annot; 
  each sublist should be a vector named by categories of 
  cooresponding annotation
## data_subcat  
- vector(), 
  including annotation names used for subsetting samples
## data_coord   
- list(), 
  of which names should be among both annotation names of 
  Coordinate-Type-columns in data_annot.
#################################################################
