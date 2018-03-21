# sciDV
##### **interactive Data Visualization for single cell RNA-Seq**  
##### 1. See `0.install.R` to install dependencies.
##### 2. Modify `global.R`, and make sure the `org` is appropriate.
##### 3. Build a customized `Rdata` file for sciDV, which should meet the following requirements.
>*Saved Rdata file should include the following variables*   
*Notice: Variable names should not be changed*

- **`data_expr`**
```R
# data.frame(),
# of which rownames are gene names 
# and colnames are sample names
> data_expr[1:5,1:5]
                 EC #01   EC #02   EC #03   EC #04 EC #05
0610005C13Rik 0.0000000 0.000000 0.000000 0.000000      0
0610007P14Rik 0.7845291 4.861415 0.000000 0.000000      0
0610009B22Rik 0.0000000 0.000000 0.000000 0.000000      0
0610009L18Rik 0.0000000 0.000000 0.000000 0.000000      0
0610009O20Rik 0.0000000 5.635482 2.344806 4.110188      0
```

- **`data_annot`**
```R
# data.frame(), 
# of which rownames are sample names (identical to colnames of data_expr) 
# and colnames are annotation names 
> data_annot[1:5,]
       Type Sample tSNE_1 tSNE_2
EC #01   EC     01    0.1    0.1
EC #02   EC     02    0.2    0.2
EC #03   EC     03    0.3    0.2
EC #04   EC     04    0.4    0.2
EC #05   EC     05    0.5    0.2
```

- **`data_color`**
```R
# Alt-1
# list(), 
# of which names should be among annotation names of data_annot; 
# each sublist should be a vector named by categories of 
# cooresponding annotation
> data_color
$Type
               EC        T1 pre-HSC
     "DarkOrange"             "red"          
# Alt-2
# vector() including sufficient color-names, RGBs.
# such as c("black","blue")
# these distinct colors will be applied to categories of each annotation in data_subcat .
```

- **`data_subcat`**
```R
# vector(), 
# including annotation names used for subsetting samples
> data_subcat
[1] "Type"
```

- **`data_coord`**
```R
# list(), 
# of which names should be among both annotation names of 
# Coordinate-Type-columns in data_annot;
# each sublist should be a vector named by annotation names
# used for 2D-dimension-plot.
> data_coord
$tSNE
[1] "tSNE_1" "tSNE_2"
```
