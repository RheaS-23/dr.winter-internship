library(dplyr)
library(Seurat)

### Understanding what a "data frame" is

  ## It's similar to a matrix: the data is organized into rows and columns
  
    # Filling a matrix with random values
    example.matrix <- matrix(data = rep(c("a", "b", "c"), times = 3), ncol = 3)
  
    example.matrix
    
  ## But unlike a matrix, data frames allow you to call columns by the name of the variable
  
    # Converting matrix into data frame
    example.dataframe <- example.matrix %>% 
      as.data.frame()
    
    example.dataframe
    
  ## By default, a dataframe will name the columns "V1", "V2", etc. But you can change those names to whatever you want
  
    # Setting new column names
    colnames(example.dataframe) <- c("column 1", "column 2", "column 3")
    colnames
    
  ## So now you can pull columns by name instead of having to figure out its index
    example.dataframe["column 2"]
  
  ## To calculate percentiles, we are going to use some functions from the "dplyr" package. 
  ## Since "dplyr" only takes data frames as input, we'll need to convert our data into a data frame.
    
### Calculating percentile for nCount_RNA
    
  ## First we want to pull the data that we want out of our Seurat object
    
    # Pulling nCount_RNA
    nCount <- training_merged_obj@meta.data$nCount_RNA
    
  ## Next we need to convert it to a data frame to be able to use our dplyr functions
    
    # Convert to data frame
    nCount <- as.data.frame(nCount)
    
  ## Now we have a data frame with one variable called "nCount"
    
  ## We can go ahead and calculate percentiles. 
  ## The ntile() does a couple things:
    # 1) Sorts the vector from smallest to largest
    # 2) Groups the data into n bins of equal size
    # 3) Returns the name of the "bin" that each original item in your list belongs to
    
  ## We will choose n = 100 bins, so that each bin will represent 1% of the all the cells in the data set
    
  ## Summarize is a helper function that applies another function to a whole vector.
  ## So running the following lines will return a vector the same length as our original nCount_RNA vector
  ## but instead of nCount values, the values are the "percentile" that each cell belongs to.
  ## So the following lines transform each nCount value into a percentile label, 
  ## and save that data under the variable "percentile" in the data frame "nCount.percent"

    # Bin nCount data into percentiles
    nCount.percent <- nCount %>%
      summarize(percentile = ntile(nCount, n = 100))

  ## Sanity check
  ## If you started with 33000 cells, and you grouped them into 100 bins of equal size, 
  ## then you should expect that there 330 cells in each bin.
  ## To make sure that this code did what we expected, lets check.
    
  ## which() will give you the index values of all the items in your vector that match your condition 
  ## (i.e. that have a value of "1")
    
  ## length() will tell us how long that list of indexes is
    
  ## Try running this next line with different values between 1-100 
    
    # Calculate many cells in the data set are in the 1% bin
    which(nCount.percent == "1") %>% length()


  ## Since the order of our percentile calculations is the same order as the original data 
  ## (i.e. the value in the 4th spot still represents to the 4th cell in the dataset),
  ## we can just input our percentile vector as a new metadata variable in our object
  
  # Insert percentile calculations as metadata variables  
  training_merged_obj@meta.data[["nCount.percentile"]] <- nCount.percent[["percentile"]]
  
  ## Note: we specified nCount.percent[["percentile"]] because we only want the raw data in that one column
  ## if we did "<- nCount.percent", then the whole dataframe (including row and column names) 
  ## would be inserted into the object
  
### Choosing a percentile threshold

  ## Lets day I want to find out which cells will be eliminated if I filtered out the bottom 10% by nCount
  ## WhichCells() gives you a list of cell names that meet a certain condition
  
    # Get names of cells in bottom 10% by nCount
    bottom10percent <- WhichCells(training_merged_obj, expression = nCount.percentile <= 10)
  
  ## You could also make a new binary variable to store this information
  ## Skip this section if you only want to use the prior method to get a list of cells
    
    # Binary variable for bottom 10% of nCount
    training_merged_obj@meta.data[["nCount.below10"]] <- ifelse(training_merged_obj@meta.data[["nCount.percentile"]] >= 10, "above.10", "below.10")

    # Set active ident of cells to binary variable so you can filter by labels
    training_merged_obj <- SetIdent(training_merged_obj, value = "nCount.below10")

    # Which cells have the "nCount.below10" label
    bottom10percent <- WhichCells(training_merged_obj, idents = "nCount.below10")
    
  ## Once you have a list of cells that would be filtered out, you can use intersect() to see the overlap between two lists
  ## i.e. if you wanted to compare which cells get filtered out using your threshold vs. using an algorithm
  
  # Get list of items that are present in both list.of.cells.1 and list.of.cells.2
  intersect(list.of.cells1, list.of.cells2)
  
  # Find how many cells overlap between cell lists
  intersect(list.of.cells1, list.of.cells2) %>% length()
  


