#############################################
# Prepare Data-                             #
# Divide into Subsets by the Two Factors     #
#############################################


#=================================================================================================================

# input:   data frame of the form (time,status, factor1, factor2)

# output:  list of :
#             - whole data set with a new variable cell, which indicates in which cell an observation is
#             - subset for each cell  (time,event,factor1, factor2, cell)

#=================================================================================================================



prepare <- function(data){

  data <- na.omit(data)
  # save the levels of the two factor variables and their length and number of observations
  levels3 <- levels(data[,3])
  a       <- length(levels3)
  levels4 <- levels(data[,4])
  b       <- length(levels4)
  n       <- dim(data)[1]

  # a vector to save which cell the observation belongs to
  cell       <- rep(0,n)
  cell_count <- 1
  data       <- cbind(data, cell)

  # assign the cell level to each individual

  for ( j in 1:a ) {

    # look trough all levels of A
    for ( k in 1:b ) {

      # look through all levels of B
      for (i in 1:n ) {

        if (data[i, 3] == levels3[j] & data[i, 4] == levels4[k]) {
          data$cell[i] <- cell_count
        }

      }
      cell_count<-cell_count+1
    }
  }

  data$cell <- as.factor(data$cell)
  result    <- list(data)

  for (i in 1:max(as.numeric(data$cell)) ) {
    datahelp <- data[data$cell == i, ]
    result   <- c(result, list(datahelp))
  }

  return (result)
}

