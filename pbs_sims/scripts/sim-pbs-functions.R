select_sims = function(path){
  # read data
  data = read.table(path, sep = " ", header = TRUE)
  # remove NAs
  data = data[!is.na(data$pbs),]
  dataFrame = c()
  # create data frame
  for (j in 1:20) {
    dataFrame = rbind(dataFrame, data[data$rep == j,])
  }
  return(dataFrame)
}
