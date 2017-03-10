transform_connectomes <- function(connectomes, transform_type){
  
  #check normality of connectome data
  normtest_no_tr <- c()
  for (i in 1:length(connectomes[1,])){
    normtest_no_tr[i] <- -log(shapiro.test(connectomes[,i])$p)
  }
  
  #transform using log or ranknormal
  if (transform_type == "log") {
    minc <- min(connectomes)
    if (minc < 0) {
      connectomes <- connectomes - minc*1.5
    }
    connectomes <- log(connectomes)
  } else if (transform_type == "ranknormal") {
    for (i in 1:length(connectomes[1,])){
      connectomes[,i] <- rntransform(connectomes[,i])
    }
  }

  #re-check normality and compare
  normtest_tr <- c()
  for (i in 1:length(connectomes[1,])){
    normtest_tr[i] <- -log(shapiro.test(connectomes[,i])$p)
  }
  
  return(list(connectomes,normtest_no_tr,normtest_tr))
}