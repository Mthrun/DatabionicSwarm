RobustNormalization = function (Data,Centered=FALSE,Capped=FALSE,na.rm=TRUE) 
{
  if (is.vector(Data)) {
    quants = quantile(Data, c(0.01, 0.5, 0.99),na.rm = na.rm)
    minX = quants[1]
    maxX = quants[3]
    Denom=maxX - minX
    if(Denom==0) Denom=1 
    Data=(Data - minX)/Denom
    #median(abs(Data-median(Data,na.rm=T)),na.rm=T)
    if(Centered){
      Data=Data-median(Data,na.rm=na.rm)
      if(Capped){
        maxX=1
        minX=-1
        Data[Data>maxX]=maxX
        Data[Data<minX]=minX
      }
    }else{
      if(Capped){
        quants = quantile(Data, c(0.01, 0.5, 0.99),na.rm = na.rm)
        minX = quants[1]
        maxX = quants[3]
        Data[Data>maxX]=maxX
        Data[Data<minX]=minX
      }
    }
    return(Data)
  }
  else if (is.matrix(Data)) {
    cols = ncol(Data)
    xtrans = Data
    for (i in 1:cols) {
      xtrans[, i] = RobustNormalization(as.vector(Data[, i]),Centered,Capped,na.rm)
    }
    return(xtrans)
  }
  else {
    tryCatch({
      warning("Data is not a vector or a matrix. Trying as.matrix")
      return(RobustNormalization(as.matrix(Data),Centered,Capped,na.rm))
    }, error = function(e) {
      stop("It did not work")
    })
  }
}