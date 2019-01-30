RobustNormalization = function (Data,Centered=FALSE,Capped=FALSE,na.rm=TRUE,WithBackTransformation=FALSE) 
{
  if(is.data.frame(Data)){
    warning('Matrix is expected but data.frame is given. Calling as.matrix().')
    Data=as.matrix(Data)
  }
if(isTRUE(na.rm)){
  # if(!is.data.frame(Data)){
    Data[!is.finite(Data)]=NaN #quantile does not accept inf,-inf
  # }else{
  #   ind=do.call(cbind, lapply(mtcars, is.finite))
  #   Data[ind]=NaN
  # }
   
}
  center=NULL
  Denom=NULL
  minX=NULL
  maxX=NULL
  if (is.vector(Data)) {
    quants = quantile(Data, c(0.01, 0.5, 0.99),na.rm = na.rm)
    minX = quants[1]
    maxX = quants[3]
    Denom=maxX - minX
    if(Denom==0) Denom=1 
    Data=(Data - minX)/Denom
    #median(abs(Data-median(Data,na.rm=T)),na.rm=T)
    if(Centered){
      center=median(Data,na.rm=na.rm)
      Data=Data-center
      if(Capped){
        maxX2=1
        minX2=-1
        Data[Data>maxX2]=maxX2
        Data[Data<minX2]=minX2
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
    if(WithBackTransformation)
      return(list(TransformedData=Data,MinX=minX,MaxX=maxX,Denom=Denom,Center=center))#
    else
      return(Data)
  }
  else if (is.matrix(Data)) {
    if(WithBackTransformation){
      cols = ncol(Data)
      xtrans = list()
      DataOut=c()
      minX=c()
      maxX=c()
      Denom=c()
      center=c()
      for (i in 1:cols) {
        xtrans = RobustNormalization(Data = as.vector(Data[, i]),Centered = Centered,Capped = Capped,na.rm = na.rm,WithBackTransformation=WithBackTransformation)
        DataOut=cbind(DataOut,as.matrix(xtrans$TransformedData))
        minX=cbind(minX,xtrans$MinX)
        maxX=cbind(maxX,xtrans$MaxX)
        Denom=cbind(Denom,xtrans$Denom)
        center=cbind(center,xtrans$Center)
      }
      
      return(list(TransformedData=DataOut,MinX=minX,MaxX=maxX,Denom=Denom,Center=center))
    }else{
    cols = ncol(Data)
    xtrans = Data
    for (i in 1:cols) {
      xtrans[, i] = RobustNormalization(as.vector(Data[, i]),Centered,Capped,na.rm,WithBackTransformation=WithBackTransformation)
    }

      return(xtrans)
    }
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