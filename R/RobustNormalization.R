RobustNormalization = function (Data,Centered=FALSE,Capped=FALSE,na.rm=TRUE,WithBackTransformation=FALSE,pmin=0.01,pmax=0.99) 
#RobustNormalization(Data,Centered=FALSE,Capped=FALSE)
#Normalizes features either between -1 to 1 (Centered=TRUE) or 0-1 (Centered=TRUE) without changing the distribution of a feature itself. For a more precise description please read [Thrun, 2018, p.17].
#INPUT
#Data							[1:n,1:d] data matrix of n cases and d features
#Centered						centered data around zero by median if TRUE
#Capped							TRUE: outliers are capped above 1 or below -1 and set to 1 or -1.
# na.rm							If TRUE, infinite vlaues are disregarded
# WithBackTransformation 		If in the case for forecasting with neural networks a backtransformation is required, this parameter can be set to 'TRUE'.
# pmin							defines outliers on the lower end of scale}
# pmax 							defines outliers on the higher end of scale}
#OUTPUT
#if WithBackTransformation=FALSE: 
# 	TransformedData				[1:n,1:d]  normalized data matrix of n cases and d features

#if WithBackTransformation=TRUE: List with
#TransformedData				[1:n,1:d]  normalized data matrix of n cases and d features
#MinX							[1:d] numerical vector used for manual back-transformation of each feature
#MaxX							[1:d] numerical vector used for manual back-transformation of each feature
#Denom							[1:d] numerical vector used for manual back-transformation of each feature
#Center							[1:d] numerical vector used for manual back-transformation of each feature

#author: MT
#[Milligan/Cooper, 1988]  Milligan, G. W., & Cooper, M. C.: A study of standardization of variables in cluster analysis, Journal of Classification, Vol. 5(2), pp. 181-204. 1988.

#[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 

{

  #if(is.data.frame(Data)){
  #  warning('Matrix is expected but data.frame is given. Calling as.matrix().')
  #  Data=as.matrix(Data)
  #}
if(isTRUE(na.rm)){
  # if(!is.data.frame(Data)){
    Data[!is.finite(Data)]=NaN #quantile does not accept inf,-inf
  # }else{
  #   ind=do.call(cbind, lapply(mtcars, is.finite))
  #   Data[ind]=NaN
  # }
   
}
  center=0
  Denom=NULL
  minX=NULL
  maxX=NULL
  if (is.vector(Data)) {
    quants = quantile(Data, c(pmin, 0.5, pmax),na.rm = na.rm)
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
        quants = quantile(Data, c(pmin, 0.5, pmax),na.rm = na.rm)
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
    #Data=checkInputDistancesOrData(Data,funname='RobustNormalization')
    if(WithBackTransformation){
      cols = ncol(Data)
      xtrans = list()
      DataOut=c()
      minX=c()
      maxX=c()
      Denom=c()
      center=rep(0,cols)
      for (i in 1:cols) {
        xtrans = RobustNormalization(Data = as.vector(Data[, i]),Centered = Centered,Capped = Capped,na.rm = na.rm,WithBackTransformation=WithBackTransformation,pmin=pmin,pmax=pmax)
        DataOut=cbind(DataOut,as.matrix(xtrans$TransformedData))
        minX[i]=xtrans$MinX
        maxX[i]=xtrans$MaxX
        Denom[i]=xtrans$Denom
        center[i]=xtrans$Center
      }
      names=colnames(Data)
      if(!is.null(names)){
        colnames(DataOut) = names
        names(minX) =  names
        names(maxX) =  names
        names(Denom)=  names
        names(center)=  names
      }
      
      return(list(TransformedData=DataOut,MinX=minX,MaxX=maxX,Denom=Denom,Center=center))
    }else{
      cols = ncol(Data)
      xtrans = Data
      for (i in 1:cols) {
        xtrans[, i] = RobustNormalization(as.vector(Data[, i]),Centered,Capped,na.rm,WithBackTransformation=WithBackTransformation,pmin=pmin,pmax=pmax)
      }
      names=colnames(Data)
      if(!is.null(names))
        colnames(xtrans)=names
      
      return(xtrans)
    }
  }
  else {
    tryCatch({
      warning("RobustNormalization:: Data is not a vector or a matrix. Trying as.matrix")
      #Data=checkInputDistancesOrData(Data,funname='RobustNormalization')
      return(RobustNormalization(as.matrix(Data),Centered,Capped,na.rm,pmin=pmin,pmax=pmax))
    }, error = function(e) {
      stop("It did not work")
    })
  }
}