RobustNorm_BackTrafo=function(TransformedData,MinX,Denom,Center=0){
  Data=TransformedData*NaN
  if(is.matrix(TransformedData)){
    for(i in 1:ncol(TransformedData))
      Data[,i]=RobustNorm_BackTrafo(TransformedData[,i],MinX[i],Denom[i],Center[i])
    return(Data)
  }else{
    return((TransformedData+Center)*Denom+MinX)
  }
}