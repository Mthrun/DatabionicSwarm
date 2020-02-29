checkInputDistancesOrData=function(Data){
  

  if(missing(Data))
    stop('Input is missing')
	
  if(!is.matrix(Data)){
    warning('Input of Data or Distances is not a matrix. Trying to transform it to a matrix')
    Data=as.matrix(Data)
  }
  
  if (!isSymmetric(unname(Data)))
    string='Data'
  else
    string='Distances'
  

  if(sum(!is.finite(Data))!=0){
    warning(paste0('Some entries in ',string,' are not finite or missing values. Computations may not work'))
  }
  if(!is.numeric(Data)){
    warning(paste0(string,' is not numeric. Trying to transform it to a numeric type.'))
    Data=as.numeric(Data)
  }
  return(Data)
}