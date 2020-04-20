checkInputDistancesOrData=function(Data,funname='checkInputDistancesOrData'){
  

  if(missing(Data))
    stop(paste0(funname,': Input of Data or Distances is missing.'))
	
  if(!is.matrix(Data)){
    warning(funname,': Input of Data or Distances is not a matrix. Trying to transform it to a matrix')
    
          if (is.data.frame(Data))
            DataOrDistance = data.matrix(Data)
          else
            DataOrDistance = as.matrix(Data)
        
  }
  
  if (!isSymmetric(unname(Data)))
    string='Data'
  else
    string='Distances'
  

  if(sum(!is.finite(Data))!=0){
    warning(paste0(funname,': Some entries in ',string,' are not finite or missing values. Computations may not work'))
  }
  if(!is.numeric(Data)){
    warning(paste0(funname,': ',string,' is not numeric. Trying to transform it to a numeric type.'))
    mode(Data)='numeric'
  }
  return(Data)
}