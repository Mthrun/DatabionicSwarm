RelativeDifference=function(X,Y,epsilon=10^-10){
  if(length(X)!=length(Y)) stop('Length of X and Y do not match.')
  if(length(X)>1) return(mapply(X,FUN = RelativeDifference,Y))
  oben=Y-X
  unten=X+Y
  if(abs(unten)<epsilon){
    warning('X and Y are too small to calcualte Relative Differences. Returning 0')
    return(0)
  }
  if(X<0) stop('Not defined for negative X values')
  if(Y<0) stop('Not defined for negative Y values')
  return(2*oben/unten)
}