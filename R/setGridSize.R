setGridSize=function(InputDistances,minp=0.01,maxp=0.99,alpha=4){
#setGridSize(InputDistances)
#setGridSize(InputDistances,minp=0.01,maxp=0.99)
# sets the size of the grid, sometimes its required to look at the distribution and choose the optional arguents manually
  #INPUT
  # Inputdistaces[n,n]      distance matrix of input pace
  # minp                    percentile of minum distanes for quantile
  # map                     percentile of maximum distanes for quantile
  # alpha										intern Parameter for Martas Swarm in Java, do not change!
	#OUTPUT
  # LC                      c(Lines,Columsn) such that L*C >4N and L*C<16N and sqrt(L^2+C^2)/1=MaxDost/MinDist
  #                         also L%%2=0 and C%%4=0 (technical requirements for grid generation)
  # author MT 03/16
  #2. Version: Okt 2016 MT & FL
  InputDistances=checkInputDistancesOrData(InputDistances)
  
  if(!isSymmetric(unname(InputDistances))) {
    warning('InputDistances are not a symmetric distance matrix')
  }
  InputDistancesTmp=InputDistances
  N = nrow(InputDistances) #InputDistances
  if(is.null(N)) stop('Distance matrix has no rows')
  InputDistances = as.numeric(InputDistances[upper.tri(InputDistances, diag = F)])
  alpha = alpha #Anzahl moeglicher Sprungpositionen der DatenBots
  p01 <- quantile(InputDistances, minp)
  p99 <- quantile(InputDistances, maxp)
  jpos = alpha * N #Verhaeltniss DataBots zu freien Gitterpositionen
  
  A = p99 / p01#Verhaeltnis der kleinsten zur groessten Distanz
  #Dieses verhaeltniss muss minimal im Gitter gewahrt bleiben
  #Verhaeltnis der kleinsten zur groessten Distanz=MaxDist/MinDist
  #C^4-A*C^2+16N^2>=0
  C = suppressWarnings(sqrt(1/2)*sqrt(A ^ 2 + sqrt(A ^ 4  - (alpha ^ 2 * N ^ 2)/4)))#Columns
  
  if(identical(C, numeric(0))) C=NaN
	if(!is.finite(C)) C=NaN
  #C=Re(sqrt(A^2/2+sqrt(A^4/4-16*N^2+0i)))#Columns
  if (is.nan(C)) {#Wenn es keine reele Loesung gibt, approximiere
    #Operator: An approximation of grid size was done.
    print('Operator: An approximation of grid size was done.')
    solutions = matrix(ncol = 2, nrow = 0)
    
    #### all candidates by first and third equation
    for (columns in 1:10000) {
      suppressWarnings({
        lines = sqrt(A ^ 2 - columns ^ 2)
      })
      if (is.na(lines)) lines=NaN
        
      ratio = columns / lines
      if (!is.nan(lines)) {
        ratio = columns / lines
        if ((ratio > 0.6) && (ratio <= 1.0)) {
          if (columns < lines)
            solutions = rbind(solutions, c(columns, lines))
          else
            solutions = rbind(solutions, c(lines, columns))
        }
      }
    }
    if(nrow(solutions)==0){
      if (minp != 0.05 & maxp != 0.95){
       # print(InputDistances)
        return(setGridSize(InputDistances=InputDistancesTmp, minp = 0.05, maxp = 0.95))
      }else{
        warning('setGridSize failed to choose the right grid, choosing standard. Maybe Change quantiles minp and/or maxp')
        if (50 * 80 >= jpos) {
          return(LC = c(50, 80))
        } else if (100 * 160 >= jpos) {
          return(LC = c(100, 160))
        } else{
          return(LC = c(200, 320))
        }
      }
        
    }
    colnames(solutions) = c('Lines', 'Columns')
    # filter everything that doesnt meet the second condition
    k = 1
    gitterausdehung = seq(from = 1, to = 500, by = 0.01)
    filteredSolutions = matrix(ncol = 2, nrow = 0)
    while (nrow(filteredSolutions) < 1) {
      solutions = solutions * gitterausdehung[k] #Gitter kann beliebig vergroessert werden
      for (i in 1:nrow(solutions)) {
        if ((solutions[i, 1] * solutions[i, 1]) >= jpos) {
          filteredSolutions = rbind(filteredSolutions, solutions[i, ])
        }
      }
      k = k + 1
      if (k > 500)
        break
      
    }
    if (nrow(filteredSolutions) < 1) {# Falls keine approximation Loesung gefunden werden konnte
      if (minp != 0.05 & maxp != 0.95)
        setGridSize(InputDistancesTmp, minp = 0.05, maxp = 0.95) #aendere quantile (R quantil funkion komisch...)
    } else{ #erfuelle nebenbedingung des hexagonalen gitters
      filteredSolutions = round(filteredSolutions)
      ind = which(filteredSolutions[, 2] %% 4 == 0)
      if (length(ind) > 0) {
        Columns = filteredSolutions[ind[1], 2]
        Lines = filteredSolutions[ind[1], 1]
        while (Lines %% 2 != 0)
          Lines = Lines + 1
      } else{
        ind = which(filteredSolutions[, 1] %% 2 == 0)
        if (length(ind) > 0) {
          Columns = filteredSolutions[ind[1], 2]
          Lines = filteredSolutions[ind[1], 1]
          while (Columns %% 4 != 0)
            Columns = Columns + 1
        } else{
          if (nrow(filteredSolutions) > 0) { #Wenn aus dem Loesungraum keine geeignete fuer Nebenbedingung da war
           # LC = filteredSolutions[1, ] #nehme erste und passe sie an die nebenbedingung an
            Columns=filteredSolutions[1, 2]
            Lines=filteredSolutions[1, 1]
            while (Columns %% 4 != 0)
              Columns = Columns + 1
            while (Lines %% 2 != 0)
              Lines = Lines + 1
          } else{#Guittergroesse ist gegenueber DBS sehr robust, Wenn alles andere fehlschlaegt, nehme standard
            warning('setGridSize failed to choose the right grid, choosing standard')
            if (50 * 80 >= jpos) {
              return(LC = c(50, 80))
            } else if (100 * 160 >= jpos) {
              return(LC = c(100, 160))
            } else{
              return(LC = c(200, 320))
            }
          }
        }
      }
    }
  } else{#  !is.nan(C), Loesung wurde gefunden
    Columns = round(C)
    Lines = round(alpha*N/C)
    if (Lines > Columns) {
      tmp = Columns
      Columns = Lines
      Lines = tmp
    }
    while ((Lines / Columns) <= 0.7) { # 2ter Wert moeglicherweise zu schraeg, passe dann an
      Lines = Lines + 2
    }
  }
#Pruefe, dass Lines<Columns und die Nebenbedingungen des hexagonalen Gitters erfuellt sind
  if (Lines > Columns) {
    tmp = Columns
    Columns = Lines
    Lines = tmp
  }
  while (Columns %% 4 != 0)
    Columns = Columns + 1
  while (Lines %% 2 != 0)
    Lines = Lines + 1
  LC = c(Lines, Columns)
  names(LC) = c('Lines', 'Columns')
  return(LC)
}
