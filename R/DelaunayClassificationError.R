DelaunayClassificationError=function(Data,ProjectedPoints,Cls,LC){
# DCE  = DelaunayClassificationError(Data,ProjectedPoints,Cls) 
# [DCE,DCEperPoint] = DelaunayClassificationError(Data,ProjectedPoints,Cls) 
# [DCE,DCEperPoint,nn,AnzNN,NNdists,HD] = DelaunayClassificationError(OutputDelaunay,InputDistances,Cls) 
# calculates DCE i.e the Delaunay classification error
#
# INPUT
# Data[1:n,1:d]					n cases d variales
# ProjectedPoints[1:n,1:2]		in xy cartesischen koordinaten					
# Cls[1:n]						classen labels, eine Ziffer je klasse, numerischer vector fuer n cases aus dem Zahlenraum von 1:p Klassen				
#
# OUPUT
# 
# DCE                            DelaunayClassificationError 
#                                NOTE the rest is just for development purposes
# DCEperPoint(1:n)               unnormalized DCE of each point: DCE = mean( DCEperPoint)
# nn                             the number of points in a relevant neghborhood: 0.5 * 85percentile(AnzNN)
# AnzNN(1:n)                     the number of points with a delaunay graph neighborhood
# NNdists(1:n,1:nn)              the distances within the relevant nehborhoot, 0 for inner cluster distances
# HD(1:nn)                       HD = HarmonicDecay(nn) i.e weight function for the NNdists: DCEperPoint = HD*NNdists

# author: MT 07/2016, Sept 2016, Ausgabeparameteranpassung an ALU
# Checked in 06/2018 by MT
## Catch Wrong Input
  if (is.list(ProjectedPoints))
    stop('ProjectedPoints is a list not a matrix')
  
  c = ncol(ProjectedPoints)
  if (c > 3 | c < 2)
    stop(paste0('Wrong number of Columns of ProjectedPoints: ', c))
  
  if (c == 3) {
    #Falls mit Key muss dieser in erster Spalte sein
    ProjectedPoints = ProjectedPoints[, 2:3]
  }
#Calculate Delaunay Graph in Outputspace  
  if (missing(LC))
    #Berechne Delaunay Graphen, eigene Funktion, falls das mal auf CRAN soll
    del = Delaunay4Points(Points=ProjectedPoints, IsToroid=F)
  else
    # Im toroiden fall gibts die ESOM definition, wo x und y vertauscht sind
    #del = Delaunay4Points(Points = ProjectedPoints, IsToroid = TRUE, Grid = LC[c(2, 1)])
    del = Delaunay4Points(Points = ProjectedPoints, IsToroid = TRUE, LC = LC)
  #######################################################################
  HarmonicDecay = function(n) {
    #  HD = HarmonicDecay(n)
    # calculates the Harmonic Decay Numbers HarmonicDecay(i,n) for i=1...n
    # HarmonicDecay(i,n) = sum(1/k) , k=i...n
    # Laws:
    # HarmonicDecay(1,n) = Harmonic number HN(n) = 1+1/2+...1/n
    #                      this approximates to:  log(n)+ 0.5772156649+ 1/(2*n)
    # HarmonicDecay(n,n) = 1/n
    #
    # INPUT
    # n      n > 0 natural number
    #
    # OUTPUT
    # HD(1:n)       HD(i) = HarmonicDecay(i,n)
    
    # MT Sept 2016
    
    if (n < 1)
      HD = NaN
    
    OneToN = c(1:n) # 1...n
    HD = OneToN * NaN
    OneOverN = 1 / OneToN    # 1/1, ... 1/n
    for (i in 1:n) {
      HD[i] = sum(OneOverN[i:n])  # 1/i+ ...+ 1/n
    }
    #plot(OneToN,HD,ylab = 'HarmonicDecay(i,n)',xlab = 'i',main = 'HarmonicDecay')
    return(HD)
  }
  ####################################################################
  
  #Gewichte mit Euklid der Delaunay kanten
  Dist = as.matrix(dist(Data)) * del
  ind = which(Dist == 0, arr.ind = T)
  Dist[ind] = NaN #Fuer spaeteres Sortieren sollten nicht direkt verbundene Punkte NaN sein, statt 0
  diag(Dist) = NaN
  
  AnzPunkte = length(Cls)
  if (length(Cls) < 3)
    stop('DelaunayClassifError: at least 3 points required')
  
  IsInterDistance = Dist * NaN     # IsInterDistance(i,j) =1 wenn die punkte i,j in verschiedene clustern
  DelaunayDists = Dist * NaN
  for (i in 1:nrow(Data)) {
    SortInd = order(Dist[, i], decreasing = F, na.last = T)
    currentClass = Cls[i]
    cc = Cls[SortInd]
    IsInterDistance[, i] = currentClass != cc   # die inneren Distanzen Null setzen
    DelaunayDists[, i]  = Dist[SortInd, i]
  }
  
  # jetzt die Anzahl der betrachten naesten Nachbarn und die Gewichtung ermitteln
  AnzNN = c()
  for (i in 1:nrow(Data)) {
    AnzNN = c(AnzNN, sum(!is.nan(Dist[, i]))) #Anzahl echter Voronoi-Zelen
  }
  NN95 = round(stats::quantile(x = AnzNN, c(0.95)), 0)
  
  # die Betrachtete nachbarschaft sind die Nachsten Nachbarn von 1...nn
  nn = round(NN95 / 2)
  
  # betrachtete Distanzen
  NNdists = IsInterDistance[1:nn, ]  # die ersten nn direkten Nachbahrn bzfl der Distanzengroesse
  
  for (i in 1:nrow(Data)) {
    #Setze 0, falls nicht direkt benachbahrt
    if (AnzNN[i] < nn) {
      NotInterInd = seq(from = AnzNN[i] + 1,
                        to = nn,
                        by = 1)
      NNdists[NotInterInd, i] = 0
    }
  }
  
  NNdists = t(NNdists)
  
  # die Gewichtungsfuktion ausrechnen
  HD = HarmonicDecay(nn)
  # Gewicht * Distanzen
  DCEperPoint = NNdists %*% HD   # vektorielle multiplikation
  
  DCE = mean(DCEperPoint, na.rm = T)  # DCE = DCEperPoint/AnzPunkte
  
  return(
    list(
      DCE = DCE,
      DCEperPoint = DCEperPoint,
      nn = nn,
      AnzNN = AnzNN,
      NNdists = NNdists,
      HD = HD,
      IsInterDistance = IsInterDistance,
      DelaunayDists = DelaunayDists
    )
  )
}    
           