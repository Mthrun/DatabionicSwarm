Delaunay4Points <- function(Points, IsToroid = TRUE,Grid=NULL,PlotIt=FALSE,Gabriel=FALSE){
# Delaunay=Delaunay4Points(BestMatches, IsToroid,Grid,PlotIt)$Delaunay
# Calculates the adjacency matrix of the delaunay graph for bestmatches in tiled form if BMs are located on a toroid grid
#
# INPUT
# BestMatches[1:n,1:3]            n by 3 matrix containing the BMKey, X and Y coordinates of the n BestMatches
#                                 BestMatches NEED NOT BE UNIQUE!
#                                 however, there is an edge in the Deaunay between duplicate points!  
#
# OPTIONAL
# Grid[2]                         A vector of length 2, containing the number of lines and columns of the Grid
# IsToroid                        logical, indicating if BM's are on a toroid grid. Default is True
# PlotIt                          Set PlotIt=TRUE, if you want to see the Plots
# OUTPUT
# DelaunayAdjazenzMatrix[1:n,1:n]  adjacency matrix of the Delaunay-Graph
#
# NOTE: Im Unterschied zu Delauany4Bestmatches hier in cartesischer Definition, Testweise auch Gabriel Graph moeglich
# 	     Die Unterfunktionen wurden fuer diesen einen Zweck noch nachoptimiert  
# authors:  MT 03/16
  if(is.list(Points))
    stop('Points is a list not a matrix')
  if (ncol(Points) > 3)
    stop('Points have wrong number of dimensions')
  if (ncol(Points) < 1)
    stop('Points have wrong number of dimensions')
  
  # if(ncol(Points)==2)
  #   if(!is.null(Grid)) stop('Points have wrong number of dimensions or LC is not NULL')
  if (ncol(Points) == 3) {
    Points = Points[, 2:3]
  }
  if (IsToroid)
    if (is.null(Grid))
      stop('Grid has to be set, if toroid=TRUE')
  if (!is.null(Grid)) {
    Xgrid = Grid[1]
    Ygrid = Grid[2]
  }
  if (IsToroid) {
    if (max(Points[, 1]) > Grid[1])
      stop('Grid[1]>max(Points[,1])')
    if (max(Points[, 2]) > Grid[2])
      stop('Grid[2]>max(Points[,2])')
  }
  #####################################################################################
  DelaunayGraphMatrix_hlp <- function(X, Y, PlotIt = FALSE) {
    # D <- DelaunayGraphMatrix(X,Y,PlotIt)
    # Calculates the Adjacency Martix for a Delaunay-Graph
    # INPUT
    # X[1:n],  Y[1:n]         point coordinates for which a planar Delaunay graph is constucted
    #                         POINTS NEED NOT BE UNIQUE!
    #                         however, there is an edge in the Delaunay between duplicate points!
    # OPTIONAL
    # PlotIt                  if TRUE Delaunay graph is plotted
    #
    # OUTPUT
    # a list of
    # Delaunay[1:n,1:n]       The adjacency matrix of the Delaunay-Graph.
    #
    
    # uses packages deldir
    
    uniquePoints_hlp <- function(data) {
      # U <- uniquePoints(data)
      # return only the unique points in data
      #
      # INPUT
      # data[1:n,1:d]   The vector/matrix with n data points of dimension d
      #				          the points are in the  rows
      #
      # OUTPUT
      # a list U containg:
      # UniqueData  <- U$unique[1:u,1:d]       the data points  without duplicate points
      # UniqueInd   <- U$sortind[1:u]		       an index vector such that unique ==  data[sortind,]
      # Uniq2DataInd<- U$mergeind[1:n] 	       an index vector such that   data ==  unique[mergeind,]
      # IsDuplicate <- U$IsDuplicate[1:n,1:n]  for i!=j IsDuplicate[i,j]== 1  iff data[i,] == data[j,]    IsDuplicate[i,i]==0
      
      eps <-
        0.0000000001         # ab dieser Distanz zwischen 2 punkten sind diese identisch
      
      if (inherits(data,"matrix")) {
        # If data is a vector.
        data <- as.matrix(data)
      }
      AnzPoints <-
        nrow(data)            # soviele punkte in den Daten
      rownames(data) <-
        c(1:AnzPoints)   # gib den Zeilen als namen ihren zeilennummer
      dists <-
        as.matrix(dist(data))     # Distance with diagonal = 0.
      IsDuplicate = (dists < eps) * 1 - diag(AnzPoints)
      dists <-
        dists - (dists * upper.tri(dists, diag = T)) + upper.tri(dists, diag = T) # ??? wozu das denn?
      
      if (length(which(dists < eps)) > 0) {
        # Duplicates found.
        ind <-
          which(dists < eps, arr.ind = TRUE) # Get indices of duplicates.
        #    rownames(ind) <- ind[,1]
        ind <-
          ind[as.character(unique(ind[, 1])), ]  # remove multiples in the duplicates, so that only their first occurance remains. Example: if 1, 3, 4, and 6 are all duplicates of the same value, ind will contain [3,1], [4,1], [4,3], [6,1], [6,3] and [6,4]. This removes all except [3,1], [4,1] and [6,1]
        uniquedata <-
          data[-as.matrix(ind)[, 1], ] #MT: Korrektur, falls genau eine Dopplung besteht
        mergeind <- c(1:AnzPoints)
        indhelp <- c(1:nrow(uniquedata))
        names(indhelp) <- rownames(uniquedata)
        mergeind[as.numeric(names(indhelp))] <- indhelp
        if (ncol(as.matrix(ind)) == 1) {
          mergeind[as.matrix(ind)[, 1]] = as.matrix(ind)[, 1]#MT: Schnellschuss Workaround
        } else{
          mergeind[ind[, 1]] <-
            mergeind[ind[, 2]] # Dataindex with marked duplicates
        }
        return(
          list(
            "unique" = as.matrix(uniquedata),
            "sortind" = as.numeric(rownames(uniquedata)),
            "mergeind" = mergeind,
            IsDuplicate = IsDuplicate
          )
        )
      } else{
        # keine duplikate gefunden
        return(
          list(
            "unique" = unique(data),
            "sortind" = c(1:AnzPoints),
            "mergeind" = c(1:AnzPoints),
            IsDuplicate = IsDuplicate
          )
        )
      }# end   if(length(which(dists<eps))>0){ # Duplicates found.
    }# end
    
    # Punkte unique machen
    unique = uniquePoints_hlp(cbind(X, Y))
    UniqXY = unique$unique
    UniqueInd   = unique$sortind
    Uniq2DataInd = unique$mergeind
    IsDuplicate = unique$IsDuplicate
    UniqX =  UniqXY[, 1]
    UniqY =  UniqXY[, 2]
    
    # Delaunay ausrechnen mit deldir
    DeldirOutput = deldir(UniqX , UniqY)  #
    PointsAndIndices = DeldirOutput$delsgs # dadrin stecken die indices des Delaunays von -> Nach
    FromInd =  PointsAndIndices$ind1     #  indices der Ausgangspunkte des Delanays
    ToInd   =  PointsAndIndices$ind2     #  indices der Endpunkte des Delanays
    
    # Adjazenzmatrix befuellen
    UniqDelaunay = matrix(0, length(UniqX), length(UniqY))  # Adjazenzmatrix initialisieren
    for (i in c(1:length(FromInd))) {
      UniqDelaunay[FromInd[i],  ToInd[i]]  <-
        1  #Only Direct neighbours A and B get an one from A to B
      UniqDelaunay[ToInd[i]  , FromInd[i]]  <-
        1  #Only Direct neighbours A and B get an one from A to B
    } # end for i neighbours
    
    # jetzt uniqe points wieder auf originale uebertragen
    Delaunay = matrix(0, length(X), length(Y))
    Delaunay = UniqDelaunay[Uniq2DataInd, Uniq2DataInd]
    Delaunay = Delaunay + IsDuplicate # noch je eine Verbindung zwischen den Doubletten eintragen
    
    # ausgabe zusammenstellen
    return(Delaunay = Delaunay)
  }# end function  DelaunayGraphMatrix
  ###############################################################################################
  ##############################################################################################
  if (IsToroid) {
    # behandlung eines pack-man (toroiden) MapSpaces
    
    TiledX = Points[, 1]
    TiledY = Points[, 2]
    TiledX = c(TiledX, TiledX + Xgrid, TiledX + Xgrid, TiledX)
    TiledY = c(TiledY, TiledY, TiledY + Ygrid, TiledY + Ygrid)
    if (Gabriel) {
      Delaunay = calcGabrielGraph2D(cbind(TiledX, TiledY), PlotIt = PlotIt)
    } else{
      Delaunay = DelaunayGraphMatrix_hlp(TiledX, TiledY, PlotIt = PlotIt)
    }
    
    
    TiledDelaunay =  Delaunay
    n = nrow(Points)
    AnzBestMatches = nrow(Points)
    
    Offset = c(0, 1, 2, 3) * AnzBestMatches
    DelaunayAdjazenzMatrix = TiledDelaunay[c(1:AnzBestMatches), c(1:AnzBestMatches)]  # initialisieren mit dem oberen viertel
    for (i in Offset) {
      for (j in Offset) {
        DelaunayAdjazenzMatrix = DelaunayAdjazenzMatrix + TiledDelaunay[c(1:AnzBestMatches) +
                                                                          i, c(1:AnzBestMatches) + j]
      }# end for j
    }# end for i
    # zurueckrechnen auf die ungekachelten werte
    DelaunayAdjazenzMatrix = (DelaunayAdjazenzMatrix > 0) * 1  # alles auf 0/1 reduzieren
    # X=TiledX
    # Y=TiledY
    
  }
  else{
    # MapSpace ist planar
    
    X = Points[, 1]
    Y = Points[, 2]
    if (Gabriel) {
      Delaunay = calcGabrielGraph2D(cbind(X, Y), PlotIt = PlotIt)
    } else{
      Delaunay <- DelaunayGraphMatrix_hlp(X, Y, PlotIt)
    }
    DelaunayAdjazenzMatrix  = Delaunay
    #TiledDelaunay = DelaunayAdjazenzMatrix
  }# end if(IsToroid)
  
  # Ausgabe zusammenstellen
  return(DelaunayAdjazenzMatrix)
}# end function