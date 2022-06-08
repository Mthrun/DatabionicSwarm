UniquePoints <- function(Datapoints,Eps) {
  # V <- UniquePoints(Datapoints)
  # return only the unique points in Datapoints
  #
  # INPUT
  # Datapoints[1:n,1:d]   [1:n,1:d] matrix of Datapoints points of dimension d
  #				                the points are in the  rows
  #
  # Eps                   Optional,scalarabove zero that defines minimum non-identical euclidean distance between two points
  # OUTPUT
  # a list V containg:
  # Unique                  [1:u,1:d]      the Datapoints points  without duplicate points
  # IsDuplicate             [1:n,1:n]      for i!=j IsDuplicate[i,j]== 1  if Datapoints[i,] == Datapoints[j,]    IsDuplicate[i,i]==0
  # UniqueInd               [1:u]		       an index vector such that Unique ==  Datapoints[UniqueInd,]
  # Uniq2DatapointsInd      [1:n] 	       an index vector such that   Datapoints ==  Unique[Uniq2DatapointsInd,]
  #MT 06/2022
  if(missing(Eps))
  Eps =0.0000000001         # ab dieser Distanz zwischen 2 punkten sind diese identisch
  
  if (inherits(Datapoints,"matrix")) {
    # If Datapoints is a vector.
    Datapoints <- as.matrix(Datapoints)
  }
  AnzPoints <-
    nrow(Datapoints)            # soviele punkte in den Daten
  rownames(Datapoints) <-
    c(1:AnzPoints)   # gib den Zeilen als namen ihren zeilennummer
  
  if(!requireNamespace('parallelDist')){
    warning("UniquePoints: package parallelDist is not installed, falling back to dist().")
    dists = as.matrix(dist(Datapoints))     # Distance with diagonal = 0.
  }else{
    dists = as.matrix(parallelDist::parDist(Datapoints, method = "euclidean"))
  }
  

  IsDuplicate = (dists < Eps) * 1 - diag(AnzPoints)
  dists <-
    dists - (dists * upper.tri(dists, diag = T)) + upper.tri(dists, diag = T) # ??? wozu das denn?
  
  if (length(which(dists < Eps)) > 0) {
    # Duplicates found.
    ind = which(dists < Eps, arr.ind = TRUE) # Get indices of duplicates.
    #    rownames(ind) <- ind[,1]

    ind = ind[as.character(unique(ind[, 1])), ,drop=FALSE] # remove multiples in the duplicates, so that only their first occurance remains. Example: if 1, 3, 4, and 6 are all duplicates of the same value, ind will contain [3,1], [4,1], [4,3], [6,1], [6,3] and [6,4]. This removes all except [3,1], [4,1] and [6,1]
   
    uniqueDatapoints = Datapoints[-as.matrix(ind)[, 1], ,drop=FALSE] #MT: Korrektur, falls genau eine Dopplung besteht
  
    mergeind <- c(1:AnzPoints)
    indhelp <- c(1:nrow(uniqueDatapoints))
    names(indhelp) <- rownames(uniqueDatapoints)
    mergeind[as.numeric(names(indhelp))] <- indhelp
    if (ncol(as.matrix(ind)) == 1) {
      mergeind[as.matrix(ind)[, 1]] = as.matrix(ind)[, 1]#MT: Schnellschuss Workaround
    } else{
      mergeind[ind[, 1]] <-
        mergeind[ind[, 2]] # Datapointsindex with marked duplicates
    }
    return(
      list(
        "Unique" = as.matrix(uniqueDatapoints),
        "UniqueInd" = as.numeric(rownames(uniqueDatapoints)),
        "Uniq2DatapointsInd" = mergeind,
        IsDuplicate = IsDuplicate
      )
    )
  } else{
    # keine duplikate gefunden
    return(
      list(
        "Unique" = Datapoints,
        "UniqueInd" = c(1:AnzPoints), #sortind
        "Uniq2DatapointsInd" = c(1:AnzPoints),
        IsDuplicate = IsDuplicate
      )
    )
  }# end   if(length(which(dists<Eps))>0){ # Duplicates found.
}# end
