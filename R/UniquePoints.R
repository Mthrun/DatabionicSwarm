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
  # a list V containing:
  # Unique                  [1:k,1:d]      the Datapoints points  without duplicate points
  # IsDuplicate             [1:n,1:n]      for i!=j IsDuplicate[i,j]== 1  if Datapoints[i,] == Datapoints[j,]    IsDuplicate[i,i]==0
  # UniqueInd               [1:k]		       an index vector 1:k such that Unique ==  Datapoints[UniqueInd,], it has k non-consecutive numbers or labels, each label defines a row number within Datapoints[1:n,1:d] of a unique data point
  # Uniq2DatapointsInd      [1:n] 	       an index vector 1:n. It has k unique index numbers representing the arbitrary labels. Each labels is mapped uniquely to a point in Unique. Logically in a way such that Datapoints ==  Unique[Uniq2DatapointsInd,] (will not work directly in R this way) 
  #
  #Description: Euclidean distance is computed and used within. Setting \code{Eps} to a very small number results in the identification of unique data points. Setting epsilon to a higher number results in the definition of mesh points within an d-dimensional R-ball grap
  #MT 06/2022
  
  # ab dieser Distanz zwischen 2 punkten sind diese identisch
  if(missing(Eps))
  Eps =0.0000000001        
  
  if (inherits(Datapoints,"matrix")) {
    # If Datapoints is a vector.
    Datapoints <- as.matrix(Datapoints)
  }
   # soviele punkte in den Daten
  AnzPoints = nrow(Datapoints)           
  rownames(Datapoints) = c(1:AnzPoints)   # gib den Zeilen als namen ihren zeilennummer
  
  if(!requireNamespace('parallelDist')){
    warning("UniquePoints: package parallelDist is not installed, falling back to dist().")
    dists = as.matrix(dist(Datapoints))     # Distance with diagonal = 0.
  }else{
    dists = as.matrix(parallelDist::parDist(Datapoints, method = "euclidean"))
  }
  
  IsDuplicate=rep(FALSE,AnzPoints)

  #in der diagonalen nicht suchen
  diag(dists)=Inf
  #nur in eine richtung suchen reicht, sonst werden beide punkte entfernt
  dists[upper.tri(dists)]=Inf

    #zeile ist die dublette
    #spalte ist der datenpunkt auf den die dublette verweist
    #kann moeglicherweise wiedereine dublette sein
    ind_all = which(dists < Eps, arr.ind = TRUE) # Get indices of duplicates.
    
  # keine duplikate gefunden
  if (length(ind_all) == 0) {
    return(
      list(
        "Unique" = Datapoints,
        "UniqueInd" = c(1:AnzPoints), #sortind
        "Uniq2DatapointsInd" = c(1:AnzPoints),
        IsDuplicate = IsDuplicate
      )
    )
  }#end ifno duplicates

# remove multiples in the duplicates, so that
# only their first occurance remains. Example: if 1, 3, 4, and 6 are all duplicate
# of the same value, ind will contain [3,1], [4,1], [4,3], [6,1], [6,3] and [6,4]. 
# This removes all except [3,1], [4,1] and [6,1]
   ind = ind_all[as.character(unique(ind_all[, 1])), ,drop=FALSE] 
   #aus irgendeinem grund ist das as.character wichtig!
 #ind = ind[unique(ind[, 1]), ,drop=FALSE] #funktioniert nicht
    IsDuplicate[ind[, 1]]=TRUE
    uniqueDatapoints = Datapoints[IsDuplicate==FALSE, ,drop=FALSE] #MT: Korrektur, falls genau eine Dopplung besteht

    #init, wir initialisieren mit den index aller,
    #damit wir fuer die dubletten mit den index von ind arbeiten koennen
     mergeind <- c(1:AnzPoints)
     #index der uniquen
     indhelp <- c(1:nrow(uniqueDatapoints))
     #index der echten punkte
     names(indhelp) <- rownames(uniqueDatapoints)
     #setze an position der echten punkte
     #den index der uniquen
     mergeind[as.numeric(names(indhelp))] <- indhelp
     #befuelle die dubletten
     #Achtung, vorher muessen wir noch sicherstellen, das auf keine
     #weiteren dubletten in der 2ten spalte von ind verwiesen wird
     #suche alle fuer die das zutrifft
     uniquerowlabels=as.numeric(rownames(uniqueDatapoints))
     #finde diejenigen verweise, welche von dublett aud dublett gehen
     notfound=which(!(ind[,2]%in%uniquerowlabels))
     #ergaenze den vektor, sodass der verweis
     #auf ein nicht dublett geht
     ind_complete=ind
     notfound_point_ind=unique(ind_complete[notfound,2])
     for(j in notfound_point_ind){
       #welches label is geeignet?
       uniquePointLabelInd=ind_complete[ind_complete[,1]==j,2]
       #koennenvielen sein
       
       #aus irgendeinem grund findet er manchmal keine ueberschneidung
       #GetuniquePointLabelInds=ind_all[ind_all[,1]==j,2]
       #wir wollen davon keine dubletten
       #uniquePointLabelInd=intersect(GetuniquePointLabelInds,uniquerowlabels)
       
       #dahin soll das label gespeichert werden:
       if(length(uniquePointLabelInd)==1){
        ind_complete[ind_complete[,2]==j,2]=uniquePointLabelInd
       }else{#sollteesnicht geben
         warning("UniquePoints: Something went wrong...")
         #message(GetuniquePointLabelInds)
         #message(uniquerowlabels)
         #ind_complete[ind_complete[,2]==j,2]=NA
       }
        
     }
     
     #der index an der stelle der duplikate ist der werte von mergeind
     # an der stelle wohin das duplikat in der distanz fuert
     #und dieses muesste ja zu einem wert von indhelp fueren
     mergeind[ind_complete[, 1]] = mergeind[ind_complete[, 2]] # Datapointsindex with marked duplicates

    return(
      list(
        "Unique" = as.matrix(uniqueDatapoints),
        "UniqueInd" = as.numeric(rownames(uniqueDatapoints)),
        "Uniq2DatapointsInd" = mergeind,
        IsDuplicate = IsDuplicate,
        ind_all=ind_all #nur fürs debugging, spaeter entfernen
      )
    )
  
}# end function
