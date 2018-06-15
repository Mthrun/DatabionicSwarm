DBSclustering=function(k,DataOrDistance,BestMatches,LC,StructureType=TRUE,PlotIt=FALSE,method='euclidean',...){
#Cls=DBSclustering(k,Data,BestMatches,LC,StructureType=TRUE,PlotIt=F,method='euclidean')
# automated Clustering approach of the DataBionicSwarm with abstact U distances
# INPUT
# k                   number of classes, how many to you see in the 3d landscape?
# DataOrDistance      Matrix of Data or Distance that will be used. One DataPoint per row
# BestMatches         Array with positions of Bestmatches=ProjectedPoints
# LC
# OPTIONAL  
# StructureType           compact structure of clusters assumed, =FALSE: connected structure of clusters assumed
# PlotIt              Plots Dendrogramm
# method              do not change 
# OUTPUT 
# Cls                 vector with selected classes of the bestmatches
#  
# author: MT 06/16 
#  
#  NOTE: das ist eine eigene Idee: Nehme die Distantz in LoetschUltsch2014 definiert und stecke sie in
#		Distanz=ShortestGraphPaths -> Cluster sind immer in sich geschlossen in 2D
#							-> Politische Karte irrelevant -> Nahteil: Falls Projektion Fehler hat (Punke des einen Clusters innerhalb des anderen Clusters)
#							-> Hat Clusterung fehler und CLusterung ist nichtmehr dichte basiert
# Alus 2016 IDEE ueber die Distanz ist in AUstarDist.R implementiert
  #requireRpackage('deldir')
  #requireRpackage('geometry')
  DataOrDistance=checkInputDistancesOrData(DataOrDistance)
  
  if (isSymmetric(DataOrDistance)) {
    InputD = DataOrDistance
    rnames=1:nrow(DataOrDistance)
  } else{
    if(!is.null(rownames(DataOrDistance)))
      rnames=rownames(DataOrDistance)
    else
      rnames=1:nrow(DataOrDistance)
    requireNamespace('parallelDist')
    InputD = as.matrix(parallelDist::parDist(DataOrDistance, method = method,...))
  }# end if(isSymmetric(DataOrDists))
  
 GabrielGraph=FALSE #gabriel graph immer schlechter...
  GOutput=Delaunay4Points(BestMatches, Grid = LC, IsToroid=T,PlotIt=F,Gabriel=GabrielGraph)
    
    Dist=ShortestGraphPathsC(GOutput,InputD)
  if(StructureType){
    pDist=as.dist(Dist)
    hc <- hclust(pDist,method="ward.D")
    m="Compact DBS clustering"
  }else{
    ind=which(GOutput==0,arr.ind=T)
    Dist2=Dist*GOutput
    Dist2[ind]=max(Dist)*2
    pDist=as.dist(Dist2)
    hc <- hclust(pDist,method="single")
    m="Connected DBS clustering"
  }
  if(PlotIt){
    x=as.dendrogram(hc)
    plot(x, main=m,xlab="No. of Data Points N", ylab="Distance",sub=" ",leaflab ="none")
    axis(1,col="black",las=1)
  }
  
  Cls=cutree(hc,k)

  counter = 0
  bool = T
  NumberOfClassesSet = k
  while (counter < 0.05 * length(Cls)) {
    counter = counter + 1
    uniqueClasses <- sort(na.last = T, unique(Cls))
    numberOfClasses <- length(uniqueClasses)
    countPerClass <- rep(0, numberOfClasses)
    for (i in 1:numberOfClasses) {
      inClassI <-
        sum(Cls == uniqueClasses[i]) # counts all occurances of uniqueClass[i] in cls
      countPerClass[i] = inClassI
    }
    
    if (max(countPerClass) > 0.95 * length(Cls)) {
      k = k + 1
      Cls = cutree(hc, k)
      uniqueClasses <- sort(na.last = T, unique(Cls))
      numberOfClasses <- length(uniqueClasses)
      countPerClass <- rep(0, numberOfClasses)
      for (i in 1:numberOfClasses) {
        inClassI <-
          sum(Cls == uniqueClasses[i]) # counts all occurances of uniqueClass[i] in cls
        countPerClass[i] = inClassI
      }
      ind = order(countPerClass, decreasing = F, na.last = T)
      n = numberOfClasses - NumberOfClassesSet
      for (l in 1:n) {
        indc = which(Cls == uniqueClasses[ind[l]])
        Cls[indc] = 99
      }
    } else{
      bool = FALSE
    }
    if (!bool)
      break
    # print(counter)
  }
  names(Cls)=rnames
  return(Cls)
}