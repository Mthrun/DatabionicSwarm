DBSclustering=function(k,DataOrDistance,BestMatches,LC,StructureType=TRUE,
                       PlotIt=FALSE,ylab,main,method='euclidean',...){
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
  DataOrDistance=checkInputDistancesOrData(DataOrDistance,funname='DBSclustering')
  if(missing(ylab)){
    ylab="Ultrametric Portion of Distance"
  }
  if(missing(main)){
    if(isTRUE(StructureType)){
      StructureTypeStr="Compact"
    }else{
      StructureTypeStr="Connected"
  }
  main=paste0(StructureTypeStr," DBS Clustering")
  }
  if(isSymmetric(unname(DataOrDistance))){
    InputD = DataOrDistance
    rnames=1:nrow(DataOrDistance)
  } else{
    if(!is.null(rownames(DataOrDistance)))
      rnames=rownames(DataOrDistance)
    else
      rnames=1:nrow(DataOrDistance)
    if(!requireNamespace('parallelDist')){
      warning("DBSclustering: package parallelDist is not installed, falling back to dist().")
      InputD = as.matrix(dist(DataOrDistance, method = method))
    }
    if(!PlotIt)
       InputD = as.matrix(parallelDist::parDist(DataOrDistance, method = method,...))
    else
      InputD = as.matrix(parallelDist::parDist(DataOrDistance, method = method))
  }# end if(isSymmetric(DataOrDists))
  
  GabrielGraph = FALSE                                                          # Gabriel graph immer schlechter...
  GOutput      = Delaunay4Points(Points   = BestMatches, LC = LC[c(2,1)],
                                 IsToroid = TRUE, PlotIt = FALSE,
                                 Gabriel = GabrielGraph)
  if(!is.matrix(GOutput)){
    GOutput = GOutput[[1]]
  }
  Dist=ShortestGraphPathsC(GOutput,InputD)
  if(StructureType){
    pDist=as.dist(Dist)
    hc <- hclust(pDist,method="ward.D")
  }else{
    ind=which(GOutput==0,arr.ind=T)
    Dist2=Dist*GOutput
    Dist2[ind]=max(Dist)*2
    pDist=as.dist(Dist2)
    hc <- hclust(pDist,method="single")
  }

  if(k>1){
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
  }else{#no clustering was performed
    Cls=rep(1,length(rnames))
  }
  names(Cls)=rnames
  
  if(PlotIt){
    x=as.dendrogram(hc)
    if(requireNamespace('dendextend',quietly = TRUE)){
      #what is the ordering of the cluster in dendrogram
      # from left to right
      Cls_tmp=Cls[order.dendrogram(x)]
      #count frequency in that ordering
      uniqueClasses <- unique(Cls_tmp)
      numberOfClasses <- length(uniqueClasses)
      
      #countPerClass <- rep(0, numberOfClasses)
      countPerClass=list()
      for (i in uniqueClasses) {
        inClassI <- sum(Cls_tmp == uniqueClasses[i])
        countPerClass[[i]] = inClassI
      }
      names(countPerClass)=uniqueClasses
      countPerClass=unlist(countPerClass)
      #get the right number of colors
      if(numberOfClasses>1)
        cols=ProjectionBasedClustering::DefaultColorSequence[1:numberOfClasses]
      else
        cols="black"
      #what would be the ordering of datra based on frequency
      data_order=order(countPerClass,decreasing = TRUE) #from highest frequency
      #what would be the orders of the branches
      unique_reordered=uniqueClasses[data_order]
      # fit that order to the colors
      cols_order = match(table = unique_reordered,uniqueClasses)
      cols=cols[cols_order]
      #branch colors with specific set of colors based on cluster frequency
      x=dendextend::set(x,"branches_k_color", k = k,cols)
    }
    plot(x, main=main,xlab="No. of Data Points N", ylab=ylab,sub=" ",leaflab ="none", ...)
    axis(1,col="black",las=1)
  }
  
  return(Cls)
}