GeneratePmatrix=function(Data,Weights,Lines,Columns,Radius=NULL,PlotIt=FALSE) {
  #  PMatrix =  GeneratePmatrix(Data,EsomNeurons,Radius,Umatrix,PlotIt)
  #INPUT
  # Data[1:n,1:d]                      A \code{[n,k]} matrix containing the data
  # EsomNeurons[1:Lines,1:Columns,1:weights]    Information stored as a List of weights in a 2D matrix,
  #                                   oder das Rlistenobjekt aus ReadEsomNeurons
  # Radius                            The radius for measuring the density within the hypersphere
  # Umatrix[1:Lines,1:Columns]        Matrix with U-Matrix Height values
  #Optional
  # PlotHist
  # Weights[1:Lines*1:Columns,1:weights]    If EsomNeurons is missing, Information stored as a List of weights
  # Lines
  # Columns
  
  #OUTPUT
  # UstarMatrix[1:Lines,1:Columns]
  #
  # MT 03/2015 aus matlab uebernommen
  # MT01/17: Rollback von vorheriger Version, Neue Arguemente nach hinten verschoben
  # FL01/17: Parameter EsomNeurons entfernt. Alte Aufrufe muessen gegebenenfalls zu Weights konvertiert werden
  
  x=as.matrix(dist(Data))
  
  if(missing(Lines)|missing(Columns))
    stop("Lines and Columns are necessary")
  if(missing(Weights))
    stop("Weights are necessary")
  
  EsomNeurons=GeneralizedUmatrix::ListAsEsomNeurons(Weights,Lines,Columns)
  
  
  if(is.null(Radius)){
  
  x=x[lower.tri(x, diag = FALSE)]
  par=quantile(x,c(0.2)) #geschaetzter paretorRadius
  xx=ABCanalysis::ABCRemoveSmallYields(x,0.5)
  x=xx$SubstantialData
  res=suppressWarnings(ABCanalysis::ABCanalysis(x))
  Radius=1/(min(x[res$Aind])/max(x[res$Cind]))*par  #Verhaeltnis vermutliche inner/Inter Clusterdistanz
  #print(min(x[res$Aind])/max(x[res$Cind]))
  #print(par)
  print(Radius)
  }
  #dimensions=dim(EsomNeurons)
  #AnzEsomNeurons=dimensions[1]
  #AnzVariablen=dimensions[2]
  
  d=dim(EsomNeurons)[3]
  UmatrixLines=dim(EsomNeurons)[1]
  UmatrixCols=dim(EsomNeurons)[2]
  
  if(is.null(d)){#EsomNeurons als liste
    stop('esom EsomNeurons has to be an array[1:Lines,1:Columns,1:Weights], use ListAsEsomNeurons')
  }
  
  AnzData=nrow(Data)
  AnzVariablen=ncol(Data)
  
  norm=2
  distances=array(0,c(UmatrixLines,UmatrixCols,AnzData))
  #EsomNeuronsAnzInKugel=vector(mode='numeric',AnzEsomNeurons)
  PMatrix <- matrix(0,UmatrixLines,UmatrixCols)
  for ( i in 1:UmatrixLines ){
    for (j in 1:UmatrixCols){
      aux <- t(Data) - EsomNeurons[i,j,] # columns
      distances[i,j,]= sqrt(colSums(aux^norm)) #quadrierte Differenzen zu  EsomNeurons(i)
  		#EsomNeuronsAnzInKugel[i]=sum(distances[i,j,]<=Radius)
      PMatrix[i,j] <- sum(distances[i,j,] <= Radius)
    }
  }
  
  ##########################################
  #ParetoRadius bestimmung alternative
  # while(!
  # median(EsomNeuronsAnzInKugel)==prctile(DistanceMatrix(Data),20))
  # {
  #   suche neuen Radius
  # }
  ##########################################
  
  #PMatrix <- matrix(0,UmatrixLines,UmatrixCols)
  # for ( i in 1:UmatrixLines ){
  #   for (j in 1:UmatrixCols){
  #     PMatrix[i,j] <- sum(distances[i,j,] <= Radius)
  #   }
  # }
  
  #PMatrix =matrix(EsomNeuronsAnzInKugel,UmatrixLines,UmatrixCols)
  # PMatrix =matrix(0,UmatrixLines,UmatrixCols)
  #
  # for( i in 1:UmatrixLines ){
  #    for (j in 1:UmatrixCols){
  #       PMatrix[i,j]=EsomNeuronsAnzInKugel[i]
  #    }
  # }
  
  #for(i in 1:AnzEsomNeurons){
  #    X = EsomNeurons[i,]
      #matlab: cumsum(((Data-ones(AnzData,1)*x).^2),2) #rowCumsums
      #x=(Data-t(matrix(x,AnzVariablen,AnzData)))^2
      #Dist2x=t(apply(x, 1, cumsum)) #rowCumsums
  
  		###folgende neue und reine Uebersetzung aus Matlab funktioniert nicht
  #    Dist2x = rowSums(((Data-X)^2)) # Summe der quadrierte Differenzen  zu  x
  #    Dist2x = sqrt(Dist2x)
  #    EsomNeuronsAnzInKugel[i] = sum(Dist2x<=Radius)
  		####
      #Dist2x = rowCumsums(((Data-t(matrix(x,AnzVariablen,AnzData)))^2)) # quadrierte Differenzen zu w(i)
      #matlab: sqrt(Dist2x(:,end))
     # Dist2x=sqrt(t(tail(t(Dist2x),1)))
      #AnzInKugelI[i] = sum(Dist2x<=Radius)
  ## P-Matrix Berechnen
  
  #B = reshape(A,m,n) returns the m-by-n matrix B whose elements are taken column-wise from A.
  #matlab: PMatrix  = (reshape(AnzInKugelI,UmatrixCols,UmatrixLines)');
  #PMatrix=t(matrix(EsomNeuronsAnzInKugel,nrow=UmatrixCols,ncol=UmatrixLines))
  
  ###
  #PMatrix=t(matrix(EsomNeuronsAnzInKugel,nrow=UmatrixCols,ncol=UmatrixLines))
  ###
  if(PlotIt){
    if(requireNamespace("ggplot2", quietly = TRUE)){
      p <- GeneralizedUmatrix::plotTopographicMap(PMatrix, Colormap=DataVisualizations::PmatrixColormap, Tiled=TRUE) +
        ggplot2::ggtitle('P-Matrix')
    }else{
      p <- GeneralizedUmatrix::plotTopographicMap(PMatrix, Colormap=DataVisualizations::PmatrixColormap, Tiled=TRUE)
    }
    print(p)
  #	optNrOfBins = OptimalNoBins(EsomNeuronsAnzInKugel)
  #	minData = min(EsomNeuronsAnzInKugel,na.rm = TRUE)
  #	maxData = max(EsomNeuronsAnzInKugel,na.rm = TRUE)
  #	i = maxData-minData
  #	optBreaks = seq(minData, maxData, i/optNrOfBins) # bins in fixed intervals
  #	hist(EsomNeuronsAnzInKugel, breaks=optBreaks,xlab='# in spheres',main='distribution of # in spheres')
  } # if
  #plotPmxTopView(PMatrix,tiling=T)
  return(PMatrix)
}
