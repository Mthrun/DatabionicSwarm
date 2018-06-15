GeneratePswarmVisualization=function(Data,ProjectedPoints,LC,PlotIt=FALSE,ComputeInR=FALSE){
# results=GeneratePswarmVisualization(Data,res$ProjectedPoints,res$LC)
# plotUmatrix(results$Umatrix,results$Bestmatches,Cls)
#Generiert die ESOMneurons(wts) und die Umatrix (umx) fuer einen Schwarm
# Ist der Spezialfall der generalisierten Umatrix, bei festem GridSize
#INPUT 
#  Data[1:n,1:d]
#  ProjectedPoints[1:n,1:2]
#  LC                                       Grid size c(Lines,Columns) of Pswarm
# OPTIONAL
# PlotIt																			=T: plots, =F doesn plot
# ComputeInR                                  =T: Rcode, =F Cpp Code
# Output
# BMUs[1:2,n]                             BestMatchingUnits
# wts                                         ESOMneurons
# umx                                         Umatrix
# GridPoints                                ProjectedPoints on Grid
# LC										Sometimes is better to choose a different grid size, e.g. to to reduce computional effort
#											contrary to SOM, here the grid size defined only the resolution of the visualizations
#											It real grid size is predfined by Pswarm, but you may choose a factor x*res$LC if you so desire.
#											Therefore, The resulting grid size is given back here.
#author MT 03/16
  Data=checkInputDistancesOrData(Data)
  # requireNamespace("Rcpp", quietly = T)
  # requireNamespace("RcppArmadillo", quietly = T)
  # WhereAmI = getSrcDirectory(function(x) {x})
  # path=sub('/R','/src/',WhereAmI)
  # tryCatch({
  #   if(Sys.info()["sysname"]=="Windows"){
  #     if(path=='.'){
  #       path=paste0(SubversionDirectory(),'PUB/dbt/DataBionicSwarm/src/')
  #     }
  #   }else{#lkonnte nur MacOs testen -_-
  #     if(path==""){
  #       path=paste0(SubversionDirectory(),'PUB/dbt/DataBionicSwarm/src/')
  #     }
  #   }
  # },error=function(ex){
  #   warning('Something is wrong with the path, assuming package instead of source files')
  #   #path=paste0(SubversionDirectory(),'PUB/dbt/DataBionicSwarm/src/')
  # })
  # tryCatch({path=suppressWarnings(normalizePath(path))
  # #print(path)
  #   if(ComputeInR){
  #     suppressWarnings(sourceCpp(paste0(path,'Delta3dWeightsC.cpp')))
  #   }else{
  #     suppressWarnings(sourceCpp(paste0(path,'trainstepC.cpp')))
  #   }
  # },error=function(ex){
  #   print('Precompiled package, sourceCPP omitted.')
  # })
  Lines=LC[1]
  Columns=LC[2]
  
if(Lines>Columns){
  warning('Lines should not be bigger then Columns. Rotating Projected Points in a 90 degree angle.')
  Linestemp=Columns
  Columns=Lines
  Lines=Linestemp
  ProjectedPoints=ProjectedPoints[,c(2,1)]
  eps=4
  
}
  

#Der Standardalgorithmus funktioniert ohne Lines und Columns, und liefer nicht
# immer exakt die Voreinstellung des Schwarmes an Lines/Columns
Points=ProjectedPoints2Grid(ProjectedPoints,Lines,Columns)

Lines=LC[1]+1
Columns=LC[2]+1
#############################################################################
##calcUmatrixToroid()
############################################################################
calcUmatrixToroid <- function(EsomNeurons){
  # Umatrix=calcUmatrix(wts)
  # Calculate the Umatrix for given EsomNeurons projection
  # INPUT
  # EsomNeurons[Lines,Columns,weights]		neuronen aus EsomNeurons
  # OPTIONAL
  # Toroid				planar=F
  # OUTPUT
  # Umatrix[Lines,Columns]
  
  #############################################
  ## Nachbarn()
  nachbarn <- function(k, EsomNeurons){
    # INPUT
    # k Gitterpunkt
    # EsomNeurons[Lines,Columns,weights]
    # Toroid
    # OUTPUT
    # nb
    M <- dim(EsomNeurons)[1]
    N <- dim(EsomNeurons)[2]
    pos1 = c(k[1]-1,k[1]-1,k[1]-1,k[1],k[1],k[1]+1,k[1]+1,k[1]+1) %% M
    pos1[which(pos1==0)] = M
    pos2 = c(k[2]-1,k[2],k[2]+1,k[2]-1,k[2]+1,k[2]-1,k[2],k[2]+1) %% N
    pos2[which(pos2==0)] = N
    nb = cbind(pos1,pos2)
    return(nb)
  }
  ############################################
  k = dim(EsomNeurons)[1]
  m = dim(EsomNeurons)[2]
  Umatrix = matrix(0,k,m)
  
  d=dim(EsomNeurons)[3]
  if(is.null(d)){#wts als liste
    stop('EsomNeurons wts has to be an array[1:Lines,1:Columns,1:Weights], use ListAsEsomNeurons')
  }
  
  for(i in 1:k){
    for(j in 1:m){
      nbs=nachbarn(c(i,j),EsomNeurons)
      wij=EsomNeurons[i,j,]
      n.nbs=dim(nbs)[1]
      for(l in 1:n.nbs){
        nij=EsomNeurons[nbs[l,1],nbs[l,2],]
        Umatrix[i,j]=Umatrix[i,j]+sqrt(sum((wij-nij)^2))
      }	
      Umatrix[i,j]=Umatrix[i,j]/n.nbs
    }
  }
  return(Umatrix)
}
#########################################################################
##end calcUmatrixToroid
#########################################################################
rr=round(max(c(Columns,Lines))/10,0)

if(rr<12){
  HeuristischerParameter=12
}else{
  HeuristischerParameter=rr
}
#Nur eine Lineare Transformation
n=nrow(Points)
c=ncol(Points)
if(c==2){
  BMUs=matrix(NaN,nrow=n,ncol=c)
  BMUs=Points[,c(2,1)]
}else if(c==3){
  BMUs=matrix(NaN,nrow=n,ncol=(c-1))
  BMUs=Points[,c(3,2)]
}else{
  stop('Error, wrong number of colums')
}
BMUs[,1]=Lines-BMUs[,1]+1
#BMUs=ProjectedPoints2Bestmatches(Points,Lines)
d=ncol(Data) #NumberOfweights

  rnd=runif(n=d*Lines*Columns, min =min(as.numeric(Data),na.rm = T), max = max(as.numeric(Data),na.rm = T)) #besser als min(data) bis max(data)
  #rnd=max(Data)
  wts<- array(rnd,c(Lines,Columns,d)) #[Lines,Columns,weights]
  print('Initializing sESOM algorithm')
 
  # BestMatches werden festgehalten
  for(i in c(1:nrow(BMUs))){
      wts[BMUs[i,1],BMUs[i,2],] = Data[i,]
  }
  #Jeder Radius sollte min. 1 Eppoche durchlaufen werden, mehr als eine Eppoche fuehrte nicht zu mehr Emergenz
  # s. auch Experimente mit iUmatrix(), wo eine Umatrix als Video pro Eppoche bei diverser Parameterwahl gezeichnet wird
  epochs=HeuristischerParameter
  AnfangsRadius=HeuristischerParameter
  #epochs=20
  #AnfangsRadius=20
vec=pmax(seq(from=AnfangsRadius-1,by=-1,length.out = HeuristischerParameter),1)
  for (i in vec){
    CurrentRadius =  i#max(AnfangsRadius-i,1) #Endradius=1
    #Algorithmus
      wts=sESOM4BMUs(BMUs,Data, wts, toroid=T, CurrentRadius,ComputeInR)
    print(paste0('Operator: getUmatrix4BMUs() at ',round(1-i/HeuristischerParameter,2)*100,'%'))
  } # end 1:epochs
 
  print('Calculating Umatrix')

  Umap=calcUmatrixToroid(wts)
LCnew=c(dim(wts)[1],dim(wts)[2])
if(PlotIt){
  requireNamespace("GeneralizedUmatrix")
  GeneralizedUmatrix::plotTopographicMap(Umap,BMUs)
}
return(list(Bestmatches=BMUs,Umatrix=Umap,WeightsOfNeurons=wts,GridPoints=Points,LC=LCnew))
}