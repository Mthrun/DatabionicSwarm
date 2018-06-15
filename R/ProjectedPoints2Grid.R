ProjectedPoints2Grid <- function(ProjectedPoints, Lines, Columns,PlotIt=F){
# ProjectedPoints2Grid(ProjectedPoints,Lines,Columns)  
# Converts numeric ProjectedPoints to integervalues corresponding to a two dimensional rectangular grid
# INPUT
# ProjectedPoints[n,1:2] or [n,1:3] matrix containing the number of n coordinates !(x,y)! of the Projection
#                                   if [n,1:3] the first column has to be the key
#
# Lines                             Default(50), widght of the retangular grid
# Optional
# Columns                           height of the rectangular grid
#                                   We are able to calculate Columns by using Lines and der Range of Data in x and y
# PlotIt                            Difference plottet in new Window, if TRUE
# 
# OUTPUT
# BestMatches[n,Key,Lines,Columns]  Integer GridPositons(1:Lines,1:Colums) corresponding to numerical !(y,x)! coordinates
#                                   can be saved with WriteBM
# author MT 07/2015
# Note: Problem: Machmal  werden verschiedener Punkte werden an die gleiche Gitterstelle transformiert
#       jetzige Loesung: Wenn es zuviele doppelte gibt, einfach Lines vergroessern!
#       mögliche spaetere Loesung: Im speateren Verlauf könnten Voronoi-Nachbahrschaften verwendet werde statt einem Gitter
#                         -> generalisierte Abstrakte Umatrix
  
  
# Behandlungen verschiedener EingabeTypen
if(!is.matrix(ProjectedPoints))
    stop('ProjectedPoints has to be a matrix')
#if(!missing(Columns))
#  if(Lines>Columns)
#    stop('Lines has to be smaller or equal Columns')

n=nrow(ProjectedPoints)
c=ncol(ProjectedPoints)
if(c>3 |c<2)
  stop(paste0('Wrong number of Columns of ProjectedPoints: ',c))

if(c==3){
   coord=ProjectedPoints[,2:3]
   coordtmp=ProjectedPoints[,2:3]
}else{
   coord=ProjectedPoints
   coordtmp=ProjectedPoints
}

# Define Range of Data
minX <- min(coord[,1])
minY <- min(coord[,2])
maxX <- max(coord[,1])
maxY <- max(coord[,2])
RangeX=maxX-minX
RangeY=maxY-minY

# Per Definition: RangeY<RangeX because Linese<Columns, if not rotate
if(RangeY>RangeX){ 
  coord=coord[,c(2,1)]
  drehen=T
  #print(RangeX)
  #print(RangeY)
}else{
  drehen=F
}
# If Columns not chosen, Calculate it respectivly to the range of Data
if(missing(Columns)){
  Columns =round(Lines/RangeY*RangeX,0)
  while(Columns*Lines<4096){ #Empirisch min 2500 Neuronen und Lines<Columns
    Columns=Columns+1
  }
  print(paste('Estimating Columnlength with',Columns))
}
# Make all coordinates positiv.
if(minX < 1){
	coord[,1] <- coord[,1]+(-minX+1) # +1 to get the minimum to 1 (not 0).
}
if(minY < 1){
	coord[,2] <- coord[,2]+(-minY+1) # +1 to get the minimum to 1 (not 0).
}

# Update min and calculate max.

maxX <- max(coord[,1])
maxY <- max(coord[,2])

minX <- min(coord[,1])
minY <- min(coord[,2])


# Subtract min to be sure that the coordinates start with 0 after scaling.
if(minX<=minY){
	rCoord <- coord-minX
} else {
	rCoord <- coord-minY
}

# Update min and calculate max.
minX <- min(rCoord[,1])
maxX <- max(rCoord[,1])
minY <- min(rCoord[,2])
maxY <- max(rCoord[,2])

# Make sure the coordinates with the biggest value is on y-axis (Columns values!)
# to use the maximum of the Grid.
# Divide by maxX or maxY so all coordinates are scaled and are between 0 and 1.
if(maxX>=maxY){
	rCoord <- rCoord/maxX

} else {
	# Make sure the coordinates with the biggest value is on x-axis (Columns values!)
	# to use the maximum of the Grid.
	rCoord <- rCoord[,c(2,1)]
	rCoord <- rCoord/maxY
}

# Update min and calculate max.
minX <- min(rCoord[,1])
maxX <- max(rCoord[,1])
minY <- min(rCoord[,2])
maxY <- max(rCoord[,2])

# Calculate the max value to multiply with (to fit the Grid).
fact1 <- min(Columns-1,Lines-1)/maxY # 
fact2 <- max(Columns-1,Lines-1)/maxX
# choose the smaller value to be sure not to extend the maximum gitter size.

fact <- min(fact1,fact2)
rCoord <- round(rCoord*fact)
# To get the min value from 0 to 1.
rCoord <- rCoord+1

# Update min and calculate max.
minX <- min(rCoord[,1])
maxX <- max(rCoord[,1])
minY <- min(rCoord[,2])
maxY <- max(rCoord[,2])

# Centering the points on the grid
if((Columns-maxX)>=2){ # Otherwise you can't move the x-coordinates.
	rCoord[,1] <- rCoord[,1]+floor((Columns-maxX)/2)
}
if((Lines-maxY)>=2){ # Otherwise you can't move the y-coordinates.
	rCoord[,2] <- rCoord[,2]+floor((Lines-maxY)/2)
}

#Fehlerabfang, das zwei Punkte auf dem selben Gitterpukt leigen, welche vorher
# nicht unique waren
AnzPoints <- nrow(coordtmp)            # soviele punkte in den Daten
dists <- as.matrix(dist(coordtmp))     # Distance with diagonal = 0.
vecdists=dists[upper.tri(dists, diag = F)]
nvorher=AnzPoints-length(which(vecdists==0)) #number of duplicates

AnzPoints <- nrow(rCoord)            # soviele punkte in den Daten
dists <- as.matrix(dist(rCoord))     # Distance with diagonal = 0.
vecdists=dists[upper.tri(dists, diag = F)]
nnachher=AnzPoints-length(which(vecdists==0))  #number of duplicates

#nvorher=nrow(uniquePoints(coordtmp)$unique)
#nnachher=nrow(uniquePoints(rCoord)$unique)
if(nvorher!=nnachher){
  print(paste('Different Points are now on the same grid position. Sum of unique points (before) (after):', nvorher,nnachher))
}

# If rotatet, rotate back
if(drehen)
  BestMatches=cbind(1:n,rCoord[,c(2:1)]) #umgekehrt Definiert wie x,y koordinaten
else
  BestMatches=cbind(1:n,rCoord[,c(1:2)]) #umgekehrt Definiert wie x,y koordina

if(PlotIt){
  #windows()
  par(mfrow = c(1,2))
  plotSwarm(coordtmp,main='ProjectedPoints')
  plotSwarm(BestMatches[,2:3],main='Corresponding BestMatches')
}
return(BestMatches)

} 