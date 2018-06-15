calcGabrielGraph2D=function(Data,PlotIt=F){
  
  requireNamespace('spdep')
  requireNamespace('sp')
  coords <- sp::coordinates(Data)
  Liste=spdep::gabrielneigh(coords)

  if(PlotIt){
    col.gab.nb<-spdep::graph2nb(Liste, sym=TRUE)
    plot(Data, border="grey")
    plot(col.gab.nb,coords,add=TRUE)
  }
  n=nrow(Data)
  Adj=matrix(0,n,n)
  for(i in seq_along(Liste$from)){
    Adj[Liste$from[i],Liste$to[i]]=1
    Adj[Liste$to[i],Liste$from[i]]=1
  }
  return(list(Delaunay=Adj))
}