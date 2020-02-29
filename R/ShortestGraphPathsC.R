ShortestGraphPathsC =function(Adj,Cost){
# res=ShortestGraphPathsC(Adj, Cost) 
# Shortest GraphPaths = geodesic distances
# #INPUT
# Adj
# 1:n,1:n        0/1 adjascency matrix, e.g. from delaunay graph or gabriel graph
# 
# Cost
# [1:n,1:n]       matrix, distances between n points (normally euclidean)
# 
# OUTPUT
# ShortestPaths[1:n,1:n]   
# vector, shortest paths (geodesic) to all other vertices including the source vertice itself
# from al vertices to all vertices, stored as a matrix
#
# author: Michael Thrun 08/2016
#  Dijkstra's SSSP (Single source shortest path) algorithm, from all points to all points
# Vertices are the points, edges have the costs defined by weights (normally a distance)
# Dijkstra, E. W.: A note on two problems in connexion with graphs, Numerische mathematik, Vol. 1(1), pp. 269-271. 1959.


  #Fehlerabfang
  if(is.list(Adj)) stop('Adj is a list and not a matrix.')
  if(is.list(Cost)) stop('Cost is a list and not a matrix.')
   
  if(!is.matrix(Adj)) {
    warning('Adj is not a matrix, as.matrix() is called.')
    Adj=as.matrix(Adj)
  }
  if(!is.matrix(Cost)) {
    warning('Cost is not a matrix, as.matrix() is called.')
    Cost=as.matrix(Cost)
  }
  
  n=nrow(Adj)
  if(n!=ncol(Adj)) stop('Adj hast not equal number of rows and colums.')
  if(nrow(Cost)!=ncol(Cost)) stop('Cost hast not equal number of rows and colums.')
  
  if(!isSymmetric(unname(Adj))) stop('Adj is not symmetric, maybe a directed instead of an undirected graph was used?')
  
  if(!isSymmetric(unname(Cost))) warning('Cost is not symmetric.')
  
  if(n!=nrow(Cost)) stop('Adj rows does not equal Cost rows.')
  
  if(!is.numeric(Adj)) {
    warning('Adj is not numeric')
  }
  if(!is.numeric(Cost)) {
    warning('Cost is not numeric.')
  }

  Dists=matrix(NaN,n,n)
  for(i in 1:nrow(Adj))
    Dists[i,]=DijkstraSSSP(Adj,Cost,i) #inspiriert durch http://ideone.com/qkmt31
  
  return(Dists)
}