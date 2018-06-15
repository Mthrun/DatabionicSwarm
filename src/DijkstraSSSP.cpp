#include <Rcpp.h>
#include <iostream>
#include <queue>
#include <vector>
#include <climits>

using namespace Rcpp;
using namespace std;
#define INF INT_MAX //Infinity
const int sz=10001; //Maximum possible number of vertices. Preallocating space for DataStructures accordingly

// [[Rcpp::export]]
NumericVector DijkstraSSSP(NumericMatrix Adj, NumericMatrix Costs, int source){
//Dijkstra's SSSP (Single source shortest path) algorithm
// INPUT
// Adj[1:n,1:n]         0/1 adjascency matrix, e.g. from delaunay graph or gabriel graph
// Costs[1:n,1:n]       matrix, distances between n points (normally euclidean)
// source               int, vertice(point) from which to calculate the geodesic distance to all other points
// OUTPUT
// ShortestPaths[1:N]   vector, shortest paths (geodesic) to all other vertices including the source vertice itself
// author: MT 08/16
// inspired by web: http://ideone.com/qkmt31
// require C++11 standard (set flag in Compiler)
//
// Description: gets the shortest path (geodesic distance) from source vertice(point) to all other vertices(points) defined by 
//              the edges of the adjasency matrix
vector<pair<int,double> > a[sz]; //Adjacency list
  double dis[sz]; //Stores shortest distance
  bool vis[sz]={0}; //Determines whether the node has been visited or not
  
  
  int n=Adj.nrow();
  int m=Adj.ncol();
  NumericVector ShortestPaths(n);
  double w;
  double cw; //the final shortest distance for this vertex
  double cwCurrent; // current shortest disane
  for(int i=0;i<n;i++) //Building Graph
    for(int j=0;j<m;j++)
    {
      if(Adj(i,j)==1){
        w=Costs(i,j);
        a[i+1].push_back(make_pair(j+1,w));
        a[j+1].push_back(make_pair(i+1,w));
      }
    }
    //for(int source =1;source<n;source++){
    for(int i=0;i<sz;i++) //Set initial distances to Infinity
      dis[i]=INF;
  
  //Custom Comparator for Determining priority for priority queue (shortest edge comes first)
  class prioritize{
    public: 
      bool operator ()(pair<int, double>&p1 ,pair<int, double>&p2){
        return p1.second>p2.second;
        }
  };
  priority_queue<pair<int,double> ,vector<pair<int,double> >, prioritize> pq; //Priority queue to store vertex,weight pairs
  pq.push(make_pair(source,dis[source]=0)); //Pushing the source with distance from itself as 0
  while(!pq.empty())
  {
    pair<int, double> curr=pq.top(); //Current vertex. The shortest distance for this has been found
    pq.pop();
    int cv=curr.first;
    cw=curr.second; //'cw' the final shortest distance for this vertex
    if(vis[cv]) //If the vertex is already visited, no point in exploring adjacent vertices
      continue;
    vis[cv]=true;
    for(unsigned int i=0;i<a[cv].size();i++) //Iterating through all adjacent vertices
      if(!vis[a[cv][i].first] && a[cv][i].second+cw<dis[a[cv][i].first]) {//If this node is not visited and the current parent node distance+distance from there to this node is shorted than the initial distace set to this node, update it
        cwCurrent=a[cv][i].second;
        pq.push(make_pair(a[cv][i].first,(dis[a[cv][i].first]=cwCurrent+cw))); //Set the new distance and add to priority queue
      }
  }
  
  for(int k=0;k<n;k++){
    dis[k+1]!=INF? ShortestPaths(k)=dis[k+1]:ShortestPaths(k)=INF;
  }
  //}// end source
  return(ShortestPaths);
}