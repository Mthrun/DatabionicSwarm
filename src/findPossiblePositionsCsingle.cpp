#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace sugar;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
ComplexVector findPossiblePositionsCsingle(NumericMatrix RadiusPositionsschablone,double jumplength,double alpha,double Lines){
  //erlaubte sprungpositionen befinden sich innerhalb eines streifens um den gewuerfelten radius
// Fuer aufrufe in PswarmCPP (einmalig) und pswarm() damit es sich nicht mit der internen CPP fkt von PswarmCPP ueberschneidet
  // author: MT 02/2016  
  // Note in R ware dies
  //erlaubePositionenaufSchablone=which(RadiusPositionsschablone<=(jumplength[i]+alpha)&RadiusPositionsschablone>(jumplength[i]-alpha),arr.ind=T)
  //Nullpunkt Verschiebung der erlaubten Indizes bezueglich eines polaren gitters von Mitte zum linken unteren Ecke
  //IndPossibleDBPosR=erlaubePositionenaufSchablone-(Lines/2+1)
    
  int n=RadiusPositionsschablone.nrow();
  int m=RadiusPositionsschablone.ncol();
  //double max =jumplength+alpha;
  //double min =jumplength-alpha;
  ComplexVector OpenPositions;
  Rcomplex tmp;
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      // Wenn Radius innerhalb von Ringbreite, speicher die #Position
      if((RadiusPositionsschablone(i,j)<=jumplength)){
      //if((RadiusPositionsschablone(i,j)<=max)&(RadiusPositionsschablone(i,j)>min)){
//Nullpunkt Verschiebung der erlaubten Indizes bezueglich eines polaren gitters von Mitte zum linken unteren Ecke
        // In R wird ab 1 gezaehlt
        //tmp.r=i-(Lines/2+1); 
        //tmp.i=j-(Lines/2+1);
        // In C++ ab 0 deswegen
        tmp.r=i-(Lines/2);  // Lines muss durch 2 teilbar sein, das wird aber in pswarm() schon abgefangen
        tmp.i=j-(Lines/2);
        OpenPositions.push_back(tmp);
      }
    }
    
  }
  return(OpenPositions);
} 
