#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace sugar;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix rDistanceToroidCsingle(NumericVector AllDataBotsPosX, NumericVector AllDataBotsPosY,NumericMatrix AllallowedDBPosR0,double Lines,double Columns, NumericVector Nullpunkt){
// Fuer aufrufe in PswarmCPP (einmalig) und pswarm() damit es sich nicht mit der internen CPP fkt von PswarmCPP ueberschneidet
  
    double DBanzahl=AllDataBotsPosX.length();
    NumericMatrix Distances(DBanzahl,DBanzahl);
    NumericVector Dx(DBanzahl);
    NumericVector Dy(DBanzahl);
    NumericVector D1(DBanzahl);
    NumericVector D2(DBanzahl);
    for(int i=0;i<DBanzahl;i++){
      Dx=abs(AllDataBotsPosX-AllDataBotsPosX[i]);
      Dy=abs(AllDataBotsPosY-AllDataBotsPosY[i]);  

      D1=Lines-Dx+1;
      D2=Columns-Dy+1;
      // -1 da im Gegensatz zu R die MAtrixenindizes bei 0 beginnen!
      Dx = pmin(Dx,D1)+Nullpunkt[0]-1; //Toroid machen und an Nullpunkt der Schablone anpassen
      Dy = pmin(Dy,D2)+Nullpunkt[1]-1; //Toroid machen und an Nullpunkt der Schablone anpassen

      for(int j=0;j<DBanzahl;j++){
        Distances(i,j)=AllallowedDBPosR0(Dx[j],Dy[j]);
      }
    }
    return Distances;
}
