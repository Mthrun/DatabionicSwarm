using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
Rcomplex modComplexC(Rcomplex x,Rcomplex y) {
  // Berechnet den Modulo zweier komplexer Zahlen
  Rcomplex res;
  //res.i=x.i-((int)x.i/y.i)*y.i;
  //res.r=x.r-((int)x.r/y.r)*y.i;
  res.i=(int)abs(x.i)%(int)y.i;
  res.r=(int)abs(x.r)%(int)y.r;
  return res;
}
// [[Rcpp::depends(RcppArmadillo)]]
Rcomplex makePolarPositionToroidC(Rcomplex DBpositionInd,Rcomplex IndPossibleDBPosR,double Lines, double Columns) {
  //INPUT
  // DBpositionInd                 complex number, Current row and column Indize of the DataBot  
  // IndPossibleDBPosR[        complex number of Indice of radius of AllallowedDBPosR0 of one possible DataBots Position in one radial jump  
  // Lines         Integer, hast to be able to be divided by 2
  // Columns       Integer, with Columns>Lines  
  // Optional
  // ToroidPositionInd       complex number of Indice of radius of all one DataBot Position in one radial jump in Packman-Universe
  // Output:
  
  // author: MT 02/2016  
  //Verschiebe alle moeglichen Positionen B bis zum DB Ort A);
  //ComplexVector ToroidPositionInds(n);
  Rcomplex ToroidPositionInd=IndPossibleDBPosR+DBpositionInd;
  
  // IndPossibleDBPosR[,1]=IndPossibleDBPosR[,1]+DBpositionInd[1]//-(Lines/2+1)# And den Ursprung des anderen Bezug-Systemes (BS) anpassen
  // IndPossibleDBPosR[,2]=IndPossibleDBPosR[,2]+DBpositionInd[2]//-(Lines/2+1)# And den Ursprung des anderen Bezug-Systemes anpassen
  //anderes BS, s.  ind=which(AllallowedDBPosR0<=Radius & AllallowedDBPosR0!=0,arr.ind=T) in getOriginPositionsByRadius()
  //  Toroides Feld anpassen
  //Trick: u.U. sind Positionen negativ und werden es mit Modulo positiv angepasst
  //cout <<db<<endl;
  //cout<<DBpositionInd(db)<<endl;
  //cout<<ToroidPositionInds(0)<<endl;
  Rcomplex eins;
  eins.i=1;
  eins.r=1;
  Rcomplex LC;
  LC.i=Columns;
  LC.r=Lines;
  
  //Der Eins Trick hat etwas mit dem Modulo bei negativen Zahlen zu tun, allerdings weis ich nichtmehr wieso
  ToroidPositionInd=modComplexC(ToroidPositionInd-eins,LC)+eins;
  
  //cout<<ToroidPositionInds<<endl;
  //IndPossibleDBPosR[,1]=(IndPossibleDBPosR[,1]-1)%Lines+1
  //IndPossibleDBPosR[,2]=(IndPossibleDBPosR[,2]-1)%Columns+1
  return ToroidPositionInd;
}
