#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
double stress4BotC(NumericVector DistInput,NumericVector DistOutput,double Radius,NumericVector Nachbahrschaftsfunktion, double DBAnzahl, double StressConstAditiv){
  // stress=stresskriterium(DistsInput,DistsOutput,Radius)
  // Stresskriterium fuer den sprung eines Databots berechnet fuer alle Databots seperat
  //
  // INPUT
  // DistsInput       vector: Pairwise distance between pairs of objects
  // DistsOutput      vector of polar projected points: Pairwise distance between pairs of objects
  // Radius          Radius der Umgebung fuer Nachbahrschaftsfunktion 
  //
  // OUTPUT
  // stress         One value
  //  
  // Autor: MT 01/2015
  // 1. Editor: CL 01/15
  // 2. Editor: MT 02/15
  
  
  //Nachbahrschaftsfunktion=1-DistOutput*DistOutput/(3.14159265*Radius*Radius);
  // double x;
  // double N=0;
  // for(int i=0;i<DBAnzahl;i++){
  //   x=1-DistOutput[i]*DistOutput[i]/(3.14159265*Radius*Radius);
  //   if(x<0){
  //     Nachbahrschaftsfunktion[i]=0;
  //   }else{
  //     Nachbahrschaftsfunktion[i]=x;
  //     N=N+x;
  //   }
  // }
  Nachbahrschaftsfunktion=1-(DistOutput*DistOutput)/(3.14159265*Radius*Radius);
  //Nachbahrschaftsfunktion[Nachbahrschaftsfunktion<0]=0;
  for(int i=0;i<DBAnzahl;i++){
    if(Nachbahrschaftsfunktion[i]<0)
      Nachbahrschaftsfunktion[i]=0;
  }
  //double N=std::accumulate(Nachbahrschaftsfunktion.begin(),Nachbahrschaftsfunktion.end());
  
  double N=sum(Nachbahrschaftsfunktion);
  
  if(N<=0.0000001){return(StressConstAditiv);}
  
  
  return(StressConstAditiv-sum(Nachbahrschaftsfunktion*DistInput)/N);
  
}

double vecminInd(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference 
  return it - x.begin();
}

double vecmaxInd(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::max_element(x.begin(), x.end());
  // we want the value so dereference 
  return it - x.begin();
}
// [[Rcpp::depends(RcppArmadillo)]]
NumericMatrix calcStressC(NumericMatrix DataDists,NumericMatrix OutputDistance,NumericMatrix OutputDistanceNeu,NumericMatrix OutputDistanceNeu2,NumericMatrix OutputDistanceNeu3,NumericMatrix OutputDistanceNeu4,double Radius,double StressConstAditiv, double DBAnzahl,NumericVector Nachbahrschaftsfunktion,NumericMatrix xx){
  
  //double DBAnzahl=DataDists.nrow();
  //NumericVector Nachbahrschaftsfunktion(DBAnzahl);
  //NumericMatrix xx(DBAnzahl, 3);
  double phi;
  //double phiNeu;
  //double phiNeu2;
  double phiTest;
  //double pos;
  for(int db =0;db<DBAnzahl;db++){  //espringen nur 15% der DataBots durch BotsJumping
    // aber es wird trotzdem fuer alle der stress neu berechnet//
    //Achtung Distanz zu sich selber auslassen!
    // phi=stress4Bot(DataDists[db,-db],OutputDistance[db,-db],Radius)  
    // phiNeu=stress4Bot(DataDists[db,-db],OutputDistanceNeu[db,-db],Radius)
    
    
    phi=stress4BotC(DataDists.row(db),OutputDistance.row(db),Radius,Nachbahrschaftsfunktion,DBAnzahl,StressConstAditiv) ; 
    NumericVector PhiNeu(4);
    PhiNeu[0]=stress4BotC(DataDists.row(db),OutputDistanceNeu.row(db),Radius,Nachbahrschaftsfunktion,DBAnzahl,StressConstAditiv);
    PhiNeu[1]=stress4BotC(DataDists.row(db),OutputDistanceNeu2.row(db),Radius,Nachbahrschaftsfunktion,DBAnzahl,StressConstAditiv);
    PhiNeu[2]=stress4BotC(DataDists.row(db),OutputDistanceNeu3.row(db),Radius,Nachbahrschaftsfunktion,DBAnzahl,StressConstAditiv);
    PhiNeu[3]=stress4BotC(DataDists.row(db),OutputDistanceNeu4.row(db),Radius,Nachbahrschaftsfunktion,DBAnzahl,StressConstAditiv);

      double ind=vecmaxInd(PhiNeu);
        // if(phiNeu<phiNeu2){
        //   phiTest=phiNeu2;
        //   pos=2;
        // }else{
        //     phiTest=phiNeu;
        // pos=1;
      // }  
      phiTest=PhiNeu[ind];
    if(phiTest>phi)
    {
      xx(db,2)=ind;
      xx(db,1)=db;
      xx(db,0)=phiTest;
      
    } //end if phiNeu<phi
    else
    {
      xx(db,2)=0;
      xx(db,1)=-1; // Edit QMS: No improvement => -1 instead of 0, so that databot index 0 implies improvement of databot 0
      xx(db,0)=phi;
      
    }
    //xx(db,2)=phi;
    //xx(db,3)=phiNeu;
  } // end for 1:DBAnzahl
  return(xx);
}
