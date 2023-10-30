#include "DataBionicSwarm.h"

using namespace std;
#include <RcppArmadilloExtensions/sample.h>

//#ifndef lmC
//#include "lmC.h"
//#endif 

//#ifndef calcStressC
//#include "calcStressC.h"
//#endif 

//#ifndef calcPolarPositionsC
//#include "calcPolarPositionsC.h"
//#endif 

//#ifndef makePolarPositionToroidC
//#include "makePolarPositionToroidC.h"
//#endif 

//#include "calcStressC.h"
//#include "calcPolarPositionsC.h"
//#include "makePolarPositionToroidC.h"

// [[Rcpp::depends(RcppArmadillo)]]
NumericVector sampleC(NumericVector x,double len) {
  bool replace=0;
  return(RcppArmadillo::sample(x,len,replace));
}

// [[Rcpp::depends(RcppArmadillo)]]
NumericMatrix rDistanceToroidC(NumericVector AllDataBotsPosX, NumericVector AllDataBotsPosY,NumericMatrix AllallowedDBPosR0,double Lines,double Columns, NumericVector Nullpunkt,double DBanzahl,NumericMatrix Distances,NumericVector Dx,NumericVector Dy,NumericVector D1,NumericVector D2){
  
  //double DBanzahl=AllDataBotsPosX.length();
  //NumericMatrix Distances(DBanzahl,DBanzahl); //
  //NumericVector Dx(DBanzahl);
  //NumericVector Dy(DBanzahl);
  //NumericVector D1(DBanzahl);
  //NumericVector D2(DBanzahl);
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
// [[Rcpp::export]]
List PswarmCurrentRadiusC2botsPositive(ComplexVector AllDataBotsPosOld,
                                       double Radius,
                                       NumericMatrix DataDists,
                                       ComplexVector IndPossibleDBPosR,
                                       NumericMatrix RadiusPositionsschablone,
                                       NumericVector pp,
                                       NumericVector Nullpunkt,
                                       double Lines,
                                       double Columns,
                                       double nBots,
                                       int limit,
                                       int steigungsverlaufind, 
                                       double StressConstAditiv, 
                                       bool debug){
  // PswarmCurrentRadiusC2botsPositive( AllDataBotsPosOld,
  //                                      Radius, DataDists,
  //                                      IndPossibleDBPosR,
  //                                      RadiusPositionsschablone,  pp,
  //                                      Nullpunkt, Lines,  Columns,
  //                                      nBots,  limit, steigungsverlaufind,  StressConstAditiv)
  //   intern function, do not use yourself
  //   Finds the weak Nash equilibirium for DataBots in one epoch(Radius), requires the setting of constants, grid, and so on in \code{\link{pswarmCpp}}
  // INPUT
  //   AllDataBotsPosOld              ComplexVector [1:n,1], DataBots position in the last Nash-Equlibriuum}
  //   Radius                         double, Radius of payoff function, neighborhood, where other DatsBots can be smelled}
  //   DataDists                      NumericMatrix, Inputdistances[1:n,1:n]}
  //   IndPossibleDBPosR              ComplexVector, see output of \code{\link{findPossiblePositionsCsingle}}}
  //   RadiusPositionsschablone       NumericMatrix, see \code{AllallowedDBPosR0} in \code{\link{setPolarGrid}}}
  //   pp                             NumericVector, number of jumping simultaneously DataBots of one eppoch (per nash-equilibirum), this vector is linearly monotonically decreasing}
  //   Nullpunkt                      NumericVector, equals \code{which(AllallowedDBPosR0==0,arr.ind=T)}, see see \code{AllallowedDBPosR0} in \code{\link{setPolarGrid}}}
  //   Lines                          double, small edge length of rectangulare grid}
  //   Columns                        double, big edge length of rectangulare grid}
  //   nBots                          double, intern constant, equals \code{round(pp[Radius]*DBAnzahl)}}
  //   limit                          int, intern constant, equals \code{ceiling(1/pp[Radius])}}
  //   steigungsverlaufind            int, intern constant}
  //   StressConstAditiv              double, intern constant, sum of payoff of all databots in random condition before the algorithm starts}
  // OUTPUT:
  //   list V of
  //   V$AllDataBotsPos           ComplexVector, indizes of DataBot Positions after a weak Nash equlibrium is found}
  //   V$stressverlauf            NumericVector, intern result, for debugging only}
  //   V$fokussiertlaufind        NumericVector, intern result, for debugging only}
  // 
  // author: Michael Thrun, 04/16
  
  ComplexVector  AllDataBotsPos=clone(AllDataBotsPosOld);
  bool Jumping=1;
  int fokussiertlaufind=0;
  int leng=AllDataBotsPos.length();
  //double DBanzahl=leng;
  NumericVector stress(leng);
  int Iteration=0;
  NumericVector slopeVec(2);
  NumericVector KeyBot(leng); //dummy
  NumericMatrix PosAndStres(leng,2);
  ComplexVector DataBotsPosNeu(leng);
  ComplexVector DataBotsPosNeu2(leng);
  ComplexVector DataBotsPosNeu3(leng);
  ComplexVector DataBotsPosNeu4(leng);
  NumericMatrix OutputDistanceNeu(leng,leng);
  NumericMatrix OutputDistanceNeu2(leng,leng);
  NumericMatrix OutputDistanceNeu3(leng,leng);
  NumericMatrix OutputDistanceNeu4(leng,leng);
  NumericMatrix OutputDistance(leng,leng);
  
  //Initialisierung calcPolarPositions und calcSressC, calcdistancetoroid
  ComplexVector ToroidPosition(nBots);
  double db=-1;
  int nBotsalsInt=nBots;
  
  double DBAnzahl=DataDists.nrow();
  NumericVector Nachbahrschaftsfunktion(DBAnzahl);
  NumericMatrix xxVergleich(DBAnzahl, 3);
  
  NumericMatrix Distances(DBAnzahl,DBAnzahl); //Alternativ muesste in toroiddistances eine deep copy machen!
  NumericMatrix Distances1(DBAnzahl,DBAnzahl);
  NumericMatrix Distances2(DBAnzahl,DBAnzahl);
  NumericMatrix Distances3(DBAnzahl,DBAnzahl);
  NumericMatrix Distances4(DBAnzahl,DBAnzahl);
  
  NumericVector Dx(DBAnzahl);
  NumericVector Dy(DBAnzahl);
  NumericVector D1(DBAnzahl);
  NumericVector D2(DBAnzahl);
  
  //NumericVector CurrentKeyBot(leng);
  NumericVector ChosenForJump(nBots);
  NumericVector AllDataBotsPosReal(DBAnzahl);
  NumericVector AllDataBotsPosImag(DBAnzahl);
  
  NumericVector PossiblePositionsIndi(nBots);
  NumericVector PossiblePositionsIndi2(nBots);
  NumericVector PossiblePositionsIndi3(nBots);
  NumericVector PossiblePositionsIndi4(nBots);
  for(int i=0;i<leng;i++){
    KeyBot(i)=i;
    //CurrentKeyBot(i)=i;
  }
  
  int KeyPossiblePositionLen=IndPossibleDBPosR.length();
  NumericVector KeyPossiblePosition(KeyPossiblePositionLen);
  ComplexVector PossiblePositions(nBots);
  ComplexVector PossiblePositions2(nBots);
  ComplexVector PossiblePositions3(nBots);
  ComplexVector PossiblePositions4(nBots);
  for(int i=0;i<KeyPossiblePositionLen;i++){
    KeyPossiblePosition(i)=i;
  }
  
  NumericVector stressverlauf;
  NumericVector KeySteigung(steigungsverlaufind);
  NumericVector stresstail(steigungsverlaufind);
  double epsilon=0.01;//Steigungsgenauigkeit, steigung nimmt kaum mehr ab
  for(int i=0;i<steigungsverlaufind;i++)
    KeySteigung(i)=i;
  
  while(Jumping){
    fokussiertlaufind=fokussiertlaufind+1;
    //NumericVector AllDataBotsPosReal(leng);
    //NumericVector AllDataBotsPosImag(leng);
    // for(int i=0;i<leng;i++){
    //    AllDataBotsPosReal(i)=AllDataBotsPos(i).r;
    //    AllDataBotsPosImag(i)=AllDataBotsPos(i).i;
    // }
    AllDataBotsPosReal=Re(AllDataBotsPos);
    AllDataBotsPosImag=Im(AllDataBotsPos);
    
    OutputDistance=rDistanceToroidC(AllDataBotsPosReal,AllDataBotsPosImag,RadiusPositionsschablone,Lines,Columns,Nullpunkt,DBAnzahl,Distances,Dx,Dy,D1,D2);
    ////// Bestimme Anzahl Position und Distanz Springender Databots   // - //
    //Zeitfaktor: Je naeher DataBots springen, desto schneller riechen sie erneut entspricht, weniger DBs springen pro Iteration
    
    //nBots=round(0.05*DBAnzahl)//s. Fast and reliable ESOM learning#
    // Sample als ziehen ohne zuruecklegen!!
    //if(CurrentKeyBot.length()<nBots){
    //NumericVector CurrentKeyBot(leng);
    //CurrentKeyBot=clone(KeyBot);
    // }
    ChosenForJump=sampleC(KeyBot,nBots); //15% Wahrscheinlichkeit das der DataBot ueberhaupt versucht zu springen
    
    //NumericVector CurrentKeyBot2=DeleteAll(CurrentKeyBot,ChosenForJump);
    //NumericVector CurrentKeyBot(CurrentKeyBot2.length());
    // CurrentKeyBot=clone(CurrentKeyBot2);
    
    //Normalerweise koennte man direkt aus de Vektor ein sample ziehen, allerdings geht die
    //sample funktion nur fuer numerischeVektoren, hier ist aber ein ComplexerVektor vorhanden
    PossiblePositionsIndi=sampleC(KeyPossiblePosition,nBots);
    PossiblePositionsIndi2=sampleC(KeyPossiblePosition,nBots);
    PossiblePositionsIndi3=sampleC(KeyPossiblePosition,nBots);
    PossiblePositionsIndi4=sampleC(KeyPossiblePosition,nBots);
    for(int i=0;i<nBots;i++){
      PossiblePositions(i)=IndPossibleDBPosR[PossiblePositionsIndi(i)];
      PossiblePositions2(i)=IndPossibleDBPosR[PossiblePositionsIndi2(i)];
      PossiblePositions3(i)=IndPossibleDBPosR[PossiblePositionsIndi3(i)];
      PossiblePositions4(i)=IndPossibleDBPosR[PossiblePositionsIndi4(i)];
    }
    //DataBotsPosNeu=calcPolarPositionsV3(AllDataBotsPos,ChosenForJump,PossiblePositions,Radius,Lines,Columns)
    // Indize abzug von 1 in ChosenForJump, da R von 1 und C++ von 0 zaehlt
    DataBotsPosNeu=calcPolarPositionsC(AllDataBotsPos,ChosenForJump,PossiblePositions,Radius,Lines,Columns, ToroidPosition, db,nBotsalsInt,DataBotsPosNeu); //DataBotsPosNeu als pointen hinten uebergen
    DataBotsPosNeu2=calcPolarPositionsC(AllDataBotsPos,ChosenForJump,PossiblePositions2,Radius,Lines,Columns, ToroidPosition, db,nBotsalsInt,DataBotsPosNeu2);
    DataBotsPosNeu3=calcPolarPositionsC(AllDataBotsPos,ChosenForJump,PossiblePositions3,Radius,Lines,Columns, ToroidPosition, db,nBotsalsInt,DataBotsPosNeu3);
    DataBotsPosNeu4=calcPolarPositionsC(AllDataBotsPos,ChosenForJump,PossiblePositions4,Radius,Lines,Columns, ToroidPosition, db,nBotsalsInt,DataBotsPosNeu4);
    
    OutputDistanceNeu=rDistanceToroidC(Re(DataBotsPosNeu),Im(DataBotsPosNeu),RadiusPositionsschablone,Lines,Columns,Nullpunkt,DBAnzahl,Distances1,Dx,Dy,D1,D2);
    OutputDistanceNeu2=rDistanceToroidC(Re(DataBotsPosNeu2),Im(DataBotsPosNeu2),RadiusPositionsschablone,Lines,Columns,Nullpunkt,DBAnzahl,Distances2,Dx,Dy,D1,D2);
    OutputDistanceNeu3=rDistanceToroidC(Re(DataBotsPosNeu3),Im(DataBotsPosNeu3),RadiusPositionsschablone,Lines,Columns,Nullpunkt,DBAnzahl,Distances3,Dx,Dy,D1,D2);
    OutputDistanceNeu4=rDistanceToroidC(Re(DataBotsPosNeu4),Im(DataBotsPosNeu4),RadiusPositionsschablone,Lines,Columns,Nullpunkt,DBAnzahl,Distances4,Dx,Dy,D1,D2);
    
    //////                                                             // - //
    ////// Vegleiche Phi                                               // - //
    // Achtung: Es werden aber als "Jumps", Alle DBs gezaehlt, welche sich in einer
    // neuen Position wohl fuehlen wuerden, selbst wenn sie keine moeglichkeit haben
    // auf diese neue Position zu springen
    PosAndStres=calcStressC(DataDists,OutputDistance,OutputDistanceNeu,OutputDistanceNeu2,OutputDistanceNeu3,OutputDistanceNeu4,Radius,StressConstAditiv,DBAnzahl,Nachbahrschaftsfunktion, xxVergleich);
    
    int j=0;
    double SumStress=0;
    for(int i=0;i<leng ;i++){
      stress(i)=PosAndStres(i,0);
      if(PosAndStres(i,1)!=-1){ // Edit QMS: Improvement must be decided with -1 versus databot index 0:NAllbots
        if(PosAndStres(i,2)==0) //           Before: improvement was decided with 0 versus databot index 0:NAllbots => databot index 0 never jumps
          AllDataBotsPos(PosAndStres(i,1))=DataBotsPosNeu(PosAndStres(i,1));
        if(PosAndStres(i,2)==1)
          AllDataBotsPos(PosAndStres(i,1))=DataBotsPosNeu2(PosAndStres(i,1));
        if(PosAndStres(i,2)==2)
          AllDataBotsPos(PosAndStres(i,1))=DataBotsPosNeu3(PosAndStres(i,1));
        if(PosAndStres(i,2)==3)
          AllDataBotsPos(PosAndStres(i,1))=DataBotsPosNeu4(PosAndStres(i,1));
        j++;
      }
      if(stress(i)!=R_PosInf){
        SumStress=SumStress+stress(i);
      }
    }
    
    stressverlauf.push_back(SumStress);///sqrt(sqrt(Nomierung)))///(pi*Radius^2)*DBAnzahl)
    Iteration=Iteration+1;
    //////                                                               // - //
    ////// Abbruchbedingung                                              // - //
    if(j==0){
      //Jumping=0; //macht den Algorithmus nur 4mal schneller sonst nichts
    }// end if 
    
    if(Iteration>limit){ //Pruefe ab 20. Iteration in einem Radius
      stresstail=tail(stressverlauf,min(steigungsverlaufind,Iteration)); //Innerhalb einer Radius oder nur letzten Iterationn?
      //slope=lm(formula=stresstail~x)$coefficients[2]
      slopeVec=lmC(KeySteigung,stresstail);
      
      if(slopeVec[1]<epsilon){ //globale Stress nicht mehr besonders zu
        Jumping=0;
        //vielleicht stattdessen random walk mit 1% der dbs?
      } // end if slope>=0
    }// end of Iteration>limit
    //////                                                               // - //
    if(debug){
      if(Iteration%100==(limit+1)){
        //CRAN limits cout,therefore this is depricated
        //std::cout<<"Steigung:"<<slopeVec[1]<<"; Iteration:"<<Iteration<<"; Stress:"<<SumStress<<";"<<std::endl;
      }
    }
  }// end While Jumping
  return Rcpp::List::create(
    Rcpp::Named("AllDataBotsPos")    = AllDataBotsPos,
    Rcpp::Named("stressverlauf")     = stressverlauf,
    Rcpp::Named("fokussiertlaufind") = fokussiertlaufind
  ) ;
}
