using namespace Rcpp;

using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
ComplexVector calcPolarPositionsC(ComplexVector DataBotsPos,NumericVector ChosenForJump,ComplexVector PossiblePositions,double Radius,double Lines,double Columns,ComplexVector ToroidPosition, double db, int n,ComplexVector DataBotsPosNeu){
  // calcPolarPositionsGauss(DataBotsPos,RadiusPositionsschablone,Radius,Lines,Columns)
  // Position random generation from the normal distribution in an two dimensional toroid grid defined by Lines and Columns
  //INPUT
  // DataBotsPos[1:AnzData,2]        complex vector of Two Indizes per Databot such that getCartesianCoordinates(DataBotPos,GridRadius,GridAngle) gets the cartesian positions on the grid
  // ChosenForJump                              numeric vector of DataBots, which where chosen to jump, s. Fast and reliable ESOM learning
  // PossiblePositions[m,2]                    complex vector of two indizes of possible jum positions
  // Radius                                        Jump radius of DataBot, around this DataBot, for each DataBots equal
  // Lines                           Integer, hast to be able to be divided by 2, sets Size of planar grid
  // Columns                         Integer, with Columns>Lines  sets Size of planar grid
  // Output:
  // DataBotsPos[1:AnzData]          random new Indizes respective to polar dataPoints positions on free grid places within a radius
  // author: MT 02/2016
  
  //ComplexVector ClosedPositions=clone(DataBotsPos) ;//Hier kommen alle Positionen rein, welche belegt sind, d.h. dorthin darf nicht gesprungen werdem
  //DataBotsPosNeu=clone(DataBotsPos); //deep copy
//ist anscheinend schneller als clone
    std::copy( DataBotsPos.begin(), DataBotsPos.end(), DataBotsPosNeu.begin() ) ;
  //IntegerVector DataBotsPosRealInt=as<IntegerVector>(DataBotsPosReal);
  //IntegerVector DataBotsPosImagInt=as<IntegerVector>(DataBotsPosImag);
  
  // deep copy, sonst diverse fehler
 // ComplexVector DataBotsPosClone = clone(DataBotsPos);
  // Muss abhaengig vom gesamtradius sein => bei kleinen Radien stehen viel weniger meogliche positionen zur verfuegung
  //int n=ChosenForJump.length();
  ///// wird nun einmalig vorher initialisiert
  //ComplexVector ToroidPosition(n);
  //double db;
  ////
  for(int i=0;i<n;i++){
    db=ChosenForJump[i];//R zaehlt um 1 anders wie C++
    Rcomplex PosOld=DataBotsPos[db];
    Rcomplex PosNew=PossiblePositions[i];
    //Nullpunkt Verschiebung der erlaubten Indizes bezueglich eines polaren gitters von Mitte zum linken unteren Ecke
    // die indizes dieser erlaubten positionen muessen nun an ein toroides gitter angepasst werden
    //ToroidPositionIndsCurrent=makePolarPositionsToroid(cbind(Re(DataBotsPos[db]),Im(DataBotsPos[db])),cbind(Re(IndPossibleDBPosR),Im(IndPossibleDBPosR)),Lines,Columns)
    ////ToroidPositionIndsCurrent=makePolarPositionsToroidC(DataBotsPos,db-1,PossiblePositions[i],Lines,Columns);
    Rcomplex ToroidPositionIndsCurrent= makePolarPositionToroidC(PosOld,PosNew,Lines,Columns);
    // besetzte Positionen muessen entfernt werden
    //OpenPositions=setdiffComplexVectors(ToroidPositionIndsCurrent,ClosedPositions)
    ;
    ToroidPosition[i]=ToroidPositionIndsCurrent;
  }
  //aehnlich zu R setdiff(a,b)
  int i;int j;
  IntegerVector Res(ToroidPosition.length());
  for(i=0; i < ToroidPosition.size();i++){
    for(j=0; j < DataBotsPos.size();j++){
      if((ToroidPosition[i].i == DataBotsPos[j].i) &&
         (ToroidPosition[i].r == DataBotsPos[j].r))
        Res[i] = 1;
    }
  }
  //Achtung es koennten theoretisch 2 DataBots auf die gleiche neue Position springen, wenn
  // a) die Gleiche Position in 2 verscheidenen samplen gewuerfelt worden ist
  // b) diese Position fuer beide DataBots einen besseren Stress hat als jeweils 3 andere Positionen
  for(int i=0;i<n;i++){
    db=ChosenForJump[i];//R zaehlt um 1 anders wie C++
   // ComplexVector tzz(1);
   // tzz=ToroidPosition[i];
    //ComplexVector OpenPositions=setdiffComplexVectorC(tzz,DataBotsPosRealInt,DataBotsPosImagInt);
    //aus den offenen positionen wird eine position gezogen und die outputvarible dort veraendert
    //int lengthDbsindtest=OpenPositions.length();
    if(Res[i]==0){
      //DataBotsPosClone[db]=OpenPositions[0];
      
      DataBotsPosNeu[db]=ToroidPosition[i];//DataBot darf weder auf alte Positions eines anderen DBs springen
      //ClosedPositions.push_back(OpenPositions[0]);
    } // end iflengthDbsindtest>0
  } // end for db in vec
  return(DataBotsPosNeu);
}
