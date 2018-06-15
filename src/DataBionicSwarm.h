#ifndef DATABIONICSWARM_H
#define DATABIONICSWARM_H
#include <RcppArmadillo.h>

using namespace Rcpp;

arma::cube Delta3DWeightsC(arma::cube x,Rcpp::NumericVector Datasample);
double stress4BotC(NumericVector DistInput,NumericVector DistOutput,double Radius,NumericVector Nachbahrschaftsfunktion, double DBAnzahl, double StressConstAditiv);
ComplexVector calcPolarPositionsC(ComplexVector DataBotsPos,NumericVector ChosenForJump,ComplexVector PossiblePositions,double Radius,double Lines,double Columns,ComplexVector ToroidPosition, double db, int n,ComplexVector DataBotsPosNeu);
NumericVector lmC(NumericVector x,NumericVector yr);
Rcomplex modComplexC(Rcomplex x,Rcomplex y);
Rcomplex makePolarPositionToroidC(Rcomplex DBpositionInd,Rcomplex IndPossibleDBPosR,double Lines, double Columns);
NumericMatrix calcStressC(NumericMatrix DataDists,NumericMatrix OutputDistance,NumericMatrix OutputDistanceNeu,NumericMatrix OutputDistanceNeu2,NumericMatrix OutputDistanceNeu3,NumericMatrix OutputDistanceNeu4,double Radius,double StressConstAditiv, double DBAnzahl,NumericVector Nachbahrschaftsfunktion,NumericMatrix xx);


#endif
