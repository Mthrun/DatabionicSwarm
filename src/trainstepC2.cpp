#include <Rcpp.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppParallel)]]
struct Delta3DWeightsC : public Worker {    // Worker for parallelization
  // inputs to read from
  //const RVector<double> esom;
  const RVector<double> DataSample;
  const int Lines;
  const int Columns;
  const int Weights;
  const RVector<double> neighmatrix;
  
  // output to write to
  RVector<double> esom;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to form the Rcpp matrix type)
  Delta3DWeightsC(const NumericVector DataSample,
                  const int Lines,
                  const int Columns,
                  const int Weights,
                  const NumericVector neighmatrix,
                  NumericVector esom):
    DataSample(DataSample),
    Lines(Lines),
    Columns(Columns),
    Weights(Weights),
    neighmatrix(neighmatrix),
    esom(esom) {}
  // function call operator that work for the specified range (begin/end)    esomwts = esomwts - (neighmatrix * inputdiff);
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t k = begin; k < end; k++){
      for(int j = 0; j < Columns; j++){
        for(int i = 0; i < Weights; i++){
          int tmpIdx1 = i * Columns * Lines + j * Lines + k;
          int tmpIdx2 = j * Lines + k;
          esom[tmpIdx1] = esom[tmpIdx1] - (neighmatrix[tmpIdx2] * (esom[tmpIdx1] - DataSample[i]));
        }
      }
    }
  }
};


// [[Rcpp::depends(RcppParallel)]]
NumericVector RcppParallelDelta3DWeights(NumericVector esom,
                                         NumericVector DataSample,
                                         NumericVector neighmatrix,
                                         int Lines, int Columns, int Weights){
  //NumericVector inputdiff(esom);
  Delta3DWeightsC delta3DWeightsC(DataSample,                   // create the worker
                                  Lines,
                                  Columns,
                                  Weights,
                                  neighmatrix,
                                  esom);
  parallelFor(0, Lines, delta3DWeightsC);                           // call it with parallelFor
  return esom;
}



// [[Rcpp::depends(RcppParallel)]]
struct ToroidDistance : public Worker {    // Worker for parallelization
  // inputs to read from
  const RVector<double> aux;
  const RMatrix<double> kmatrix;
  const RMatrix<double> mmatrix;
  const RMatrix<double> bm1;
  const RMatrix<double> bm2;
  const int Lines;
  const int Columns;
  const int LCS;
  
  // output to write to
  RMatrix<double> OutputDistances;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to form the Rcpp matrix type)
  ToroidDistance(const NumericVector aux,
                 const NumericMatrix kmatrix,
                 const NumericMatrix mmatrix,
                 const NumericMatrix bm1,
                 const NumericMatrix bm2,
                 const int Lines,
                 const int Columns,
                 const int LCS,
                 NumericMatrix OutputDistances):
    aux(aux),
    kmatrix(kmatrix),
    mmatrix(mmatrix),
    bm1(bm1),
    bm2(bm2),
    Lines(Lines),
    Columns(Columns),
    LCS(LCS),
    OutputDistances(OutputDistances) {}
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      for(int j = 0; j < Columns; j++){
        int auxIdx1 = j*Lines + i;
        int auxIdx2 = LCS + j*Lines + i;
        double FirstPart = 0.5*sqrt(pow(kmatrix(i,j) - abs(2 * abs(aux[auxIdx1] - bm1(i,j)) - kmatrix(i,j)), 2));
        double SecondPart = 0.5*sqrt(pow(mmatrix(i,j) - abs(2 * abs(aux[auxIdx2] - bm2(i,j)) - mmatrix(i,j)), 2));
        OutputDistances(i,j) = FirstPart + SecondPart;
      }
    }
  }
};


// [[Rcpp::depends(RcppParallel)]]
NumericMatrix RcppParallelToroidDistance(NumericVector aux,
                                         NumericMatrix kmatrix,
                                         NumericMatrix mmatrix,
                                         NumericMatrix bm1,
                                         NumericMatrix bm2,
                                         int Lines,
                                         int Columns,
                                         int LCS,
                                         NumericMatrix OutputDistances){
  //NumericVector inputdiff(esom);
  ToroidDistance toroidDistance(aux,                   // create the worker
                                kmatrix,
                                mmatrix,
                                bm1,
                                bm2,
                                Lines,
                                Columns,
                                LCS,
                                OutputDistances);
  parallelFor(0, Lines, toroidDistance);                           // call it with parallelFor
  return OutputDistances;
}

// [[Rcpp::depends(RcppParallel)]]
struct NonToroidDistance : public Worker {    // Worker for parallelization
  // inputs to read from
  const RVector<double> aux;
  const RMatrix<double> bm1;
  const RMatrix<double> bm2;
  const int Lines;
  const int Columns;
  const int LCS;
  
  // output to write to
  RMatrix<double> OutputDistances;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to form the Rcpp matrix type)
  NonToroidDistance(const NumericVector aux,
                    const NumericMatrix bm1,
                    const NumericMatrix bm2,
                    const int Lines,
                    const int Columns,
                    const int LCS,
                    NumericMatrix OutputDistances):
    aux(aux),
    bm1(bm1),
    bm2(bm2),
    Lines(Lines),
    Columns(Columns),
    LCS(LCS),
    OutputDistances(OutputDistances) {}
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      for(int j = 0; j < Columns; j++){
        // sqrt(pow(aux.slice(0)-bm1,2) + pow(aux.slice(1)-bm2,2));
        int auxIdx1 = j*Lines + i;
        int auxIdx2 = LCS + j*Lines + i;
        OutputDistances(i,j) = sqrt(pow(aux[auxIdx1] - bm1(i,j), 2) + pow(aux[auxIdx2] - bm2(i,j), 2));
      }
    }
  }
};


// [[Rcpp::depends(RcppParallel)]]
NumericMatrix RcppParallelNonToroidDistance(NumericVector aux,
                                            NumericMatrix bm1,
                                            NumericMatrix bm2,
                                            int Lines,
                                            int Columns,
                                            int LCS,
                                            NumericMatrix OutputDistances){
  //NumericVector inputdiff(esom);
  NonToroidDistance nonToroidDistance(aux,                   // create the worker
                                      bm1,
                                      bm2,
                                      Lines,
                                      Columns,
                                      LCS,
                                      OutputDistances);
  parallelFor(0, Lines, nonToroidDistance);                           // call it with parallelFor
  return OutputDistances;
}




// [[Rcpp::depends(RcppParallel)]]
struct NeighborMatrix : public Worker {    // Worker for parallelization
  // inputs to read from
  const RMatrix<double> OutputDistances;
  const double Radius;
  const double Lines;
  const double Columns;
  
  // output to write to
  RMatrix<double> neighmatrix;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to form the Rcpp matrix type)
  NeighborMatrix(const NumericMatrix OutputDistances,
                 const double Radius,
                 const double Lines,
                 const double Columns,
                 NumericMatrix neighmatrix):
    OutputDistances(OutputDistances),
    Radius(Radius),
    Lines(Lines),
    Columns(Columns),
    neighmatrix(neighmatrix) {}
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      for(int j = 0; j < Columns; j++){
        double tmpVal = 1 - (pow(OutputDistances(i,j),2) / (3.14159265*pow(Radius,2)));
        if(tmpVal < 0){
          tmpVal = 0;
        }
        neighmatrix(i,j) = tmpVal;
      }
    }
  }
};
//neighmatrix = 1 - (OutputDistances % OutputDistances)/(3.14159265*Radius*Radius);


// [[Rcpp::depends(RcppParallel)]]
NumericMatrix RcppParallelNeighborMatrix(NumericMatrix OutputDistances,
                                         double Radius,
                                         double Lines,
                                         double Columns,
                                         NumericMatrix neighmatrix){
  //NumericVector inputdiff(esom);
  NeighborMatrix neighborMatrix(OutputDistances,                   // create the worker
                                Radius,
                                Lines,
                                Columns,
                                neighmatrix);
  parallelFor(0, Lines, neighborMatrix);                           // call it with parallelFor
  return neighmatrix;
}


// [[Rcpp::export]]
NumericVector trainstepC2(NumericVector esomwts,
                          NumericVector aux,
                          NumericMatrix DataSampled,
                          NumericMatrix BMUsampled,
                          double Lines,
                          double Columns,
                          double Weights,
                          double Radius,
                          bool toroid) {

  NumericVector DataSample(DataSampled.rows());
  NumericVector bmpos(BMUsampled.rows());
  
  NumericMatrix OutputDistances(Lines, Columns);
  NumericMatrix neighmatrix(Lines, Columns);
  NumericMatrix kmatrix(Lines, Columns);
  NumericMatrix mmatrix(Lines, Columns);
  NumericMatrix bm1(Lines, Columns);
  NumericMatrix bm2(Lines, Columns);
  
  int LCS = Lines * Columns;

  int xsize = Lines * Columns;
  for (int i = 0; i < xsize; i++) { // Fill with value
    kmatrix[i] = Lines-1;
  }
  for (int i = 0; i < xsize; i++) { // Fill with value
    mmatrix[i] = Columns-1;
  }

  int NumberOfDataSamples = DataSampled.nrow();
  
  for(int p = 0; p < NumberOfDataSamples; p++){
    DataSample = DataSampled.row(p);
    bmpos      = BMUsampled.row(p);
    for (int i = 0; i < xsize; i++) { // Fill with value
      bm1[i] = bmpos(0);
    }
    for (int i = 0; i < xsize; i++) { // Fill with value
      bm2[i] = bmpos(1);
    }
    if(toroid){
      OutputDistances = RcppParallelToroidDistance(aux,
                                                   kmatrix, mmatrix,
                                                   bm1, bm2,
                                                   Lines, Columns, LCS,
                                                   OutputDistances);
    }else{
      OutputDistances = RcppParallelNonToroidDistance(aux,
                                                      bm1, bm2,
                                                      Lines, Columns, LCS,
                                                      OutputDistances);
    }
    neighmatrix = RcppParallelNeighborMatrix(OutputDistances,
                                             Radius,
                                             Lines,
                                             Columns,
                                             neighmatrix);
    esomwts = RcppParallelDelta3DWeights(esomwts,
                                         DataSample,
                                         neighmatrix,
                                         Lines,
                                         Columns,
                                         Weights);
  }
  return(esomwts);
}
