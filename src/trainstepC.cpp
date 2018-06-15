#include "DataBionicSwarm.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube trainstepC(Rcpp::NumericVector vx,Rcpp::NumericVector vy,Rcpp::NumericMatrix DataSampled,Rcpp::NumericMatrix BMUsampled,double Lines, double Columns,double Radius,bool toroid) {
  
  Rcpp::IntegerVector x_dims = vx.attr("dim");
  Rcpp::IntegerVector y_dims = vy.attr("dim");
  int NumberOfweights=x_dims[2];
  double k=Lines;
  double m=Columns;
  
  arma::cube esomwts(vx.begin(), x_dims[0], x_dims[1], x_dims[2], false);
  arma::cube aux(vy.begin(), y_dims[0], y_dims[1], y_dims[2], false);
  
  Rcpp::NumericVector DataSample(DataSampled.rows());
  Rcpp::NumericVector bmpos(BMUsampled.rows());
  arma::mat OutputDistances(k,m);
  arma::mat neighmatrix(k,m);
  
  arma::cube neigharray(k, m, NumberOfweights); 
  arma::cube inputdiff(k, m, NumberOfweights); 
  
  arma::mat kmatrix(x_dims[0], x_dims[1]);
  kmatrix.fill(k-1);
  arma::mat mmatrix(x_dims[0], x_dims[1]);
  mmatrix.fill(m-1);
  arma::mat bm1(x_dims[0], x_dims[1]);
  arma::mat bm2(x_dims[0], x_dims[1]);
  
  int NumberOfDataSamples=DataSampled.nrow();
  for(int p=0;p<NumberOfDataSamples;p++){
    // std::cout<<"anfang:p"<<p<<std::endl;
    //// Begin One Learnstep for one inputvector (1 Datenzeile)
    DataSample=DataSampled.row(p);
    bmpos=BMUsampled.row(p);
    // toroid map: different distances
    // example: on a toroid 5 x 4-grid the distance between the points (3,4) and (1,1) is sqrt(8)
    //dass geht so in arma nicht, dass muss ich elementweise machen
    bm1.fill(bmpos(0));
    bm2.fill(bmpos(1));
    //  if(p==1)
    //  std::cout<<kmatrix*2<<" "<<sqrt(half)<<std::endl;
    //std::cout<<kmatrix-56<<" "<<abs(kmatrix-56)<<std::endl;
    //std::cout<<kmatrix<<" "<<pow(kmatrix,2)<<endl;
    //   std::cout<<"slice1p"<<p<<std::endl;
    if (toroid){
      //   OutputDistances <- 0.5*sqrt(  (k-abs( 2*abs(aux[,,1]-bmpos[1])-k ))^2 + (m-abs( 2*abs(aux[,,2]-bmpos[2])-m ))^2)
      OutputDistances = 0.5*sqrt(  pow(kmatrix-abs( 2*abs(aux.slice(0)-bm1)-kmatrix ),2) + pow(mmatrix-abs( 2*abs(aux.slice(1)-bm2)-mmatrix ),2));
      //OutputDistances =OutputDistances%0.5;
    }else{
      OutputDistances =  sqrt(pow(aux.slice(0)-bm1,2) + pow(aux.slice(1)-bm2,2));
    } //end if toroid
    //Bestimme Nachbahrschaftsfunktion innerhalb Radius
    //  neighborfunction='kreis'
    neighmatrix=1-(OutputDistances%OutputDistances)/(3.14159265*Radius*Radius);
    //neighmatrix[neighmatrix<0]=0
    for(unsigned int i=0;i<neighmatrix.n_rows;i++){
      for(unsigned int j=0;j<neighmatrix.n_cols;j++){
        if(neighmatrix(i,j)<0)
          neighmatrix(i,j)=0;
      }
    }
    // neigharray=array(neighmatrix,c(k,m,NumberOfweights)) //Bringe auf gleiche Dimension wie esom und inputdiff
    //   std::cout<<"slice2p "<<p<<std::endl;
    for(int i=0;i<NumberOfweights;i++){
      neigharray.slice(i)=neighmatrix;
    }
    //Das muesste auch vektoriesierbar sein, momentan klappts aber nur als for schleife 
    //Injeder Dimension der inputdiff matrize wirdein Wert des Datenvekotr abgezogen
    //dadurch wird defakto von jedem gewicht der datenvektor komplett abgezogen
    //das dunktioniert, weil R automatisch den Wert auf die korrekte Matrizengroesse vergroessert
    
    // inputdiff <- esom
    // for (w in 1:NumberOfweights) {
    //   //neigharray[,,w] <- neighmatrix //Nachbahrschaftsmatrix ist fuer alle Gewichtsvektoren konstant
    //   inputdiff[,,w] <- inputdiff[,,w]-DataSample[w]
    // } //end for 1:NumberOfweights
    inputdiff=Delta3DWeightsC(esomwts,DataSample); //MT: Ohne Multiplikation mit 1 ist esom==inputdiff, kA wieso
    //CurrentLearnrate is always 1
    esomwts = esomwts-neigharray%inputdiff; // element-wise cube multiplicatio with %
    //// End One Learnstep for one inputvector (1 Datenzeilen)
    // Hold BMUs
    //std::cout<<"hold p"<<p<<std::endl;
    for(int i=0;i<NumberOfDataSamples;i++){
      for(int j=0;j<NumberOfweights;j++){
        
        // std::cout<<i<<" "<<j<<" "<<BMUsampled(i,0)<<" "<< BMUsampled(i,1)<<std::endl;
        esomwts(BMUsampled(i,0),BMUsampled(i,1),j) = DataSampled(i,j);
      }
    } // end for hold bmus
    //std::cout<<"end p"<<p<<std::endl;
  } //end for 1:NumberOfDataSamples
  return(esomwts);
}