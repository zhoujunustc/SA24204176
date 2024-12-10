#include <Rcpp.h>
using namespace Rcpp;

double FDP(NumericVector M, double t) {
  int num = 0, den = 0;
  for (int j = 0; j < M.size(); j++){
    if (M[j] <= -t) num++;
    if (M[j] >=  t) den++;
  }
  den = den >= 1 ? den : 1;
  return 1.0 * num / den;
}


//' @title Selection features
//' @description function to perform selection giving Mirror statistic and FDR level
//' @param M mirror statistic of all features
//' @param q the nominal FDR level
//' @return the selected features that are considered significant 
//' @export
// [[Rcpp::export]]
IntegerVector Selection(NumericVector M, double q) {
  if (q <= 0 || q >= 1) {
    stop("Error: Invalid q value");
  }
  
  double tao = 0.0; double interval = max(M) / 1000;
  while(FDP(M, tao) > q)  tao += interval;
  
  IntegerVector features;
  
  for(int j = 0; j < M.size(); j++) {
    if (M[j] >= tao) {
      features.push_back(j + 1);
    } 
  }
  return features;
}