#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcpp_hamunyela(NumericVector num) {
  const int window = 2;
  const int n = num.size();
  double threshold = -0.01;
  NumericVector dat = clone(num);

  int before_d = 0;
  int after_d = 0;

  for (int j=1; j < (n-1); ++j) {
    before_d = (dat[j] - dat[j-1]) < threshold*dat[j-1];
    after_d = (dat[j] - dat[j+1]) < threshold*dat[j+1];

    if (before_d && after_d) {
      dat[j] = (dat[j-1]+dat[j+1])/2;
    }
  }
  return dat;
}
