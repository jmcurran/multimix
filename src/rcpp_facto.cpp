#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
  int facto_(int *n);
}

//' Test of calling a Fortran Function
 //'
 //' @param n an integer
 //' @return n factorial (n! = n * (n -1 ) * ... * 1)
 //' @examples
 //' test_function(9)
 // [[Rcpp::export]]
int test_function(int n) {
   int result = facto_(&n);
   return result;
 }
