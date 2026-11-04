#ifndef NOC_ORIGINAL_BASIS_DECLARATION_V01
#define NOC_ORIGINAL_BASIS_DECLARATION_V01
#include <RcppArmadillo.h>
Rcpp::List noc_basis_cpp(const arma::mat& B,const arma::vec& a,bool project);
#endif
