#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP do_rowMax(SEXP arg_matrix) {
    // returns the vector of maximum elements in each row for a matrix
    Matrix<REALSXP> m(arg_matrix);
    NumericVector max(m.nrow());
    int i, j;
    for(i = 0; i < m.nrow(); i++) {
	max(i) = m(i, 0);
	for(j = 1; j < m.ncol(); j++) {
	    if(max(i) < m(i, j)) {
		max(i) = m(i, j);
	    }
	}
    }
    return max;
}
