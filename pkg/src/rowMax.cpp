#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP do_rowMax(SEXP arg_matrix, SEXP arg_na_rm) {
    // returns the vector of maximum elements in each row for a matrix
    Matrix<REALSXP> m(arg_matrix);
    LogicalVector na_rm(arg_na_rm);
    bool narm = na_rm(0);
    NumericVector max(m.nrow());
    int i, j;
    for(i = 0; i < m.nrow(); i++) {
	bool valid = false;
	// max is not valid (yet)
	for(j = 0; j < m.ncol(); j++) {
	    if(ISNA(m(i, j))) {
		if(!narm) {
		    // NA in the row -> NA answer
		    max(i) = NA_REAL;
		    break;
		}
		// otherwise we ignore it
	    }
	    else {
		// was not NA
		if(!valid) {
		    max(i) = m(i, j);
		    valid = true;
		}
		else {
		    if(max(i) < m(i, j)) {
			max(i) = m(i, j);
		    }
		}
	    }
	    if(!valid)
		max(i) = NA_REAL;
	}
    }
    return max;
}
