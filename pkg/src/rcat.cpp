//#define PRINT_DEBUG 
#include <iostream>
#include <Rcpp.h>

#ifdef PRINT_DEBUG
#define DEB(m) std::cout << m
#else
#define DEB(m)
#endif

using namespace Rcpp;

RcppExport SEXP do_rcat(SEXP arg_prob, SEXP arg_na_rm) {
// Generates a vector of categorical random variables where probabilities
// may differ for each draw
//
// arg_prob: matrix of case-specific probabilities.  These should sum to unity
// by rows
    Matrix<REALSXP> prob(arg_prob);
    LogicalVector na_rm(arg_na_rm);
    bool narm = na_rm(0);
    int i, j;
    RNGScope scope;
    // needed to initialize/release RNG
    NumericVector rnd(runif(prob.nrow()));
    // uniform random numbers
    IntegerVector r(prob.nrow());
    // final state based on rows
    for(i = 0; i < prob.nrow(); i++) {
	double s = 0;
	r(i) = NA_INTEGER;
	// if all the probs are NA for this row, we return NA as well
	for(j = 0; j < prob.ncol(); j++) {
	    if(ISNA(prob(i,j))) {
		if(!narm)
		    Rf_error("NA in argument 'prob'");
		}
	    else
		s += prob(i, j);
	    if(rnd(i) < s) {
		r(i) = j + 1;
		// +1 -- we start counting from 1
		break;
	    }
	}
    }
    return r;
}
