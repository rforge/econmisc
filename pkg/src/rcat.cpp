//#define PRINT_DEBUG 
#include <iostream>
#include <Rcpp.h>

#ifdef PRINT_DEBUG
#define DEB(m) std::cout << m
#else
#define DEB(m)
#endif

using namespace Rcpp;

RcppExport SEXP do_rcat(SEXP arg_prob) {
// Generates a vector of categorical random variables where probabilities
// may differ for each draw
//
// arg_prob: matrix of case-specific probabilities.  These should sum to unity
// by rows
    Matrix<REALSXP> prob(arg_prob);
    int i, j;
    RNGScope scope;
    // needed to initialize/release RNG
    NumericVector rnd(runif(prob.nrow()));
    // uniform random numbers
    IntegerVector r(prob.nrow());
    // final state based on rows
    for(i = 0; i < prob.nrow(); i++) {
	double s = 0;
	for(j = 0; j < prob.ncol(); j++) {
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
