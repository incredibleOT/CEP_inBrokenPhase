#ifndef CEP_WITHFULLBOSDET_H
#define CEP_WITHFULLBOSDET_H

/*
This computes the CEP when the gaussian part contains as much as possible as it can be found in chapter 10 of the notes.
Eventually one might also to include the possibility to minimize the potential wrt the zero mode
So far, only computing the potential for a set of parameters
*/
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <utility>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_min.h"

class CEP_withFullBosDet
{
	private:
	//physical parameters
	int L0, L1, L2, L3; //only L3 can be chosen antiperiodic
	bool antiperiodicBC_L3; //default false, only for fermions
	
	//physical parameters of model
	double m0Squared;
	double yukawa_t;
	double yukawa_b;
	double lambda;
	double lambda_6;
	int N_f; //default 1
	
	bool ignore_goldstone_modes; //simply all P_G are set to zero
	
	//numerical parameters
	double rho; //default 1
	double one_ov_two_rho;
	double r;   //default 0.5
	
	//Arrays for faster computation of fermionic contribution and propagatorsums
	std::complex< double > *eigenvalues_of_overlap;
	double *factor_for_eigenvalues;
	double *four_times_sum_of_sinSquared; //store 4*\sum( sin(p_mu/2)^2)
	double *factor_for_sum_of_sinSquared;
	
	int numberOfDistingtMomenta_fermions;
	int numberOfDistingtMomenta_bosons;
	
	//stuff for minimizer
	int max_numberOfIterations; //default 100
	double relative_Accuracy;   //default 1e-7
	double absolute_Accuracy;   //default 1e-7
	int minimizationAlgorithm; // 1 - gsl_min_fminimizer_goldensection
	                           // 2 - gsl_min_fminimizer_brent
	                           // 3 - gsl_min_fminimizer_quad_golden
	                           //default 1
	
	//stuff for gsl
	gsl_min_fminimizer *minimizer;
	gsl_function functionHandler;
	const gsl_min_fminimizer_type *algorithmForMinimization;
	bool minimizerInitialized;
	int iterator_status;//stores the output of the iterator (some gsl enum)
	
	
	public:
	
	//constructor and destructor
	CEP_withFullBosDet(int l0, int l1, int l2, int l3, bool anti);
	~CEP_withFullBosDet();
	
	//set parameters
	void set_m0Squared( double new_m0Squared );
	void set_yukawas( double y_t, double y_b );
	void set_lambda( double new_lambda );
	void set_lambda_6( double new_lambda_6 );
	void set_N_f( int new_N_f);
	
	void set_ignore_goldstone_modes(bool new_ignore);
	
	void set_rho( double new_rho );
	void set_r( double new_r);
	
	
	double get_m0Squared();
	double get_yukawa_t();
	double get_yukawa_b();
	double get_lambda();
	double get_lambda_6();
	int get_N_f();
	
	double get_rho();
	double get_r();
	
	
	//returns eigenvalue nu^+(p) of the free overlap
	std::complex< double > computeAnalyticalEigenvalue(double p0, double p1, double p2, double p3);
	
	//fills the arrays for faster computation
	void fill_eigenvaluesAndFactors();
	void fill_sinSquaredAndFactors();
	
	double compute_propagatorSum( double massSquared ); //computes sum_p( 1/(\hat p^2 + massSquared)) NOTE excludes zero mode
	
	//returns the CEP(\hat v), note, the argument v is relatet to the magnetization via v=sqrt(2*kappa)*mag
	double compute_CEP_withFullBosDet( double value );
	//parts for the CEP
	double compute_treeLevel( double value );
	double compute_fermionicContribution( double value );
	double compute_BosDetContribution( double value );
	double compute_firstOrderInLambdas( double value );
	
	//bool load_fermionicContribution( const std::string &fileName );
	
	
// 	double compute_CEP_withFullBosDet_secondDerivative( double value );
// 	//parts for CEP_secondDerivative
// 	double compute_treeLevel_secondDerivative( double value );
// 	double compute_fermionicContribution_secondDerivative( double value );
// 	double compute_BosDetContribution_secondDerivative( double value );
// 	double compute_firstOrderInLambdas_secondDerivative( double value );
	
	//stuff for the minimizer
	//set parameter
	void set_max_numberOfIterations(int new_max);
	void set_relative_Accuracy( double new_rel_acc);
	void set_absolute_Accuracy( double new_abs_acc);
	void set_minimizationAlgorithm( int new_algorithm );
	//get parameters
	int get_max_numberOfIterations();
	double get_relative_Accuracy();
	double get_absolute_Accuracy();
	int get_minimizationAlgorithm(); 
	
	bool reInitialize_minimizer( double minimum, double lower, double upper );
	bool initialize_minimizer( double minimum, double lower, double upper );
	
	//will return true if no -inf occurs during scan, otherwise returns false
	//if there is no global minimum, but for all values the function evaluated correctly, out* are all set to the same value
	bool determine_startingPoints( double inLower, double inUpper, double inStep, double &outMinimum, double &outLower,double &outUpper);
	
	int iterate_minimizer();
	double get_actual_minimum();
	double get_potentialAtMinimum();
	void get_actual_Interval( double &outMinimum, double &outLower,double &outUpper );
	bool check_convergence();//calls the gsl routine for checking the convergence
	int iterate_minimizer_until_convergence(); //returns number of iteration needed,
	                                           //returns -1 if maxIter is reached without convergence, 
	                                           //returns -2 if iterator returns GSL_EBADFUNC meaning nan or inf occured
	
};

//wrapper for gsl since member function does not work
double wrapper_compute_CEP_withFullBosDet_gsl(double value, void *params);
#endif
