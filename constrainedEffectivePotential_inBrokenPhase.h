#ifndef CONSTRAINEDEFFECTIVEPOTENTIAL_INBROKENPHASE_H
#define CONSTRAINEDEFFECTIVEPOTENTIAL_INBROKENPHASE_H

/*
This Program should determin the "phase structure" from the CEP, while still assuming to be in the broken phase.

The prrinciple is the following:

 - we use the continuum notation (m_0^2, lambda, y, lambda_6), but keep all lattice expressions like sums etc
 - The potential we use is the one for the determination of the lower bound
 - we set all parameters, also m_0^2
 - now the program iterates
   {
 	- determine minimum
 	- determin mass in the minimum
 	- set mass in propagator sum
   }
   until the minimum is stable
*/
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <utility>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_min.h"

class constrainedEffectivePotential_inBrokenPhase
{
	private:
	//physical parameters
	int L0, L1, L2, L3; //only L3 can be chosen antiperiodic
	bool antiperiodicBC_L3;
	
	//physical parameters of model
	double m0Squared;
	double yukawa_t;
	double yukawa_b;
	double lambda;
	double lambda_6;
	int N_f; //default 1
	
	//numerical parameters
	double rho; //default 1
	double one_ov_two_rho;
	double r;   //default 0.5
	
	//Arrays for faster computation of fermionic contribution and propagatorsums
	std::complex< double > *eigenvalues_of_overlap;
	double *factor_for_eigenvalues;
	double *four_times_sum_of_sinSquared; //store 4*\sum( sin(p_mu/2)^2)
	double *factor_for_sum_of_sinSquared;
	
	double propagatorSum_withZeroMass;
	double propagatorSum_withHiggsMass;
	double actual_HiggsMassSquared;
	
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
	constrainedEffectivePotential_inBrokenPhase(int l0, int l1, int l2, int l3, bool anti);
	~constrainedEffectivePotential_inBrokenPhase();
	
	//set parameters
	void set_m0Squared( double new_m0Squared );
	void set_yukawas( double y_t, double y_b );
	void set_lambda( double new_lambda );
	void set_lambda_6( double new_lambda_6 );
	void set_N_f( int new_N_f);
	
	void set_rho( double new_rho );
	void set_r( double new_r);
	
	void init_HiggsMassSquared(); //set higgs mass to m0^2 if m0^2 is non-negative. elsewise it is set to zero, also sets PropSum (if lambdas are != 0)
	void set_HigsMassSquared( double mHSquared ); //also sets the propagatorSum
	
	void reinitialize();
	
	double get_moSquared();
	double get_yukawa_t();
	double get_yukawa_b();
	double get_lambda();
	double get_lambda_6();
	int get_N_f();
	
	double get_rho();
	double get_r();
	
	double get_actual_HiggsMassSquared();
	
	
	
	//returns eigenvalue nu^+(p) of the free overlap
	std::complex< double > computeAnalyticalEigenvalue(double p0, double p1, double p2, double p3);
	
	//fills the arrays for faster computation
	void fill_eigenvaluesAndFactors();
	void fill_sinSquaredAndFactors();
	
	double compute_propagatorSum( double massSquared );
	
	//returns the CEP(\hat v), note, the argument v is relatet to the magnetization via v=sqrt(2*kappa)*mag
	double compute_CEP_inBrokenPhase( double value );
	//parts for the CEP
	double compute_treeLevel( double value );
	double compute_fermionicContribution( double value );
	double compute_firstOrderInLambdas( double value );
	
	double compute_CEP_inBrokenPhase_secondDerivative( double value );
	//parts for CEP_secondDerivative
	double compute_treeLevel_secondDerivative( double value );
	double compute_fermionicContribution_secondDerivative( double value );
	double compute_firstOrderInLambdas_secondDerivative( double value );
	
	
	//stuff for the minimizer
	void set_minimizationAlgorithm( int new_algorithm );
	bool reInitialize_minimizer( double minimum, double lower, double upper );
	bool initialize_minimizer( double minimum, double lower, double upper );
// 	void determine_startingAndInitialize( double lower, double upper, double step ); //scans CEP and initializes with best starting point
	void determine_startingPoints( double inLower, double inUpper, double inStep, double &outMinimum, double &outLower,double &outUpper);
	
	int iterate_minimizer();
	double get_actual_minimum();
	void get_actual_Interval( double &outMinimum, double &outLower,double &outUpper );
	bool check_convergence();//calls the gsl routine for checking the convergence
	int iterate_minimizer_until_convergence(); //returns number of iteration needed,
	                                           //returns -1 if maxIter is reached without convergence, 
	                                           //returns -2 if iterator returns GSL_EBADFUNC meaning nan or inf occured
	                                           
};

//wrapper for gsl since member function does not work
double wrapper_compute_CEP_inBrokenPhase_gsl(double value, void *params);

#endif