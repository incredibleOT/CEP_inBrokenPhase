#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "CEP_withFullBosDet.h"


using std::cout;
using std::cerr;
using std::endl;

int main()
{
	
	int L0(32), L1(32), L2(32), L3(64);
	bool anti(false);
	
	
	double m0Squared=0.15;
	double lambda_6=0.001;
	double lambda=-0.008;
// 	double lambda=-0.08;
	double y(175.0/246.0);
	
	//parameters for minimizer
	int max_numberOfIterations(100); //default 100
	double relative_Accuracy(1.0e-7);   //default 1e-7
	double absolute_Accuracy(1.0e-7);   //default 1e-7
	int minimizationAlgorithm(2); // 1 - gsl_min_fminimizer_goldensection
	
	//trial range for start of minimizer
	
	double v_min=0.01;
	double v_max=1.00;
	double v_step=0.05;
	
	
	cout <<"initialize CEP_withFullBosDet:" <<endl;
	CEP_withFullBosDet CEP(L0,L1,L2,L3,anti);
	cout <<"done" <<endl <<endl;
	
	cout <<"set parameters" <<endl;
	CEP.set_lambda( lambda );
	cout <<"lambda = " <<CEP.get_lambda() <<endl;
	CEP.set_lambda_6( lambda_6 );
	cout <<"lambda_6 = " <<CEP.get_lambda_6() <<endl;
	CEP.set_m0Squared( m0Squared );
	cout <<"m0Squared = " <<CEP.get_m0Squared() <<endl;
	CEP.set_yukawas( y,y );
	cout <<"y_t = " <<CEP.get_yukawa_t() <<endl;
	cout <<"y_b = " <<CEP.get_yukawa_b() <<endl;
	cout <<"done" <<endl <<endl;
	
	
	cout <<"Set parametes for minimizer" <<endl;
	CEP.set_max_numberOfIterations(max_numberOfIterations);
	CEP.set_relative_Accuracy(relative_Accuracy);
	CEP.set_absolute_Accuracy(absolute_Accuracy);
	CEP.set_minimizationAlgorithm(minimizationAlgorithm);
	cout <<"done" <<endl <<endl;
	
	/*
	double plus_inf=1.0/log(1.0);
	double minus_inf=-1.0/log(1.0);
	double number=3.0;
	cout <<"plus_inf:  " <<plus_inf <<endl;
	cout <<"minus_inf: " <<minus_inf <<endl;
	cout <<plus_inf <<" < " <<number <<" = " <<(plus_inf < number) <<endl;
	cout <<minus_inf <<" < " <<number <<" = " <<(minus_inf < number) <<endl;
	*/
	
	//check potential for initial guess of minimizer
	double min_test(0.0), lower_test(0.0), upper_test(0.0);
	//function scans roughly a given region (v_min...v_max) and puts the minimum in this range into min_test. The left and right values are
	// stored in lower_test and upper_test
	// if global minimum is at the border, all three values are set to the same value 
	bool return_from_initial_scan( CEP.determine_startingPoints( v_min, v_max, v_step, min_test, lower_test, upper_test) );
	
	cout <<"lower_test = " <<lower_test <<"    min_test = " <<min_test <<"    upper_test = " <<upper_test <<endl; 
	
	if( !return_from_initial_scan )
	{
		cerr <<"Error, negative divergence occured during initial scan" <<endl;
		exit(EXIT_FAILURE);
	}
	
	if(lower_test == upper_test)
	{
		cerr <<"Error, initial testrange did not yield global maximum within the range" <<endl;
		exit(EXIT_FAILURE);
	}
	
	cout <<"initialize minimizer" <<endl;
	if( !CEP.initialize_minimizer( min_test, lower_test, upper_test ) )
	{
		cerr <<"Error initializing the minimizer" <<endl;
		exit(EXIT_FAILURE);
	}
	cout <<"done" <<endl <<endl;
	
	cout <<"begin minimizetaion process" <<endl;
	int num_of_iterations( CEP.iterate_minimizer_until_convergence() );
	if( num_of_iterations == -1 )
	{
		cerr <<"Error, minimizer did not converge in max_numberOfIterations=" <<CEP.get_max_numberOfIterations() <<" iterations" <<endl;
		exit(EXIT_FAILURE);
	}
	else if( num_of_iterations == -2 )
	{
		cerr <<"Error, nan of inf occured during minimization process" <<endl;
		exit(EXIT_FAILURE);
	}
	else if( num_of_iterations == -3 )
	{
		cerr <<"Error, minimizer could not improve without reaching convergence" <<endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		cout <<"Minimizer converged after " <<num_of_iterations <<" iterations" <<endl;
	}
	cout.precision(12);
	cout <<"minimum found at v = " <<CEP.get_actual_minimum() <<endl;
	cout <<"potential at minimmum U(v) = " <<CEP.get_potentialAtMinimum() <<endl;
	

	
	//NOTE continue here!!
	//NOTE continue here!!
	//NOTE continue here!!
	//NOTE continue here!!
	//NOTE continue here!!
	
	
	
	
	
// 	double v_min=atof( arg[1] );
// 	double v_max=atof( arg[2] );
// 	double v_step=atof( arg[3] );
// 	if( v_max - v_min < 0.0 || v_step <=0.0 )
// 	{
// 		cerr << "Error, negative range" <<endl;
// 		exit(EXIT_FAILURE);
// 	}
// 	int counterMax=100000;
// 	if( (v_max-v_min)/v_step >= static_cast< double >(counterMax) )
// 	{
// 		cerr<<"To many steps" <<endl;
// 		exit(EXIT_FAILURE);
// 	}
// 	
// 	std::ofstream outputFile( arg[4] );
// 	if(!outputFile.good())
// 	{
// 		cerr <<"Error opening output file" <<endl <<arg[4] <<endl;
// 		exit(EXIT_FAILURE);
// 	}
// 	outputFile <<"# L0 = " <<L0 <<"  L1 = " <<L1 <<"  L2 = " <<L2 <<"  L3 = " <<L3 <<endl;
// 	outputFile <<"# lambda = " <<CEP.get_lambda() <<endl;
// 	outputFile <<"# lambda_6 = " <<CEP.get_lambda_6() <<endl;
// 	outputFile <<"# m0Squared = " <<CEP.get_m0Squared() <<endl;
// 	outputFile <<"# y_t = " <<CEP.get_yukawa_t() <<endl;
// 	outputFile <<"# y_b = " <<CEP.get_yukawa_b() <<endl;
// 	cout <<"# v     U_tree    U_ferm    U_det    U_1st   U_total" <<endl;
// 	
// 	double v=v_min;
// 	double U_tree(0.0),U_ferm(0.0),U_det(0.0),U_1st(0.0),U_total(0.0);
// 	cout.precision(10);
// 	outputFile.precision(10);
// 	for(int counter=0; counter <=counterMax; ++counter)
// 	{
// 		U_tree=CEP.compute_treeLevel(v);
// 		U_ferm=CEP.compute_fermionicContribution(v);
// 		U_det=CEP.compute_BosDetContribution(v);
// 		U_1st=CEP.compute_firstOrderInLambdas(v);
// 		U_total=CEP.compute_CEP_withFullBosDet(v);
// 			
// 	
// 		outputFile <<v <<"  " <<U_tree <<"  " <<U_ferm <<"  " <<U_det <<"  " <<U_1st <<"  " <<U_total <<endl;
// 		cout <<v <<"  " <<U_tree <<"  " <<U_ferm <<"  " <<U_det <<"  " <<U_1st <<"  " <<U_total <<endl;
// 		v+=v_step;
// 		if(v>v_max){ break; }
// 	}
// 	
// 	if(!outputFile.good())
// 	{
// 		cerr <<"WARNING!!! something went wrong during output" <<endl;
// 	}
// 	outputFile.close();
	
	
	
	
	
	
	
	
	
	
	
	
	
}
