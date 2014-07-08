#include <iostream>

#include "CEP_withFullBosDet.h"

using std::cout;
using std::cerr;
using std::endl;

int main()
{
// 	int L0(4), L1(4), L2(4), L3(4);
// 	int L0(8), L1(8), L2(8), L3(8);
int L0(24), L1(24), L2(24), L3(32);
// 	int L0(32), L1(32), L2(32), L3(64);
	bool anti(false);
// 	double y_t(175.0/246.0);
	double y_t=0.711382;
	double m0Squared(0.15);
	double lambda(-0.008);
	double lambda_6(0.001);
	
	CEP_withFullBosDet CEP(L0,L1,L2,L3,anti);
	
	CEP.set_yukawas(y_t,y_t);
	CEP.set_m0Squared(m0Squared);
	CEP.set_lambda(lambda);
	CEP.set_lambda_6(lambda_6);
	
	double test_low(1.0e-7), test_high(3.0), test_step(0.05);
	double start_low(0.0), start_high(0.0), start_min(0.0);
	CEP.determine_startingPoints(test_low, test_high, test_step, start_min, start_low, start_high);
	cout <<"start_min=" <<start_min <<"   start_low=" <<start_low <<"   start_high=" <<start_high <<endl;
	if(! CEP.initialize_minimizer( start_min, start_low, start_high ) )
	{ 
		cerr <<"Errror with initialization" <<endl;
		exit(EXIT_FAILURE); 
		
	}
	
	int num_of_iter(CEP.iterate_minimizer_until_convergence());
	cout <<"minimizer converged after " <<num_of_iter <<" iteration" <<endl;
	double minimum=CEP.get_actual_minimum();
	double potential=CEP.get_potentialAtMinimum();
	cout.precision(15);
	cout <<"Minimum: " <<minimum <<"   potential: " <<potential <<endl;
	
	cout <<"contributions:" <<endl;
	
	cout <<"U_tree = " <<CEP.compute_treeLevel(minimum) <<endl;
	cout <<"U_ferm = " <<CEP.compute_fermionicContribution(minimum) <<endl;
	cout <<"U_bosDet = " <<CEP.compute_BosDetContribution(minimum) <<endl;
	cout <<"U_1st = " <<CEP.compute_firstOrderInLambdas(minimum) <<endl;
	cout <<"U_total = " <<CEP.compute_CEP_withFullBosDet(minimum) <<endl;
	cout <<endl;
	
	cout <<"contributions to 2nd derivative:" <<endl;
	cout <<"U''_tree = " <<CEP.compute_treeLevel_secondDerivative(minimum)   <<endl;
	cout <<"U''_ferm = " <<CEP.compute_fermionicContribution_secondDerivative(minimum)   <<endl;
	cout <<"U''_bosDet = " <<CEP.compute_BosDetContribution_secondDerivative(minimum)   <<endl;
	cout <<"U''_1st = " <<CEP.compute_firstOrderInLambdas_secondDerivative(minimum)   <<endl;
	cout <<"U''_total = " <<CEP.compute_CEP_withFullBosDet_secondDerivative(minimum)   <<endl;
	
	
	
}