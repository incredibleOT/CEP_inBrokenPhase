#include <iostream>

#include "constrainedEffectivePotential_inBrokenPhase.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int narg,char **arg)
{
// 	int L0(6), L1(6), L2(6), L3(6);
	int L0(32), L1(32), L2(32), L3(32);
	bool anti(false);
// 	double y_t(175.0/246.0);
	double y_t(0.75);
	double m0Squared(-0.941952694810136);
	double lambda(0.1);
	double lambda_6(0.1);
// 	double HiggsMassSquared(0.556140879127285);
	
	constrainedEffectivePotential_inBrokenPhase CEP(L0,L1,L2,L3,anti);
	
	CEP.set_yukawas(y_t,y_t);
	CEP.set_m0Squared(m0Squared);
	CEP.set_lambda(lambda);
	CEP.set_lambda_6(lambda_6);
	
	CEP.init_HiggsMassSquared();
	
// 	CEP.set_HigsMassSquared(HiggsMassSquared);
	cout.precision(15);
	cout <<"actual mHSquared: " <<CEP.get_actual_HiggsMassSquared() <<endl;
	cout <<"propsum(0.0): " <<CEP.compute_propagatorSum(0.0) <<endl;
	cout <<"propsum(mHSquared): " <<CEP.compute_propagatorSum(CEP.get_actual_HiggsMassSquared()) <<endl;
	double mag=0.1;
	cout <<"tree: " <<CEP.compute_treeLevel(mag) <<endl;
	cout <<"ferm: " <<CEP.compute_fermionicContribution(mag) <<endl;
	cout <<"1st:  " <<CEP.compute_firstOrderInLambdas(mag) <<endl;
	cout <<"full: " <<CEP.compute_CEP_inBrokenPhase(mag) <<endl;
	cout <<"full second derivative at: " <<CEP.compute_CEP_inBrokenPhase_secondDerivative(mag) <<endl;
	
	double inLower(0.01), inUpper(5.01), inStep(0.1), minimum(0.0), lower(0.0), upper(0.0);
	CEP.determine_startingPoints(inLower, inUpper, inStep, minimum, lower, upper);
	cout <<"result of first scan: minimum at: " <<minimum <<"   lower end: " <<lower <<"   upper end: " <<upper <<endl;
	CEP.set_minimizationAlgorithm(2);
	if( CEP.initialize_minimizer(minimum, lower, upper) )
	{
		cout <<"minimizer initiallized successfully" <<endl;
	}
// 	for(int i=1; i<150; ++i)
// 	{
// 		int iterOutput=CEP.iterate_minimizer();
// 		double min(0.0), low(0.0), up(0.0);
// 		CEP.get_actual_Interval(min,low,up);
// 		cout <<"iteration " <<i <<"   actual minimum at: " <<min <<"   lower: " <<low <<"   upper: " <<up <<"    status: " <<iterOutput 
// 		     <<"   converged: " <<CEP.check_convergence() <<endl;
// 	}
	
	
	bool toContinue(true);
	double limitForHiggsMass(1.0e-5);
	int counter(0), maxCounter(100);
	while(toContinue)
	{
		++counter;
		double oldHiggsMassSquared=CEP.get_actual_HiggsMassSquared();
		int numberOfIteration(0);
// 		cout <<"HiggsMassSquared before iteration: " <<oldHiggsMassSquared <<endl;
		numberOfIteration=CEP.iterate_minimizer_until_convergence();
		cout <<"Convergence reached after " <<numberOfIteration <<" iterations.   Minimum at " <<CEP.get_actual_minimum() <<endl;
		double newHiggsMassSquared=CEP.compute_CEP_inBrokenPhase_secondDerivative(CEP.get_actual_minimum());
		cout <<"HiggsMassSquared after iteration: " <<newHiggsMassSquared <<endl;
		if( std::abs( 2.0*(newHiggsMassSquared - oldHiggsMassSquared)/(oldHiggsMassSquared + newHiggsMassSquared)) <limitForHiggsMass )
		{
			toContinue=false;
			cout <<"Higgs mass determination converged after " <<counter << " iterations" <<endl;
			cout <<"Final values:" <<endl;
			cout <<"minimum at : " <<CEP.get_actual_minimum() <<"   corresponding to cutoff of " <<246.0/CEP.get_actual_minimum() <<" GeV" <<endl;
			cout <<"HiggsMassSquared: " <<newHiggsMassSquared <<"   corresponding to mass of "
			     <<sqrt(newHiggsMassSquared)*246.0/CEP.get_actual_minimum() <<" GeV" <<endl;
		}
		if( counter== maxCounter )
		{
			cout <<"No convergence after " <<maxCounter <<" iterations" <<endl;
			break;
		}
		CEP.set_HigsMassSquared(newHiggsMassSquared);
		CEP.determine_startingPoints(inLower, inUpper, inStep, minimum, lower, upper);
		CEP.reInitialize_minimizer(minimum, lower, upper);
	}
	
	
}