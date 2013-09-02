#include <iostream>

#include "constrainedEffectivePotential_inBrokenPhase.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int narg,char **arg)
{
// 	int L0(6), L1(6), L2(6), L3(6);
	int L0(192), L1(192), L2(192), L3(384);
	bool anti(false);
	double y_t(175.0/246.0);
// 	double y_t(175);
	double m0Squared(0.537);
	double lambda(-0.202);
	double lambda_6(0.05);
// 	double HiggsMassSquared(0.556140879127285);
	
	constrainedEffectivePotential_inBrokenPhase CEP(L0,L1,L2,L3,anti);
	
	CEP.set_yukawas(y_t,y_t);
	CEP.set_m0Squared(m0Squared);
	CEP.set_lambda(lambda);
	CEP.set_lambda_6(lambda_6);
	
	CEP.init_HiggsMassSquared();
	
	int maxCount=100;
	double testvalue=0.2;
	double result(0.0);
	for(int i=0; i<maxCount; ++i)
	{
		result=CEP.compute_fermionicContribution( testvalue );
	}
	cout.precision(15);
	cout <<"result = " <<result <<endl;
	
	
	
// 	double value_1(0.000030);
// 	
// 	double value_2(value_1+ 1e-7);
// 	double step(0.3e-7);
// 	double value_aim(value_1 + step);
	
// 	double diff(1.0e-6);
// 	double step(0.9e-6);
// 	double x_l(0.000032);
// 	double x_ll(x_l-diff);
// 	double x_u(x_l+diff);
// 	double x_uu(x_l+diff+diff);
// 	
// 	double x_aim(x_l+step);
// 	
// 	long double res_l(CEP.compute_fermionicContribution( x_l ));
// 	long double res_ll(CEP.compute_fermionicContribution( x_ll ));
// 	long double res_u(CEP.compute_fermionicContribution( x_u ));
// 	long double res_uu(CEP.compute_fermionicContribution( x_uu ));
// 	
// 	long double res_aim(CEP.compute_fermionicContribution( x_aim ));
// 	
// 	long double res_aim_lin = res_l + step*( res_u-res_l)/diff;
// 	long double res_aim_quad1 = res_l + step*( res_u-res_l)/diff + 0.5*step*step*(res_u + res_ll - 2.0*res_l)/(diff*diff);
// 	long double res_aim_quad2 = res_u - (diff-step)*( res_u-res_l)/diff + 0.5*(diff-step)*(diff-step)*(res_uu + res_l - 2.0*res_u)/(diff*diff);
// 	
// 	cout <<"x_ll: " <<x_ll <<"   x_l: " <<x_l <<"   x_u: " <<x_u <<"   x_uu: " <<x_uu << endl;
// 	cout <<"x_aim: "<<x_aim <<endl;
// 	cout <<"res_ll:       " << res_ll <<endl;
// 	cout <<"res_l:        " << res_l <<endl;
// 	cout <<"res_u:        " << res_u <<endl;
// 	cout <<"res_uu:       " << res_uu <<endl;
// 	cout <<"res_u-res_l:  " << res_u-res_l <<endl;
// 	cout <<"res_aim:      " << res_aim <<endl;
// 	cout <<"res_aim_lin:  " << res_aim_lin <<"    diff: " <<res_aim-res_aim_lin <<endl;
// 	cout <<"res_aim_quad1:" << res_aim_quad1 <<"    diff: " <<res_aim-res_aim_quad1 <<endl;
// 	cout <<"res_aim_quad2:" << res_aim_quad2 <<"    diff: " <<res_aim-res_aim_quad2 <<endl;
	
	
	
	
}