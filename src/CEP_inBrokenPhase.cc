#include "CEP_inBrokenPhase.h"


//constructor
CEP_inBrokenPhase::CEP_inBrokenPhase(int l0, int l1, int l2, int l3, bool anti):
L0(l0), L1(l1), L2(l2), L3(l3), antiperiodicBC_L3(anti),
m0Squared(0.0), yukawa_t(0.0), yukawa_b(0.0), lambda(0.0), lambda_6(0.0), N_f(1),ignore_goldstone_modes(false),
rho(1.0), one_ov_two_rho(0.5), r(0.5),
eigenvalues_of_overlap(0), factor_for_eigenvalues(0), four_times_sum_of_sinSquared(0), factor_for_sum_of_sinSquared(0),
propagatorSum_withZeroMass(-1.0), propagatorSum_withHiggsMass(-1.0), actual_HiggsMassSquared(0.0),
numberOfDistingtMomenta_fermions(0), numberOfDistingtMomenta_bosons(0),
max_numberOfIterations(100), relative_Accuracy(1.0e-7), absolute_Accuracy(1.0e-7),minimizationAlgorithm(1),
minimizer(NULL), algorithmForMinimization(gsl_min_fminimizer_goldensection), minimizerInitialized(false), iterator_status(0),
fermionicContributions_are_valid(false)
{
	if( L0%2 + L1%2 + L2%2 + L3%2 != 0 )
	{
		std::cerr <<"Error, only even lattice extends supported!" <<std::endl;
		exit(EXIT_FAILURE);
	}
	fill_sinSquaredAndFactors();
	propagatorSum_withZeroMass=compute_propagatorSum(0.0);
}


//destructor
CEP_inBrokenPhase::~CEP_inBrokenPhase()
{
	delete [] eigenvalues_of_overlap;
	delete [] factor_for_eigenvalues;
	delete [] four_times_sum_of_sinSquared;
	delete [] factor_for_sum_of_sinSquared;
	gsl_min_fminimizer_free(minimizer);
}


//setting parameters
void CEP_inBrokenPhase::set_m0Squared( double new_m0Squared )
{ 
	if( m0Squared != new_m0Squared ){ m0Squared=new_m0Squared; reinitialize(); } 
}
void CEP_inBrokenPhase::set_yukawas( double y_t, double y_b )
{ 
	if( yukawa_t != y_t || yukawa_b != y_b )
	{ 
		yukawa_t=y_t; yukawa_b=y_b; if(eigenvalues_of_overlap!=0){fill_eigenvaluesAndFactors(); }
		//
		fermionicContributions.clear();
		fermionicContributions_are_valid=false;
	} 
		
	
	reinitialize();
}
void CEP_inBrokenPhase::set_lambda( double new_lambda ){ if(lambda != new_lambda){ lambda = new_lambda; reinitialize(); } }
void CEP_inBrokenPhase::set_lambda_6( double new_lambda_6 ){ if(lambda_6 != new_lambda_6 ){ lambda_6=new_lambda_6; reinitialize(); } }
void CEP_inBrokenPhase::set_N_f( int new_N_f){ if(N_f != new_N_f){ N_f=new_N_f; reinitialize(); } }

void CEP_inBrokenPhase::set_ignore_goldstone_modes(bool new_ignore){ ignore_goldstone_modes=new_ignore; }

void CEP_inBrokenPhase::set_rho( double new_rho )
{
	if(rho != new_rho){ rho=new_rho; one_ov_two_rho=0.5/rho; if(eigenvalues_of_overlap!=0){fill_eigenvaluesAndFactors(); } reinitialize(); }
}
void CEP_inBrokenPhase::set_r( double new_r)
{
	if(r != new_r){ r=new_r; if(eigenvalues_of_overlap!=0){fill_eigenvaluesAndFactors(); } reinitialize(); }
}

void CEP_inBrokenPhase::init_HiggsMassSquared()
{
	actual_HiggsMassSquared=(m0Squared>=0)? m0Squared : 0.0;
	if(actual_HiggsMassSquared==0.0){ propagatorSum_withHiggsMass=propagatorSum_withZeroMass; }
	else if( lambda!=0.0 || lambda_6!=0.0 )
	{
		propagatorSum_withHiggsMass=compute_propagatorSum( actual_HiggsMassSquared );
	}
}

void CEP_inBrokenPhase::set_HigsMassSquared( double mHSquared )
{
	if( mHSquared!=actual_HiggsMassSquared )
	{
		actual_HiggsMassSquared=mHSquared;
		propagatorSum_withHiggsMass=compute_propagatorSum( actual_HiggsMassSquared );
	}
}
	
	
	
void CEP_inBrokenPhase::reinitialize(){}


//getting parameters
double CEP_inBrokenPhase::get_m0Squared(){ return m0Squared; }
double CEP_inBrokenPhase::get_yukawa_t(){ return yukawa_t; }
double CEP_inBrokenPhase::get_yukawa_b(){ return yukawa_b; }
double CEP_inBrokenPhase::get_lambda(){ return lambda; }
double CEP_inBrokenPhase::get_lambda_6(){ return lambda_6;}
int CEP_inBrokenPhase::get_N_f(){ return N_f; }

double CEP_inBrokenPhase::get_rho(){ return rho; }
double CEP_inBrokenPhase::get_r(){ return r; }

double CEP_inBrokenPhase::get_actual_HiggsMassSquared(){ return actual_HiggsMassSquared; }

std::complex< double > CEP_inBrokenPhase::computeAnalyticalEigenvalue(double p0, double p1, double p2, double p3)
{
	//computes \nu^{+} from philipp's thesis (eq 3.9)
	// \nu(p) = \rho/a + \rho/a * [i\sqrt{\tilde{p}^2} + a*r *\hat{p}^2 - \rho/a]/[\srqt{\tilde{p}^2 + (a*r *\hat{p}^2 - \rho/a)^2}]
	//a is set to 1
	//\hat{p}^2 = 1/a^2 \sum_{\mu} 4 \sin{a p_\mu/2}^2
	//\tilde{p}^2 = 1/a^2 \sum_{\mu} \sin{a p_\mu}^2
	//NOTE: in philipp's thesis it's r/2 NOTE May be an ambiguity in the definition of r...
	double dummy=sin(p0); double p_tilde_sq=dummy*dummy;
	dummy=sin(p1); p_tilde_sq+=dummy*dummy;
	dummy=sin(p2); p_tilde_sq+=dummy*dummy;
	dummy=sin(p3); p_tilde_sq+=dummy*dummy;
	
	dummy=sin(0.5*p0); double p_hat_sq = dummy*dummy;
	dummy=sin(0.5*p1); p_hat_sq += dummy*dummy;
	dummy=sin(0.5*p2); p_hat_sq += dummy*dummy;
	dummy=sin(0.5*p3); p_hat_sq += dummy*dummy;
	p_hat_sq*=4.0;
	
	double one_ov_denom=1.0/sqrt(p_tilde_sq + (r*p_hat_sq - rho)*(r*p_hat_sq - rho) );
	
	return std::complex< double >(rho + rho*one_ov_denom*(r*p_hat_sq - rho) , rho*one_ov_denom*sqrt(p_tilde_sq) );
}



void CEP_inBrokenPhase::fill_eigenvaluesAndFactors()
{
	//fills the arrays with the eigenvalues of the overlap operator and a corresponding verctor with the
	//frequency of its occurence
	delete [] eigenvalues_of_overlap;
	delete [] factor_for_eigenvalues;
	
	const double PI(atan(1) * 4.0);
	double toAddForL3 = antiperiodicBC_L3 ? PI/static_cast< double >( L3 ) : 0.0;
	
	if(!(L0==L1 && L0==L2))
	{
		int L0_half=L0/2;
		int L1_half=L1/2;
		int L2_half=L2/2;
		int L3_half=L3/2;
		int dummyVolume(0);
		double two_PI_ov_L0=2.0*PI/static_cast< double >( L0 );
		double two_PI_ov_L1=2.0*PI/static_cast< double >( L1 );
		double two_PI_ov_L2=2.0*PI/static_cast< double >( L2 );
		double two_PI_ov_L3=2.0*PI/static_cast< double >( L3 );
		if(antiperiodicBC_L3)
		{
			dummyVolume=(L0_half+1)*(L1_half+1)*(L2_half+1)*(L3_half);
		}
		else
		{
			dummyVolume=(L0_half+1)*(L1_half+1)*(L2_half+1)*(L3_half+1);
		}
		numberOfDistingtMomenta_fermions=dummyVolume;
		std::cout <<"numberOfDistingtMomenta_fermions: "<<numberOfDistingtMomenta_fermions <<std::endl;
		eigenvalues_of_overlap = new std::complex< double > [numberOfDistingtMomenta_fermions];
		factor_for_eigenvalues = new double [numberOfDistingtMomenta_fermions];
		double fac_l0,fac_l1,fac_l2,fac_l3;
		size_t counter=0;
		for(int l0=0; l0<=L0_half; ++l0)
		{
			double p0 = two_PI_ov_L0 * l0;
			if(l0==0 || l0==L0_half){ fac_l0=1.0; }else{ fac_l0=2.0; }
			for(int l1=0; l1<=L1_half; ++l1)
			{
				double p1 = two_PI_ov_L1 * l1;
				if(l1==0 || l1==L1_half){ fac_l1=1.0; }else{ fac_l1=2.0; }
				for(int l2=0; l2<=L2_half; ++l2)
				{
					double p2 = two_PI_ov_L2 * l2;
					if(l2==0 || l2==L2_half){ fac_l2=1.0; }else{ fac_l2=2.0; }
					for(int l3=0; l3<=L3_half; ++l3)
					{
						double p3 = two_PI_ov_L3 * l3 + toAddForL3;
						if(!antiperiodicBC_L3 &&( l3==0 || l3==L3_half) ){ fac_l3=1.0; }
						else if(!antiperiodicBC_L3 && !( l3==0 || l3==L3_half) ){ fac_l3=2.0; }
						else if(antiperiodicBC_L3 && l3==L3_half){ continue; }
						else{ fac_l3=2.0; }
						
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( p0,p1,p2,p3);
						factor_for_eigenvalues[counter] = fac_l0*fac_l1*fac_l2*fac_l3;
						++counter;
					}
				}
			}
		}
	}
	else
	{
		int L_half=L0/2;
		int LT_half=L3/2;
		
		double two_PI_ov_L=2.0*PI/static_cast< double >( L2 );
		double two_PI_ov_LT=2.0*PI/static_cast< double >( L3 );
		size_t counter=0;
		//do it twice to count
		for(int l3=1; l3<LT_half; ++l3)
		{
			for(int l0=1; l0<L_half; ++l0)
			{
				for(int l1=l0+1; l1<L_half; ++l1)
				{
					for(int l2=l1+1; l2<L_half; ++l2)
					{
						//p,q,r 48
						l2=l2; //no compiler warning
						++counter;
					}
					{
						int l2=l1;
						l2=l2; 
						//p,q,q, 24
						++counter;
						l2=l0;
						//p,p,q, 24
						++counter;
						l2=L_half;
						//p,q,L/2 24
						++counter;
						l2=0;
						//0,p,q, 24
						++counter;
					}
				}
				{
					int l1=l0;
					{
						int l2=l1;
						l2=l2;//no compiler warning
						//p,p,p, 8
						++counter;
						l2=L_half;
						//p,p,L/2 12
						++counter;
						l2=0;
						//0,p,p, 12
						++counter;
					}
					l1=L_half;
					{
						int l2=L_half;
						l2=l2;//no compiler warning
						//p,L/2,L/2, 6
						++counter;
						l2=0;
						//0,p,L/2, 12
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,p  6
						++counter;
					}
				}
			}
			{
				int l0=L_half;
				{
					int l1=L_half;
					l1=l1;//no compiler warning
					{
						int l2=L_half;
						l2=l2;//no compiler warning
						//L/2,L/2,L/2, 1
						++counter;
						l2=0;
						//0,L/2,L/2, 3
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						++counter;
					}
				}
				l0=0;
				l0=l0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,0, 1
						++counter;
					}
				}
			}
		}
		//now do the same thing for l3=0 and l3=Lt_half
		for(int l3=0; l3<=LT_half; l3+=(antiperiodicBC_L3 ? LT_half*2 :LT_half))
		{
			for(int l0=1; l0<L_half; ++l0)
			{
				for(int l1=l0+1; l1<L_half; ++l1)
				{
					for(int l2=l1+1; l2<L_half; ++l2)
					{
						//p,q,r 48
						++counter;
					}
					{
						int l2=l1;
						l2=l2;//no compiler warning
						//p,q,q, 24
						++counter;
						l2=l0;
						//p,p,q, 24
						++counter;
						l2=L_half;
						//p,q,L/2 24
						++counter;
						l2=0;
						//0,p,q, 24
						++counter;
					}
				}
				{
					int l1=l0;
					l1=l1;//no compiler warning
					{
						int l2=l1;
						l2=l2;//no compiler warning
						//p,p,p, 8
						++counter;
						l2=L_half;
						//p,p,L/2 12
						++counter;
						l2=0;
						//0,p,p, 12
						++counter;
					}
					l1=L_half;
					{
						int l2=L_half;
						l2=l2;//no compiler warning
						//p,L/2,L/2, 6
						++counter;
						l2=0;
						//0,p,L/2, 12
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2;//no compiler warning
						l2=l2; //no compiler warning
						//0,0,p  6
						++counter;
					}
				}
			}
			{
				int l0=L_half;
				l0=l0;//no compiler warning
				{
					int l1=L_half;
					l1=l1;//no compiler warning
					{
						int l2=L_half;
						l2=l2;//no compiler warning
						//L/2,L/2,L/2, 1
						++counter;
						l2=0;
						//0,L/2,L/2, 3
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2;//no compiler warning
						//0,0,0, 1
						++counter;
					}
				}
			}
		}
		numberOfDistingtMomenta_fermions=counter;
		std::cout <<"numberOfDistingtMomenta_fermions: "<<numberOfDistingtMomenta_fermions <<std::endl;
		eigenvalues_of_overlap = new std::complex< double > [numberOfDistingtMomenta_fermions];
		factor_for_eigenvalues = new double [numberOfDistingtMomenta_fermions];
		//factor has to be doubled due to l3
		counter=0;
		for(int l3=1; l3<LT_half; ++l3)
		{
			for(int l0=1; l0<L_half; ++l0)
			{
				for(int l1=l0+1; l1<L_half; ++l1)
				{
					for(int l2=l1+1; l2<L_half; ++l2)
					{
						//p,q,r 48
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 96.0;
						++counter;
					}
					{
						int l2=l1;
						//p,q,q, 24
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 48.0;
						++counter;
						l2=l0;
						//p,p,q, 24
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 48.0;
						++counter;
						l2=L_half;
						//p,q,L/2 24
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 48.0;
						++counter;
						l2=0;
						//0,p,q, 24
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 48.0;
						++counter;
					}
				}
				{
					int l1=l0;
					{
						int l2=l1;
						//p,p,p, 8
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 16.0;
						++counter;
						l2=L_half;
						//p,p,L/2 12
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 24.0;
						++counter;
						l2=0;
						//0,p,p, 12
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 24.0;
						++counter;
					}
					l1=L_half;
					{
						int l2=L_half;
						//p,L/2,L/2, 6
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 12.0;
						++counter;
						l2=0;
						//0,p,L/2, 12
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 24.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						//0,0,p  6
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 12.0;
						++counter;
					}
				}
			}
			{
				int l0=L_half;
				{
					int l1=L_half;
					{
						int l2=L_half;
						//L/2,L/2,L/2, 1
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 2.0;
						++counter;
						l2=0;
						//0,L/2,L/2, 3
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 6.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						//0,0,L/2, 3
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 6.0;
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					{
						int l2=0;
						//0,0,0, 1
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = 2.0;
						++counter;
					}
				}
			}
		}
		double BC_factor=(antiperiodicBC_L3 ? 2.0 :1.0);
		//now do the same thing for l3=0 and l3=Lt_half
		for(int l3=0; l3<=LT_half; l3+=(antiperiodicBC_L3 ? LT_half*2 :LT_half))
		{
			for(int l0=1; l0<L_half; ++l0)
			{
				for(int l1=l0+1; l1<L_half; ++l1)
				{
					for(int l2=l1+1; l2<L_half; ++l2)
					{
						//p,q,r 48
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*48.0;
						++counter;
					}
					{
						int l2=l1;
						//p,q,q, 24
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*24.0;
						++counter;
						l2=l0;
						//p,p,q, 24
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*24.0;
						++counter;
						l2=L_half;
						//p,q,L/2 24
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*24.0;
						++counter;
						l2=0;
						//0,p,q, 24
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*24.0;
						++counter;
					}
				}
				{
					int l1=l0;
					{
						int l2=l1;
						//p,p,p, 8
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*8.0;
						++counter;
						l2=L_half;
						//p,p,L/2 12
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*12.0;
						++counter;
						l2=0;
						//0,p,p, 12
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*12.0;
						++counter;
					}
					l1=L_half;
					{
						int l2=L_half;
						//p,L/2,L/2, 6
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*6.0;
						++counter;
						l2=0;
						//0,p,L/2, 12
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*12.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						//0,0,p  6
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*6.0;
						++counter;
					}
				}
			}
			{
				int l0=L_half;
				{
					int l1=L_half;
					{
						int l2=L_half;
						//L/2,L/2,L/2, 1
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*1.0;
						++counter;
						l2=0;
						//0,L/2,L/2, 3
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*3.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						//0,0,L/2, 3
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*3.0;
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					{
						int l2=0;
						//0,0,0, 1
						eigenvalues_of_overlap[counter]=computeAnalyticalEigenvalue( two_PI_ov_L*l0, two_PI_ov_L*l1, two_PI_ov_L*l2, two_PI_ov_LT*l3+toAddForL3 );
						factor_for_eigenvalues[counter] = BC_factor*1.0;
						++counter;
					}
				}
			}
		}
	}
}



void CEP_inBrokenPhase::fill_sinSquaredAndFactors()
{
	//fills the arrays 4*sum_mu(sin(p_mu/2)^2) and the corresponding
	//frequency of its occurence
	delete [] four_times_sum_of_sinSquared;
	delete [] factor_for_sum_of_sinSquared;
	
	const double PI(atan(1) * 4.0);
	
	if(  !(L0==L1 && L0==L2))
	{
		
		int L0_half=L0/2, L1_half=L1/2, L2_half=L2/2, L3_half=L3/2;
		double PI_ov_L0=PI/static_cast< double >( L0 );
		double PI_ov_L1=PI/static_cast< double >( L1 );
		double PI_ov_L2=PI/static_cast< double >( L2 );
		double PI_ov_L3=PI/static_cast< double >( L3 );
		double *sinSquared_p0_half, *sinSquared_p1_half, *sinSquared_p2_half, *sinSquared_p3_half;
		sinSquared_p0_half= new double [L0_half+1];
		sinSquared_p1_half= new double [L1_half+1];
		sinSquared_p2_half= new double [L2_half+1];
		sinSquared_p3_half= new double [L3_half+1];
		
		//collect the sin()^2 first
		for(int i=0; i< L0_half+1 ; ++i){ sinSquared_p0_half[i]=sin( PI_ov_L0*i ); sinSquared_p0_half[i]*=sinSquared_p0_half[i]; }
		for(int i=0; i< L1_half+1 ; ++i){ sinSquared_p1_half[i]=sin( PI_ov_L1*i ); sinSquared_p1_half[i]*=sinSquared_p1_half[i]; }
		for(int i=0; i< L2_half+1 ; ++i){ sinSquared_p2_half[i]=sin( PI_ov_L2*i ); sinSquared_p2_half[i]*=sinSquared_p2_half[i]; }
		for(int i=0; i< L3_half+1 ; ++i){ sinSquared_p3_half[i]=sin( PI_ov_L3*i ); sinSquared_p3_half[i]*=sinSquared_p3_half[i]; }
		int dummyVolume=(L0_half+1)*(L1_half+1)*(L2_half+1)*(L3_half+1)-1;
		numberOfDistingtMomenta_bosons=dummyVolume;
		std::cout <<"numberOfDistingtMomenta_bosons: "<<numberOfDistingtMomenta_bosons <<std::endl;
		four_times_sum_of_sinSquared = new double [dummyVolume];
		factor_for_sum_of_sinSquared = new double [dummyVolume];
		double fac_l0,fac_l1,fac_l2,fac_l3;
		
			
		size_t counter=0;
		for(int l0=0; l0<=L0_half; ++l0)
		{
			if(l0==0 || l0==L0_half){ fac_l0=1.0; }else{ fac_l0=2.0; }
			for(int l1=0; l1<=L1_half; ++l1)
			{
				if(l1==0 || l1==L1_half){ fac_l1=1.0; }else{ fac_l1=2.0; }
				for(int l2=0; l2<=L2_half; ++l2)
				{
					if(l2==0 || l2==L2_half){ fac_l2=1.0; }else{ fac_l2=2.0; }
					for(int l3=((l0+l1+l2)?(0):(1)); l3<=L3_half; ++l3)
					{
						if(l3==0 || l3==L3_half){ fac_l3=1.0; }else{ fac_l3=2.0; }
						
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p0_half[l0] + sinSquared_p0_half[l1] + sinSquared_p0_half[l2] + sinSquared_p0_half[l3] );
						factor_for_sum_of_sinSquared[counter] = fac_l0*fac_l1*fac_l2*fac_l3;
						++counter;
					}
				}
			}
		}
		delete [] sinSquared_p0_half; delete [] sinSquared_p1_half; delete [] sinSquared_p2_half; delete [] sinSquared_p3_half;
	}
	else
	{
		int L_half=L0/2;
		int LT_half=L3/2;
		double PI_ov_L=PI/static_cast< double >( L0 );
		double PI_ov_LT=PI/static_cast< double >( L3 );
		double *sinSquared_p_half, *sinSquared_pT_half;
		sinSquared_p_half = new double [L_half +1];
		sinSquared_pT_half = new double [LT_half +1];
		
		//collect the sin()^2 first
		for(int i=0; i< L_half+1 ; ++i){ sinSquared_p_half[i]=sin( PI_ov_L*i ); sinSquared_p_half[i]*=sinSquared_p_half[i]; }
		for(int i=0; i< LT_half+1 ; ++i){ sinSquared_pT_half[i]=sin( PI_ov_LT*i ); sinSquared_pT_half[i]*=sinSquared_pT_half[i]; }
		
		size_t counter=0;
		//do it twice to count
		for(int l3=1; l3<LT_half; ++l3)
		{
			for(int l0=1; l0<L_half; ++l0)
			{
				for(int l1=l0+1; l1<L_half; ++l1)
				{
					for(int l2=l1+1; l2<L_half; ++l2)
					{
						//p,q,r 48
						l2=l2; //no compiler warning
						++counter;
					}
					{
						int l2=l1;
						l2=l2;//no compiler warning
						//p,q,q, 24
						++counter;
						l2=l0;
						//p,p,q, 24
						++counter;
						l2=L_half;
						//p,q,L/2 24
						++counter;
						l2=0;
						//0,p,q, 24
						++counter;
					}
				}
				{
					int l1=l0;
					{
						int l2=l1;
						l2=l2;//no compiler warning
						//p,p,p, 8
						++counter;
						l2=L_half;
						//p,p,L/2 12
						++counter;
						l2=0;
						//0,p,p, 12
						++counter;
					}
					l1=L_half;
					{
						int l2=L_half;
						l2=l2;//no compiler warning
						//p,L/2,L/2, 6
						++counter;
						l2=0;
						//0,p,L/2, 12
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,p  6
						++counter;
					}
				}
			}
			{
				int l0=L_half;
				l0=l0;//no compiler warning
				{
					int l1=L_half;
					l1=l1;//no compiler warning
					{
						int l2=L_half;
						l2=l2;//no compiler warning
						//L/2,L/2,L/2, 1
						++counter;
						l2=0;
						//0,L/2,L/2, 3
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,0, 1
						++counter;
					}
				}
			}
		}
		//now do the same thing for l3=0 and l3=Lt_half
		for(int l3=0; l3<=LT_half; l3+=LT_half )
		{
			for(int l0=1; l0<L_half; ++l0)
			{
				for(int l1=l0+1; l1<L_half; ++l1)
				{
					for(int l2=l1+1; l2<L_half; ++l2)
					{
						//p,q,r 48
						++counter;
					}
					{
						int l2=l1;
						l2=l2;//no compiler warning
						//p,q,q, 24
						++counter;
						l2=l0;
						//p,p,q, 24
						++counter;
						l2=L_half;
						//p,q,L/2 24
						++counter;
						l2=0;
						//0,p,q, 24
						++counter;
					}
				}
				{
					int l1=l0;
					{
						int l2=l1;
						l2=l2;//no compiler warning
						//p,p,p, 8
						++counter;
						l2=L_half;
						//p,p,L/2 12
						++counter;
						l2=0;
						//0,p,p, 12
						++counter;
					}
					l1=L_half;
					{
						int l2=L_half;
						l2=l2;//no compiler warning
						//p,L/2,L/2, 6
						++counter;
						l2=0;
						//0,p,L/2, 12
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2;//no compiler warning
						l2=l2; //no compiler warning
						//0,0,p  6
						++counter;
					}
				}
			}
			{
				int l0=L_half;
				l0=l0;//no compiler warning
				{
					int l1=L_half;
					l1=l1;//no compiler warning
					{
						int l2=L_half;
						l2=l2;//no compiler warning
						//L/2,L/2,L/2, 1
						{++counter; }
						l2=0;
						//0,L/2,L/2, 3
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2;//no compiler warning
						//0,0,0, 1
						if(l3 != 0){++counter;}
					}
				}
			}
		}
		numberOfDistingtMomenta_bosons=counter;
		std::cout <<"numberOfDistingtMomenta_bosons: "<<numberOfDistingtMomenta_bosons <<std::endl;
		four_times_sum_of_sinSquared = new double [numberOfDistingtMomenta_bosons];
		factor_for_sum_of_sinSquared = new double [numberOfDistingtMomenta_bosons];
		counter=0;
		for(int l3=1; l3<LT_half; ++l3)
		{
			for(int l0=1; l0<L_half; ++l0)
			{
				for(int l1=l0+1; l1<L_half; ++l1)
				{
					for(int l2=l1+1; l2<L_half; ++l2)
					{
						//p,q,r 48
						l2=l2; //no compiler warning
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 96.0;
						++counter;
					}
					{
						int l2=l1;
						//p,q,q, 24
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 48.0;
						++counter;
						l2=l0;
						//p,p,q, 24
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 48.0;
						++counter;
						l2=L_half;
						//p,q,L/2 24
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 48.0;
						++counter;
						l2=0;
						//0,p,q, 24
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 48.0;
						++counter;
					}
				}
				{
					int l1=l0;
					{
						int l2=l1;
						//p,p,p, 8
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 16.0;
						++counter;
						l2=L_half;
						//p,p,L/2 12
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 24.0;
						++counter;
						l2=0;
						//0,p,p, 12
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 24.0;
						++counter;
					}
					l1=L_half;
					{
						int l2=L_half;
						//p,L/2,L/2, 6
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 12.0;
						++counter;
						l2=0;
						//0,p,L/2, 12
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 24.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,p  6
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 12.0;
						++counter;
					}
				}
			}
			{
				int l0=L_half;
				{
					int l1=L_half;
					{
						int l2=L_half;
						//L/2,L/2,L/2, 1
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 2.0;
						++counter;
						l2=0;
						//0,L/2,L/2, 3
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 6.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 6.0;
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,0, 1
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 2.0;
						++counter;
					}
				}
			}
		}
		//now do the same thing for l3=0 and l3=Lt_half
		for(int l3=0; l3<=LT_half; l3+=LT_half )
		{
			for(int l0=1; l0<L_half; ++l0)
			{
				for(int l1=l0+1; l1<L_half; ++l1)
				{
					for(int l2=l1+1; l2<L_half; ++l2)
					{
						//p,q,r 48
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 48.0;
						++counter;
					}
					{
						int l2=l1;
						//p,q,q, 24
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 24.0;
						++counter;
						l2=l0;
						//p,p,q, 24
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 24.0;
						++counter;
						l2=L_half;
						//p,q,L/2 24
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 24.0;
						++counter;
						l2=0;
						//0,p,q, 24
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 24.0;
						++counter;
					}
				}
				{
					int l1=l0;
					{
						int l2=l1;
						//p,p,p, 8
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 8.0;
						++counter;
						l2=L_half;
						//p,p,L/2 12
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 12.0;
						++counter;
						l2=0;
						//0,p,p, 12
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 12.0;
						++counter;
					}
					l1=L_half;
					{
						int l2=L_half;
						//p,L/2,L/2, 6
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 6.0;
						++counter;
						l2=0;
						//0,p,L/2, 12
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 12.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,p  6
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 6.0;
						++counter;
					}
				}
			}
			{
				int l0=L_half;
				{
					int l1=L_half;
					{
						int l2=L_half;
						//L/2,L/2,L/2, 1
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 1.0;
						++counter; 
						l2=0;
						//0,L/2,L/2, 3
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 3.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
						factor_for_sum_of_sinSquared[counter] = 3.0;
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2;//no compiler warning
						//0,0,0, 1
						if(l3 != 0)
						{
							four_times_sum_of_sinSquared[counter]=4.0*( sinSquared_p_half[l0] + sinSquared_p_half[l1] + sinSquared_p_half[l2] + sinSquared_pT_half[l3] );
							factor_for_sum_of_sinSquared[counter] = 1.0;
							++counter;
						}
					}
				}
			}
		}
		delete [] sinSquared_p_half; delete [] sinSquared_pT_half;
	}
}



double CEP_inBrokenPhase::compute_propagatorSum( double massSquared )
{
	if( four_times_sum_of_sinSquared == 0 ){ fill_sinSquaredAndFactors() ; }; //should not happen...
	double dummy(0.0);
	for(int counter=0; counter<numberOfDistingtMomenta_bosons; ++counter)
	{
		dummy+=factor_for_sum_of_sinSquared[counter]/(four_times_sum_of_sinSquared[counter] + massSquared);
	}
	dummy/=static_cast< double >(L0); dummy/=static_cast< double >(L1); dummy/=static_cast< double >(L2); dummy/=static_cast< double >(L3); 
	return dummy;
}

double CEP_inBrokenPhase::compute_CEP_inBrokenPhase( double value )
{
	//this function computes the constrained effective potential in the broken Phase:
	// U(v) = U_tree + U_ferm + U_1st
	// U_tree(v) = 1/2 * m0Squared * v^2  + lambda * v^4  + lambda_6 * v^6
	// U_ferm(v) = -2*N_f/V * sum_p[ log ( z_t * z_t^* )  +  log ( z_b * z_b^* ) ]
	//              z_t/b = nu  +  y_t/b * v * ( 1 - 1/(2 rho)*nu)  with nu being the eigenvalue of the overlap 
	// U_1st(v) = lambda * 6 * v^2 * ( P_G + P_H )  +  lambda_6 * ( v^2 * ( 45*P_G^2  +  54*P_G*P_H  +  45*P_H^2 )  +  v^4 * ( 9*P_G + 15*P_H ) )
	
	double result( compute_treeLevel( value ) );
	if(yukawa_t != 0.0 || yukawa_b != 0.0){ result += compute_fermionicContribution( value ); }
	if(lambda != 0.0 || lambda_6 != 0.0){ result += compute_firstOrderInLambdas( value ); }
	return result;
}


double CEP_inBrokenPhase::compute_treeLevel( double value )
{
	return 0.5*m0Squared*value*value + lambda*value*value*value*value + lambda_6*value*value*value*value*value*value;
}

double CEP_inBrokenPhase::compute_fermionicContribution( double value )
{
	// U_ferm(v) = -2*N_f/V * sum_p[ log ( z_t * z_t^* )  +  log ( z_b * z_b^* ) ]
	//              z_t/b = nu  +  y_t/b * v * ( 1 - 1/(2 rho)*nu)  with nu being the eigenvalue of the overlap 
	
	if(eigenvalues_of_overlap==0){ fill_eigenvaluesAndFactors(); }
	
	if(fermionicContributions_are_valid)
	{
		//NOTE lower is a bad name, lower is equal or lager than value!
		std::map< double, double >::const_iterator lower=fermionicContributions.lower_bound( value );
		if(lower->first==value){ return lower->second; }
		else if(lower != fermionicContributions.end() && lower != fermionicContributions.begin())
		{
			std::map< double, double >::const_iterator lower_min_one=lower;
			--lower_min_one;
			double deriv( (lower->second - lower_min_one->second)/(lower->first - lower_min_one->first) );
			return (lower_min_one->second + (value - lower_min_one->first)*deriv);
		}
	}
	
	long double dummy(0.0);
// 	double dummy(0.0);
	double fac_t=yukawa_t*value;
	double fac_b=yukawa_b*value;
	std::complex< double > z_t, z_b;
	for(int counter=0; counter <numberOfDistingtMomenta_fermions; ++counter)
	{
		z_t=eigenvalues_of_overlap[counter] + fac_t*(1.0 - eigenvalues_of_overlap[counter]*one_ov_two_rho);
		z_b=eigenvalues_of_overlap[counter] + fac_b*(1.0 - eigenvalues_of_overlap[counter]*one_ov_two_rho);
		dummy += static_cast< long double >( factor_for_eigenvalues[counter]*( log( real( z_t * conj(z_t) ) * real( z_b * conj(z_b) ) ) ) );
// 		dummy += ( factor_for_eigenvalues[counter]*( log( real( z_t * conj(z_t) ) ) + log( real( z_b * conj(z_b) ) ) ) );
	}
	dummy*=-2.0*static_cast< double >(N_f);
	dummy/= static_cast< double >(L0); dummy/= static_cast< double >(L1); dummy/= static_cast< double >(L2); dummy/= static_cast< double >(L3);
	return static_cast< double >(dummy);
// 	return (dummy);
}

double CEP_inBrokenPhase::compute_firstOrderInLambdas( double value )
{
	// U_1st(v) = lambda * 6 * v^2 * ( P_G + P_H )  +  lambda_6 * ( v^2 * ( 45*P_G^2  +  54*P_G*P_H  +  45*P_H^2 )  +  v^4 * ( 9*P_G + 15*P_H ) )
	if(  ignore_goldstone_modes  )
	{
			return 6.0 * lambda * value*value * (propagatorSum_withHiggsMass ) 
	                             + lambda_6 * ( value*value * 45.0*propagatorSum_withHiggsMass*propagatorSum_withHiggsMass
			              + value*value*value*value * 15.0*propagatorSum_withHiggsMass );
	}
	else
	{
		return 6.0 * lambda * value*value * (propagatorSum_withHiggsMass + propagatorSum_withZeroMass) 
	       + lambda_6 * ( value*value * ( 45.0*propagatorSum_withZeroMass*propagatorSum_withZeroMass 
	                                    + 54.0*propagatorSum_withZeroMass*propagatorSum_withHiggsMass
	                                    + 45.0*propagatorSum_withHiggsMass*propagatorSum_withHiggsMass)
			              + value*value*value*value * ( 9.0*propagatorSum_withZeroMass + 15.0*propagatorSum_withHiggsMass ) );
	}
}

double CEP_inBrokenPhase::compute_CEP_inBrokenPhase_secondDerivative( double value )
{
	//computes U''(v) = U_tree'' + U_ferm'' + U_1st''
	// U_tree''(v) =  m0Squared  + 12*lambda * v^2  + 30*lambda_6 * v^4
	// U_ferm''(v) = U_fpp=sum_p{-2*y_t^2*Re[w^2/z_t^2]} + , sum_p{-2*y_b^2*Re[w^2/z_b^2]} with w = 1 - 1/(2 rho)*nu, z_t/b=(nu +  y_t/b * v *w)
	// U_1st''(v) = lambda * 12 ( P_G + P_H )  +  lambda_6 * ( ( 90*P_G^2  +  108*P_G*P_H  +  90*P_H^2 )  +  v^2 * ( 108*P_G + 180*P_H ) )
	double result( compute_treeLevel_secondDerivative( value ) );
	if(yukawa_t != 0.0 || yukawa_b != 0.0){ result += compute_fermionicContribution_secondDerivative( value ); }
	if(lambda != 0.0 || lambda_6 != 0.0){ result += compute_firstOrderInLambdas_secondDerivative( value ); }
	return result;
}
	
double CEP_inBrokenPhase::compute_treeLevel_secondDerivative( double value )
{
	// U_tree''(v) =  m0Squared  + 12*lambda * v^2  + 30*lambda_6 * v^4
	return m0Squared + 12.0*lambda*value*value + 30.0*lambda_6*value*value*value*value;
}
	
double CEP_inBrokenPhase::compute_fermionicContribution_secondDerivative( double value )
{
	// U_ferm''(v) = U_fpp=sum_p{-2*y_t^2*Re[w^2/z_t^2]} + , sum_p{-2*y_b^2*Re[w^2/z_b^2]} with w = 1 - 1/(2 rho)*nu, z_t/b=(nu +  y_t/b * v *w)
	if(eigenvalues_of_overlap==0){ fill_eigenvaluesAndFactors(); }
	double dummy(0.0);
	double fac_t=yukawa_t*value;
	double fac_b=yukawa_b*value;
	double y_t_Squared(yukawa_t*yukawa_t);
	double y_b_Squared(yukawa_b*yukawa_b);
	std::complex< double > z_t, z_b,w;
	for(int counter=0; counter <numberOfDistingtMomenta_fermions; ++counter)
	{
		w=(1.0 - eigenvalues_of_overlap[counter]*one_ov_two_rho);
		z_t=eigenvalues_of_overlap[counter] + fac_t*w;
		z_b=eigenvalues_of_overlap[counter] + fac_b*w;
		dummy += -2.0*factor_for_eigenvalues[counter]*( y_t_Squared * real( w*w/(z_t*z_t) ) + y_b_Squared * real( w*w/(z_b*z_b) ) );
	}
	dummy*=-2.0*static_cast< double >(N_f);
	dummy/= static_cast< double >(L0); dummy/= static_cast< double >(L1); dummy/= static_cast< double >(L2); dummy/= static_cast< double >(L3);
	return dummy;
}

double CEP_inBrokenPhase::compute_firstOrderInLambdas_secondDerivative( double value )
{
	// U_1st''(v) = lambda * 12 ( P_G + P_H )  +  lambda_6 * ( ( 90*P_G^2  +  108*P_G*P_H  +  90*P_H^2 )  +  v^2 * ( 108*P_G + 180*P_H ) )
	if( ignore_goldstone_modes )
	{
		return 12.0 * lambda * (propagatorSum_withHiggsMass ) 
	                     + lambda_6 * ( 90.0*propagatorSum_withHiggsMass*propagatorSum_withHiggsMass
			                + value*value * ( 180.0*propagatorSum_withHiggsMass ) );
	}
	else
	{
		return 12.0 * lambda * (propagatorSum_withHiggsMass + propagatorSum_withZeroMass) 
	                     + lambda_6 * ( ( 90.0*propagatorSum_withZeroMass*propagatorSum_withZeroMass 
								 + 108.0*propagatorSum_withZeroMass*propagatorSum_withHiggsMass
								 + 90.0*propagatorSum_withHiggsMass*propagatorSum_withHiggsMass)
			                + value*value * ( 108.0*propagatorSum_withZeroMass + 180.0*propagatorSum_withHiggsMass ) );
	}
}

bool CEP_inBrokenPhase::load_fermionicContribution( const std::string &fileName )
{
	fermionicContributions.clear();
	fermionicContributions_are_valid=false;
	double tolForCheck=1.0e-14;
	if(yukawa_t==0.0 && yukawa_b==0.0)
	{
		std::cerr <<"Error, no yukawa couplings set" <<std::endl;
		return false;
	}
	std::ifstream inputFile(fileName.c_str());
	std::string line,word;
	double dummy1(0.0),dummy2(0.0);
	if(!(inputFile.good()))
	{
		std::cerr <<"Error opening inputfile " <<fileName <<std::endl;
		return false;
	}
	std::map< double, double >::iterator hintIter=fermionicContributions.begin();
	while(inputFile)
	{
		getline(inputFile, line);
		if(line.size()==0 || line[0] =='#' || line.find_first_not_of(' ') == std::string::npos){ continue; }
		std::istringstream strm(line);
		if(!(strm >> dummy1 >> dummy2))
		{
			std::cerr <<"Error reading list of fermionic contributions from " <<fileName <<std::endl;
			fermionicContributions.clear();
			return false;
		}
		hintIter=fermionicContributions.insert( hintIter, std::make_pair(dummy1, dummy2) );
		if(inputFile.eof()){ break; }
	}
	inputFile.close(); inputFile.clear();
	
	//check first and last value for agreement...
	double resBegin=compute_fermionicContribution(fermionicContributions.begin()->first);
	double resEnd  =compute_fermionicContribution(fermionicContributions.rbegin()->first);
	//debug
	/*std::cout.precision(15);
	std::cout <<"resBegin = " <<resBegin <<std::endl;
	std::cout <<"from List= " <<fermionicContributions.begin()->second <<std::endl;
	std::cout <<"rel dev  = " <<(0.5*std::abs( resBegin - fermionicContributions.begin()->second  )/(resBegin+fermionicContributions.begin()->second) ) <<std::endl; */
	
	if( 0.5*std::abs( (resBegin - fermionicContributions.begin()->second )/(resBegin+fermionicContributions.begin()->second) ) > tolForCheck || 0.5*std::abs( (resEnd- fermionicContributions.rbegin()->second )/(resEnd+fermionicContributions.rbegin()->second) ) > tolForCheck )
	{
		std::cerr <<"Error, loaded list does not agree within accuracy" <<std::endl;
		fermionicContributions.clear();
		return false;
	}
	fermionicContributions_are_valid=true;
	return true;
}//load_fermionicContribution




void CEP_inBrokenPhase::set_max_numberOfIterations(int new_max){ max_numberOfIterations=new_max; }
void CEP_inBrokenPhase::set_relative_Accuracy( double new_rel_acc){ relative_Accuracy=new_rel_acc; }
void CEP_inBrokenPhase::set_absolute_Accuracy( double new_abs_acc){ absolute_Accuracy=new_abs_acc; }

void CEP_inBrokenPhase::set_minimizationAlgorithm( int new_algorithm )
{
	// 1 - gsl_min_fminimizer_goldensection
	// 2 - gsl_min_fminimizer_brent
	// 3 - gsl_min_fminimizer_quad_golden
	switch(new_algorithm)
	{
		case 1: 
			minimizationAlgorithm=new_algorithm;
			algorithmForMinimization = (gsl_min_fminimizer_goldensection);
			break;
		case 2: 
			minimizationAlgorithm=new_algorithm;
			algorithmForMinimization = (gsl_min_fminimizer_brent);
			break;
		case 3: 
			minimizationAlgorithm=new_algorithm;
			algorithmForMinimization = (gsl_min_fminimizer_quad_golden);
			break;
		default:
			std::cerr <<"Error, algorithm " <<new_algorithm <<" is not valid in CEP_inBrokenPhase::set_minimizationAlgorithm(int new_alg) " <<std::endl;
			std::cerr <<"Use one of the following:" <<std::endl;
			std::cerr <<"1 = gsl_min_fminimizer_goldensection" <<std::endl;
			std::cerr <<"2 = gsl_min_fminimizer_brent" <<std::endl;
			std::cerr <<"3 = gsl_min_fminimizer_quad_golden" <<std::endl;
			exit(EXIT_FAILURE);
	}
}

bool CEP_inBrokenPhase::reInitialize_minimizer( double minimum, double lower, double upper )
{
	if(!minimizerInitialized){ return false; }
	int dummy( gsl_min_fminimizer_set( minimizer, &functionHandler, minimum, lower, upper ) );
	if(dummy==GSL_EINVAL){ return false; }
	iterator_status=0;
	return true;
}

bool CEP_inBrokenPhase::initialize_minimizer( double minimum, double lower, double upper )
{
	if(minimizerInitialized){ return reInitialize_minimizer( minimum, lower, upper ); }
	minimizer = gsl_min_fminimizer_alloc( algorithmForMinimization );
	functionHandler.function = &wrapper_compute_CEP_inBrokenPhase_gsl;
	functionHandler.params= (void *) this;
	int dummy=gsl_min_fminimizer_set( minimizer, &functionHandler, minimum, lower, upper );
	if(dummy==GSL_EINVAL){ return false; }
	minimizerInitialized=true;
	return true;
}




// int CEP_inBrokenPhase::determine_startingAndInitialize( double lower, double upper, double step )
void CEP_inBrokenPhase::determine_startingPoints( double inLower, double inUpper, double inStep, double &outMinimum, double &outLower,double &outUpper)
{
	const int MAXSCAN=1000;
	if( (inUpper - inLower <= 0.0 || inStep <=0.0) )
	{
		std::cerr <<"Error, inconsitent scan range in CEP_inBrokenPhase::determine_startingAndInitialize" <<std::endl;
		exit(EXIT_FAILURE);
	}
	int numberOfEntries=static_cast< int >( (inUpper-inLower)/inStep + 1.5);
	if( numberOfEntries > MAXSCAN )
	{
		std::cerr << "Error, to many testValues in CEP_inBrokenPhase::determine_startingAndInitialize." <<std::endl;
		std::cerr << "Change MAXSCAN if neccessary" <<std::endl;
	}
	std::set< double > trialPoints;
	for( int i=0; i<numberOfEntries; ++i)
	{
		trialPoints.insert( inLower + static_cast< double >(i) *inStep);
	}
	std::map< double, double > trial_functionAndPoints;
	for( std::set< double >::const_iterator iter=trialPoints.begin(); iter!=trialPoints.end(); ++iter )
	{
		trial_functionAndPoints.insert( std::make_pair( compute_CEP_inBrokenPhase( *iter ), *iter ) );
	}
	//debug
	/*for( std::map< double, double >::const_iterator iter=trial_functionAndPoints.begin(); iter !=trial_functionAndPoints.end(); ++iter)
	{
		std::cout <<"trial: value= " <<iter->second <<"   pot: " <<iter->first <<std::endl;
	}*/
	double lowest=trial_functionAndPoints.begin()->second;
	if( lowest==*trialPoints.begin() || lowest==*trialPoints.rbegin() )
	{
		outMinimum=lowest; outLower=lowest; outUpper=lowest;
	}
	else
	{
		outMinimum=lowest; outLower=lowest-inStep; outUpper=lowest+inStep;
	}
}

int CEP_inBrokenPhase::iterate_minimizer()
{
	if(!minimizerInitialized)
	{
		std::cerr <<"Error, minimizer not initialized in CEP_inBrokenPhase::iterate_Minimizer()" <<std::endl;
		exit(EXIT_FAILURE);
	}
	iterator_status = gsl_min_fminimizer_iterate(minimizer);
	return iterator_status;
}

double CEP_inBrokenPhase::get_actual_minimum()
{
	if(!minimizerInitialized)
	{
		std::cerr <<"Error, minimizer not initialized in CEP_inBrokenPhase::get_actual_minimum()" <<std::endl;
		exit(EXIT_FAILURE);
	}
	return gsl_min_fminimizer_x_minimum(minimizer);
}

double CEP_inBrokenPhase::get_potentialAtMinimum()
{
	if(!minimizerInitialized)
	{
		std::cerr <<"Error, minimizer not initialized in CEP_inBrokenPhase::get_actual_minimum()" <<std::endl;
		exit(EXIT_FAILURE);
	}
	return gsl_min_fminimizer_f_minimum(minimizer);
}

void CEP_inBrokenPhase::get_actual_Interval( double &outMinimum, double &outLower,double &outUpper )
{
	//returns current best interval of minimizer
	if(!minimizerInitialized)
	{
		std::cerr <<"Error, minimizer not initialized in CEP_inBrokenPhase::get_actual_Interval" <<std::endl;
		exit(EXIT_FAILURE);
	}
	outMinimum = gsl_min_fminimizer_x_minimum(minimizer);
	outLower   = gsl_min_fminimizer_x_lower(minimizer);
	outUpper   = gsl_min_fminimizer_x_upper(minimizer);
}

bool CEP_inBrokenPhase::check_convergence()
{
	if( gsl_min_test_interval(gsl_min_fminimizer_x_lower(minimizer), gsl_min_fminimizer_x_upper(minimizer), absolute_Accuracy, relative_Accuracy) == GSL_SUCCESS )
	{
		return true;
	}
	return false;
}

int CEP_inBrokenPhase::iterate_minimizer_until_convergence()
{
	int counter=0;
	for(int i=0; i<max_numberOfIterations; ++i)
	{
		iterate_minimizer();
		++counter;
		if( check_convergence() )
		{
			return counter;
		}
		if(iterator_status==GSL_EBADFUNC)
		{
			std::cerr <<"Inf or NaN occured in itartion" <<std::endl;
			return -2;
		}
		if(iterator_status==GSL_FAILURE)
		{
			std::cerr <<"Iterator could not improved without reaching convergence" <<std::endl;
			return -3;
		}
	}
	return -1;
}




//wrapper for gsl since member function does not work
double wrapper_compute_CEP_inBrokenPhase_gsl(double value, void *params)
{
	CEP_inBrokenPhase *CEP = (CEP_inBrokenPhase *)params;
	return CEP->compute_CEP_inBrokenPhase(value);
}

