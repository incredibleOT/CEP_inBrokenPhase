#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>


#include "CEP_inBrokenPhase.h"


using std::cout;
using std::cerr;
using std::endl;

int main(int narg,char **arg)
{
//generates a list of fermionic contributrion
//input is    Ls   Lt    y_t    y_ratio     v_min    v_max    v_step     output
	if( narg != 9 )
	{
		cerr <<"Error, start with the folowing arguments:" <<endl;
		cerr <<"L_s   L_t   yukawa_t    yukawa_ratio    v_min   v_max   v_step   output_file" <<endl;
		exit(EXIT_FAILURE);
	}
	
	int Ls=atoi( arg[1] );
	int Lt=atoi( arg[2] );
	double y_t=atof( arg[3] );
	double y_b=y_t*atof( arg[4] );
	long double v_min=atof( arg[5] );
	long double v_max=atof( arg[6] );
	long double v_step=atof( arg[7] );
	std::string output_file_name(arg[8]);
	
	bool anti(false);
	
	if(v_min >= v_max || v_step <= 0.0)
	{
		cerr <<"Error, inconsistent range" <<endl;
	
	}
	CEP_inBrokenPhase CEP(Ls,Ls,Ls,Lt,anti);
	
	CEP.set_yukawas(y_t,y_b);
	
	std::map< double, double > results;
	
	long double actual_v(v_min);
	
	cout.precision(8);
	while(actual_v <= v_max+0.5*v_step)
	{
// 		cout <<"v=" <<actual_v <<endl;
		results.insert( std::make_pair( actual_v, CEP.compute_fermionicContribution( actual_v ) ) );
		actual_v+=v_step;
	}
	
		
	//output
	cout <<"write output to: " <<endl <<output_file_name <<endl;
	std::ofstream outputFile( output_file_name.c_str() );
	if(!(outputFile.good()) )
	{
		cerr <<"Error opening " <<endl <<output_file_name <<endl;
		exit(EXIT_FAILURE);
	}
	
	outputFile.precision(15);
	outputFile <<"# fermionic contribiution for the CEP in the broken phase" <<endl;
	outputFile <<"# Ls=" <<Ls  <<"   Lt=" <<Lt <<"   antiperiodic BC in time=" <<anti <<endl;
	outputFile <<"# y_t=" <<y_t <<"   y_b=" <<y_b <<endl;
	outputFile <<"# v_min=" <<v_min <<"   v_max=" <<v_max <<"   v_step=" <<v_step <<endl;
	outputFile <<"# format:  v    U_f" <<endl;
	for(std::map< double, double >::const_iterator iter=results.begin(); iter!=results.end(); ++iter)
	{
		outputFile <<iter->first <<"  " <<iter->second <<endl;
	}
	
	
	
}