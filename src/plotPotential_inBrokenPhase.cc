//This program will plot the CEP in the broken phase for a given range of fieldvalues

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "constrainedEffectivePotential_inBrokenPhase.h"
#include "plotPotential_inBrokenPhase_helper.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int narg,char **arg)
{
	std::map< std::string, double > parametersDouble;
	std::map< std::string, int > parametersInt;
	std::map< std::string, std::string > parametersString;
	std::map< std::string, bool > parametersIsSet;
	
	plotPotential_inBrokenPhase_helper::prepareParameterMaps(parametersDouble, parametersInt, parametersString, parametersIsSet);
	size_t numberOfParameters=parametersDouble.size()+parametersInt.size()+parametersString.size();
	
	if(narg!=2)
	{
		cerr <<"Error, start program with:" <<endl;
		cerr <<arg[0] <<"   inputFile" <<endl;
		exit(EXIT_FAILURE);
	}
	
	if( !plotPotential_inBrokenPhase_helper::loadParameterMapsFromFile(parametersDouble, parametersInt, parametersString, parametersIsSet, arg[1]) )
	{
		cerr <<"Error loading input file" <<endl <<arg[1] <<endl;
		exit(EXIT_FAILURE);
	}
	
	cout <<endl <<"Parameters loaded:" <<endl;
	plotPotential_inBrokenPhase_helper::streamSetParameterMaps(parametersDouble, parametersInt, parametersString, parametersIsSet,cout);
	
	cout <<endl <<"Check consistency of parameters" <<endl;
	if(!(plotPotential_inBrokenPhase_helper::checkConsistencyOfParameters(parametersDouble, parametersInt, parametersString, parametersIsSet)))
	{
		cerr <<"Failed! Some parameters are inconsistent!" <<endl;
		exit(EXIT_FAILURE);
	}
	cout <<"passed" <<endl;
	
	
	int Ls(parametersInt["L_s"]), Lt(parametersInt["L_t"]);
	constrainedEffectivePotential_inBrokenPhase CEP(Ls,Ls,Ls,Lt,parametersInt["antiperiodic_L_t"]);
	
	
	
	//set N_f, rho, r if set an different
	if(parametersIsSet["N_f"] && parametersInt["N_f"]!=CEP.get_N_f()){ CEP.set_N_f(parametersInt["N_f"]); }
	if(parametersIsSet["rho"] && parametersDouble["rho"]!=CEP.get_rho()){ CEP.set_rho(parametersDouble["rho"]); }
	if(parametersIsSet["r"] && parametersDouble["r"]!=CEP.get_r()){ CEP.set_rho(parametersDouble["r"]); }
	
	
	
	//set tolerances and max iterations
	if(parametersIsSet["absolut_tolerance_for_minimization"])
	{
		CEP.set_absolute_Accuracy( parametersDouble["absolut_tolerance_for_minimization"] ); 
	}
	if(parametersIsSet["relative_tolerance_for_minimization"])
	{
		CEP.set_relative_Accuracy( parametersDouble["relative_tolerance_for_minimization"] ); 
	}
	if(parametersIsSet["max_numer_of_iterations_minimizer"])
	{
		CEP.set_max_numberOfIterations( parametersInt["max_numer_of_iterations_minimizer"] );
	}
	//minimizer
	if(parametersIsSet["minimization_algorithm"])
	{
		CEP.set_minimizationAlgorithm( parametersInt["minimization_algorithm"] );
	}
	
	//set parameters
	CEP.set_yukawas(parametersDouble["yukawa_t"], parametersDouble["yukawa_t"] * parametersDouble["yukawa_ratio"]);
	CEP.set_lambda(parametersDouble["lambda"]);
	CEP.set_lambda_6(parametersDouble["lambda_6"]);
	//load list of fermionic contributions, if given
	if(parametersIsSet["use_listOfFermContr"] && parametersInt["use_listOfFermContr"])
	{
		cout <<"Load list of fermionic contributions from" <<endl << parametersString["listOfFermContr"] <<endl;
		if( ! (CEP.load_fermionicContribution( parametersString["listOfFermContr"] ) ) )
		{
			cerr <<"Error loading list of fermionic contributions!" <<endl;
			exit(EXIT_FAILURE);
		}
	}
	if(parametersIsSet["use_kappa"] && parametersInt["use_kappa"])
	{
		if(parametersDouble["kappa"]==0.0)
		{
			cerr <<"Error, cannot deal with zero kappa!  ignore value" <<endl;
			exit(EXIT_FAILURE);
		}
		//m0^2=(1-8*Nf*lambda*kappa^2-8*kappa)/(kappa)
		CEP.set_m0Squared( ( 1.0 - 8.0*CEP.get_N_f()*CEP.get_lambda()*parametersDouble["kappa"]*parametersDouble["kappa"] - 8.0*parametersDouble["kappa"])/(parametersDouble["kappa"]) );
	}
	else
	{
		CEP.set_m0Squared( parametersDouble["m0Squared"] );
	}
	
	//determine higgsmass
	int counter=0;
	bool toContinue(true);
	bool skipValue(false);
	double old_HiggsMassSquared(0.0), new_HiggsMassSquared(0.0);
	while(toContinue)
	{
		++counter;
		if(counter==1){ CEP.set_HigsMassSquared(0.000001); }
		else{ CEP.set_HigsMassSquared( new_HiggsMassSquared ); }
		double min(0.0), lower(0.0), upper(0.0);
		CEP.determine_startingPoints(parametersDouble["testvalue_min"], parametersDouble["testvalue_max"], parametersDouble["testvalue_step"], min, lower, upper);
		if(lower==upper)
		{
			cerr <<"Error, no minimum found in testvalue interval." <<endl;
			
			if(lower==0.0 || lower== parametersDouble["testvalue_min"])
			{
				cerr <<"assume true minimum between first two values" <<endl;
				bool success(false);
				double newLower=parametersDouble["testvalue_min"];
				double newStep=parametersDouble["testvalue_step"];
				double newUpper(0.0);
				for(int trial=0; trial <10; ++trial)
				{
					newLower=parametersDouble["testvalue_min"];
					newStep/=10.0;
					newUpper=newLower+10.0*newStep;
					CEP.determine_startingPoints(newLower, newUpper, newStep, min, lower, upper);
					if(lower!=upper){ success=true; break; }
				}
				if(success){ cerr <<"success" <<endl; }
				else
				{
					skipValue=true;
				}
			}
			else
			{
				skipValue=true;
			}
		}
		if(lower==0.0){lower=parametersDouble["testvalue_min"]; }
		if(skipValue)
		{
			cerr <<"Problems occured, during initializing of minimization process" <<endl;
			exit(EXIT_FAILURE);
		}
// 							cout <<"start with interval " <<lower <<"<" <<min <<"<" <<upper <<endl;
		CEP.initialize_minimizer(min, lower, upper);
		int n_of_iter=CEP.iterate_minimizer_until_convergence();
		old_HiggsMassSquared=CEP.get_actual_HiggsMassSquared();
		if(n_of_iter<=0)
		{
			cerr <<"Error, minimizer did not converge properly, skip value" <<endl;
			exit(EXIT_FAILURE);
			
		}
		else
		{
			new_HiggsMassSquared=CEP.compute_CEP_inBrokenPhase_secondDerivative(CEP.get_actual_minimum());
			cout <<"minimum: " <<CEP.get_actual_minimum() <<" (" <<n_of_iter <<" iter.)";
			cout <<"  newMassSquared: " <<new_HiggsMassSquared <<endl;
			if( std::abs( 2.0*(new_HiggsMassSquared - old_HiggsMassSquared)/(old_HiggsMassSquared + new_HiggsMassSquared)) <parametersDouble["tolerance_for_HiggsMassSquared"] )
			{
				cout <<"Mass determination converged after " <<counter <<" iterations" ;
				cout <<"   minimum: " <<CEP.get_actual_minimum() <<"  mHSquared: " <<new_HiggsMassSquared <<endl;
				toContinue=false;
			}
			if(counter==parametersInt["max_numer_of_iterations_HiggsMassSquared"])
			{
				cout <<"Mass determination did not converge after " <<parametersInt["max_numer_of_iterations_HiggsMassSquared"];
				cout <<", will use last values" <<endl;
				toContinue=false;
// 				skipValue=true;
			}
		}
	}
	if(skipValue)
	{
		cerr <<"Problems occured during mass determination" <<endl;
		exit(EXIT_FAILURE);
	}
	
	std::set< double > field;
	plotPotential_inBrokenPhase_helper::fillSetWithRange( parametersDouble["field_min"], parametersDouble["field_max"], parametersDouble["field_step"], field);
	
	//stores <field, potential>
	std::map< double, double > result;
	for( std::set< double >::const_iterator iter=field.begin(); iter!=field.end(); ++iter )
	{
		double potential=CEP.compute_CEP_inBrokenPhase(*iter);
		result.insert( std::make_pair( *iter, potential ) );
		cout <<"CEP( " <<*iter <<" ) = " <<potential <<endl;
	}
	
	if(parametersIsSet["outputfile"])
	{
		std::string outputFileName( plotPotential_inBrokenPhase_helper::generate_outputFileName( parametersString["outputfile"], parametersDouble, parametersInt, parametersIsSet) );
		
		cout <<"print output to: " <<outputFileName <<endl;
		
		std::ofstream outputFile( outputFileName.c_str() );
		if(!outputFile.good())
		{
			cerr <<"Error opening output file" <<endl <<outputFileName <<endl;
			cerr <<"send ouput to cout!" <<endl;
		}
		else
		{
			outputFile.precision(14);
			outputFile <<"# Output of plotPotential_inBrokenPhase" <<endl;
			outputFile <<"# parameters set:" <<endl;
			plotPotential_inBrokenPhase_helper::streamSetParameterMaps( parametersDouble, parametersInt, parametersString, parametersIsSet, outputFile, "#");
			outputFile <<"# Location of minimum: " <<CEP.get_actual_minimum() <<endl;
			outputFile <<"# curvature at minimum: " <<CEP.get_actual_HiggsMassSquared() <<endl;
			outputFile <<"# Output format is:" <<endl;
			outputFile <<"# fieldvalue      potential" <<endl;
			
			for(std::map< double, double >::const_iterator iter=result.begin(); iter!=result.end(); ++iter)
			{
				outputFile <<iter->first <<"   " <<iter->second <<endl;
			}
			outputFile.close();
		}
	}
	
	
	if(parametersDouble.size()+parametersInt.size()+parametersString.size() != numberOfParameters || parametersIsSet.size()!=numberOfParameters )
	{
		cerr <<"Error, number of parameters changed!" <<endl;
		cerr <<"parametersDouble.size()=" <<parametersDouble.size() <<endl;
		cerr <<"parametersInt.size()=" <<parametersInt.size() <<endl;
		cerr <<"parametersString.size()=" <<parametersString.size() <<endl;
		cerr <<"parametersIsSet.size()=" <<parametersIsSet.size() <<endl;
		exit(EXIT_FAILURE);
	}
}