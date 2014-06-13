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

#include "CEP_withFullBosDet.h"
#include "CEPscan_helper.h"

using std::cout;
using std::cerr;
using std::endl;


struct result_type
{
	double field; 
	double U_tree;
	double U_ferm;
	double U_det;
	double U_1st; 
// 		result_type( const &result_type o): field(o.field), U_tree(o.U_tree), U_ferm(o.U_ferm), U_det(o.U_det), U_1st(o.U_1st)
// 		{}
};


int main(int narg,char **arg)
{
	std::map< std::string, double > parametersDouble;
	std::map< std::string, int > parametersInt;
	std::map< std::string, std::string > parametersString;
	std::map< std::string, bool > parametersIsSet;
	
	CEPscan_helper::prepareParameterMaps_plotPotential_withFullBosDet(parametersDouble, parametersInt, parametersString, parametersIsSet);
	size_t numberOfParameters=parametersDouble.size()+parametersInt.size()+parametersString.size();
	
	if(narg!=2)
	{
		cerr <<"Error, start program with:" <<endl;
		cerr <<arg[0] <<"   inputFile" <<endl;
		exit(EXIT_FAILURE);
	}
	
	if( !CEPscan_helper::loadParameterMapsFromFile(parametersDouble, parametersInt, parametersString, parametersIsSet, arg[1]) )
	{
		cerr <<"Error loading input file" <<endl <<arg[1] <<endl;
		exit(EXIT_FAILURE);
	}
	
	cout <<endl <<"Parameters loaded:" <<endl;
	CEPscan_helper::streamSetParameterMaps(parametersDouble, parametersInt, parametersString, parametersIsSet,cout);
	
	cout <<endl <<"Check consistency of parameters" <<endl;
	if(!(CEPscan_helper::checkConsistencyOfParameters_plotPotential_withFullBosDet(parametersDouble, parametersInt, parametersString, parametersIsSet)))
	{
		cerr <<"Failed! Some parameters are inconsistent!" <<endl;
		exit(EXIT_FAILURE);
	}
	cout <<"passed" <<endl;
	
	
	int Ls(parametersInt["L_s"]), Lt(parametersInt["L_t"]);
	CEP_withFullBosDet CEP(Ls,Ls,Ls,Lt,parametersInt["antiperiodic_L_t"]);
	
	
	
	//set N_f, rho, r if set an different
	if(parametersIsSet["N_f"] && parametersInt["N_f"]!=CEP.get_N_f()){ CEP.set_N_f(parametersInt["N_f"]); }
	if(parametersIsSet["rho"] && parametersDouble["rho"]!=CEP.get_rho()){ CEP.set_rho(parametersDouble["rho"]); }
	if(parametersIsSet["r"] && parametersDouble["r"]!=CEP.get_r()){ CEP.set_rho(parametersDouble["r"]); }
	
	
	
	//set parameters
	CEP.set_yukawas(parametersDouble["yukawa_t"], parametersDouble["yukawa_t"] * parametersDouble["yukawa_ratio"]);
	CEP.set_lambda(parametersDouble["lambda"]);
	CEP.set_lambda_6(parametersDouble["lambda_6"]);
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
	
	
	
	
	std::set< double > field;
	CEPscan_helper::fillSetWithRange( parametersDouble["field_min"], parametersDouble["field_max"], parametersDouble["field_step"], field);
	
	//stores <field, potential>
	
	

	
	std::vector< result_type > result;
	
	result_type dummy;
	
	for( std::set< double >::const_iterator iter=field.begin(); iter!=field.end(); ++iter )
	{
		dummy.field=*iter;
		dummy.U_tree=CEP.compute_treeLevel(*iter);
		dummy.U_ferm=CEP.compute_fermionicContribution(*iter);
		dummy.U_det=CEP.compute_BosDetContribution(*iter);
		dummy.U_1st=CEP.compute_firstOrderInLambdas(*iter);
		result.push_back( dummy );
		cout <<"CEP( " <<*iter <<" ) = " <<dummy.U_tree + dummy.U_ferm + dummy.U_det + dummy.U_1st <<endl;
	}

	if(parametersIsSet["outputfile"])
	{
		std::string outputFileName( CEPscan_helper::generate_outputFileName_plotPotential_inBrokenPhase( parametersString["outputfile"], parametersDouble, parametersInt, parametersIsSet) );
		
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
			outputFile <<"# Output of plotPotential_withFullBosDet" <<endl;
			outputFile <<"# parameters set:" <<endl;
			CEPscan_helper::streamSetParameterMaps( parametersDouble, parametersInt, parametersString, parametersIsSet, outputFile, "#");
			outputFile <<"# Output format is:" <<endl;
			outputFile <<"# fieldvalue      U_tree   U_ferm   U_det   U_1st   U_total" <<endl;
			
			for(std::vector< result_type >::const_iterator iter=result.begin(); iter!=result.end(); ++iter)
			{
				outputFile <<iter->field <<"   " <<iter->U_tree <<"   " <<iter->U_ferm <<"   " <<iter->U_det <<"   " <<iter->U_1st <<"   " 
				           <<iter->U_tree + iter->U_ferm + iter->U_det + iter->U_1st <<endl;
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