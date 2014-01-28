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
#include "CEPscan_inBrokenPhase_helper.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int narg,char **arg)
{
	std::map< std::string, double > parametersDouble;
	std::map< std::string, int > parametersInt;
	std::map< std::string, std::string > parametersString;
	std::map< std::string, bool > parametersIsSet;
	
	CEPscan_inBrokenPhase_helper::prepareParameterMaps(parametersDouble, parametersInt, parametersString, parametersIsSet);
	size_t numberOfParameters=parametersDouble.size()+parametersInt.size()+parametersString.size();
	
	if(narg!=2)
	{
		cerr <<"Error, start program with:" <<endl;
		cerr <<arg[0] <<"   inputFile" <<endl;
		exit(EXIT_FAILURE);
	}
	
	if( !CEPscan_inBrokenPhase_helper::loadParameterMapsFromFile(parametersDouble, parametersInt, parametersString, parametersIsSet, arg[1]) )
	{
		cerr <<"Error loading input file" <<endl <<arg[1] <<endl;
		exit(EXIT_FAILURE);
	}
	
	cout <<endl <<"Parameters loaded:" <<endl;
	CEPscan_inBrokenPhase_helper::streamSetParameterMaps(parametersDouble, parametersInt, parametersString, parametersIsSet,cout);
	
	cout <<endl <<"Check consistency of parameters" <<endl;
	if(!(CEPscan_inBrokenPhase_helper::checkConsistencyOfParameters(parametersDouble, parametersInt, parametersString, parametersIsSet)))
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
	
	//set whether goldstones should be excluded
	CEP.set_ignore_goldstone_modes(parametersInt["exclude_goldstones"]);
	
	//set tolerances and max iterations
	CEP.set_absolute_Accuracy( parametersDouble["absolut_tolerance_for_minimization"] ); 
	CEP.set_relative_Accuracy( parametersDouble["relative_tolerance_for_minimization"] ); 
	
	if(parametersIsSet["max_numer_of_iterations_minimizer"])
	{
		CEP.set_max_numberOfIterations( parametersInt["max_numer_of_iterations_minimizer"] );
	}
	//minimizer
	if(parametersIsSet["minimization_algorithm"])
	{
		CEP.set_minimizationAlgorithm( parametersInt["minimization_algorithm"] );
	}
	
	
	
	//preparing the ranges for scanning
	std::set< double > m0Squared_or_kappa_values, lambda_values, lambda_6_values, yukawa_t_values, yukawa_ratio_values;
	
	//m0Squared_or_kappa
	if(parametersIsSet["use_kappa"] && parametersInt["use_kappa"])
	{
		if(parametersIsSet["scan_kappa"] && parametersInt["scan_kappa"] )
		{
			CEPscan_inBrokenPhase_helper::fillSetWithRange( parametersDouble["kappa_min"], parametersDouble["kappa_max"], parametersDouble["kappa_step"], m0Squared_or_kappa_values );
		}
		else{ m0Squared_or_kappa_values.insert(parametersDouble["kappa"]); }
	}
	else
	{
		if(parametersIsSet["scan_m0Squared"] && parametersInt["scan_m0Squared"] )
		{
			CEPscan_inBrokenPhase_helper::fillSetWithRange( parametersDouble["m0Squared_min"], parametersDouble["m0Squared_max"], parametersDouble["m0Squared_step"], m0Squared_or_kappa_values );
		}
		else{ m0Squared_or_kappa_values.insert(parametersDouble["m0Squared"]); }
	}
	//lambda
	if(parametersIsSet["scan_lambda"] && parametersInt["scan_lambda"] )
	{
		CEPscan_inBrokenPhase_helper::fillSetWithRange( parametersDouble["lambda_min"], parametersDouble["lambda_max"], parametersDouble["lambda_step"], lambda_values );
	}
	else{ lambda_values.insert(parametersDouble["lambda"]); }
	//lambda_6
	if(parametersIsSet["scan_lambda_6"] && parametersInt["scan_lambda_6"] )
	{
		CEPscan_inBrokenPhase_helper::fillSetWithRange( parametersDouble["lambda_6_min"], parametersDouble["lambda_6_max"], parametersDouble["lambda_6_step"], lambda_6_values );
	}
	else{ lambda_6_values.insert(parametersDouble["lambda_6"]); }
	//yukawa_t
	if(parametersIsSet["scan_yukawa_t"] && parametersInt["scan_yukawa_t"] )
	{
		CEPscan_inBrokenPhase_helper::fillSetWithRange( parametersDouble["yukawa_t_min"], parametersDouble["yukawa_t_max"], parametersDouble["yukawa_t_step"], yukawa_t_values );
	}
	else{ yukawa_t_values.insert(parametersDouble["yukawa_t"]); }
	//yukawa_ratio
	if(parametersIsSet["scan_yukawa_ratio"] && parametersInt["scan_yukawa_ratio"] )
	{
		CEPscan_inBrokenPhase_helper::fillSetWithRange( parametersDouble["yukawa_ratio_min"], parametersDouble["yukawa_ratio_max"], parametersDouble["yukawa_ratio_step"], yukawa_ratio_values );
	}
	else{ yukawa_ratio_values.insert(parametersDouble["yukawa_ratio"]); }

	
	
	
	
	cout <<"start scanning" <<endl;
	cout.precision(12);
	std::vector< CEPscan_inBrokenPhase_helper::resultForOutput > results;
	double lastHiggsMassSquared(0.0000001);
	//now iterate
	for(std::set< double >::const_iterator yukawa_t=yukawa_t_values.begin(); yukawa_t!=yukawa_t_values.end(); ++yukawa_t)
	{
		for(std::set< double >::const_iterator yukawa_ratio=yukawa_ratio_values.begin(); yukawa_ratio!=yukawa_ratio_values.end(); ++yukawa_ratio)
		{
			CEP.set_yukawas(*yukawa_t, *yukawa_t * (*yukawa_ratio));
			cout <<"y_t=" <<CEP.get_yukawa_t() <<"  y_b=" <<CEP.get_yukawa_b() <<endl;
			if(parametersIsSet["use_listOfFermContr"] && parametersInt["use_listOfFermContr"])
			{
				cout <<"Load list of fermionic contributions from" <<endl << parametersString["listOfFermContr"] <<endl;
				if( ! (CEP.load_fermionicContribution( parametersString["listOfFermContr"] ) ) )
				{
					cerr <<"Error loading list of fermionic contributions!" <<endl;
					exit(EXIT_FAILURE);
				}
			}	
			for(std::set< double >::const_iterator lambda_6=lambda_6_values.begin(); lambda_6!=lambda_6_values.end(); ++lambda_6)
			{
				CEP.set_lambda_6(*lambda_6);
				cout <<"     lambda_6=" <<CEP.get_lambda_6() <<endl;
				for(std::set< double >::const_iterator lambda=lambda_values.begin(); lambda!=lambda_values.end(); ++lambda)
				{
					CEP.set_lambda(*lambda);
					cout <<"          lambda=" <<CEP.get_lambda() <<endl;
					for(std::set< double >::const_iterator m0Squared_or_kappa=m0Squared_or_kappa_values.begin(); m0Squared_or_kappa!=m0Squared_or_kappa_values.end(); ++m0Squared_or_kappa)
					{
						if(parametersIsSet["use_kappa"] && parametersInt["use_kappa"])
						{
							if(*m0Squared_or_kappa==0.0)
							{
								cerr <<"Error, cannot deal with zero kappa!  ignore value" <<endl;
								continue;
							}
							//m0^2=(1-8*Nf*lambda*kappa^2-8*kappa)/(kappa)
							CEP.set_m0Squared( (1.0 - 8.0*CEP.get_N_f()*CEP.get_lambda()*(*m0Squared_or_kappa)*(*m0Squared_or_kappa) - 8.0*(*m0Squared_or_kappa))/(*m0Squared_or_kappa) );
						}
						else
						{
							CEP.set_m0Squared( *m0Squared_or_kappa );
						}
						cout <<"               m0Squared=" <<CEP.get_m0Squared();
						if(parametersIsSet["use_kappa"] && parametersInt["use_kappa"])
						{
							cout <<" (kappa=" <<*m0Squared_or_kappa <<")";
						}
						cout <<endl;
						
						
						
						//do the magic
						int counter=0;
						bool toContinue(true);
						bool skipValue(false);
						double old_HiggsMassSquared(0.0), new_HiggsMassSquared(0.0);
						
						//the following is for the case of periodic solutions (i.e. arriving at a mass during the iteration where we were already)
						std::map< int, double > HiggsMassesSquared; //<counter, HiggsMasses>, to test for periodic solutions
						std::map< int, double > minima; //<counter, minimum>, to test for periodic solutions
						bool periodic_solution(false);
						const double accuracy_for_periodicity=1.0e-12;//when two masses are assumed to be equal
						const double maximal_spread_for_periodicity=5.0e-3; //if spread for periodicity is too large, its not an artifact but maybe a fluctuation between two minima
						double av_massSquared=0.0;
						double av_minimum=0.0;
// 						bool secondIteration(false); //set to true, if periodic solution is found
						while(toContinue)
						{
// 							if(secondIteration)
// 							{
// 								CEP.set_absolute_Accuracy( parametersDouble["absolut_tolerance_for_minimization"]*0.1 ); 
// 								CEP.set_relative_Accuracy( parametersDouble["relative_tolerance_for_minimization"]*0.1);
// 							}
// 							else
// 							{
// 								CEP.set_absolute_Accuracy( parametersDouble["absolut_tolerance_for_minimization"] ); 
// 								CEP.set_relative_Accuracy( parametersDouble["relative_tolerance_for_minimization"] );
// 							}
							++counter;
							if(counter==1){ CEP.set_HigsMassSquared(lastHiggsMassSquared); }
							else{ CEP.set_HigsMassSquared( new_HiggsMassSquared ); }
							HiggsMassesSquared.insert( std::make_pair( counter, CEP.get_actual_HiggsMassSquared() ) );
							double min(0.0), lower(0.0), upper(0.0);
							CEP.determine_startingPoints(parametersDouble["testvalue_min"], parametersDouble["testvalue_max"], parametersDouble["testvalue_step"], min, lower, upper);
							if(lower==upper)
							{
// 								cerr <<"Error, no minimum found in testvalue interval." <<endl;
								
								if(lower==0.0 || lower== parametersDouble["testvalue_min"])
								{
// 									cerr <<"assume true minimum between first two values" <<endl;
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
// 									if(success){ cerr <<"success" <<endl; }
									if(!success)
									{
										cerr <<"Error, no minimum found in testvalue interval." <<endl;
										skipValue=true;
										break;
									}
								}
								else
								{
									cerr <<"Error, no minimum found in testvalue interval." <<endl;
									skipValue=true;
									break;
								}
// 								skipValue=true;
// 								break;
							}
							if(lower==0.0){lower=parametersDouble["testvalue_min"]; }
							
							
// 							cout <<"start with interval " <<lower <<"<" <<min <<"<" <<upper <<endl;
							CEP.initialize_minimizer(min, lower, upper);
							int n_of_iter=CEP.iterate_minimizer_until_convergence();
							old_HiggsMassSquared=CEP.get_actual_HiggsMassSquared();
							if(n_of_iter<=0)
							{
								cerr <<"Error, minimizer did not converge properly, skip value" <<endl;
								skipValue=true;
								break;
							}
							else
							{
								new_HiggsMassSquared=CEP.compute_CEP_inBrokenPhase_secondDerivative(CEP.get_actual_minimum());
// 								cout <<"minimum: " <<CEP.get_actual_minimum() <<" (" <<n_of_iter <<" iter.)";
// 								cout <<"  newMassSquared: " <<new_HiggsMassSquared <<endl;
								cout <<counter <<":  min: " <<CEP.get_actual_minimum() <<"  mHSquared: " <<new_HiggsMassSquared <<endl;
								minima.insert( std::make_pair( counter, CEP.get_actual_minimum() ) );
								if( std::abs( 2.0*(new_HiggsMassSquared - old_HiggsMassSquared)/(old_HiggsMassSquared + new_HiggsMassSquared)) <parametersDouble["tolerance_for_HiggsMassSquared"] )
								{
									//mass converged
									cout <<"Mass determination converged after " <<counter <<" iterations" ;
									cout <<"   minimum: " <<CEP.get_actual_minimum() <<"  mHSquared: " <<CEP.get_actual_HiggsMassSquared() <<endl;
									toContinue=false;
									break; //redundant....
								}
								//check for periodic solution
								std::map< int, double >::iterator closest=CEPscan_inBrokenPhase_helper::findClosestMass( HiggsMassesSquared,  new_HiggsMassSquared );
								if( std::abs( 2.0*( closest->second - new_HiggsMassSquared)/(closest->second + new_HiggsMassSquared)) < accuracy_for_periodicity )
								{
									int period=counter+1-closest->first;
									cout <<"Periodic solution found! Period=" <<period <<endl;
									//for determin spread of mSquared and minimum
									double largest_mSquared=closest->second;
									double smallest_mSquared=closest->second;
									double largest_minimum= (minima.find(closest->first)!=minima.end()) ? ((minima.find(closest->first)->second)):-1.0;
									double smallest_minimum=largest_minimum;
									
									av_massSquared=0.0;
									av_minimum=0.0;
									while(closest != HiggsMassesSquared.end())
									{
										av_massSquared+=closest->second;
										std::map< int, double >::iterator minIter=minima.find(closest->first);
										if(minIter==minima.end())
										{
											cerr<<"Error, there is no corresponding minimum. Skip value" <<endl;
											skipValue=true;
											toContinue=false;
											break;
										}
										av_minimum+=minIter->second;
										if(closest->second < smallest_mSquared){ smallest_mSquared=closest->second; }
										if(closest->second > largest_mSquared){ largest_mSquared=closest->second; }
										if(minIter->second < smallest_minimum){ smallest_minimum=minIter->second; }
										if(minIter->second > largest_minimum){ largest_minimum=minIter->second; }
										++closest;
// 										cout <<"av=" <<av <<endl;
									}
									//if spread is too large, ignore values
									if(   (0.5*(largest_mSquared-smallest_mSquared)/(largest_mSquared+smallest_mSquared) > maximal_spread_for_periodicity)   ||    (0.5*(largest_minimum-smallest_minimum)/(largest_minimum+smallest_minimum) > maximal_spread_for_periodicity)    )
									{
										cout <<"Spread of values during period too large, skip value" <<endl;
										skipValue=true;
										toContinue=false;
									}
									else
									{
										cout <<"take average of values" <<endl;
										av_massSquared/=static_cast< double >(period);
										av_minimum/=static_cast< double >(period);
										toContinue=false;
										periodic_solution=true;
										cout <<"Mass determination converged after " <<counter <<" iterations" ;
										cout <<"   minimum: " <<av_minimum <<"  mHSquared: " <<av_massSquared <<endl;
									}
									
								}
								if(counter==parametersInt["max_numer_of_iterations_HiggsMassSquared"])
								{
									cout <<"Mass determination did not converge after " <<parametersInt["max_numer_of_iterations_HiggsMassSquared"];
									cout <<", skip value" <<endl;
									toContinue=false;
									skipValue=true;
								}
							}
						}
						if(skipValue){ continue; }
						CEPscan_inBrokenPhase_helper::resultForOutput dummy;
						dummy.m0Squared = CEP.get_m0Squared();
						dummy.lambda    = CEP.get_lambda();
						dummy.lambda_6  = CEP.get_lambda_6();
						dummy.yukawa_t  = CEP.get_yukawa_t();
						dummy.yukawa_b  = CEP.get_yukawa_b();
						dummy.minimum   = (periodic_solution)?(av_minimum):(CEP.get_actual_minimum());
						dummy.mHSquared = (periodic_solution)?(av_massSquared):(CEP.get_actual_HiggsMassSquared());
						dummy.potential = CEP.get_potentialAtMinimum();
						results.push_back(dummy);
						lastHiggsMassSquared=dummy.mHSquared;
					}//m0Squared_or_kappa
				}//lambda
			}//lambda_6
		}//yukawa_ratio
	}//yukawa_t


	
	//output to cout
	{
		cout <<"Results (m0Squared,   lambda,   lambda_6,   yukawa_t,   yukawa_b,   minimum,   mHSquared,   potential:" <<endl;
		CEPscan_inBrokenPhase_helper::printResultsVectorToStream( results, cout );
	}
	
	
	if(parametersIsSet["outputfile"])
	{
		std::string outputFileName( CEPscan_inBrokenPhase_helper::generate_outputFileName( parametersString["outputfile"], parametersDouble, parametersInt, parametersIsSet) );
		
		cout <<"print output to: " <<outputFileName <<endl;
		
		std::ofstream outputFile( outputFileName.c_str() );
		if(!outputFile.good())
		{
			cerr <<"Error opening output file" <<endl <<outputFileName <<endl;
			cerr <<"send ouput to cout!" <<endl;
		}
		else
		{
			outputFile <<"# Output of CEPscan_inBrokenPhase" <<endl;
			outputFile <<"# parameters set:" <<endl;
			CEPscan_inBrokenPhase_helper::streamSetParameterMaps( parametersDouble, parametersInt, parametersString, parametersIsSet, outputFile, "#");
			outputFile <<"# Output format is:" <<endl;
			outputFile <<"# m0Squared   lambda   lambda_6   yukawa_t   yukawa_b   minimum   mHSquared   potential" <<endl;
			CEPscan_inBrokenPhase_helper::printResultsVectorToStream(results, outputFile);
			outputFile.close();
		}
	}
	
	
	if(parametersDouble.size()+parametersInt.size()+parametersString.size() != numberOfParameters || parametersIsSet.size()!=numberOfParameters )
	{
		cerr <<"Error, number of parameters changed! It was " <<numberOfParameters <<"in the beginning." <<endl;
		cerr <<"parametersDouble.size()=" <<parametersDouble.size() <<endl;
		cerr <<"parametersInt.size()=" <<parametersInt.size() <<endl;
		cerr <<"parametersString.size()=" <<parametersString.size() <<endl;
		cerr <<"parametersIsSet.size()=" <<parametersIsSet.size() <<endl;
		CEPscan_inBrokenPhase_helper::streamParameterMaps(parametersDouble, parametersInt, parametersString, cout);
		exit(EXIT_FAILURE);
	}
	
}