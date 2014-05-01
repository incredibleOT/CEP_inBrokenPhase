#ifndef CEPSCAN_HELPER_H
#define CEPSCAN_HELPER_H

#include <cstdlib>
#include <cmath>
#include <istream>
#include <iostream>
#include <fstream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace CEPscan_helper
{
	
	//prepares the parameterMaps for CEPscan_inBrokenPhase
	void prepareParameterMaps_inBrokenPhase( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet );
	
	//prepares the parameterMaps for plotPotential_inBrokenPhase
	void prepareParameterMaps_plotPotential_inBrokenPhase( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet );
	
	
	bool loadParameterMapsFromFile( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, const std::string &fileName );
	
	void streamSetParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, std::ostream &output, const std::string &prefix = "");
	
	void streamParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::ostream &output, const std::string &prefix = "");
	
	bool checkConsistencyOfParameters_inBrokenPhase( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet );
	
	bool checkConsistencyOfParameters_plotPotential_inBrokenPhase( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet );
	
	
	void fillSetWithRange( const double min, const double max, const double step, std::set< double > &toFill );
	
	struct resultForOutput
	{
		double m0Squared;
		double lambda;
		double lambda_6;
		double yukawa_t;
		double yukawa_b;
		double minimum;
		double mHSquared;
		double potential;
		
	};
	
	bool printResultsVectorToStream(const std::vector< resultForOutput > &results, std::ostream &output);
	
	std::map< int, double >::iterator findClosestMass( std::map< int, double > &HiggsMassesSquared, double value );
	
	std::string generate_outputFileName_inBrokenPhase(const std::string &baseName, std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, bool > &paraIsSet );
	
	std::string generate_outputFileName_plotPotential_inBrokenPhase(const std::string &baseName, std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, bool > &paraIsSet );
	
}


#endif
