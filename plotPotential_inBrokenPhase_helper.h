#ifndef PLOTPOTENTIAL_INBROKENPHASE_HELPER_H
#define PLOTPOTENTIAL_INBROKENPHASE_HELPER_H

#include <cstdlib>
#include <istream>
#include <iostream>
#include <fstream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace plotPotential_inBrokenPhase_helper
{
	void prepareParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet );
	
	bool loadParameterMapsFromFile( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, const std::string &fileName );
	
	void streamSetParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, std::ostream &output, const std::string &prefix = "");
	
	bool checkConsistencyOfParameters( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet );
	
	void fillSetWithRange( const double min, const double max, const double step, std::set< double > &toFill );
	
	std::string generate_outputFileName(const std::string &baseName, std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, bool > &paraIsSet );
}

#endif