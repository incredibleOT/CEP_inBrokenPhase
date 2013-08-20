#include "CEPscan_inBrokenPhase_helper.h"


void CEPscan_inBrokenPhase_helper::prepareParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet )
{
	paraD.clear(); paraI.clear(); paraS.clear(); paraIsSet.clear();
	
	paraI["L_s"]=-1;
	paraI["L_t"]=-1;
	paraI["antiperiodic_L_t"]=-1;
	
	paraI["scan_m0Squared"]   = 0;
	paraD["m0Squared"]        = 0.0;
	paraD["m0Squared_min"]    = 0.0;
	paraD["m0Squared_max"]    = 0.0;
	paraD["m0Squared_step"]   = 0.0;

	paraI["use_kappa"]        = 0;
	paraI["scan_kappa"]       = 0.0;
	paraD["kappa"]            = 0.0;
	paraD["kappa_min"]        = 0.0;
	paraD["kappa_max"]        = 0.0;
	paraD["kappa_step"]       = 0.0;

	paraI["scan_lambda"]      = 0;  
	paraD["lambda"]           = 0.0;
	paraD["lambda_min"]       = 0.0;
	paraD["lambda_max"]       = 0.0;
	paraD["lambda_step"]      = 0.0;

	paraI["scan_lambda_6"]    = 0;
	paraD["lambda_6"]         = 0.0;
	paraD["lambda_6_min"]     = 0.0;
	paraD["lambda_6_max"]     = 0.0;
	paraD["lambda_6_step"]    = 0.0;

	paraI["scan_yukawa_t"]    = 0;  
	paraD["yukawa_t"]         = 0.0;
	paraD["yukawa_t_min"]     = 0.0;
	paraD["yukawa_t_max"]     = 0.0;
	paraD["yukawa_t_step"]    = 0.0;

	paraI["scan_yukawa_ratio"]= 0;  
	paraD["yukawa_ratio"]     = 0.0;
	paraD["yukawa_ratio_min"] = 0.0;        
	paraD["yukawa_ratio_max"] = 0.0;
	paraD["yukawa_ratio_step"]= 0.0;


	// default: N_f=1, rho=1, r=0.5
	paraI["N_f"]              = 1;
	paraD["rho"]              = 1.0;
	paraD["r"]                = 0.5;

	// default: absolut_tolerance_for_minimization=1.0e-1, relative_tolerance_for_minimization=1.0e-7, tolerance_for_HiggsMassSquared=1.0e-5
	paraD["absolut_tolerance_for_minimization"]  =1.0e-7;
	paraD["relative_tolerance_for_minimization"] =1.0e-7;
	paraD["tolerance_for_HiggsMassSquared"]      =1.0e-5;

	//default 100 each
	paraI["max_numer_of_iterations_minimizer"]        = 100;
	paraI["max_numer_of_iterations_HiggsMassSquared"] = 100;
		
	// 1 gsl_min_fminimizer_goldensection, 2 gsl_min_fminimizer_brent, 3 gsl_min_fminimizer_quad_golden
	paraI["minimization_algorithm"] = 0;

	paraD["testvalue_min"]              = 0.0; 
	paraD["testvalue_max"]              = 0.0; 
	paraD["testvalue_step"]             = 0.0;

	paraS["outputfile"]             ="";
	
	for( std::map< std::string, double >::const_iterator iter=paraD.begin(); iter!=paraD.end(); ++iter )
	{
		paraIsSet[iter->first]=false;
	}
	for( std::map< std::string, int >::const_iterator iter=paraI.begin(); iter!=paraI.end(); ++iter )
	{
		paraIsSet[iter->first]=false;
	}
	for( std::map< std::string, std::string >::const_iterator iter=paraS.begin(); iter!=paraS.end(); ++iter )
	{
		paraIsSet[iter->first]=false;
	}
}//prepareParameterMaps



bool CEPscan_inBrokenPhase_helper::loadParameterMapsFromFile( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, const std::string &fileName )
{
	std::ifstream inputFile(fileName.c_str());
	std::string line,word;
	if(!(inputFile.good()))
	{
		std::cerr <<"Error opening inputfile " <<fileName <<std::endl;
		return false;
	}
	while(inputFile)
	{
		getline(inputFile, line);
		if(line.size()==0 || line[0] =='#' || line.find_first_not_of(' ') == std::string::npos){ continue; }
		std::istringstream strm(line);
		if(!(strm >> word)){ continue; }
		if( paraI.find(word) != paraI.end() )
		{
			if(paraIsSet[word]){ std::cerr<<"Error, parameter " <<word <<" already set" <<std::endl; return false; } 
			if(!(strm >> paraI[word]))
			{
				std::cerr <<"Error loading value of " <<word <<std::endl; 
				std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; 
			}
			else{ paraIsSet[word]=true;}
		}
		else if( paraD.find(word) != paraD.end() )
		{ 
			if(paraIsSet[word]){ std::cerr<<"Error, parameter " <<word <<" already set" <<std::endl; return false; }
			if(!(strm >> paraD[word]))
			{
				std::cerr <<"Error loading value of " <<word <<std::endl; 
				std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; 
			}
			else{ paraIsSet[word]=true;}
		}
		else if( paraS.find(word) != paraS.end() )
		{ 
			if(paraIsSet[word]){ std::cerr<<"Error, parameter " <<word <<" already set" <<std::endl; return false; }
			if(!(strm >> paraS[word]))
			{
				std::cerr <<"Error loading value of " <<word <<std::endl; 
				std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; 
			}
			else{ paraIsSet[word]=true;}
		}
		else{ std::cerr <<"Warning, there is no parameter called \"" <<word <<"\""<<std::endl; return false;}
		if(inputFile.eof()){ break; }
	}
	inputFile.close(); inputFile.clear();
	return true;
}//loadParameterMapsFromFile




void CEPscan_inBrokenPhase_helper::streamSetParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, std::ostream &output, const std::string &prefix)
{
	for( std::map< std::string, double >::const_iterator iter=paraD.begin(); iter!=paraD.end(); ++iter )
	{
		if(paraIsSet[iter->first])
		{
			if(!(prefix=="")){output <<prefix <<" ";}  
			output <<iter->first <<"     " <<iter->second <<std::endl;
		}
	}
	for( std::map< std::string, int >::const_iterator iter=paraI.begin(); iter!=paraI.end(); ++iter )
	{
		if(paraIsSet[iter->first])
		{
			if(!(prefix=="")){output <<prefix <<" ";}  
			output <<iter->first <<"     " <<iter->second <<std::endl;
		}
	}
	for( std::map< std::string, std::string >::const_iterator iter=paraS.begin(); iter!=paraS.end(); ++iter )
	{
		if(paraIsSet[iter->first])
		{
			if(!(prefix=="")){output <<prefix <<" ";}  
			output <<iter->first <<"     " <<iter->second <<std::endl;
		}
	}
}//streamSetParameterMaps




bool CEPscan_inBrokenPhase_helper::checkConsistencyOfParameters( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet )
{
	//extend related
	if(!(paraIsSet["L_s"] && paraIsSet["L_t"]) || paraI["L_s"]<=0 || paraI["L_t"]<=0)
	{
		std::cerr <<"Error, no or non-positive lattice extends given" <<std::endl;
		return false;
	}
	if(paraI["L_s"]%2!=0 || paraI["L_t"]%2!=0)
	{
		std::cerr <<"Error, only even lattice extend allowed" <<std::endl;
		return false;
	}
	if(!paraIsSet["antiperiodic_L_t"])
	{
		std::cerr <<"Error, antiperiodic_L_t is not specified" <<std::endl;
		return false;
	}
	//m0Squared (only if not use_kappa is set)
	if(!paraI["use_kappa"])
	if( (!paraIsSet["scan_m0Squared"] || paraI["scan_m0Squared"]==0) && !paraIsSet["m0Squared"])
	{
		std::cerr <<"Error, no m0Squared given" <<std::endl;
		return false;
	}
	if( paraIsSet["scan_m0Squared"] && paraI["scan_m0Squared"] )
	{
		if( !( paraIsSet["m0Squared_min"] && paraIsSet["m0Squared_max"] && paraIsSet["m0Squared_step"] ) )
		{
			std::cerr <<"Error, no scan range in m0Squared given" <<std::endl;
			return false;
		}
		if( paraD["m0Squared_max"] < paraD["m0Squared_min"] || paraD["m0Squared_step"] <=0.0 )
		{
			std::cerr <<"Error, inconsistent scan range in m0Squared" <<std::endl;
			return false;
		}
	}
	//
	if(paraI["use_kappa"])
	if( (!paraIsSet["scan_kappa"] || paraI["scan_kappa"]==0) && !paraIsSet["kappa"])
	{
		std::cerr <<"Error, no kappa given" <<std::endl;
		return false;
	}
	if( paraIsSet["scan_kappa"] && paraI["scan_kappa"] )
	{
		if( !( paraIsSet["kappa_min"] && paraIsSet["kappa_max"] && paraIsSet["kappa_step"] ) )
		{
			std::cerr <<"Error, no scan range in kappa given" <<std::endl;
			return false;
		}
		if( paraD["kappa_max"] < paraD["kappa_min"] || paraD["kappa_step"] <=0.0 )
		{
			std::cerr <<"Error, inconsistent scan range in kappa" <<std::endl;
			return false;
		}
	}
	//lambda
	if( (!paraIsSet["scan_lambda"] || paraI["scan_lambda"]==0) && !paraIsSet["lambda"])
	{
		std::cerr <<"Error, no lambda given" <<std::endl;
		return false;
	}
	if( paraIsSet["scan_lambda"] && paraI["scan_lambda"] )
	{
		if( !( paraIsSet["lambda_min"] && paraIsSet["lambda_max"] && paraIsSet["lambda_step"] ) )
		{
			std::cerr <<"Error, no scan range in lambda given" <<std::endl;
			return false;
		}
		if( paraD["lambda_max"] < paraD["lambda_min"] || paraD["lambda_step"] <=0.0 )
		{
			std::cerr <<"Error, inconsistent scan range in lambda" <<std::endl;
			return false;
		}
	}
	//lambda_6
	if( (!paraIsSet["scan_lambda_6"] || paraI["scan_lambda_6"]==0) && !paraIsSet["lambda_6"])
	{
		std::cerr <<"Error, no lambda_6 given" <<std::endl;
		return false;
	}
	if( paraIsSet["scan_lambda_6"] && paraI["scan_lambda_6"] )
	{
		if( !( paraIsSet["lambda_6_min"] && paraIsSet["lambda_6_max"] && paraIsSet["lambda_6_step"] ) )
		{
			std::cerr <<"Error, no scan range in lambda_6 given" <<std::endl;
			return false;
		}
		if( paraD["lambda_6_max"] < paraD["lambda_6_min"] || paraD["lambda_6_step"] <=0.0 )
		{
			std::cerr <<"Error, inconsistent scan range in lambda_6" <<std::endl;
			return false;
		}
	}
	//yukawa_t
	if( (!paraIsSet["scan_yukawa_t"] || paraI["scan_yukawa_t"]==0) && !paraIsSet["yukawa_t"])
	{
		std::cerr <<"Error, no yukawa_t given" <<std::endl;
		return false;
	}
	if( paraIsSet["scan_yukawa_t"] && paraI["scan_yukawa_t"] )
	{
		if( !( paraIsSet["yukawa_t_min"] && paraIsSet["yukawa_t_max"] && paraIsSet["yukawa_t_step"] ) )
		{
			std::cerr <<"Error, no scan range in yukawa_t given" <<std::endl;
			return false;
		}
		if( paraD["yukawa_t_max"] < paraD["yukawa_t_min"] || paraD["yukawa_t_step"] <=0.0 )
		{
			std::cerr <<"Error, inconsistent scan range in yukawa_t" <<std::endl;
			return false;
		}
	}
	//yukawa_ratio
	if( (!paraIsSet["scan_yukawa_ratio"] || paraI["scan_yukawa_ratio"]==0) && !paraIsSet["yukawa_ratio"])
	{
		std::cerr <<"Error, no yukawa_ratio given" <<std::endl;
		return false;
	}
	if( paraIsSet["scan_yukawa_ratio"] && paraI["scan_yukawa_ratio"] )
	{
		if( !( paraIsSet["yukawa_ratio_min"] && paraIsSet["yukawa_ratio_max"] && paraIsSet["yukawa_ratio_step"] ) )
		{
			std::cerr <<"Error, no scan range in yukawa_ratio given" <<std::endl;
			return false;
		}
		if( paraD["yukawa_ratio_max"] < paraD["yukawa_ratio_min"] || paraD["yukawa_ratio_step"] <=0.0 )
		{
			std::cerr <<"Error, inconsistent scan range in yukawa_ratio" <<std::endl;
			return false;
		}
	}
	//testvalue
	if(!(paraIsSet["testvalue_min"] && paraIsSet["testvalue_max"] && paraIsSet["testvalue_step"]))
	{
		std::cerr <<"Error, no scan range for testvalue given" <<std::endl;
		return false;
	}
	if( paraD["testvalue_max"] < paraD["testvalue_min"] || paraD["testvalue_step"] <=0.0 )
	{
		std::cerr <<"Error, inconsistent scan range in testvalue" <<std::endl;
		return false;
	}
	
	return true;
}//checkConsistencyOfParameters




void CEPscan_inBrokenPhase_helper::fillSetWithRange( const double min, const double max, const double step, std::set< double > &toFill)
{
	toFill.clear();
	
	if( (max-min)/step < 0 )
	{
		std::cerr <<"Error, unreasonable scanrange: min: " <<min <<", max: " <<max <<", step: " <<step <<std::endl;
		exit(EXIT_FAILURE);
	}
	
	int numberOfEntries=static_cast< int >( (max-min)/step + 1.5);
	
	for( int i=0; i<numberOfEntries; ++i)
	{
		toFill.insert( min + static_cast< double >(i) *step);
	}
}



bool CEPscan_inBrokenPhase_helper::printResultsVectorToStream(const std::vector< resultForOutput > &results, std::ostream &output)
{
	//increase precision, but store the old one
	std::streamsize oldPrec=output.precision();
	output.precision(12);
	for(std::vector< resultForOutput >::const_iterator iter=results.begin(); iter!=results.end(); ++iter)
	{
		if( !( output <<iter->m0Squared <<" " <<iter->lambda <<" " <<iter->lambda_6 <<" " <<iter->yukawa_t <<" " <<iter->yukawa_b <<" " <<iter->minimum <<" " <<iter->mHSquared <<" " <<iter->potential <<" "   <<std::endl ) )
		{
			std::cerr <<"Error during output of results" <<std::endl;
			return false;
		}
	}
	output.precision(oldPrec);
	return true;
}


std::string CEPscan_inBrokenPhase_helper::generate_outputFileName(const std::string &baseName, std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, bool > &paraIsSet )
{
	//replaces placeholders in the filename
	// [Ls]->Lxx, [Lt]->Txx, [m0Sq_k]->m0Sq_xxx_xxx or k_xxx_xxx, [l]->l_xxx_xxx ->[l6]->l6_xxx_xxx
	// [yt]->yt_xxx_xxx
	std::string outputFileName(baseName);
	if( outputFileName.find("[Ls]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"L" <<paraI["L_s"];
		outputFileName.replace(outputFileName.find("[Ls]"),4, ss.str() );
	}
	if( outputFileName.find("[Lt]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"T" <<paraI["L_t"];
		outputFileName.replace(outputFileName.find("[Lt]"),4, ss.str() );
	}
	if( outputFileName.find("[m0Sq_k]")!=std::string::npos )
	{
		if(paraIsSet["use_kappa"] && paraI["use_kappa"])
		{
			std::ostringstream ss;
			ss <<"k_";
			ss.precision(7);
			if(paraI["scan_kappa"]){ss <<paraD["kappa_min"] <<"_"<<paraD["kappa_max"]; }
			else {ss <<paraD["kappa"]; }
			outputFileName.replace(outputFileName.find("[m0Sq_k]"),8, ss.str() );
		}
		else
		{
			std::ostringstream ss;
			ss <<"m0Sq_";
			ss.precision(7);
			if(paraI["scan_kappa"]){ss <<paraD["m0Squared_min"] <<"_"<<paraD["m0Squared_max"]; }
			else {ss <<paraD["m0Squared"]; }
			outputFileName.replace(outputFileName.find("[m0Sq_k]"),8, ss.str() );
		}
	}
	if( outputFileName.find("[l]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"l_";
		ss.precision(7);
		if(paraI["scan_lambda"]){ss <<paraD["lambda_min"] <<"_"<<paraD["lambda_max"]; }
		else {ss <<paraD["lambda"]; }
		outputFileName.replace(outputFileName.find("[l]"),3, ss.str() );
	}
	if( outputFileName.find("[l6]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"l6_";
		ss.precision(7);
		if(paraI["scan_lambda_6"]){ss <<paraD["lambda_6_min"] <<"_"<<paraD["lambda_6_max"]; }
		else {ss <<paraD["lambda_6"]; }
		outputFileName.replace(outputFileName.find("[l6]"),4, ss.str() );
	}
	if( outputFileName.find("[yt]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"yt_";
		ss.precision(7);
		if(paraI["scan_yukawa_t"]){ss <<paraD["yukawa_t_min"] <<"_"<<paraD["yukawa_t_max"]; }
		else {ss <<paraD["yukawa_t"]; }
		outputFileName.replace(outputFileName.find("[yt]"),4, ss.str() );
	}
	if( outputFileName.find("[yr]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"yr_";
		ss.precision(7);
		if(paraI["scan_yukawa_ratio"]){ss <<paraD["yukawa_ratio_min"] <<"_"<<paraD["yukawa_ratio_max"]; }
		else {ss <<paraD["yukawa_ratio"]; }
		outputFileName.replace(outputFileName.find("[yr]"),4, ss.str() );
	}
	return outputFileName;
}
	
	
	
	
	


