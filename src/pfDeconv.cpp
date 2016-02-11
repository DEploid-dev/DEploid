#include <iostream>
#include <vector>
#include <string>
#include <sstream>      // std::stringstream
#include<stdlib.h>     /* strtol, strtod */
#include<fstream>
//#include <boost/math/special_functions/gamma.hpp>
#include<math.h>
#include "mcmc.hpp"
//#include<iostream>
//#include<vector>
//#include<cassert>       // assert
#include<stdexcept>      // std::invalid_argument
//#include<stdlib.h>     /* strtol, strtod */
//#include<climits> // INT_MAX


using namespace std;
//#include "param.hpp"


double calcLLK( double ref, double alt, double unadjustedWsaf, double err = 0.01, double fac=100 ) {
    double adjustedWsaf = unadjustedWsaf+err*(1-2*unadjustedWsaf);
    //double llk = lbeta(alt+adjustedWsaf*fac, ref+(1-adjustedWsaf)*fac)-lbeta(adjustedWsaf*fac,(1-adjustedWsaf)*fac);
    double llk = lgamma(fac*adjustedWsaf+alt)+lgamma(fac*(1-adjustedWsaf)+ref)-lgamma(fac*adjustedWsaf)-lgamma(fac*(1-adjustedWsaf));
    return llk;
}

void buildPanel(){

}

void readFileLines(const char inchar[], vector <double> & out_vec){
    ifstream in_file( inchar );
    string tmp_line;
    if ( in_file.good() ){
        getline ( in_file, tmp_line ); // skip the first line, which is the header
        getline ( in_file, tmp_line );
        while ( tmp_line.size() > 0 ){
            size_t field_start = 0;
            size_t field_end = 0;
            size_t field_index = 0;

            while ( field_end < tmp_line.size() ){
                field_end = min ( tmp_line.find('\t',field_start), tmp_line.find('\n', field_start) );
                string tmp_str = tmp_line.substr( field_start, field_end - field_start );
                field_end = min ( tmp_line.find('\t',field_start), tmp_line.find('\n', field_start) );

                if ( field_index == 2 ){
                    out_vec.push_back( strtod(tmp_str.c_str(), NULL) );
                }
                field_start = field_end+1;
                field_index++;
          }
          getline ( in_file, tmp_line );
      }
    } else {
        throw std::invalid_argument("Invalid input file. " + string (inchar) );

    }
    in_file.close();
}


int main(){
    // Read in input
    vector <double> plaf;
    vector <double> refCount;
    vector <double> altCount;


    (void)readFileLines( "tests/PG0390_first100ref.txt", refCount);
    dout <<" refCount.size() = "<< refCount.size()<<endl;
    (void)readFileLines( "tests/PG0390_first100alt.txt", altCount);
    dout <<" altCount.size() = "<< altCount.size()<<endl;
    (void)readFileLines( "tests/labStrains_first100_PLAF.txt", plaf);
    dout <<" plaf.size() = "<< plaf.size()<<endl;


    //size_t nLoci = plaf.size();

    // Initilize mcmc
    McmcSample mcmcSamples;
    mcmcSamples.runMcmcChain();

    // Export log
    cout << "place holder!" << std::endl;
}
