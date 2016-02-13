#include <iostream> // std::cout
//#include <string>
//#include <boost/math/special_functions/gamma.hpp>
//#include <math.h>
#include "mcmc.hpp"
#include "panel.hpp"
#include "param.hpp"
//#include<iostream>
//#include<vector>
//#include<cassert>       // assert
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


int main(){

    Input input( "tests/labStrains_first100_PLAF.txt",
           "tests/PG0390_first100ref.txt",
           "tests/PG0390_first100alt.txt",
           (size_t)5);

    Panel panel("tests/lab_first100_Panel.txt");
    //size_t nLoci = plaf.size();

    // Initilize mcmc
    //McmcSample mcmcSamples(&input);
    //mcmcSamples.runMcmcChain();

    // Export log

    cout << "place holder!" << std::endl;
}
