#include <iostream> // std::cout
//#include <string>
//#include <boost/math/special_functions/gamma.hpp>
//#include <math.h>
#include "mcmc.hpp"
#include "panel.hpp"
#include "param.hpp"
//#include<stdlib.h>     /* strtol, strtod */
//#include<climits> // INT_MAX


using namespace std;
//#include "param.hpp"



int main(){

    Input input( "tests/labStrains_first100_PLAF.txt",
           "tests/PG0390_first100ref.txt",
           "tests/PG0390_first100alt.txt",
           (size_t)5);

    Panel panel("tests/lab_first100_Panel.txt");

    // Initilize mcmc
    McmcMachinery McmcMachinerys(&input, &panel);
    McmcMachinerys.runMcmcChain();

    // Export log

    cout << "place holder!" << std::endl;
}
