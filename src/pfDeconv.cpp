/*
 * pfDeconv is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016, Sha (Joe) Zhu, Jacob Almagro and Prof. Gil McVean
 *
 * This file is part of pfDeconv.
 *
 * scrm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <iostream> // std::cout
#include <stdio.h>
//#include <string>
//#include <boost/math/special_functions/gamma.hpp>
//#include <math.h>
#include "mcmc.hpp"
#include "panel.hpp"
#include "param.hpp"
//#include "utility.hpp"
//#include<stdlib.h>     /* strtol, strtod */
//#include<climits> // INT_MAX


using namespace std;
//#include "param.hpp"



int main(){

    remove( "tmp.llk" );
    remove( "tmp.prop" );
    remove( "tmp.hap" );

    //Input input( "tests/labStrains_first100_PLAF.txt",
           //"tests/PG0390_first100ref.txt",
           //"tests/PG0390_first100alt.txt",
           //(size_t)5);

    //Panel panel("tests/lab_first100_Panel.txt");

    Input input( "tests/labStrains_samples_PLAF.txt",
           "tests/PG0394_ref.txt",
           "tests/PG0394_alt.txt",
           (size_t)5);
    Panel panel("tests/clonalPanel.csv");


    McmcSample * mcmcSample = new McmcSample();

    bool repeat = true;
    size_t repeat_counter = 0;
    while ( repeat && repeat_counter < 100 ){
        mcmcSample->clear();
        cout << "repeat_counter = "<< repeat_counter<<endl;
        // Initilize mcmc
        McmcMachinery McmcMachinerys(&input, &panel, mcmcSample, 1000, 5, 1);
        McmcMachinerys.runMcmcChain();

        //if ( mcmcSample->llkRange() < 600.0 ){
            repeat = false;
        //}
        repeat_counter++;
    }
    // Export log
    mcmcSample->output();

    delete mcmcSample;
    cout << "place holder!" << std::endl;
}
