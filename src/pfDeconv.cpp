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
