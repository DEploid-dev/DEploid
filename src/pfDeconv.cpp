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
#include "mcmc.hpp"
#include "panel.hpp"
#include "pfDeconvIO.hpp"

using namespace std;

int main( int argc, char *argv[] ){
    try {

        PfDeconvIO pfDeconvIO( argc, argv );
        if ( pfDeconvIO.help() ){
            pfDeconvIO.printHelp();
            return EXIT_SUCCESS;
        }

        Panel panel(pfDeconvIO.panelFileName_.c_str());

        McmcSample * mcmcSample = new McmcSample();

        bool repeat = true;
        size_t repeat_counter = 0;
        while ( repeat && repeat_counter < 100 ){
            mcmcSample->clear();
            cout << "repeat_counter = "<< repeat_counter<<endl;
            // Initilize mcmc
            McmcMachinery McmcMachinerys(&pfDeconvIO, &panel, mcmcSample, 1000, 5 );
            McmcMachinerys.runMcmcChain();

            //if ( mcmcSample->llkRange() < 1000.0 ){
                repeat = false;
            //}

            repeat_counter++;
        }
        // Export log
        mcmcSample->output();

        delete mcmcSample;
    }
    catch (const exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}
