/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016, Sha (Joe) Zhu, Jacob Almagro and Prof. Gil McVean
 *
 * This file is part of dEploid.
 *
 * dEploid is free software: you can redistribute it and/or modify
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
 *
 */

#include <iostream> // std::cout
#include "mcmc.hpp"
#include "panel.hpp"
#include "dEploidIO.hpp"

using namespace std;

int main( int argc, char *argv[] ){
    try {

        DEploidIO dEploidIO(argc, argv);
        std::ostream *output = &std::cout;

        if ( dEploidIO.version() ){
            dEploidIO.printVersion(*output);
            return EXIT_SUCCESS;
        }

        if ( dEploidIO.help() ){
            dEploidIO.printHelp(*output);
            return EXIT_SUCCESS;
        }

        if ( dEploidIO.doPainting() ){
            dEploidIO.chromPainting();
        } else{
            McmcSample * mcmcSample = new McmcSample();
            MersenneTwister rg(dEploidIO.randomSeed());

            McmcMachinery mcmcMachinery(&dEploidIO, mcmcSample, &rg);
            mcmcMachinery.runMcmcChain(true);
            delete mcmcSample;
        }
        // Finishing, write log
        dEploidIO.wrapUp();
    }
    catch (const exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}
