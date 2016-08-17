/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016, Sha (Joe) Zhu, Jacob Almagro and Prof. Gil McVean
 *
 * This file is part of dEploid.
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
#include "dEploidIO.hpp"

using namespace std;

int main( int argc, char *argv[] ){
    try {

        DEploidIO dEploidIO;
        (void)dEploidIO.core( argc, argv );

        if ( dEploidIO.version() ){
            dEploidIO.printVersion();
            return EXIT_SUCCESS;
        }

        if ( dEploidIO.help() ){
            dEploidIO.printHelp();
            return EXIT_SUCCESS;
        }

        Panel *panel = NULL; // Move panel to dEploidIO

        if ( dEploidIO.usePanel() ){
            panel = new Panel();
            panel->readFromFile(dEploidIO.panelFileName_.c_str());
            if ( dEploidIO.excludeSites() ){
                panel->findAndKeepMarkers( dEploidIO.excludedMarkers );
            }

            panel->computeRecombProbs( dEploidIO.averageCentimorganDistance(), dEploidIO.Ne(), dEploidIO.useConstRecomb(), dEploidIO.constRecombProb(), dEploidIO.forbidCopyFromSame() );
            panel->checkForExceptions( dEploidIO.nLoci(), dEploidIO.panelFileName_ );
        }

        McmcSample * mcmcSample = new McmcSample();

        McmcMachinery mcmcMachinery(&dEploidIO, panel, mcmcSample);
        mcmcMachinery.runMcmcChain();

        dEploidIO.write(mcmcSample, panel);

        if ( panel ){
            delete panel;
        }
        delete mcmcSample;
    }
    catch (const exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}
