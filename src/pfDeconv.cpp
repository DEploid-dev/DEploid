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

        Panel *panel = NULL;

        if ( pfDeconvIO.usePanel() ){
            panel = new Panel(pfDeconvIO.panelFileName_.c_str());
            if ( pfDeconvIO.exclude_sites_ ){
                panel->removeMarkers( pfDeconvIO.excludedMarkers );
            }

            panel->computeRecombProbs( pfDeconvIO.averageCentimorganDistance(), pfDeconvIO.Ne(), pfDeconvIO.useConstRecomb(), pfDeconvIO.constRecombProb() );
            panel->checkForExceptions( pfDeconvIO.nLoci(), pfDeconvIO.panelFileName_ );
        }

        McmcSample * mcmcSample = new McmcSample();

        McmcMachinery mcmcMachinery(&pfDeconvIO, panel, mcmcSample);
        mcmcMachinery.runMcmcChain();

        pfDeconvIO.write(mcmcSample, panel);

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
