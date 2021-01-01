/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
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

#include <iostream>  // std::cout
#include "dEploidIO.hpp"

int main(int argc, char *argv[]) {
    try {
        DEploidIO dEploidIO(argc, argv);
        std::ostream *output = &std::cout;

        if ( dEploidIO.version() ) {
            dEploidIO.operation_printVersion(*output);
            return EXIT_SUCCESS;
        }

        if ( dEploidIO.help() ) {
            dEploidIO.operation_printHelp(*output);
            return EXIT_SUCCESS;
        }

        if (dEploidIO.doComputeLLK()) {
            dEploidIO.computeLLKfromInitialHap();
        } else if (dEploidIO.doLsPainting()) {
            dEploidIO.operation_chromPainting();
        } else if (dEploidIO.doIbdPainting()) {
            dEploidIO.operation_paintIBD();
        } else if (dEploidIO.doIbdViterbiPainting()) {
            dEploidIO.operation_paintIBDviterbi();
        } else if (dEploidIO.useLasso()) {  // DEploid-Lasso
            dEploidIO.workflow_lasso();
        } else if (dEploidIO.useBestPractice()) {  // best practice
            dEploidIO.workflow_best();
        } else {  // classic version, and DEploid-IBD
            dEploidIO.workflow_ibd();
        }
        // Finishing, write log
        dEploidIO.wrapUp();
    }
    catch (const exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
