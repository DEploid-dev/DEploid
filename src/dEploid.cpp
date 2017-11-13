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

#include <iostream> // std::cout
#include "mcmc.hpp"
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

            if (dEploidIO.useIBD()){ // ibd
                McmcSample * ibdMcmcSample = new McmcSample();
                MersenneTwister ibdRg(dEploidIO.randomSeed());

                McmcMachinery ibdMcmcMachinery(&dEploidIO, ibdMcmcSample, &ibdRg, true);
                ibdMcmcMachinery.runMcmcChain(true, // show progress
                                              true);  // use IBD
                delete ibdMcmcSample;
            }
            McmcSample * mcmcSample = new McmcSample();
            MersenneTwister rg(dEploidIO.randomSeed());

            McmcMachinery mcmcMachinery(&dEploidIO, mcmcSample, &rg,
                                        false); // use IBD
            mcmcMachinery.runMcmcChain(true, // show progress
                                       false); // use IBD

            dEploidIO.paintIBD();
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



void DEploidIO::paintIBD(){
    vector <double> goodProp;
    vector <size_t> goodStrainIdx;
    for ( size_t i = 0; i < this->finalProp.size(); i++){
        if (this->finalProp[i] > 0.01){
            goodProp.push_back(this->finalProp[i]);
            goodStrainIdx.push_back(i);
        }
    }

    if (goodProp.size() == 1){
        return;
    }

    DEploidIO tmpDEploidIO; // (*this);
    tmpDEploidIO.setKstrain(goodProp.size());
    tmpDEploidIO.setInitialPropWasGiven(true);
    tmpDEploidIO.initialProp = goodProp;
    tmpDEploidIO.finalProp = goodProp;
    tmpDEploidIO.refCount_ = this->refCount_;
    tmpDEploidIO.altCount_ = this->altCount_;
    tmpDEploidIO.plaf_ = this->plaf_;
    tmpDEploidIO.nLoci_= this->nLoci();
    tmpDEploidIO.position_ = this->position_;
    tmpDEploidIO.chrom_ = this->chrom_;
    //tmpDEploidIO.writeLog (&std::cout);

    MersenneTwister tmpRg(this->randomSeed());
    McmcSample * tmpMcmcSample = new McmcSample();
    McmcMachinery tmpIBDmcmc(&tmpDEploidIO, tmpMcmcSample, &tmpRg, true);
    tmpIBDmcmc.buildPathProbabilityForPainting();

    vector < vector <double> > reshapedProbs = tmpIBDmcmc.reshapeFm(tmpIBDmcmc.hprior.stateIdx);
    this->ibdProbsHeader = tmpIBDmcmc.getIBDprobsHeader();
    this->ibdProbsIntegrated = tmpIBDmcmc.getIBDprobsIntegrated(reshapedProbs);
    this->writeIBDpostProb(reshapedProbs, this->ibdProbsHeader);

    delete tmpMcmcSample;
}

