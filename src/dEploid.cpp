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
#include "mcmc.hpp"
#include "dEploidIO.hpp"
#include "chooseK.hpp"

int main(int argc, char *argv[]) {
    try {
        DEploidIO dEploidIO(argc, argv);
        std::ostream *output = &std::cout;

        if ( dEploidIO.version() ) {
            dEploidIO.printVersion(*output);
            return EXIT_SUCCESS;
        }

        if ( dEploidIO.help() ) {
            dEploidIO.printHelp(*output);
            return EXIT_SUCCESS;
        }

        if (dEploidIO.doComputeLLK()) {
            dEploidIO.computeLLKfromInitialHap();
        } else if (dEploidIO.doLsPainting()) {
            dEploidIO.chromPainting();
        } else if (dEploidIO.doIbdPainting()) {
            dEploidIO.paintIBD();
        } else if (dEploidIO.useLasso()) {
            dEploidIO.dEploidLasso();
            MersenneTwister lassoRg(dEploidIO.randomSeed());
            DEploidIO tmpIO(dEploidIO);
            vector < vector <double> > hap;
            for (size_t chromi = 0;
                 chromi < dEploidIO.indexOfChromStarts_.size();
                 chromi++ ) {
                tmpIO.position_.clear();
                tmpIO.position_.push_back(dEploidIO.position_.at(chromi));
                tmpIO.indexOfChromStarts_.clear();
                tmpIO.indexOfChromStarts_.push_back(0);
                string job = string("DEploid-Lasso learning chromosome ");
                job.append(dEploidIO.chrom_[chromi]).append(" haplotypes");
                McmcSample * lassoMcmcSample = new McmcSample();
                McmcMachinery lassoMcmcMachinery(
                                            &dEploidIO.lassoPlafs.at(chromi),
                                            &dEploidIO.lassoRefCount.at(chromi),
                                            &dEploidIO.lassoAltCount.at(chromi),
                                            dEploidIO.lassoPanels.at(chromi),
                                            &tmpIO,
                                            job,
                                            "lasso",
                                            lassoMcmcSample,
                                            &lassoRg,
                                            false);
                lassoMcmcMachinery.runMcmcChain(true,   // show progress
                                                false);  // use IBD
                for (size_t snpi = 0;
                     snpi < lassoMcmcSample->hap.size(); snpi++) {
                    hap.push_back(vector <double> (
                                        lassoMcmcSample->hap[snpi].begin(),
                                        lassoMcmcSample->hap[snpi].end()));
                }
                delete lassoMcmcSample;
            }
            dEploidIO.writeHap(hap, "lasso");
        } else if (dEploidIO.useBestPractice()) {  // best practice
            MersenneTwister rg(dEploidIO.randomSeed());
            if (dEploidIO.usePanel()) {
                // #############################################################
                // ################# DEploid-LASSO to learn K ##################
                // #############################################################

                ChooseK chooseK;

                DEploidIO toLearnK(dEploidIO);
                toLearnK.dEploidLassoTrimfirst();
                DEploidIO toLearnKtmp(dEploidIO);
                for (size_t chromi = 0;
                     chromi < toLearnK.indexOfChromStarts_.size();
                     chromi++ ) {
                    toLearnKtmp.position_.clear();
                    toLearnKtmp.position_.push_back(
                                                toLearnK.position_.at(chromi));
                    toLearnKtmp.indexOfChromStarts_.clear();
                    toLearnKtmp.indexOfChromStarts_.push_back(0);

                    McmcSample * lassoMcmcSample = new McmcSample();
                    string job = string("DEploid-Lasso learning K ");
                    job.append(toLearnK.chrom_[chromi]).append(" leanrning K");
                    McmcMachinery lassoMcmcMachinery(
                                            &toLearnK.lassoPlafs.at(chromi),
                                            &toLearnK.lassoRefCount.at(chromi),
                                            &toLearnK.lassoAltCount.at(chromi),
                                            toLearnK.lassoPanels.at(chromi),
                                            &toLearnKtmp,
                                            job,
                                            job,
                                            lassoMcmcSample,
                                            &rg,
                                            false);
                    lassoMcmcMachinery.runMcmcChain(true,   // show progress
                                                    false);  // use IBD
                    toLearnKtmp.initialProp = toLearnKtmp.finalProp;
                    toLearnKtmp.setInitialPropWasGiven(true);
                    chooseK.appendProportions(toLearnKtmp.finalProp);
                    delete lassoMcmcSample;
                }
                // chooseK.findKmode();

                dEploidIO.initialProp.clear();
                vector <double> initialP = chooseK.chosenP();
                for (auto const& value : initialP) {
                    if (value > 0.01) {
                        dEploidIO.initialProp.push_back(value);
                    }
                }
                dEploidIO.setKstrain(dEploidIO.initialProp.size());
                dEploidIO.setInitialPropWasGiven(true);
            }
            if (dEploidIO.initialProp.size() > 1) {
                // #############################################################
                // ###################### DEploid-IBD   ########################
                // #############################################################
                McmcSample * ibdMcmcSample = new McmcSample();
                DEploidIO tmpIO2(dEploidIO);
                tmpIO2.ibdTrimming();
                McmcMachinery ibdMcmcMachinery(&tmpIO2.plaf_,
                                      &tmpIO2.refCount_,
                                      &tmpIO2.altCount_,
                                      tmpIO2.panel,
                                      &tmpIO2,
                                      string("DEploid-IBD learning proportion"),
                                      string("ibd"),
                                      ibdMcmcSample,
                                      &rg,
                                      true);
                ibdMcmcMachinery.runMcmcChain(true,   // show progress
                                              true);  // use IBD
                // if (dEploidIO.useIbdOnly()) {
                    // tmpIO.paintIBD();
                    // dEploidIO.finalProp = tmpIO.initialProp;
                // }

                tmpIO2.paintIBD();

                vector <double> initialP;
                for (auto const& value : tmpIO2.finalProp) {
                    if (value > 0.01) {
                        initialP.push_back(value);
                    }
                }
                cout << endl;
                dEploidIO.initialProp = initialP;
                dEploidIO.finalProp = initialP;
                dEploidIO.setKstrain(initialP.size());
                dEploidIO.setInitialPropWasGiven(true);
                dEploidIO.setDoUpdateProp(false);
                delete ibdMcmcSample;
            }
            // #################################################################
            // ###################### DEploid-LASSO ############################
            // #################################################################
            // *** Frist split the reference panel etc
            dEploidIO.dEploidLasso();

            DEploidIO dEploidLassoIO(dEploidIO);
            dEploidLassoIO.initialProp = dEploidIO.initialProp;
            dEploidLassoIO.setDoUpdateProp(false);
            dEploidLassoIO.setInitialPropWasGiven(true);
            dEploidLassoIO.setKstrain(dEploidIO.kStrain());
            vector < vector <double> > hap;
            for (size_t chromi = 0;
                 chromi < dEploidIO.indexOfChromStarts_.size();
                 chromi++ ) {
                dEploidLassoIO.position_.clear();
                dEploidLassoIO.position_.push_back(
                                            dEploidIO.position_.at(chromi));
                dEploidLassoIO.indexOfChromStarts_.clear();
                dEploidLassoIO.indexOfChromStarts_.push_back(0);

                McmcSample * lassoMcmcSample = new McmcSample();
                string job = string("DEploid-Lasso learning chromosome ");
                job.append(dEploidIO.chrom_[chromi]).append(" haplotypes");
                McmcMachinery lassoMcmcMachinery(
                                            &dEploidIO.lassoPlafs.at(chromi),
                                            &dEploidIO.lassoRefCount.at(chromi),
                                            &dEploidIO.lassoAltCount.at(chromi),
                                            dEploidIO.lassoPanels.at(chromi),
                                            &dEploidLassoIO,
                                            job,
                                            job,
                                            lassoMcmcSample,
                                            &rg,
                                            false);
                lassoMcmcMachinery.runMcmcChain(true,   // show progress
                                                false);  // use IBD
                // *** Append haplotypes to the end
                for (size_t snpi = 0;
                     snpi < lassoMcmcSample->hap.size(); snpi++) {
                    hap.push_back(vector <double> (
                                        lassoMcmcSample->hap[snpi].begin(),
                                        lassoMcmcSample->hap[snpi].end()));
                }
                delete lassoMcmcSample;
            }

            dEploidIO.paintIBD();
            dEploidIO.writeHap(hap, "final");
        } else { // classic version, with IBD flag
            if (dEploidIO.useIBD()) {  // ibd
                McmcSample * ibdMcmcSample = new McmcSample();
                MersenneTwister ibdRg(dEploidIO.randomSeed());
                McmcMachinery ibdMcmcMachinery(&dEploidIO.plaf_,
                                               &dEploidIO.refCount_,
                                               &dEploidIO.altCount_,
                                               dEploidIO.panel,
                                               &dEploidIO,
                                               "DEploid-IBD",
                                               "DEploid-IBD",
                                               ibdMcmcSample,
                                               &ibdRg,
                                               true);
            }

            McmcSample * mcmcSample = new McmcSample();
            MersenneTwister rg(dEploidIO.randomSeed());
            McmcMachinery mcmcMachinery(&dEploidIO.plaf_,
                                &dEploidIO.refCount_,
                                &dEploidIO.altCount_,
                                dEploidIO.panel,
                                &dEploidIO,
                                "DEploid classic version",
                                "classic", // brief
                                mcmcSample,
                                &rg,
                                false);  // use IBD
            mcmcMachinery.runMcmcChain(true,     // show progress
                               false);   // use IBD
            dEploidIO.paintIBD();
            dEploidIO.writeHap(mcmcSample->hap, "final");
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
