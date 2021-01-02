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

void DEploidIO::workflow_lasso() {
    this->dEploidLasso();
    MersenneTwister lassoRg(this->randomSeed());
    DEploidIO tmpIO(*this);
    vector < vector <double> > hap;
    for (size_t chromi = 0;
         chromi < this->indexOfChromStarts_.size();
         chromi++ ) {
        tmpIO.position_.clear();
        tmpIO.position_.push_back(this->position_.at(chromi));
        tmpIO.indexOfChromStarts_.clear();
        tmpIO.indexOfChromStarts_.push_back(0);
        string job = string("DEploid-Lasso learning chromosome ");
        job.append(this->chrom_[chromi]).append(" haplotypes");
        McmcSample * lassoMcmcSample = new McmcSample();
        McmcMachinery lassoMcmcMachinery(
                                    &this->lassoPlafs.at(chromi),
                                    &this->lassoRefCount.at(chromi),
                                    &this->lassoAltCount.at(chromi),
                                    this->lassoPanels.at(chromi),
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
    this->writeHap(hap, "lasso");
}


void DEploidIO::workflow_ibd() {
    if (this->useIBD()) {  // ibd
        McmcSample * ibdMcmcSample = new McmcSample();
        MersenneTwister ibdRg(this->randomSeed());
        McmcMachinery ibdMcmcMachinery(&this->plaf_,
                                       &this->refCount_,
                                       &this->altCount_,
                                       this->panel,
                                       this,
                                       "DEploid-IBD",
                                       "ibd",
                                       ibdMcmcSample,
                                       &ibdRg,
                                       true);
        ibdMcmcMachinery.runMcmcChain(true,   // show progress
                                      true);  // use IBD
        delete ibdMcmcSample;
    }

    McmcSample * mcmcSample = new McmcSample();
    MersenneTwister rg(this->randomSeed());
    McmcMachinery mcmcMachinery(&this->plaf_,
                        &this->refCount_,
                        &this->altCount_,
                        this->panel,
                        this,
                        "DEploid classic version",
                        "classic",  // brief
                        mcmcSample,
                        &rg,
                        false);  // use IBD
    mcmcMachinery.runMcmcChain(true,     // show progress
                       false);   // use IBD
    this->operation_paintIBD();
    this->writeHap(mcmcSample->hap, "final");
    delete mcmcSample;
}


void DEploidIO::workflow_best() {
    MersenneTwister rg(this->randomSeed());

    // #############################################################
    // ################# DEploid-LASSO to learn K ##################
    // #############################################################

    vector < vector <double> > chooseKhap;

    DEploidIO toLearnK(*this);
    toLearnK.dEploidLassoTrimfirst();
    DEploidIO toLearnKtmp(*this);
    for (size_t chromi = 0;
         chromi < toLearnK.indexOfChromStarts_.size();
         chromi++ ) {
        bool tryAgain = true;
        while (tryAgain) {
          toLearnKtmp.position_.clear();
          toLearnKtmp.position_.push_back(toLearnK.position_.at(chromi));
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
                                          false,  // use IBD
                                          true,   // notInR
                                          false);  // averageP
          toLearnKtmp.initialProp = toLearnKtmp.finalProp;
          toLearnKtmp.setInitialPropWasGiven(true);
          if (toLearnKtmp.acceptRatio() != 0) {
            this->chooseK.appendProportions(toLearnKtmp.finalProp);
            tryAgain = false;
          } else {
            toLearnKtmp.setInitialPropWasGiven(false);
            tryAgain = true;
          }

          if (!tryAgain) {
            for (size_t snpi = 0;
                snpi < lassoMcmcSample->hap.size(); snpi++) {
                chooseKhap.push_back(vector <double> (
                                    lassoMcmcSample->hap[snpi].begin(),
                                    lassoMcmcSample->hap[snpi].end()));
            }
          }
          delete lassoMcmcSample;
        }
    }
    // Gathering haplotype information
    toLearnK.writeHap(chooseKhap, "chooseK");

    // Gathering proportion information
    this->initialProp.clear();
    vector <double> initialP = this->chooseK.chosenP();
    for (auto const& value : initialP) {
        if (value > 0.01) {
            this->initialProp.push_back(value);
        }
    }
    normalizeBySum(this->initialProp);
    this->setKstrain(this->initialProp.size());
    this->setInitialPropWasGiven(true);

    if (this->inferBestPracticeP()) {
      if (this->initialProp.size() > 1) {
        // #############################################################
        // ###################### DEploid-IBD   ########################
        // #############################################################
        McmcSample * ibdMcmcSample = new McmcSample();
        DEploidIO tmpIO2(*this);
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
        // if (this->useIbdOnly()) {
            // tmpIO.paintIBD();
            // this->finalProp = tmpIO.initialProp;
        // }

        tmpIO2.operation_paintIBD();

        vector <double> initialP;
        for (auto const& value : tmpIO2.finalProp) {
            if (value > 0.01) {
                initialP.push_back(value);
            }
        }
        cout << endl;
        this->initialProp = initialP;
        this->finalProp = initialP;
        this->setKstrain(initialP.size());
        this->setInitialPropWasGiven(true);
        this->setDoUpdateProp(false);
        delete ibdMcmcSample;
      } else {
        this->finalProp.clear();
        this->finalProp.push_back(1.0);
        this->setKstrain(1);
      }
    }

    if (this->inferBestPracticeHap()) {
    // #################################################################
    // ###################### DEploid-LASSO ############################
    // #################################################################
    // *** Frist split the reference panel etc
      this->dEploidLasso();

      DEploidIO dEploidLassoIO(*this);
      dEploidLassoIO.initialProp = this->initialProp;
      dEploidLassoIO.setDoUpdateProp(false);
      dEploidLassoIO.setInitialPropWasGiven(true);
      dEploidLassoIO.setKstrain(this->kStrain());
      vector < vector <double> > hap;
      for (size_t chromi = 0;
           chromi < this->indexOfChromStarts_.size();
           chromi++ ) {
          dEploidLassoIO.position_.clear();
          dEploidLassoIO.position_.push_back(
                                      this->position_.at(chromi));
          dEploidLassoIO.indexOfChromStarts_.clear();
          dEploidLassoIO.indexOfChromStarts_.push_back(0);

          McmcSample * lassoMcmcSample = new McmcSample();
          string job = string("DEploid-Lasso learning chromosome ");
          job.append(this->chrom_[chromi]).append(" haplotypes");
          McmcMachinery lassoMcmcMachinery(
                                    &this->lassoPlafs.at(chromi),
                                    &this->lassoRefCount.at(chromi),
                                    &this->lassoAltCount.at(chromi),
                                    this->lassoPanels.at(chromi),
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
      this->writeHap(hap, "final");
      this->writeVcf(hap, dEploidLassoIO.initialProp, "final");
    }

    if (this->inferBestPracticeP() & (this->initialProp.size() > 1)) {
        this->operation_paintIBD();
    }
}

