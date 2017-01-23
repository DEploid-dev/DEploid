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

#include <vector>
#include <iostream>
#include <iomanip>      // std::setw
#include "mersenne_twister.hpp"
#include "dEploidIO.hpp"
#include "panel.hpp"
#include "utility.hpp"
#include "global.h"
#include "randomSample.hpp"

#ifndef MCMC
#define MCMC

using namespace std;


class McmcSample {
#ifdef UNITTEST
  friend class TestMcmcMachinery;
#endif
  friend class McmcMachinery;
  friend class DEploidIO;
  friend class RMcmcSample;
  public:
    McmcSample();
    ~McmcSample();
    void clear(){
        proportion.clear();
        sumLLKs.clear();
        moves.clear();
    }

    vector < vector <double> > proportion;
    vector < vector <double> > hap;
    vector < double > sumLLKs;
  private:
    vector < int > moves;
};


class McmcMachinery {
#ifdef UNITTEST
  friend class TestMcmcMachinery;
#endif
  public:
    //McmcMachinery();
    McmcMachinery( DEploidIO* dEploidIO, Panel *panel, McmcSample *mcmcSample, RandomGenerator* rg_ );
    ~McmcMachinery();
    void runMcmcChain( bool showProgress = true );

  private:
    McmcSample* mcmcSample_;
  /* Variables */
    DEploidIO* dEploidIO_;
    Panel* panel_;
    size_t kStrain_;
    size_t nLoci_;

    double burnIn_;
    size_t maxIteration_;
    size_t mcmcThresh_;
    size_t McmcMachineryRate_;
    int eventInt_;

    size_t strainIndex1_;
    size_t strainIndex2_;
    size_t strainIndex3_;

    size_t seed_;
    RandomGenerator* hapRg_;
    RandomGenerator* mcmcEventRg_;
    RandomGenerator* propRg_;
    RandomGenerator* initialHapRg_;

    //std::normal_distribution<double>* initialTitre_normal_distribution_;// (MN_LOG_TITRE, SD_LOG_TITRE);
    //std::normal_distribution<double>* deltaX_normal_distribution_;// (0, 1/PROP_SCALE);
    StandNormalRandomSample* stdNorm_;
    double initialTitreNormalVariable(){ return this->stdNorm_->genReal() * SD_LOG_TITRE + MN_LOG_TITRE; }
    double deltaXnormalVariable(){ return this->stdNorm_->genReal() * 1.0/PROP_SCALE + MN_LOG_TITRE; }
    double MN_LOG_TITRE;
    double SD_LOG_TITRE;
    double PROP_SCALE;

    size_t currentMcmcIteration_;
    vector <double> currentTitre_;
    double currentLogPriorTitre_;
    vector <double> currentProp_;
    vector <double> currentLLks_;
    vector < vector <double> > currentHap_;
    vector < double > currentExpectedWsaf_;

  /* Methods */
    void calcMaxIteration( size_t nSample, size_t McmcMachineryRate, double burnIn );
   /* Initialize */
    void initializeMcmcChain();
    void initializeProp();
    void initializeTitre();
    void initializeHap();
    void initializellk();
    void initializeExpectedWsaf();

    vector <double> calcExpectedWsaf(vector <double> &proportion );
    vector <double> titre2prop(vector <double> & tmpTitre);

    double calcLogPriorTitre( vector <double> &tmpTitre);
    double rBernoulli(double p);

    void printArray ( vector <double> array ){
        for (auto const& value: array){
            cout << value << " ";
        }
        cout << endl;
    }

    void sampleMcmcEvent();
    void recordMcmcMachinery();
    bool recordingMcmcBool_;
    void writeLastFwdProb();

  /* Moves */
    void updateProportion();
     vector <double> calcTmpTitre();
     double deltaLLKs ( vector <double> &newLLKs );

    void updateSingleHap();
     void findUpdatingStrainSingle( );

    void updatePairHaps();
     //vector <size_t> sampleNoReplace(MersenneTwister* rg, vector <double> & proportion, size_t nSample );
     void findUpdatingStrainPair( );

    void updateThreeHaps();
     void findUpdatingStrainThree( );

  /* Debug */
    bool doutProp();
    bool doutLLK();
};

#endif

