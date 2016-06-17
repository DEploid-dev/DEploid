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

#include <vector>
#include <iostream>
#include <iomanip>      // std::setw
#include "mersenne_twister.hpp"
#include "pfDeconvIO.hpp"
#include "panel.hpp"
#include "utility.hpp"
#include "global.h"

#ifndef MCMC
#define MCMC

using namespace std;


class McmcSample {
#ifdef UNITTEST
  friend class TestMcmcMachinery;
#endif
  friend class McmcMachinery;
  friend class PfDeconvIO;
  public:
    McmcSample();
    ~McmcSample();
    void clear(){
        proportion.clear();
        sumLLKs.clear();
        moves.clear();
    }


  private:
    vector < int > moves;
    vector < vector <double> > proportion;
    vector < vector <double> > hap;
    vector < double > sumLLKs;
};


class McmcMachinery {
#ifdef UNITTEST
  friend class TestMcmcMachinery;
#endif
  public:
    //McmcMachinery();
    McmcMachinery( PfDeconvIO* pdfDeconfIO, Panel *panel, McmcSample *mcmcSample );
    ~McmcMachinery();
    void runMcmcChain( );

  private:
    McmcSample* mcmcSample_;
  /* Variables */
    PfDeconvIO* pfDeconvIO_;
    Panel* panel_;
    size_t kStrain_;
    size_t nLoci_;

    double burnIn_;
    size_t maxIteration_;
    size_t mcmcThresh_;
    size_t McmcMachineryRate_;
    int eventInt_;

    size_t strainIndex_;
    size_t strainIndex1_;
    size_t strainIndex2_;

    size_t seed_;
    MersenneTwister* hapRg_;
    MersenneTwister* mcmcEventRg_;
    MersenneTwister* propRg_;
    MersenneTwister* initialHapRg_;

    std::default_random_engine* std_generator_;// (this->seed_);
    std::normal_distribution<double>* initialTitre_normal_distribution_;// (MN_LOG_TITRE, SD_LOG_TITRE);
    std::normal_distribution<double>* deltaX_normal_distribution_;// (0, 1/PROP_SCALE);
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

  /* Debug */
    bool doutProp();
    bool doutLLK();
};

#endif

