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

//enum EventType {PROPORTION, SINGLE, PAIR};

class McmcSample {
  friend class McmcMachinery;
  public:
    McmcSample(){}
    ~McmcSample(){}
    double llkRange(){
        return (maxOfVec(this->sumLLKs) - minOfVec(this->sumLLKs) );
    }
    void clear(){
        proportion.clear();
        sumLLKs.clear();
        moves.clear();
    }

    void output(){
        ofstream prop_file( "tmp.prop", ios::out | ios::app | ios::binary );
        for ( size_t i = 0; i < this->proportion.size(); i++){
            for ( size_t ii = 0; ii < this->proportion[i].size(); ii++){
                prop_file << setw(10) << this->proportion[i][ii];
                prop_file << ((ii < (this->proportion[i].size()-1)) ? "\t" : "\n") ;
            }
        }
        prop_file.close();


        ofstream llk_file( "tmp.llk", ios::out | ios::app | ios::binary );
        for ( size_t i = 0; i < this->sumLLKs.size(); i++){
            //if ( this->moves[i] != 2)
            llk_file << this->moves[i] << "\t" << this->sumLLKs[i] << endl;
        }
        llk_file.close();

        ofstream hap_file( "tmp.hap", ios::out | ios::app | ios::binary );
        for ( size_t i = 0; i < this->hap.size(); i++ ){
            for ( size_t ii = 0; ii < this->hap[i].size(); ii++){
                hap_file << this->hap[i][ii];
                hap_file << ((ii < (this->hap[i].size()-1)) ? "\t" : "\n") ;
            }
        }
        hap_file.close();
    }


  private:
    vector < int > moves;
    vector < vector <double> > proportion;
    vector < vector <double> > hap;
    vector < double > sumLLKs;
};


class McmcMachinery {
  public:
    McmcMachinery( PfDeconvIO* input, Panel *panel, McmcSample *mcmcSample,
                size_t nSample = 100, size_t McmcMachineryRate = 5 );
    ~McmcMachinery();
    void runMcmcChain( );

  private:
    McmcSample *mcmcSample_;
  /* Variables */
    PfDeconvIO* input_;
    Panel* panel_;
    size_t kStrain_;
    size_t nLoci_;

    double burnIn_;
    size_t maxIteration_;
    size_t mcmcThresh_;
    size_t McmcMachineryRate_;
    int eventInt_;

    size_t seed_;
    MersenneTwister* rg_;
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
    void calcMaxIteration( size_t nSample, size_t McmcMachineryRate );
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

  /* Moves */
    void updateProportion();
     vector <double> calcTmpTitre();
     double deltaLLKs ( vector <double> &newLLKs );

    void updateSingleHap();
    void updatePairHaps();

};

#endif

