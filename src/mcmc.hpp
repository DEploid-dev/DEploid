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

#include <vector>
#include <iostream>
#include <iomanip>      // std::setw
#include "mersenne_twister.hpp"
#include "dEploidIO.hpp"
#include "panel.hpp"
#include "randomSample.hpp"   // src/codeCogs/randomSample.hpp
#include "ibd.hpp"

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

    vector <double> IBDpathChangeAt;
    vector <double> siteOfTwoSwitchOne;
    vector <double> siteOfTwoMissCopyOne;
    vector <double> siteOfTwoSwitchTwo;
    vector <double> siteOfTwoMissCopyTwo;
    vector <double> siteOfOneSwitchOne;
    vector <double> siteOfOneMissCopyOne;

    vector <double> currentIBDpathChangeAt;
    vector <double> currentsiteOfTwoSwitchOne;
    vector <double> currentsiteOfTwoMissCopyOne;
    vector <double> currentsiteOfTwoSwitchTwo;
    vector <double> currentsiteOfTwoMissCopyTwo;
    vector <double> currentsiteOfOneSwitchOne;
    vector <double> currentsiteOfOneMissCopyOne;
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
  friend class DEploidIO;
  public:
    //McmcMachinery();
    McmcMachinery( DEploidIO* dEplioidIO, McmcSample *mcmcSample, RandomGenerator* rg_, bool useIBD = false );
    ~McmcMachinery();
    void runMcmcChain( bool showProgress = true, bool useIBD = false, bool notInR = true );

  private:
    McmcSample* mcmcSample_;
  /* Variables */
    DEploidIO* dEploidIO_;
    Panel* panel_;
    IBDrecombProbs ibdRecombProbs;
    size_t kStrain_;
    void setKstrain ( const size_t setTo ){ this->kStrain_ = setTo;}
    size_t kStrain() const { return this->kStrain_;}

    size_t nLoci_;
    void setNLoci ( const size_t setTo ){ this->nLoci_ = setTo;}
    size_t nLoci() const { return this->nLoci_; }

    double burnIn_;
    size_t maxIteration_;
    size_t mcmcThresh_;
    size_t McmcMachineryRate_;
    int eventInt_;

    size_t strainIndex_;
    size_t strainIndex1_;
    size_t strainIndex2_;

    size_t seed_;
    RandomGenerator* hapRg_;
    RandomGenerator* mcmcEventRg_;
    RandomGenerator* propRg_;
    RandomGenerator* initialHapRg_;
    RandomGenerator* ibdRg_;

    //std::normal_distribution<double>* initialTitre_normal_distribution_;// (MN_LOG_TITRE, SD_LOG_TITRE);
    //std::normal_distribution<double>* deltaX_normal_distribution_;// (0, 1/PROP_SCALE);
    StandNormalRandomSample* stdNorm_;
    double initialTitreNormalVariable(){ return this->stdNorm_->genReal() * SD_LOG_TITRE + MN_LOG_TITRE; }
    //double deltaXnormalVariable(){ return this->stdNorm_->genReal() * 1.0/PROP_SCALE + MN_LOG_TITRE; }
    double deltaXnormalVariable(){ return this->stdNorm_->genReal() * SD_LOG_TITRE* 1.0/PROP_SCALE + MN_LOG_TITRE; }
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
    vector < double > cumExpectedWsaf_;

  /* Methods */
    void calcMaxIteration( size_t nSample, size_t McmcMachineryRate, double burnIn );
   /* Initialize */
    void initializeMcmcChain(bool useIBD);
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

    void sampleMcmcEvent(bool useIBD = false);
    void recordMcmcMachinery();
    bool recordingMcmcBool_;
    void writeLastFwdProb(bool useIBD);
    void updateReferencePanel(size_t inbreedingPanelSizeSetTo, size_t excludedStrain);
    void initializeUpdateReferencePanel(size_t inbreedingPanelSizeSetTo);
    void computeDiagnostics();

   /* IBD */
    double theta_;
    void setTheta(const double setTo) {this->theta_ = setTo;}
    double theta() const {return this->theta_;}

    double fSum;
    Hprior hprior;
    vector < vector <double> > llkSurf;
    vector <int> uniqueEffectiveKCount;
    vector < vector<double> > ibdTransProbs;
    vector < vector <double> > fm;
    vector <double> fSumState;
    vector <size_t> ibdPath;

    vector <double> computeStatePrior(double theta);
    vector <size_t> findWhichIsSomething(vector <size_t> tmpOp, size_t something);
    vector <double> computeLlkOfStatesAtSiteI( size_t siteI, double err = 0.01);
    vector <double> computeLlkAtAllSites(double err = 0.01);
    vector <double> averageProportion(vector < vector <double> > &proportion );

    void ibdInitializeEssentials();
    void makeLlkSurf(vector <double> altCount,
                     vector <double> refCount,
                     double scalingConst = 100.0,
                     double err = 0.01,
                     size_t gridSize=99);
    void ibdSampleMcmcEventStep();
    void makeIbdTransProbs();
    void initializePropIBD();
    void computeIbdPathFwdProb(vector <double> statePrior);
    void ibdSamplePath(vector <double> statePrior);
    void computeUniqueEffectiveKCount();
    void updateFmAtSiteI(vector <double> & prior,
                         vector <double> & llk);
    void ibdUpdateHaplotypesFromPrior();
    void ibdUpdateProportionGivenHap(vector <double> &llkAtAllSites);
    void computeAndUpdateTheta();
    vector < vector <double> > reshapeFm(vector <size_t> stateIdx);
    vector <string> getIBDprobsHeader();
    vector <double> getIBDprobsIntegrated(vector < vector <double> > &prob);
    // For painting IBD
    void buildPathProbabilityForPainting();
    void computeIbdPathBwdProb();
    void combineFwdBwd();


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
    int acceptUpdate;
};

#endif

