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
#include <exceptions.hpp>
#include <sstream>
#include "utility.hpp"
#include "mersenne_twister.hpp"
#include "dEploidIO.hpp"

#ifndef IBD
#define IBD

using namespace std;

int nchoose2(int n);
bool twoVectorsAreSame(vector<int> vec1, vector<int> vec2);
vector < vector<int> > unique( vector < vector<int> > &mat );
vector<int> convertIntToBinary(int x, size_t len);
vector < vector <int> > enumerateBinaryMatrixOfK(size_t k);

// The IBDconfiguration is used for index, which should be non-negative, use int, any thing below zero should throw.
class IBDconfiguration{
#ifdef UNITTEST
  friend class TestIBDconfig;
  friend class TestHprior;
#endif
  friend class Hprior;
  friend class McmcMachinery;
    IBDconfiguration();
    ~IBDconfiguration();

    void buildIBDconfiguration(size_t k = 5);
    size_t kStrain_;
    void setKstrain(const size_t setTo) {this->kStrain_ = setTo;}
    size_t kStrain() const {return this->kStrain_;}

    vector < vector<int> > op;
    vector < vector<int> > pairToEmission;
    vector < vector<size_t> > pairList;
    vector < vector<int> > states;
    vector < size_t > effectiveK;

    size_t stateSize() const { return this->states.size(); }
    void enumerateOp();
    void makePairList();
    void makePairToEmission();
    void findUniqueState();
    void findEffectiveK();

    vector <int> makeTmpRow();
    vector <size_t> findWhichIsOne(vector <int> tmpOp);
    bool twoVectorsAreSame(vector<int> vec1, vector<int> vec2);
    vector <string> getIBDconfigureHeader();
};


class Hprior{
#ifdef UNITTEST
  friend class TestHprior;
  friend class TestMcmcMachinery;
  friend class TestIBDpath;
#endif
  friend class IBDpath;
  friend class McmcMachinery;
  friend class DEploidIO;
    Hprior();
    ~Hprior();

    void buildHprior(size_t kStrain, vector <double> &plaf);

    IBDconfiguration ibdConfig;
    size_t kStrain_;
    void setKstrain(const size_t setTo) {this->kStrain_ = setTo;}
    size_t kStrain() const {return this->kStrain_;}

    size_t nLoci_;
    void setnLoci(const size_t setTo) {this->nLoci_ = setTo;}
    size_t nLoci() const {return this->nLoci_;}

    vector <double> plaf_;
    vector < vector <double> > priorProb; // size: nState x nLoci
    vector < vector <double> > priorProbTrans; // size: nLoci x nState
    void transposePriorProbs();

    vector <size_t> stateIdx; // size: nState
    vector <size_t> stateIdxFreq;

    vector <vector <int> > hSet; // size: nState x kStrain
    size_t nState_;
    size_t nState() const {return this->nState_;}

    vector < size_t > effectiveK;
    size_t nPattern() const {return this->effectiveK.size();}
    vector <string> getIBDconfigureHeader();
};


class IBDpath{
#ifdef UNITTEST
  friend class TestIBDpath;
#endif
  friend class McmcMachinery;
  friend class DEploidIO;
    RandomGenerator* ibdRg_;

    double fSum;
    Hprior hprior;
    IBDrecombProbs ibdRecombProbs;
    vector < vector<double> > ibdTransProbs;
    vector < vector <double> > fm;
    vector <double> fSumState;
    vector <size_t> ibdConfigurePath;

    vector < vector <double> > bwd;
    vector < vector <double> > fwdbwd;

    IBDpath();

    ~IBDpath();

    size_t kStrain_;
    void setKstrain ( const size_t setTo ){ this->kStrain_ = setTo;}
    size_t kStrain() const { return this->kStrain_;}

    size_t nLoci_;
    void setNLoci ( const size_t setTo ){ this->nLoci_ = setTo;}
    size_t nLoci() const { return this->nLoci_; }

    double theta_;
    void setTheta(const double setTo) {this->theta_ = setTo;}
    double theta() const {return this->theta_;}

    vector <double> currentIBDpathChangeAt;

    vector < vector <double> > llkSurf;
    vector <int> uniqueEffectiveKCount;

    vector <double> IBDpathChangeAt;
    // Methods
    void computeAndUpdateTheta();
    void updateFmAtSiteI(vector <double> & prior,
                         vector <double> & llk);
    void ibdSamplePath(vector <double> statePrior);
    void makeIbdTransProbs();
    vector <double> computeEffectiveKPrior(double theta);
    vector <double> computeStatePrior(vector <double> effectiveKPrior);
    void makeLlkSurf(vector <double> altCount,
                     vector <double> refCount,
                     double scalingConst = 100.0,
                     double err = 0.01,
                     size_t gridSize=99);
    void computeUniqueEffectiveKCount();
    vector <double> computeLlkOfStatesAtSiteI(vector<double> proportion, size_t siteI, double err = 0.01);
    vector <size_t> findWhichIsSomething(vector <size_t> tmpOp, size_t something);

    // For painting IBD
    void buildPathProbabilityForPainting(vector <double> proportion);
    void computeIbdPathFwdProb(vector <double> proportion, vector <double> statePrior);
    void computeIbdPathBwdProb(vector <double> proportion, vector <double> effectiveKPrior, vector <double> statePrior);
    void combineFwdBwd(vector < vector <double>> &reshapedFwd, vector < vector <double>> &reshapedBwd);
    vector < vector <double> > reshapeProbs(vector < vector <double> >& probs);
    double bestPath(vector <double> proportion, double err = 0.01);

public:
    vector <string> getIBDprobsHeader();
    void init(DEploidIO &dEploidIO, RandomGenerator* rg);
};

#endif
