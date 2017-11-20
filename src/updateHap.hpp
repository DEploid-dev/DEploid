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
#include "utility.hpp"
#include "panel.hpp"

#ifndef HAP
#define HAP

using namespace std;


class UpdateHap{
#ifdef UNITTEST
  friend class TestUpdatePairHap;
  friend class TestUpdateSingleHap;
  friend class TestUpdateHap;
#endif
  friend class McmcMachinery;
  friend class UpdateSingleHap;
  friend class UpdatePairHap;
  friend class DEploidIO;

  public:
    size_t nPanel() const { return this->nPanel_; }
  private:
    UpdateHap();
    UpdateHap( vector <double> &refCount,
               vector <double> &altCount,
               vector <double> &plaf,
               vector <double> &expectedWsaf,
               vector <double> &proportion,
               vector < vector <double> > &haplotypes,
               RandomGenerator* rg,
               size_t segmentStartIndex,
               size_t nLoci,
               Panel* panel,
               double missCopyProb,
               double scalingFactor);
    virtual ~UpdateHap();

    Panel* panel_;
    double missCopyProb_;
    RandomGenerator* recombRg_;
    RandomGenerator* recombLevel2Rg_;
    RandomGenerator* missCopyRg_;

    size_t kStrain_;
    size_t nPanel_;
    void setPanelSize ( const size_t setTo ){ this->nPanel_ = setTo; }

    vector <double> newLLK;

    size_t segmentStartIndex_;
    size_t nLoci_;

    vector < vector <double> > emission_;
    double scalingFactor_;
    double scalingFactor() const {return this->scalingFactor_; }
    void setScalingFactor ( const double setTo ){ this->scalingFactor_ = setTo; }

    // Methods
    virtual void core(vector <double> &refCount,
                           vector <double> &altCount,
                           vector <double> &plaf,
                           vector <double> &expectedWsaf,
                           vector <double> &proportion,
                           vector < vector <double> > &haplotypes );
    virtual void calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes);
    virtual void calcHapLLKs( vector <double> &refCount, vector <double> &altCount);
    virtual void buildEmission( double missCopyProb );
    // calcFwdProbs() differ for class UpdateSingleHap and UpdatePairHap
    //virtual void calcFwdProbs();
    virtual void samplePaths();
    virtual void addMissCopying( double missCopyProb );
    virtual void updateLLK();
    virtual void sampleHapIndependently(vector <double> &plaf);
};


class UpdateSingleHap : public UpdateHap{
#ifdef UNITTEST
  friend class TestUpdateSingleHap;
#endif
 friend class McmcMachinery;
 friend class DEploidIO;
  //public:
  private:
    UpdateSingleHap ();
    UpdateSingleHap( vector <double> &refCount,
                      vector <double> &altCount,
                      vector <double> &plaf,
                      vector <double> &expectedWsaf,
                      vector <double> &proportion,
                      vector < vector <double> > &haplotypes,
                      RandomGenerator* rg,
                      size_t segmentStartIndex,
                      size_t nLoci,
                      Panel* panel, double missCopyProb,
                      double scalingFactor,
                      size_t strainIndex );
    ~UpdateSingleHap();

    vector <double> siteOfOneSwitchOne;
    vector <double> siteOfOneMissCopyOne;
    vector < vector <double> > fwdProbs_;
    vector < vector < double > > bwdProbs_;
    vector < vector <double> > fwdBwdProbs_;

    size_t strainIndex_;
    vector <double> expectedWsaf0_;
    vector <double> expectedWsaf1_;
    vector <double> llk0_;
    vector <double> llk1_;

    vector <double> path_;
    vector <double> hap_;

    // Methods
    void core( vector <double> &refCount,
               vector <double> &altCount,
               vector <double> &plaf,
               vector <double> &expectedWsaf,
               vector <double> &proportion,
               vector < vector <double> > &haplotypes );
    void painting( vector <double> &refCount,
                   vector <double> &altCount,
                   vector <double> &expectedWsaf,
                   vector <double> &proportion,
                   vector < vector <double> > &haplotypes );
    void calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes);
    void calcHapLLKs( vector <double> &refCount, vector <double> &altCount);
    void buildEmission( double missCopyProb );
    void buildEmissionBasicVersion( double missCopyProb );
    void calcFwdProbs();
    void calcBwdProbs();
    void calcFwdBwdProbs();
    void samplePaths();
    void addMissCopying( double missCopyProb );
    void sampleHapIndependently(vector <double> &plaf);
    void updateLLK();
};


class UpdatePairHap : public UpdateHap{
#ifdef UNITTEST
 friend class TestUpdatePairHap;
#endif
 friend class McmcMachinery;
 friend class DEploidIO;
  public:
     UpdatePairHap();
     UpdatePairHap( vector <double> &refCount,
                      vector <double> &altCount,
                      vector <double> &plaf,
                      vector <double> &expectedWsaf,
                      vector <double> &proportion,
                      vector < vector <double> > &haplotypes,
                      RandomGenerator* rg,
                      size_t segmentStartIndex,
                      size_t nLoci,
                      Panel* panel, double missCopyProb,
                      double scalingFactor, bool forbidCopyFromSame,
                      size_t strainIndex1,
                      size_t strainIndex2 );
    ~UpdatePairHap();

  private:
    vector <double> siteOfTwoSwitchOne;
    vector <double> siteOfTwoMissCopyOne;
    vector <double> siteOfTwoSwitchTwo;
    vector <double> siteOfTwoMissCopyTwo;
    vector< vector < vector <double> > > fwdProbs_;

    size_t strainIndex1_;
    size_t strainIndex2_;
    bool forbidCopyFromSame_;

    vector <double> expectedWsaf00_;
    vector <double> expectedWsaf01_;
    vector <double> expectedWsaf10_;
    vector <double> expectedWsaf11_;
    vector <double> llk00_;
    vector <double> llk01_;
    vector <double> llk10_;
    vector <double> llk11_;
    vector <double> path1_;
    vector <double> path2_;
    vector <double> hap1_;
    vector <double> hap2_;

    // Methods
    void core(vector <double> &refCount,
                           vector <double> &altCount,
                           vector <double> &plaf,
                           vector <double> &expectedWsaf,
                           vector <double> &proportion,
                           vector < vector <double> > &haplotypes );

    void calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes);
    void calcHapLLKs( vector <double> &refCount, vector <double> &altCount);
    void buildEmission( double missCopyProb );
    void calcFwdProbs( bool forbidCopyFromSame );
    void samplePaths();
    void addMissCopying( double missCopyProb );
    void sampleHapIndependently(vector <double> &plaf);
    void updateLLK();

    // Own methods
    vector <double> computeRowMarginalDist( vector < vector < double > > & probDist );
    vector <double> computeColMarginalDist( vector < vector < double > > & probDist );
    vector <size_t> sampleMatrixIndex( vector < vector < double > > &probDist );
};

#endif
