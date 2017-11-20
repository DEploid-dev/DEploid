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

#ifndef PANEL
#define PANEL

#include "txtReader.hpp"
#include "exceptions.hpp"

class Panel: public TxtReader{
#ifdef UNITTEST
 friend class TestPanel;
 friend class TestInitialHaplotypes;
 friend class TestUpdateHap;
 friend class TestUpdatePairHap;
 friend class TestUpdateSingleHap;
#endif
 friend class McmcMachinery;
 friend class UpdateSingleHap;
 friend class UpdatePairHap;
 friend class UpdateHap;
 friend class DEploidIO;
 friend class InitialHaplotypes;
  private:
    // Members
    vector < double > pRec_;
    // Used in update single haplotype
    vector < double > pRecEachHap_; // = pRec / nPanel_;
    vector < double > pNoRec_; // = 1.0 - pRec;
    // Used in update pair of haplotypes
    vector < double > pRecRec_; // pRecEachHap * pRecEachHap;
    vector < double > pRecNoRec_; // pRecEachHap * pNoRec;
    vector < double > pNoRecNoRec_; // pNoRec * pNoRec;

    size_t truePanelSize_;
    void setTruePanelSize ( const size_t setTo ){ this->truePanelSize_ = setTo; }

    size_t inbreedingPanelSize_;
    void setInbreedingPanelSize ( const size_t setTo ){ this->inbreedingPanelSize_ = setTo; }

    size_t inbreedingPanelSize() const { return this->inbreedingPanelSize_; }
    size_t truePanelSize() const { return this->truePanelSize_; }
    Panel();
    //Panel(const char inchar[] );
    virtual ~Panel(){};

    // Methods
    void readFromFile( const char inchar[] );
    void computeRecombProbs( double averageCentimorganDistance, double Ne, bool useConstRecomb, double constRecombProb, bool forbidCopyFromSame );
    void checkForExceptions( size_t nLoci, string panelFileName );
    void initializeUpdatePanel( size_t inbreedingPanelSizeSetTo);
    void updatePanelWithHaps( size_t inbreedingPanelSizeSetTo, size_t excludedStrain, vector < vector<double> > & haps);

    void print();
    void buildExamplePanelContent();
    void buildExamplePanel1();
    void buildExamplePanel2();
};


class InitialHaplotypes: public Panel{
#ifdef UNITTEST
 friend class TestInitialHaplotypes;
#endif
 friend class DEploidIO;
    InitialHaplotypes():Panel(){}
    ~InitialHaplotypes(){}
};


class IBDrecombProbs: public VariantIndex{
 friend class IBDpath;

#ifdef UNITTEST
 friend class TestIBDpath;
#endif

  private:
    vector < double > pRec_;
    vector < double > pNoRec_; // = 1.0 - pRec;

    void computeRecombProbs( double averageCentimorganDistance, double Ne, bool useConstRecomb, double constRecombProb );

  public:
    IBDrecombProbs():VariantIndex(){};
    IBDrecombProbs(vector < vector < int> > position, size_t nLoci){
        this->position_ = position;
        this->nLoci_ = nLoci;
    }
    ~IBDrecombProbs(){}
};

#endif
