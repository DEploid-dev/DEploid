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

#ifndef VARIANTINDEX
#define VARIANTINDEX


#include <vector>
#include <string>
#include <cassert>
#include "global.h"


using namespace std;

class ExcludeMarker;

class VariantIndex {
#ifdef UNITTEST
 friend class TestPanel;
 friend class TestTxtReader;
 friend class TestInitialHaplotypes;
#endif
 friend class DEploidIO;
 friend class TxtReader;
 friend class ExcludeMarker;
 friend class Panel;
 friend class IBDrecombProbs;
 friend class VcfReader;

  private:
    // Members
    vector <string> chrom_;
    vector < size_t > indexOfChromStarts_;
    vector < vector < int> > position_;
    vector < vector < int> > keptPosition_;
    size_t nLoci_;

    // For removing markers and positions
    void findWhoToBeKept (ExcludeMarker* excludedMarkers );
    virtual void removeMarkers ();

    /* Index of content/info will be kept */
    vector < size_t > indexOfContentToBeKept;
    /* Index of positions entry to be kept, this will have the same size as this->chrom_, */
    vector < vector < size_t > > indexOfPosToBeKept;

    bool doneGetIndexOfChromStarts_;
    bool doneGetIndexOfChromStarts() const { return doneGetIndexOfChromStarts_; }
    void setDoneGetIndexOfChromStarts ( const bool setTo ){ this->doneGetIndexOfChromStarts_ = setTo; }

    // Methods
    void init();
    void getIndexOfChromStarts();
    void removePositions();

  public:
    VariantIndex();
    virtual ~VariantIndex(){};
    void findAndKeepMarkers( ExcludeMarker* excludedMarkers );
};


#endif
