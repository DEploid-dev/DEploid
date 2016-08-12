/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016, Sha (Joe) Zhu, Jacob Almagro and Prof. Gil McVean
 *
 * This file is part of dEploid.
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

#ifndef VARIANTINDEX
#define VARIANTINDEX


#include <vector>
#include <string>
#include <cassert>
#include "global.h"

using namespace std;

class ExcludeMarker;

class VariantIndex {
 friend class DEploidIO;
 friend class TxtReader;
 friend class ExcludeMarker;
 friend class Panel;
 friend class VcfReader;

  private:
    // Members
    vector <string> chrom_;
    vector < size_t > indexOfChromStarts_;
    vector < vector < double> > position_;
    size_t nLoci_;
    vector < size_t > findWhoToBeRemoved (ExcludeMarker* excludedMarkers );
    virtual void removeMarkers ( vector < size_t > & index ){};

  public:
    VariantIndex(){};
    virtual ~VariantIndex(){};
    void findAndRemoveMarkers( ExcludeMarker* excludedMarkers ){
        vector < size_t > indexOfContentToBeRemoved = this->findWhoToBeRemoved (excludedMarkers );
        this->removeMarkers (indexOfContentToBeRemoved);
    }
};


#endif
