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

#ifndef PANEL
#define PANEL

#include "atMarker.hpp"
#include "exceptions.hpp"

class Panel: public AtMarker{
 friend class TestPanel;
 friend class UpdateSingleHap;
 friend class UpdatePairHap;
 friend class UpdateHap;

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

    size_t nPanel_;

    void computeRecombProbs( double averageCentimorganDistance = 15000.0, double Ne = 10.0 );

  public:
    Panel(const char inchar[], size_t nLociForChecking ):AtMarker(inchar){
        if ( this->content_.size() != nLociForChecking ){
            throw LociNumberUnequal( string(inchar) );
        }
        this->nPanel_ = this->nInfoLines_;
        this->computeRecombProbs();
    };
    ~Panel(){};

    // Methods
    void print();
};

#endif
