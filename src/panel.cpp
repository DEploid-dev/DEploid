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

#include "panel.hpp"
#include <math.h>
#include <iostream>

Panel::Panel():InputMarker(){
    this->nPanel_ = this->nInfoLines_;
};

Panel::Panel(const char inchar[] ):InputMarker(inchar){
    this->nPanel_ = this->nInfoLines_;
};

Panel::~Panel(){};

void Panel::checkForExceptions( size_t nLoci, string panelFileName ){
    if ( this->content_.size() != nLoci ){
        throw LociNumberUnequal( panelFileName );
    }

    if ( this->pRec_.size() != nLoci ){
        throw LociNumberUnequal( panelFileName );
    }

}

void Panel::computeRecombProbs( double averageCentimorganDistance, double Ne, bool useConstRecomb, double constRecombProb, bool forbidCopyFromSame ){
    assert(pRec_.size() == 0 );
    assert(pRecEachHap_.size() == 0 );
    assert(pNoRec_.size() == 0 );
    assert(pRecRec_.size() == 0 );
    assert(pRecNoRec_.size() == 0 );
    assert(pNoRecNoRec_.size() == 0 );

    double averageMorganDistance = averageCentimorganDistance * 100;
    double geneticDistance;
    double rho;
    double nPanelDouble = (double)this->nPanel_;
    double nPanlelMinus1 = nPanelDouble - 1.0;
    for ( size_t i = 0; i < this->position_.size(); i++){
        for ( size_t j = 1; j < this->position_[i].size(); j++){
            geneticDistance = (this->position_[i][j] - this->position_[i][j-1])/averageMorganDistance ;
            rho = geneticDistance * 2 * Ne;

            double pRecTmp = ( useConstRecomb ) ? constRecombProb : 1.0 - exp(-rho);
            this->pRec_.push_back( pRecTmp );

            double pRecEachHapTmp = pRecTmp / nPanelDouble;
            this->pRecEachHap_.push_back( pRecTmp / nPanelDouble );

            double pNoRecTmp = 1.0 - pRecTmp;
            this->pNoRec_.push_back( pNoRecTmp );

            double secondPRecEachHapTmp = ( forbidCopyFromSame ) ? (pRecTmp / nPanlelMinus1) : pRecEachHapTmp; // allowing copy from the same

            this->pRecRec_.push_back ( pRecEachHapTmp * secondPRecEachHapTmp );
            this->pRecNoRec_.push_back ( secondPRecEachHapTmp * pNoRecTmp );
            this->pNoRecNoRec_.push_back ( pNoRecTmp * pNoRecTmp );
        }
        this->pRec_.push_back(1.0);
        this->pRecEachHap_.push_back( 1.0 / nPanelDouble );
        this->pNoRec_.push_back(0.0);
        this->pRecRec_.push_back ( ( ( forbidCopyFromSame ) ? (1.0 / nPanelDouble / nPanlelMinus1) : (1.0 / nPanelDouble / nPanelDouble) ) );
        this->pRecNoRec_.push_back ( 0.0 );
        this->pNoRecNoRec_.push_back ( 0.0 );
    }
    assert(pRec_.size() == this->nLoci_ );
    assert(pRecEachHap_.size() == this->nLoci_ );
    assert(pNoRec_.size() == this->nLoci_ );
    assert(pRecRec_.size() == this->nLoci_ );
    assert(pRecNoRec_.size() == this->nLoci_ );
    assert(pNoRecNoRec_.size() == this->nLoci_ );
}


//vector<vector<double>> outtrans(out[0].size(),
                                    //vector<double>(out.size()));
    //for (size_t i = 0; i < out.size(); ++i)
        //for (size_t j = 0; j < out[0].size(); ++j)
            //outtrans[j][i] = out[i][j];
