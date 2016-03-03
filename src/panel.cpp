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

#include <fstream>
#include <iostream>
#include <cassert>
#include <math.h>
#include "panel.hpp"
#include "exceptions.hpp"

Panel::Panel(const char inchar[]){
    chromInex_ = -1;

    //cout << "Build Panel" << endl;
    ifstream in_file( inchar );
    string tmp_line;
    if ( in_file.good() ){
        getline ( in_file, tmp_line ); // skip the first line, which is the header
        getline ( in_file, tmp_line );
        while ( tmp_line.size() > 0 ){
            size_t field_start = 0;
            size_t field_end = 0;
            size_t field_index = 0;
            vector <double> contentRow;
            while ( field_end < tmp_line.size() ){
                field_end = min ( tmp_line.find(',',field_start), tmp_line.find('\n', field_start) );
                string tmp_str = tmp_line.substr( field_start, field_end - field_start );

                if ( field_index > 1 ){
                    contentRow.push_back( strtod(tmp_str.c_str(), NULL) );
                } else if ( field_index == 0 ){
                    this->extractChrom( tmp_str );
                } else if ( field_index == 1 ){
                    this->extractPOS ( tmp_str );
                }

                field_start = field_end+1;
                field_index++;
            }
            this->content_.push_back(contentRow);
            getline ( in_file, tmp_line );
        }
    } else {
        throw InvalidInputFile( string (inchar) );

    }
    in_file.close();

    this->position_.push_back( this->tmpPosition_ );

    this->nLoci_ = this->content_.size();
    this->nPanel_ = this->content_.back().size();
    //pRec_ = vector <double> (this->nLoci_, 0.0);
    //pRec_ = vector <double> (this->nLoci_, 0.01); // use constant recombination probability
    this->computeRecombProbs();

    assert ( chromInex_ > -1 );
}


void Panel::extractChrom( string & tmp_str ){
    if ( chromInex_ >= 0 ){
        if ( tmp_str != this->chrom_.back() ){
            chromInex_++;
            // save current positions
            this->position_.push_back(this->tmpPosition_);

            // start new chrom
            this->tmpPosition_.clear();
            this->chrom_.push_back(tmp_str);
        }
    } else {
        chromInex_++;
        assert (this->chrom_.size() == 0);
        this->chrom_.push_back( tmp_str );
        assert ( this->tmpPosition_.size() == 0 );
        assert ( this->position_.size() == 0);
    }
}

void Panel::extractPOS( string & tmp_str ){
    this->tmpPosition_.push_back(strtod(tmp_str.c_str(), NULL));
}


void Panel::print(){
    for ( size_t i = 0; i < this->content_.size(); i++){
        for (size_t j = 0; j < this->content_[i].size(); j++){
            cout <<this->content_[i][j]<<" ";
        }
        cout<<endl;
    }
}


void Panel::computeRecombProbs( double averageCentimorganDistance, double Ne){
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
    for ( size_t i = 0; i < this->position_.size(); i++){
        for ( size_t j = 1; j < this->position_[i].size(); j++){
            geneticDistance = (this->position_[i][j] - this->position_[i][j-1])/averageMorganDistance ;
            rho = geneticDistance * 2 * Ne;

            double pRecTmp = 1.0 - exp(-rho);
            this->pRec_.push_back( pRecTmp );

            double pRecEachHapTmp = pRecTmp / nPanelDouble;
            this->pRecEachHap_.push_back( pRecTmp / nPanelDouble );

            double pNoRecTmp = 1.0 - pRecTmp;
            this->pNoRec_.push_back( pNoRecTmp );

            this->pRecRec_.push_back ( pRecEachHapTmp * pRecEachHapTmp );
            this->pRecNoRec_.push_back ( pRecEachHapTmp * pNoRecTmp );
            this->pNoRecNoRec_.push_back ( pNoRecTmp * pNoRecTmp );
        }
        this->pRec_.push_back(1.0);
        this->pRecEachHap_.push_back( 1.0 / nPanelDouble );
        this->pNoRec_.push_back(0.0);
        this->pRecRec_.push_back ( 1.0 / nPanelDouble / nPanelDouble );
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


