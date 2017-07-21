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

#include <fstream>
#include <iostream>
#include <iterator>     // std::distance
#include "exceptions.hpp"
#include "txtReader.hpp"

void TxtReader::readFromFileBase(const char inchar[]){
    fileName = string (inchar);
    tmpChromInex_ = -1;

    ifstream in_file(inchar);
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
                field_end = min ( min ( tmp_line.find(',',field_start),
                                        tmp_line.find('\t',field_start) ),
                                  tmp_line.find('\n', field_start) );

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
        throw InvalidInputFile(fileName);

    }
    in_file.close();

    this->position_.push_back( this->tmpPosition_ );

    this->nLoci_ = this->content_.size();
    this->nInfoLines_ = this->content_.back().size();

    if ( this->nInfoLines_ == 1 ){
        this->reshapeContentToInfo();
    }

    this->getIndexOfChromStarts();
    assert ( tmpChromInex_ > -1 );
    assert ( chrom_.size() == position_.size() );
    assert(this->doneGetIndexOfChromStarts_ == true);
}


void TxtReader::extractChrom( string & tmp_str ){
    if ( tmpChromInex_ >= 0 ){
        if ( tmp_str != this->chrom_.back() ){
            tmpChromInex_++;
            // save current positions
            this->position_.push_back(this->tmpPosition_);

            // start new chrom
            this->tmpPosition_.clear();
            this->chrom_.push_back(tmp_str);
        }
    } else {
        tmpChromInex_++;
        assert (this->chrom_.size() == 0);
        this->chrom_.push_back( tmp_str );
        assert ( this->tmpPosition_.size() == 0 );
        assert ( this->position_.size() == 0);
    }
}


void TxtReader::extractPOS( string & tmp_str ){
    int ret;
    try {
        ret = stoi(tmp_str.c_str(), NULL);
    } catch ( const std::exception &e ){
        throw BadConversion(tmp_str, fileName);
    }
    this->tmpPosition_.push_back(ret);
}


void TxtReader::reshapeContentToInfo(){
    assert ( this->info_.size() == 0 );
    for ( size_t i = 0; i < this->content_.size(); i++){
        this->info_.push_back( this->content_[i][0] );
    }
}


void TxtReader::removeMarkers ( ){
    assert( this->keptContent_.size() == (size_t)0 );
    for ( auto const &value: this->indexOfContentToBeKept){
        this->keptContent_.push_back(this->content_[value] );
    }

    this->content_.clear();
    assert(this->content_.size() == (size_t)0);
    this->content_ = this->keptContent_;
    this->keptContent_.clear();

    if ( this->nInfoLines_ == 1 ){
        this->info_.clear();
        this->reshapeContentToInfo();
    }
    this->nLoci_ = this->content_.size();
}
