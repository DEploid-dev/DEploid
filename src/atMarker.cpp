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

#include <algorithm> // find
#include <fstream>
#include <iostream>
#include <iterator>     // std::distance

#include "exceptions.hpp"
#include "atMarker.hpp"


void AtMarker::extractChrom( string & tmp_str ){
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


void AtMarker::extractPOS( string & tmp_str ){
    this->tmpPosition_.push_back(strtod(tmp_str.c_str(), NULL));
}



AtMarker::AtMarker(const char inchar[]){
    tmpChromInex_ = -1;

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
        throw InvalidInputFile( string (inchar) );

    }
    in_file.close();

    this->position_.push_back( this->tmpPosition_ );

    this->nLoci_ = this->content_.size();
    this->nInfoLines_ = this->content_.back().size();

    if ( this->nInfoLines_ == 1 ){
        this->reshapeContentToInfo();
    }

    assert ( tmpChromInex_ > -1 );
    assert ( chrom_.size() == position_.size() );
    this->getIndexOfChromStarts();
}


void AtMarker::reshapeContentToInfo(){
    assert ( this->info_.size() == 0 );
    for ( size_t i = 0; i < this->nLoci_; i++){
        this->info_.push_back( this->content_[i][0] );
    }
}

void AtMarker::getIndexOfChromStarts(){
    this->indexOfChromStarts_.clear();
    assert( indexOfChromStarts_.size() == 0 );
    indexOfChromStarts_.push_back( (size_t) 0);
    for ( size_t tmpChrom = 0 ; indexOfChromStarts_.size() < this->chrom_.size(); tmpChrom++ ){
        indexOfChromStarts_.push_back(indexOfChromStarts_.back()+this->position_[tmpChrom].size());
    }
    assert( indexOfChromStarts_.size() == this->chrom_.size() );
}


void InputMarker::removeMarkers( ExcludeMarker* excludedMarkers ){

    vector < size_t > indexOfHapRemovals;
    vector < vector < size_t > > indexOfPosRemovals;

    for ( size_t chromI = 0; chromI < this->chrom_.size(); chromI++){
        vector < size_t > tmpIndexOfPosRemovals;

        // detemine if something needs to be removed from the current chrom.
        vector<string>::iterator chromIt;

        chromIt = find ( excludedMarkers->chrom_.begin(),  excludedMarkers->chrom_.end(), this->chrom_[chromI]);
        if ( chromIt == excludedMarkers->chrom_.end() ) {
            dout << " Remove markers skipping " << chrom_[chromI] << endl;
            indexOfPosRemovals.push_back(tmpIndexOfPosRemovals);
            continue;
        }

        size_t hapIndex = indexOfChromStarts_[ std::distance(excludedMarkers->chrom_.begin(), chromIt) ];

        for ( size_t posI = 0; posI < this->position_[chromI].size(); posI++){
            //double currentPos = this->position_[chromI][posI];
            if ( std::find(excludedMarkers->position_[std::distance(excludedMarkers->chrom_.begin(), chromIt)].begin(), excludedMarkers->position_[std::distance(excludedMarkers->chrom_.begin(), chromIt)].end(), this->position_[chromI][posI]) != excludedMarkers->position_[std::distance(excludedMarkers->chrom_.begin(), chromIt)].end() ){
                indexOfHapRemovals.push_back(hapIndex);
                tmpIndexOfPosRemovals.push_back(posI);
            }
            hapIndex++;
        }
        reverse( tmpIndexOfPosRemovals.begin(), tmpIndexOfPosRemovals.end() );
        indexOfPosRemovals.push_back(tmpIndexOfPosRemovals);
    }

    reverse( indexOfHapRemovals.begin(), indexOfHapRemovals.end() );

    for ( size_t chromI = 0; chromI < this->chrom_.size(); chromI++){
        for ( auto const &value: indexOfPosRemovals[chromI]){// std::vector<size_t>::iterator it=indexOfHapRemovals.begin(); it!=indexOfHapRemovals.end(); it++ ){
            this->position_[chromI].erase ( this->position_[chromI].begin() + value );
        }
    }

    for ( auto const &value: indexOfHapRemovals){// std::vector<size_t>::iterator it=indexOfHapRemovals.begin(); it!=indexOfHapRemovals.end(); it++ ){
        this->content_.erase(this->content_.begin() + value );
    }

    if ( this->nInfoLines_ == 1 ){
        for ( auto const &value: indexOfHapRemovals){// std::vector<size_t>::iterator it=indexOfHapRemovals.begin(); it!=indexOfHapRemovals.end(); it++ ){
            this->info_.erase(this->info_.begin() + value );
        }
    }
    this->getIndexOfChromStarts();
    this->nLoci_ = this->content_.size();
}
