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


#include "txtReader.hpp"
#include <algorithm> // find
#include <iostream>

vector < size_t > VariantIndex::findWhoToBeRemoved( ExcludeMarker* excludedMarkers ){

    /* Index of content/info will be removed */
    vector < size_t > indexOfContentToBeRemoved;
    /* Index of positions entry to be removed, this will have the same size as this->chrom_, */
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

        size_t hapIndex = indexOfChromStarts_[ chromI ];
        size_t chromIndexInExclude = std::distance(excludedMarkers->chrom_.begin(), chromIt);

        for ( size_t posI = 0; posI < this->position_[chromI].size(); posI++){
            if ( std::find(excludedMarkers->position_[chromIndexInExclude].begin(), excludedMarkers->position_[chromIndexInExclude].end(), this->position_[chromI][posI]) != excludedMarkers->position_[chromIndexInExclude].end() ){
                indexOfContentToBeRemoved.push_back(hapIndex);
                tmpIndexOfPosRemovals.push_back(posI);
            }
            hapIndex++;
        }
        reverse( tmpIndexOfPosRemovals.begin(), tmpIndexOfPosRemovals.end() );
        indexOfPosRemovals.push_back(tmpIndexOfPosRemovals);
    }
    assert ( indexOfPosRemovals.size() == this->chrom_.size() );

    reverse( indexOfContentToBeRemoved.begin(), indexOfContentToBeRemoved.end() );

    for ( size_t chromI = 0; chromI < this->chrom_.size(); chromI++){
        for ( auto const &value: indexOfPosRemovals[chromI]){
            this->position_[chromI].erase ( this->position_[chromI].begin() + value );
        }
    }

    return indexOfContentToBeRemoved;

}

