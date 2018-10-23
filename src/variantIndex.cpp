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

#include <algorithm>  // find
#include <iostream>
#include "exceptions.hpp"
#include "txtReader.hpp"
#include "variantIndex.hpp"

using std::min;
using std::endl;


VariantIndex::VariantIndex() {
    this->init();
}

void VariantIndex::findWhoToBeKept(ExcludeMarker* excludedMarkers) {
    dout << " Starts findWhoToBeKept " << endl;
    assert(this->indexOfContentToBeKept.size() == 0);
    assert(this->indexOfPosToBeKept.size() == 0);

    for (size_t chromI = 0; chromI < this->chrom_.size(); chromI++) {
        dout << "   Going through chrom "<< chrom_[chromI];
        vector < size_t > tmpindexOfPosToBeKept;

        // detemine if something needs to be removed from the current chrom.
        vector<string>::iterator chromIt = find(excludedMarkers->chrom_.begin(),
                                                excludedMarkers->chrom_.end(),
                                                this->chrom_[chromI]);

        size_t hapIndex = indexOfChromStarts_[chromI];
        size_t chromIndexInExclude = std::distance(
                                        excludedMarkers->chrom_.begin(),
                                        chromIt);
        for (size_t posI = 0; posI < this->position_[chromI].size(); posI++) {
            if (chromIt == excludedMarkers->chrom_.end()) {
                indexOfContentToBeKept.push_back(hapIndex);
                tmpindexOfPosToBeKept.push_back(posI);
            } else if (
              std::find(excludedMarkers->position_[chromIndexInExclude].begin(),
                        excludedMarkers->position_[chromIndexInExclude].end(),
                        this->position_[chromI][posI]) ==
                        excludedMarkers->position_[chromIndexInExclude].end()) {
                indexOfContentToBeKept.push_back(hapIndex);
                tmpindexOfPosToBeKept.push_back(posI);
            }
            hapIndex++;
        }
        indexOfPosToBeKept.push_back(tmpindexOfPosToBeKept);

        dout << " keeping " << tmpindexOfPosToBeKept.size() << endl;
    }
    assert(indexOfPosToBeKept.size() == this->chrom_.size());

    dout << indexOfContentToBeKept.size() << " sites need to be Kept " << endl;
}


void VariantIndex::findAndKeepMarkers(ExcludeMarker* excludedMarkers) {
    this->setDoneGetIndexOfChromStarts(false);
    dout << " findAndKeepMarkers called" <<endl;
    this->findWhoToBeKept(excludedMarkers);
    this->removePositions();
    this->getIndexOfChromStarts();
    this->removeMarkers();
}


void VariantIndex::findWhoToBeKeptGivenIndex(
         const vector <size_t> & givenIndex) {
    dout << " Starts findWhoToBeKeptGivenIndex " << endl;
    assert(this->indexOfContentToBeKept.size() == 0);
    indexOfContentToBeKept = vector <size_t> (givenIndex.begin(),
                                              givenIndex.end());
    //assert(this->indexOfPosToBeKept.size() == 0);
    vector <string> oldChrom = vector <string> (chrom_.begin(), chrom_.end());
    this->chrom_.clear();

    vector < vector < int > > oldposition = this->position_;
    this->position_.clear();

    for (size_t chromI = 0; chromI < oldChrom.size(); chromI++) {
        dout << "   Going through chrom "<< oldChrom[chromI] << endl;
        size_t hapIndex = indexOfChromStarts_[chromI];
        vector <int> newTrimmedPos;
        for (size_t posI = 0; posI < oldposition[chromI].size(); posI++) {
            if (std::find(givenIndex.begin(),givenIndex.end(), hapIndex)
                    != givenIndex.end()){
                if (newTrimmedPos.size() == 0) {
                    this->chrom_.push_back(oldChrom[chromI]);
                }
                newTrimmedPos.push_back(oldposition[chromI][posI]);
            }

            hapIndex++;
        }
        this->position_.push_back(newTrimmedPos);
    }

    //assert(indexOfPosToBeKept.size() == this->chrom_.size());

    dout << indexOfContentToBeKept.size() << " sites need to be Kept " << endl;
}


void VariantIndex::findAndKeepMarkersGivenIndex(
        const vector <size_t> & givenIndex) {
    this->setDoneGetIndexOfChromStarts(false);
    dout << " findAndKeepMarkersGivenIndex called " << givenIndex.size() <<endl;
    this->findWhoToBeKeptGivenIndex(givenIndex);
    //this->removePositions();
    this->getIndexOfChromStarts();
    this->removeMarkers();
}


void VariantIndex::removePositions() {
    assert(this->keptPosition_.size() == (size_t)0);
    for (size_t chromI = 0; chromI < this->chrom_.size(); chromI++) {
        vector <int> tmpKeptPosition_;
        for (size_t i = 0; i < this->indexOfPosToBeKept[chromI].size(); i++) {
            tmpKeptPosition_.push_back(
                this->position_[chromI][this->indexOfPosToBeKept[chromI][i]]);
        }
        this->keptPosition_.push_back(tmpKeptPosition_);
    }
    this->position_.clear();
    this->position_ = this->keptPosition_;
    this->keptPosition_.clear();
}


void VariantIndex::getIndexOfChromStarts() {
    assert(this->doneGetIndexOfChromStarts_ == false);
    this->indexOfChromStarts_.clear();
    assert(indexOfChromStarts_.size() == 0);
    this->indexOfChromStarts_.push_back((size_t)0);
    for (size_t tmpChrom = 0;
            indexOfChromStarts_.size() < this->chrom_.size(); tmpChrom++ ) {
        indexOfChromStarts_.push_back(
            indexOfChromStarts_.back()+this->position_[tmpChrom].size());
    }
    assert(indexOfChromStarts_.size() == this->chrom_.size());
    this->setDoneGetIndexOfChromStarts(true);
}


void VariantIndex::init() {
    this->setDoneGetIndexOfChromStarts(false);
}

void VariantIndex::removeMarkers() { throw VirtualFunctionShouldNotBeCalled();}


void VariantIndex::checkSortedPositions(string fileName) {
    for (size_t chromI = 0; chromI < this->chrom_.size(); chromI++) {
        int previousPosition_ = 0;
        for (auto const &value : this->position_[chromI]) {
            if (value < previousPosition_) {
                throw PositionUnsorted(fileName);
            }
            previousPosition_ = value;
        }
    }
}
