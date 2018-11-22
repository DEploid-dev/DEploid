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

#include <cstddef>
#include <iostream>
#include "global.hpp"
#include "chooseK.hpp"

using std::cout;
using std::endl;

void ChooseK::appendProportions(const vector <double> &p) {
    this->proportions_.push_back(vector<double>(p.begin(), p.end()));
}


void ChooseK::gatherKs() {
    assert(this->ks.size() == 0);
    for (auto const& vec : this->proportions_) {
        size_t k = 0;
        for (auto const& v : vec) {
            if (v > 0.01) {
                k++;
            }
        }
        ks.push_back(k);
    }
    assert(this->ks.size() == this->proportions_.size());
}


void ChooseK::findKmode() {
    this->gatherKs();

    vector <size_t> uniqueK;
    for (size_t i = 0; i < this->proportions_[0].size(); i++) {
        uniqueK.push_back(i+1);
    }
    size_t maxCount = 0;
    this->max_at_ = 0;
    vector <size_t> kCount(uniqueK.size(), 0);
    for (auto const& k : ks) {
        kCount[k-1]++;
        dout << k << " " << kCount[k-1] << endl;
        if (kCount[k-1] > maxCount) {
            max_at_ = k-1;
            maxCount = kCount[max_at_];
        }
    }

    dout << "k most frequent: " << max_at_+1
         << " with " << maxCount << "count" << endl;
    this->confidence_ = static_cast<double>(maxCount) /
                        static_cast<double>(ks.size());
}


vector <double> ChooseK::chosenP() {
    this->findKmode();
    size_t idx = max_at_;
    for (size_t i = idx; i < this->ks.size(); i++) {
        if (ks[i] == ks[idx]) {
            idx = i;
        }
    }
    this->haveChosenP_ = true;
    return this->proportions_[idx];
}

