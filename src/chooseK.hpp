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

#include <cassert>
#include <vector>
using std::vector;

#ifndef CHOOSEK
#define CHOOSEK

class ChooseK {
  friend class DEploidIO;
 public:
    void appendProportions(const vector <double> &p);
    vector <double> chosenP();

 private:
    ChooseK() {this->haveChosenP_ = false;}
    ~ChooseK() {}
    size_t chosenK() const {
        assert(this->haveChosenP_ == true);
        return this->max_at_+1; }

    void findKmode();
    vector < vector<double> > proportions_;
    vector < size_t > ks;
    void gatherKs();
    size_t max_at_;
    bool haveChosenP_;
    double confidence_;
    double confidence() const { return this->confidence_; }
};

#endif
