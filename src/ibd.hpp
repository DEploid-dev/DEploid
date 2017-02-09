/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016, Sha (Joe) Zhu, Jacob Almagro and Prof. Gil McVean
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

#include <vector>
#include <iostream>
#include <exceptions.hpp>


#ifndef IBD
#define IBD

using namespace std;


// The IBDconfiguration is used for index, which should be non-negative, use int, any thing below zero should throw.
class IBDconfiguration{
#ifdef UNITTEST
  friend class TestIBDconfig;
#endif
    IBDconfiguration();
    ~IBDconfiguration();
    IBDconfiguration(int k = 5);
    int kStrain_;
    void setKstrain(const int setTo) {this->kStrain_ = setTo;}
    int kStrain() const {return this->kStrain_;}

    vector < vector<int> > op;
    vector < vector<int> > pairToEmission;
    vector < vector<size_t> > pairList;
    vector < vector<int> > states;
    vector < size_t > effectiveK;

    size_t stateSize() const { return this->states.size(); }
    void enumerateOp();
    void makePairList();
    void makePairToEmission();
    void findUniqueState();
    void findEffectiveK();

    vector <int> makeTmpRow();
    vector <size_t> findWhichIsOne(vector <int> tmpOp);
    bool twoVectorsAreSame(vector<int> vec1, vector<int> vec2);
};


class Hprior{
#ifdef UNITTEST
  friend class TestHprior;
#endif
    Hprior();
    Hprior(IBDconfiguration ibdConfig, size_t nLoci);
    ~Hprior();

    int kStrain_;
    void setKstrain(const int setTo) {this->kStrain_ = setTo;}
    int kStrain() const {return this->kStrain_;}

    size_t nLoci_;
    void setKstrain(const size_t setTo) {this->nLoci_ = setTo;}
    size_t nLoci() const {return this->nLoci_;}

};


#endif
