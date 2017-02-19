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
#include "utility.hpp"

#ifndef IBD
#define IBD

using namespace std;

vector < vector<int> > unique( vector < vector<int> > &mat );
bool twoVectorsAreSame(vector<int> vec1, vector<int> vec2);


vector<int> convertIntToBinary(int x, size_t len);
vector < vector <int> > enumerateBinaryMatrixOfK(size_t k);

struct OutOfVectorSize : std::exception{

  explicit OutOfVectorSize(){ }
  virtual ~OutOfVectorSize() throw() {}
  virtual const char* what () const noexcept {
      return string("Out of vector size!").c_str();
  }
};

int nchoose2(int n);

// The IBDconfiguration is used for index, which should be non-negative, use int, any thing below zero should throw.
class IBDconfiguration{
#ifdef UNITTEST
  friend class TestIBDconfig;
  friend class TestHprior;
#endif
  friend class Hprior;
    IBDconfiguration();
    ~IBDconfiguration();

    void buildIBDconfiguration(size_t k = 5);
    size_t kStrain_;
    void setKstrain(const size_t setTo) {this->kStrain_ = setTo;}
    size_t kStrain() const {return this->kStrain_;}

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
  friend class McmcMachinery;
    Hprior();
    ~Hprior();

    void buildHprior(size_t kStrain, vector <double> &plaf);

    IBDconfiguration ibdConfig;
    size_t kStrain_;
    void setKstrain(const size_t setTo) {this->kStrain_ = setTo;}
    size_t kStrain() const {return this->kStrain_;}

    size_t nLoci_;
    void setnLoci(const size_t setTo) {this->nLoci_ = setTo;}
    size_t nLoci() const {return this->nLoci_;}

    vector <double> plaf_;
    vector < vector <double> > priorProb; // size: nState x nLoci
    vector < vector <double> > priorProbTrans; // size: nLoci x nState

    vector <size_t> stateIdx; // size: nState

    vector <vector <int> > hSet; // size: nState x kStrain
    size_t nState_;
    size_t nState() const {return this->nState_;}
};




#endif
