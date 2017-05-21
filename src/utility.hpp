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

#include <vector>
#include <iostream>
#include <algorithm>    /* min_element, max_element */

#include "mersenne_twister.hpp"
#include "global.h"

#ifndef UTILITY
#define UTILITY

using namespace std;

template <typename T>
vector <T> vecDiff ( vector<T> &vecA, vector<T> &vecB ){
    assert(vecA.size() == vecB.size());
    vector <T> difference (vecA.size(), (T)0);
    for ( size_t i = 0; i < vecA.size(); i++ ){
        difference[i] = vecA[i] - vecB[i];
    }
    return difference;
}


template <typename T>
vector <T> vecSum ( vector<T> &vecA, vector<T> &vecB ){
    assert(vecA.size() == vecB.size());
    vector <T> tmpSum (vecA.size(), (T)0);
    for ( size_t i = 0; i < vecA.size(); i++ ){
        tmpSum[i] = vecA[i] + vecB[i];
    }
    return tmpSum;
}


template <typename T>
vector <T> vecProd ( vector<T> &vecA, vector<T> &vecB ){
    assert(vecA.size() == vecB.size());
    vector <T> tmpProd (vecA.size(), (T)0);
    for ( size_t i = 0; i < vecA.size(); i++ ){
        tmpProd[i] = vecA[i] * vecB[i];
    }
    return tmpProd;
}


template <typename T>
T sumOfVec( vector <T>& array ){
    T tmp = 0;
    for (auto const& value: array){
        tmp += value;
    }
    return tmp;
}

/*! \brief Compute factorial of a \return double a! */
template < class T > T factorial ( T a ){
    if (a > 1) return (a * factorial (a-1));
    else       return (1);
}

/*! \brief Compute a permutations of n \return double */
template < class T > T n_permu_a ( T n, T a ){
    if   ( a > 1 ) return (n*n_permu_a(n-1,a-1));
    else if (a==1) return (n);
    else           return (1);
}

/*! \brief Compute n choose k \return double */
template < class T > T n_choose_k ( T n, T k ){
    if ( k < ( n/2 ) ) return (n_choose_k(n,n-k));
    else               return (n_permu_a(n,k)/factorial(k));
}


double normal_pdf(double x, double m, double s);
double min_value( vector <double> x);
double max_value( vector <double> x);
vector <double> computeCdf ( vector <double> & dist );
double sumOfMat( vector <vector <double> > & matrix );
void normalizeBySum ( vector <double> & array );
void normalizeByMax ( vector <double> & array );
void normalizeBySumMat ( vector <vector <double> > & matrix );
vector <double> calcLLKs( vector <double> &refCount, vector <double> &altCount, vector <double> &expectedWsaf, size_t firstIndex, size_t length, double fac, double err = 0.01);
double calcLLK( double ref, double alt, double unadjustedWsaf, double err, double fac) ;
size_t sampleIndexGivenProp ( RandomGenerator* rg, vector <double> proportion );
vector <double> reshapeMatToVec ( vector < vector <double> > &Mat );
double betaPdf(double x, double a, double b);
double logBetaPdf(double x, double a, double b);
double binomialPdf(int s, int n, double p);
//double betaDistConst( double a , double b);
double rBeta(double alpha, double beta, RandomGenerator* rg);
#endif
