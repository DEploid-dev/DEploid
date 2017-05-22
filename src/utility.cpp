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

#include "logbeta.h"
#include "utility.hpp"
#include <iterator>     // std::distance
#include <algorithm> // find
#include "loggammasum.h" // which includes log_gamma.h
#include "gamma.h"



double normal_pdf(double x, double m, double s) { // See http://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-(0.5) * a * a);
}


double min_value( vector <double> x){
    assert( x.size() > 0 );
    auto tmpMaxIt = std::min_element(std::begin(x), std::end(x));
    return *tmpMaxIt;
}


double max_value( vector <double> x){
    assert( x.size() > 0 );
    auto tmpMaxIt = std::max_element(std::begin(x), std::end(x));
    return *tmpMaxIt;
}


vector <double> computeCdf ( vector <double> & dist ){
    vector <double> cdf;
    double cumsum = 0;
    for (double p : dist){
        cumsum += p;
        cdf.push_back( cumsum );
    }
    assert( cdf.size() == dist.size() );
    //assert( cdf.back() == 1.0 );
    return cdf;
}


double sumOfMat( vector <vector <double> > & matrix ){
    double tmp = 0.0;
    for (auto const& array: matrix){
        for (auto const& value: array){
            tmp += value;
        }
    }
    return tmp;
}


void normalizeBySum ( vector <double> & array ){
    double sumOfArray = sumOfVec(array);
    for( vector<double>::iterator it = array.begin(); it != array.end(); ++it) {
        *it /= sumOfArray;
    }
}


void normalizeByMax ( vector <double> & array ){
    double maxOfArray = max_value(array);
    for( vector<double>::iterator it = array.begin(); it != array.end(); ++it) {
        *it /= maxOfArray;
    }
}


void normalizeBySumMat ( vector <vector <double> > & matrix ){
    double tmpsum = sumOfMat(matrix);
    for( size_t i = 0; i < matrix.size(); i++ ){
        for( vector<double>::iterator it = matrix[i].begin(); it != matrix[i].end(); ++it) {
            *it /= tmpsum;
        }
    }
}


vector <double> calcLLKs( vector <double> &refCount, vector <double> &altCount, vector <double> &expectedWsaf, size_t firstIndex, size_t length, double fac, double err){
    assert ( expectedWsaf.size() == length );
    vector <double> tmpLLKs (length, 0.0);
    size_t index = firstIndex;
    for ( size_t i = 0; i < length; i++ ){
        assert (expectedWsaf[i] >= 0);
        //assert (expectedWsaf[i] <= 1);
        tmpLLKs[i] = calcLLK( refCount[index],
                              altCount[index],
                              expectedWsaf[i],
                              err,
                              fac);
        index++;
    }
    return tmpLLKs;
}


double calcLLK( double ref, double alt, double unadjustedWsaf, double err, double fac ) {
    double adjustedWsaf = unadjustedWsaf+err*(1-2*unadjustedWsaf);
    double llk = Maths::Special::Gamma::logBeta(alt+adjustedWsaf*fac, ref+(1-adjustedWsaf)*fac) - Maths::Special::Gamma::logBeta(adjustedWsaf*fac,(1-adjustedWsaf)*fac);
    return llk;
}


size_t sampleIndexGivenProp ( RandomGenerator* rg, vector <double> proportion ){
    #ifndef NDEBUG
        auto biggest = std::max_element(std::begin(proportion), std::end(proportion));
        return std::distance(proportion.begin(), biggest);
    #else
        vector <size_t> tmpIndex;
        for ( size_t i = 0; i < proportion.size(); i++ ){
            tmpIndex.push_back(i);
        }
        vector <double> tmpCdf = computeCdf(proportion);

        double u = rg->sample();
        size_t i = 0;
        for ( ; i < tmpCdf.size() ; i++){
            if ( u < tmpCdf[i] ){
                break;
            }
        }
        return i;
    #endif
}


vector <double> reshapeMatToVec ( vector < vector <double> > &Mat ){
    vector <double> tmp;
    for (auto const& array: Mat){
        for (auto const& value: array){
            tmp.push_back(value);
        }
    }
    return tmp;
}


double betaPdf(double x, double a, double b){
    assert(x>=0 && x<=1);
    assert(a>=0);
    assert(b>=0);
    double p = Maths::Special::Gamma::gamma(a + b) / (Maths::Special::Gamma::gamma(a) * Maths::Special::Gamma::gamma(b));
    double q = pow(1 - x, b - 1) * pow(x, a - 1);
    return p * q;
}


double logBetaPdf(double x, double a, double b){
    assert(x>=0 && x<=1);
    assert(a>=0);
    assert(b>=0);
    double ret = Maths::Special::Gamma::log_gamma(a+b) -
                 Maths::Special::Gamma::log_gamma(a) -
                 Maths::Special::Gamma::log_gamma(b) +
                 (b-1) * log(1-x) + (a-1) * log(x);
    return ret;
}


double binomialPdf(int s, int n, double p){
    assert(p>=0 && p<=1);
    double ret=n_choose_k(n, s);
    ret *= pow(p, (double)s);
    ret *= pow((1.0-p), (double)(n-s));
    return ret;
}


//double betaDistConst( double a , double b){
    //double ret = Maths::Special::Gamma::gamma(a + b) / (Maths::Special::Gamma::gamma(a) * Maths::Special::Gamma::gamma(b));
    //return ret;
//}


double rBeta(double alpha, double beta, RandomGenerator* rg){
    double mxAt = (alpha-1.0)/(alpha+beta-2.0);
    double mx = betaPdf(mxAt, alpha, beta);
    double y, u;
    do {
        u = rg->sample();
        y = rg->sample();
    } while ( u > (betaPdf(y, alpha, beta)/mx));

    return y;
}
