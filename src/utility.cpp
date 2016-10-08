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

#include "logbeta.h"
#include "utility.hpp"
#include <iterator>     // std::distance
#include <algorithm> // find


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
    for (size_t i = 0; i < dist.size(); i++){
        cumsum += dist[i];
        cdf.push_back( cumsum );
    }
    assert( cdf.size() == dist.size() );
    //assert( cdf.back() == 1.0 );
    return cdf;
}


double sumOfVec( vector <double>& array ){
    double tmp = 0.0;
    for (auto const& value: array){
        tmp += value;
    }
    return tmp;
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

void normalizeBySumMat ( vector <vector <double> > & matrix ){
    double tmpsum = sumOfMat(matrix);
    for( size_t i = 0; i < matrix.size(); i++ ){
        for( vector<double>::iterator it = matrix[i].begin(); it != matrix[i].end(); ++it) {
            *it /= tmpsum;
        }
    }
}



vector <double> calcLLKs( vector <double> &refCount, vector <double> &altCount, vector <double> &expectedWsaf, size_t firstIndex, size_t length ){
    assert ( expectedWsaf.size() == length );
    vector <double> tmpLLKs (length, 0.0);
    size_t index = firstIndex;
    for ( size_t i = 0; i < length; i++ ){
        assert (expectedWsaf[i] >= 0);
        //assert (expectedWsaf[i] <= 1);
        tmpLLKs[i] = calcLLK( refCount[index],
                              altCount[index],
                              expectedWsaf[i]);
        index++;
    }
    return tmpLLKs;
}


double calcLLK( double ref, double alt, double unadjustedWsaf, double err, double fac ) {
    double adjustedWsaf = unadjustedWsaf+err*(1-2*unadjustedWsaf);
    double llk = Maths::Special::Gamma::logBeta(alt+adjustedWsaf*fac, ref+(1-adjustedWsaf)*fac) - Maths::Special::Gamma::logBeta(adjustedWsaf*fac,(1-adjustedWsaf)*fac);
    return llk;
}


size_t sampleIndexGivenProp ( MersenneTwister* rg, vector <double> proportion ){
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
