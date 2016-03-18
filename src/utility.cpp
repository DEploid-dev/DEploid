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

#include "logbeta.h"
#include "utility.hpp"
#include <iterator>     // std::distance
#include <algorithm> // find


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
        //if ( expectedWsaf[i] < 0 || expectedWsaf[i] > 1){
            //cout << "i = "<<i <<", "<< expectedWsaf[i] << endl;
        //}
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
    //cout << "unadjustedWsaf = " << unadjustedWsaf<<endl;
    double adjustedWsaf = unadjustedWsaf+err*(1-2*unadjustedWsaf);
    //cout << "adjustedWsaf = " << adjustedWsaf<<endl;
    //cout << "alt+adjustedWsaf*fac = "<< alt+adjustedWsaf*fac<<endl;
    //cout << "ref+(1-adjustedWsaf)*fac = "<< ref+(1-adjustedWsaf)*fac<<endl;

    //cout << "adjustedWsaf*fac = "<< adjustedWsaf*fac << endl;
    //cout << "(1-adjustedWsaf)*fac = " << (1-adjustedWsaf)*fac << endl;
    double llk = Maths::Special::Gamma::logBeta(alt+adjustedWsaf*fac, ref+(1-adjustedWsaf)*fac) - Maths::Special::Gamma::logBeta(adjustedWsaf*fac,(1-adjustedWsaf)*fac);
    //cout << "llk = "<<llk<<endl;
    //double llk = log(boost::math::beta(alt+adjustedWsaf*fac, ref+(1-adjustedWsaf)*fac)) - log(boost::math::beta(adjustedWsaf*fac,(1-adjustedWsaf)*fac));
//    llk<-lgamma(fac*f.samp+cov.alt)+lgamma(fac*(1-f.samp)+cov.ref)-lgamma(fac*f.samp)-lgamma(fac*(1-f.samp));
    //double llk = lgamma(fac*adjustedWsaf+alt)+lgamma(fac*(1-adjustedWsaf)+ref)-lgamma(fac*adjustedWsaf)-lgamma(fac*(1-adjustedWsaf));
    return llk;
}


//vector <size_t> sampleNoReplace( MersenneTwister* rg, vector <double> & proportion, size_t nSample ){
    //vector <size_t> indexReturn;
    //assert( indexReturn.size() == 0 );
    //vector <double> tmpDist(proportion) ;
    //vector <size_t> tmpIndex;
    //for ( size_t i = 0; i < proportion.size(); i++ ){
        //tmpIndex.push_back(i);
    //}
    //for ( size_t nSampleRemaining = nSample; nSampleRemaining > 0; nSampleRemaining-- ){
        //// Compute cdf of tmpDist
        //vector <double> tmpCdf = computeCdf(tmpDist);
        //double u = rg->sample();
        //size_t i = 0;
        //for ( ; i < tmpCdf.size() ; i++){
            //if ( u < tmpCdf[i] ){
                //indexReturn.push_back(tmpIndex[i]);
                //break;
            //}
        //}
        //// Reduce tmpDist and tmpIndex
        //tmpDist.erase(tmpDist.begin()+i);
        //(void)normalizeBySum(tmpDist);
        //tmpIndex.erase(tmpIndex.begin()+i);
    //}
    //return indexReturn;
//}


size_t sampleIndexGivenProp ( MersenneTwister* rg, vector <double> proportion ){
    #ifndef NDEBUG
        double tmpMax = maxOfVec(proportion);
        vector<double>::iterator maxIt = std::find ( proportion.begin(),  proportion.end(), tmpMax);
        return std::distance(proportion.begin(), maxIt);
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
