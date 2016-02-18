#include <vector>
#include <iostream>
#include "mersenne_twister.hpp"
//#include "logbeta.h"

#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

#ifndef UTILITY
#define UTILITY

using namespace std;


template <typename T> // See http://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
T normal_pdf(T x, T m, T s)
{
    static const T inv_sqrt_2pi = 0.3989422804014327;
    T a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-T(0.5) * a * a);
}


template <typename T>
T minOfVec( vector <T> x){
    assert( x.size() > 0 );
    T tmpMin = x[0];
    for ( auto const& value: x ){
        tmpMin = ( value < tmpMin ) ? value : tmpMin;
    }
    return tmpMin;
}


template <typename T>
T maxOfVec( vector <T> x){
    assert( x.size() > 0 );
    T tmpMax = x[0];
    for ( auto const& value: x ){
        tmpMax = ( value > tmpMax ) ? value : tmpMax;
    }
    return tmpMax;
}


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


vector <double> computeCdf ( vector <double> & dist );
double sumOfVec( vector <double>& array );
void normalizeBySum ( vector <double> & array );
vector <size_t> sampleNoReplace( vector <double> proportion, MersenneTwister* rg, size_t nSample = 1);
size_t sampleIndexGivenProp ( vector <double> proportion, MersenneTwister* rg );
vector <double> calcLLKs( vector <double> &refCount, vector <double> &altCount, vector <double> &expectedWsaf );
double calcLLK( double ref, double alt, double unadjustedWsaf, double err = 0.01, double fac=100 ) ;

#endif
