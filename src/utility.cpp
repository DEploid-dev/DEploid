#include "utility.hpp"

vector <size_t> sampleNoReplace( vector <double> proportion, MersenneTwister* rg, size_t nSample ){
    vector <size_t> indexReturn;
    assert( indexReturn.size() == 0 );
    vector <double> tmpDist(proportion) ;
    vector <size_t> tmpIndex;
    for ( size_t i = 0; i < proportion.size(); i++ ){
        tmpIndex.push_back(i);
    }
    for ( size_t nSampleRemaining = nSample; nSampleRemaining > 0; nSampleRemaining-- ){
        // Compute cdf of tmpDist
        vector <double> tmpCdf = computeCdf(tmpDist);
        double u = rg->sample();
        size_t i = 0;
        for ( ; i < tmpCdf.size() ; i++){
            if ( u < tmpCdf[i] ){
                indexReturn.push_back(tmpIndex[i]);
                break;
            }
        }
        // Reduce tmpDist and tmpIndex
        tmpDist.erase(tmpDist.begin()+i);
        (void)normalizeBySum(tmpDist);
        tmpIndex.erase(tmpIndex.begin()+i);
    }
    return indexReturn;
}


size_t sampleIndexGivenProp ( vector <double> proportion, MersenneTwister* rg ){
    vector <size_t> strainIndex = sampleNoReplace( proportion, rg);
    assert( strainIndex.size() == 1);
    return strainIndex[0];
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


void normalizeBySum ( vector <double> & array ){
    double sumOfArray = sumOfVec(array);
    for( vector<double>::iterator it = array.begin(); it != array.end(); ++it) {
        *it /= sumOfArray;
    }
}


vector <double> calcLLKs( vector <double> &refCount, vector <double> &altCount, vector <double> &expectedWsaf ){
    vector <double> tmpLLKs (expectedWsaf.size(), 0.0);
    for ( size_t i = 0; i < tmpLLKs.size(); i++ ){
        tmpLLKs[i] = calcLLK( refCount[i],
                              altCount[i],
                              expectedWsaf[i]);
    }
    return tmpLLKs;
}


double calcLLK( double ref, double alt, double unadjustedWsaf, double err, double fac ) {
    double adjustedWsaf = unadjustedWsaf+err*(1-2*unadjustedWsaf);
    //double llk = logBeta(alt+adjustedWsaf*fac, ref+(1-adjustedWsaf)*fac)-logBeta(adjustedWsaf*fac,(1-adjustedWsaf)*fac);
//    llk<-lgamma(fac*f.samp+cov.alt)+lgamma(fac*(1-f.samp)+cov.ref)-lgamma(fac*f.samp)-lgamma(fac*(1-f.samp));
    double llk = lgamma(fac*adjustedWsaf+alt)+lgamma(fac*(1-adjustedWsaf)+ref)-lgamma(fac*adjustedWsaf)-lgamma(fac*(1-adjustedWsaf));
    return llk;
}
