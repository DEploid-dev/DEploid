#include "updateHap.hpp"

UpdateHap::UpdateHap( vector <double> &refCount,
                                  vector <double> &altCount,
                                  vector <double> &expectedWsaf,
                                  vector <double> &proportion,
                                  vector < vector <double> > &haplotypes, MersenneTwister* rg){
    this->missCopyProb = 0.01;

    this->kStrain_ = proportion.size();
    this->nLoci_ = expectedWsaf.size();
    this->rg_ = rg;

    this->findUpdatingStrain();
    this->calcExpectedWsaf( expectedWsaf, proportion, haplotypes);
    calcHapLLKs(refCount, altCount);
}


UpdateSingleHap::UpdateSingleHap( vector <double> &refCount,
                                  vector <double> &altCount,
                                  vector <double> &expectedWsaf,
                                  vector <double> &proportion,
                                  vector < vector <double> > &haplotypes, MersenneTwister* rg):
                UpdateHap(refCount, altCount, expectedWsaf, proportion, haplotypes, rg){}


void UpdateSingleHap::findUpdatingStrain(){
    this->strainIndex_ = (size_t)this->rg_->sampleInt(this->kStrain_);
    //vector <size_t> strainIndex = sampleNoReplace( (size_t)1, this->currentProp_, this->rg_ );
    //assert( strainIndex.size() == 1);
    //size_t ws = strainIndex[0];

}

void UpdateSingleHap::calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes ){
    //expected.WSAF.0 <- bundle$expected.WSAF - (bundle$prop[ws] * bundle$h[,ws]);
    this->expectedWsaf0_ = expectedWsaf;
    for ( size_t i = 0; i < expectedWsaf0_.size(); i++ ){
        expectedWsaf0_[i] -= proportion[strainIndex_] * haplotypes[i][strainIndex_];
    }

    //expected.WSAF.1 <- expected.WSAF.0 + bundle$prop[ws] ;
    this->expectedWsaf1_ = expectedWsaf0_;
    for ( size_t i = 0; i < expectedWsaf1_.size(); i++ ){
        expectedWsaf1_[i] += haplotypes[i][strainIndex_];
    }
}


void UpdateSingleHap::buildEmission(){

    //vector <double> t1omu = llk0_ , vector <double> (llklog(1-miss.copy.rate)
    //t2omu = logemiss[,2]+log(1-miss.copy.rate)
    //t1u = logemiss[,1]+log(miss.copy.rate)
    //t2u = logemiss[,2]+log(miss.copy.rate)
    //tmax = apply(cbind(t1omu, t2omu, t1u, t2u), 1, max)
    //emiss = cbind( exp(t1omu-tmax) + exp(t2u-tmax),
                   //exp(t2omu-tmax) + exp(t1u-tmax))

}


void UpdateSingleHap::calcHapLLKs( vector <double> &refCount,
                                vector <double> &altCount){
    llk0_ = calcLLKs( refCount, altCount, expectedWsaf0_ );
    llk1_ = calcLLKs( refCount, altCount, expectedWsaf1_ );
}
