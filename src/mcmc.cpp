#include "mcmc.hpp"
#include <math.h>       /* ceil */


McmcSample::McmcSample( size_t nSample, size_t mcmcSampleRate, size_t kStrain, size_t seed ){
    this->kStrain_ = kStrain;
    this->seed_ = seed;
    this->rg_ = new MersenneTwister(this->seed_);
    //this->currentMcmcIteration_ = 0;

    // Initialization
    this->initializeProp();
    this->initializeTitre();

    this->burnIn_ = 0.5;
    this->calcMaxIteration(nSample, mcmcSampleRate);
}


McmcSample::~McmcSample(){
    delete rg_;
}


void McmcSample::calcMaxIteration( size_t nSample, size_t mcmcSampleRate ){
    this->mcmcSampleRate_ = mcmcSampleRate;
    this->maxIteration_ = (size_t)ceil((double)nSample*(double)mcmcSampleRate/(1.0-this->burnIn_));
    this->mcmcThresh_ = (size_t)ceil((double)nSample*(double)mcmcSampleRate*this->burnIn_/(1.0-this->burnIn_));
}


void McmcSample::runMcmcChain( ){
    // Run mcmc
    //size_t mcmcSampleRate = 5;

    for ( this->currentMcmcIteration_ = 0 ; currentMcmcIteration_ < this->maxIteration_ ; currentMcmcIteration_++){
        this->sampleMcmcEvent();
    }
}


void McmcSample::sampleMcmcEvent( ){
    this->eventInt_ = this->rg_->sampleInt(3);
    if ( this->eventInt_ == 0 ){
        this->updateProportion();
    } else if ( this->eventInt_ == 1 ){
        this->updateSingleHap();
    } else if ( this->eventInt_ == 2 ){
        this->updatePairHaps();
    } else {
        dout << "eventInt = " << this->eventInt_ << endl;
        throw(" should never reach here!!!");
    }
    if ( currentMcmcIteration_ > this->mcmcThresh_ && currentMcmcIteration_ % this->mcmcSampleRate_ == 0){
        this->recordMcmcSample();
    }
}


void McmcSample::initializeProp( ){
    //initial.propotion[1:initial.k]<-fun.sim_prop(initial.k, TRUE);
    //initial.propotion<-initial.propotion/sum(initial.propotion);
}

void McmcSample::initializeTitre(){

}


void McmcSample::recordMcmcSample(){
    dout << "Record mcmc sample " <<endl;
}

void McmcSample::updateProportion(){
    dout << "update Proportion "<<endl;
}

void McmcSample::updateSingleHap(){
    dout << "update Single Hap "<<endl;
}

void McmcSample::updatePairHaps(){
    dout << "update Pair Hap "<<endl;
}
