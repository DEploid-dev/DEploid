#include "mcmc.hpp"
#include <math.h>       /* ceil */
#include <random>

McmcMachinery::McmcMachinery(Input* input,
                       size_t nSample, size_t McmcMachineryRate, size_t seed ){ // initialiseMCMCmachinery
    this->seed_ = seed;
    this->rg_ = new MersenneTwister(this->seed_);
    this->burnIn_ = 0.5;
    this->calcMaxIteration(nSample, McmcMachineryRate);


    this->MN_LOG_TITRE = 0.0;
    this->SD_LOG_TITRE = 3.0;

    std_generator_ = new std::default_random_engine(this->seed_);
    std_normal_distribution_ = new std::normal_distribution<double>(MN_LOG_TITRE, SD_LOG_TITRE);

    this->input_ = input;
    this->kStrain_ = this->input_->kStrain_;
    this->nLoci_ = this->input_->plaf.size();
    this->initializeMcmcChain( );

    //this->normalizeBySum(this->input_->refCount );
    //for (auto const& value: this->input_->refCount){
        //cout << value << endl;
    //}
}


McmcMachinery::~McmcMachinery(){
    this->rg_->clearFastFunc();
    delete std_generator_;
    delete std_normal_distribution_;
    delete rg_;
}


void McmcMachinery::calcMaxIteration( size_t nSample, size_t McmcMachineryRate ){
    this->McmcMachineryRate_ = McmcMachineryRate;
    this->maxIteration_ = (size_t)ceil((double)nSample*(double)McmcMachineryRate/(1.0-this->burnIn_));
    this->mcmcThresh_ = (size_t)ceil((double)nSample*(double)McmcMachineryRate*this->burnIn_/(1.0-this->burnIn_));
}


void McmcMachinery::initializeMcmcChain( ){
    // Initialization
    this->initializeTitre();
    this->initializeHap();
    this->initializeExpectedWsaf();
    this->calcExpectedWsaf();
    this->initializellk();
    this->calcCurrentLLKs();
}


void McmcMachinery::initializeHap(){
    assert( currentHap_.size() == 0);
    for ( size_t i = 0; i < this->input_->plaf.size(); i++ ){
        double currentPlaf = this->input_->plaf[i];
        vector <double> tmpVec;
        for ( size_t k = 0; k < this->kStrain_; k++){
            tmpVec.push_back( this->rBernoulli(currentPlaf) );
        }
        this->currentHap_.push_back(tmpVec);
    }
}


double McmcMachinery::rBernoulli(double p){
    double u = this->rg_->sample();
    return ( u < p ) ? 1 : 0;
}


void McmcMachinery::initializeExpectedWsaf(){
    assert( this->currentExpectedWsaf_.size() == 0);
    for ( size_t i = 0; i < this->input_->plaf.size(); i++ ){
        this->currentExpectedWsaf_.push_back(0.0);
    }
}


void McmcMachinery::initializellk(){
    assert( this->currentLLks_.size() == 0);
    for ( size_t i = 0; i < this->input_->plaf.size(); i++ ){
        this->currentLLks_.push_back(0.0);
    }
}


void McmcMachinery::initializeProp( ){
    /* fun.titre2prop<-function(titre) {
      //tmp<-exp(titre);
      //return(tmp/sum(tmp));
    }*/
    vector <double> tmp;
    for ( auto const& value: this->currentTitre_ ){
        tmp.push_back( exp(value) );
    }
    double tmpSum = sum(tmp);
    assert ( currentProp_.size() == 0 );
    for ( auto const& value: tmp ){
        currentProp_.push_back( value/tmpSum );
    }
}


void McmcMachinery::calcLogPriorTitre(){
    //sum(dnorm(titre, MN_LOG_TITRE, SD_LOG_TITRE, log=TRUE));
    vector <double> tmp;
    for ( auto const& value: this->currentTitre_ ){
        tmp.push_back( log(normal_pdf(value, MN_LOG_TITRE, SD_LOG_TITRE)) );
    }
    this->currentLogPriorTitre_ = this->sum(tmp);
}


void McmcMachinery::initializeTitre(){
    /*   titre<-rnorm(initial.k, MN_LOG_TITRE, SD_LOG_TITRE); */
    assert( currentTitre_.size() == 0);
    for ( size_t k = 0; k < this->kStrain_; k++){
        double tmp = (*this->std_normal_distribution_)((*std_generator_));
        //cout << tmp << endl;
        currentTitre_.push_back( tmp );
    }

    this->calcLogPriorTitre();
    this->initializeProp();
}


void McmcMachinery::runMcmcChain( ){
    for ( this->currentMcmcIteration_ = 0 ; currentMcmcIteration_ < this->maxIteration_ ; currentMcmcIteration_++){
        this->sampleMcmcEvent();
    }
}


void McmcMachinery::sampleMcmcEvent( ){
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
    if ( currentMcmcIteration_ > this->mcmcThresh_ && currentMcmcIteration_ % this->McmcMachineryRate_ == 0){
        this->recordMcmcMachinery();
    }
}


void McmcMachinery::calcExpectedWsaf(){
    for ( size_t i = 0; i < currentHap_.size(); i++ ){
        assert( kStrain_ == currentHap_[i].size() );
        double tmp = 0.0;
        for ( size_t k = 0; k < kStrain_; k++){
            tmp += currentHap_[i][k] * currentProp_[k];
        }
        currentExpectedWsaf_[i] = tmp;
    }
}


void McmcMachinery::calcCurrentLLKs(){
    for ( size_t i = 0; i < this->currentLLks_.size(); i++ ){
        this->currentLLks_[i] = calcLLK( this->input_->refCount[i],
                                         this->input_->altCount[i],
                                         currentExpectedWsaf_[i]);
        //cout << this->currentLLks_[i] << endl;
    }
}


double McmcMachinery::calcLLK( double ref, double alt, double unadjustedWsaf, double err, double fac ) {
    double adjustedWsaf = unadjustedWsaf+err*(1-2*unadjustedWsaf);
    //double llk = lbeta(alt+adjustedWsaf*fac, ref+(1-adjustedWsaf)*fac)-lbeta(adjustedWsaf*fac,(1-adjustedWsaf)*fac);
    double llk = lgamma(fac*adjustedWsaf+alt)+lgamma(fac*(1-adjustedWsaf)+ref)-lgamma(fac*adjustedWsaf)-lgamma(fac*(1-adjustedWsaf));
    return llk;
}


void McmcMachinery::recordMcmcMachinery(){
    dout << "Record mcmc sample " <<endl;
}


void McmcMachinery::updateProportion(){
    dout << "update Proportion "<<endl;
}


void McmcMachinery::updateSingleHap(){
    dout << "update Single Hap "<<endl;
}


void McmcMachinery::updatePairHaps(){
    dout << "update Pair Hap "<<endl;
}
