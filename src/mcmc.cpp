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
    this->PROP_SCALE = 40.0;

    std_generator_ = new std::default_random_engine(this->seed_);
    initialTitre_normal_distribution_ = new std::normal_distribution<double>(MN_LOG_TITRE, SD_LOG_TITRE);
    deltaX_normal_distribution_ = new std::normal_distribution<double>(MN_LOG_TITRE, 1.0/PROP_SCALE);

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
    delete initialTitre_normal_distribution_;
    delete deltaX_normal_distribution_;
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
    this->currentLogPriorTitre_ = this->calcLogPriorTitre(this->currentTitre_);

    this->initializeHap();
    this->initializeProp();

    this->initializeExpectedWsaf();

    this->currentLLks_ = calcLLKs( this->currentExpectedWsaf_ );
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
    this->currentExpectedWsaf_ = vector <double> (this->input_->plaf.size(), 0.0);
    assert( this->currentExpectedWsaf_.size() == this->nLoci_ );
    this->calcExpectedWsaf( this->currentProp_ );
}


void McmcMachinery::initializellk(){
    assert( this->currentLLks_.size() == (size_t)0);
    //for ( size_t i = 0; i < this->input_->plaf.size(); i++ ){
        //this->currentLLks_.push_back(0.0);
    //}
    this->currentLLks_ = vector <double> (this->nLoci_, 0.0);
    assert( this->currentLLks_.size() == this->nLoci_);
}


void McmcMachinery::initializeProp( ){
    assert( this->currentProp_.size() == (size_t)0 );
    this->currentProp_ = this->titre2prop( this->currentTitre_ );
}


double McmcMachinery::calcLogPriorTitre( vector <double> &tmpTitre ){
    //sum(dnorm(titre, MN_LOG_TITRE, SD_LOG_TITRE, log=TRUE));
    vector <double> tmp;
    for ( auto const& value: tmpTitre ){
        tmp.push_back( log(normal_pdf(value, MN_LOG_TITRE, SD_LOG_TITRE)) );
    }
    return sum(tmp);
}


void McmcMachinery::initializeTitre(){
    /*   titre<-rnorm(initial.k, MN_LOG_TITRE, SD_LOG_TITRE); */
    assert( currentTitre_.size() == 0);
    currentTitre_ = vector <double> (this->kStrain_, 0.0);
    for ( size_t k = 0; k < this->kStrain_; k++){
        double tmp = (*this->initialTitre_normal_distribution_)((*std_generator_));
        this->currentTitre_[k] = tmp;
    }
    assert( currentTitre_.size() == this->kStrain_);
}

vector <double> McmcMachinery::titre2prop(vector <double> & tmpTitre){
    vector <double> tmpExpTitre;
    for ( auto const& value: tmpTitre ){
        tmpExpTitre.push_back( exp(value) );
    }
    double tmpSum = sum(tmpExpTitre);

    vector <double> tmpProp;
    for ( auto const& value: tmpExpTitre ){
        tmpProp.push_back( value/tmpSum );
    }
    return tmpProp;
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


vector <double> McmcMachinery::calcExpectedWsaf( vector <double> &proportion ){
    vector <double> expectedWsaf (this->nLoci_, 0.0);
    for ( size_t i = 0; i < currentHap_.size(); i++ ){
        assert( kStrain_ == currentHap_[i].size() );
        double tmp = 0.0;
        for ( size_t k = 0; k < kStrain_; k++){
            tmp += currentHap_[i][k] * proportion[k];
        }
        //cout << tmp << endl;
        expectedWsaf[i] = tmp;
    }
    return expectedWsaf;
}


vector <double> McmcMachinery::calcLLKs( vector <double> &expectedWsaf ){
    vector <double> tmpLLKs (this->nLoci_, 0.0);
    for ( size_t i = 0; i < this->currentLLks_.size(); i++ ){
        tmpLLKs[i] = calcLLK( this->input_->refCount[i],
                                         this->input_->altCount[i],
                                         expectedWsaf[i]);
    }
    return tmpLLKs;
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
    if ( this->kStrain_ < 2 ) return;

    // calculate dt
    vector <double> tmpTitre = calcTmpTitre();
    vector <double> tmpProp = titre2prop(tmpTitre);
    if ( minOfVec(tmpProp) < 0 || maxOfVec(tmpProp) > 1 ) return;

    vector <double> tmpExpecedWsaf = calcExpectedWsaf(tmpProp);
    vector <double> tmpLLKs = calcLLKs (tmpExpecedWsaf);
    double diffLLKs = this->deltaLLKs(tmpLLKs);
    double tmpLogPriorTitre = calcLogPriorTitre( tmpTitre );
    double priorPropRatio = exp(tmpLogPriorTitre - this->currentLogPriorTitre_ );
    double hastingsRatio = 1.0;

    //runif(1)<prior.prop.ratio*hastings.ratio*exp(del.llk))
    if ( this->rg_->sample() > priorPropRatio*hastingsRatio*exp(diffLLKs) ) return;

    dout << "update Proportion "<<endl;
    this->currentExpectedWsaf_ = tmpExpecedWsaf;
    this->currentLLks_ = tmpLLKs;
    this->currentLogPriorTitre_ = tmpLogPriorTitre;
    this->currentTitre_ = tmpTitre;

}

double McmcMachinery::deltaLLKs ( vector <double> &newLLKs ){
    double tmp = 0.0;
    for ( size_t i = 0; i < newLLKs.size(); i++){
        tmp += ( newLLKs[i] - this->currentLLks_[i] );
    }
    return tmp;
}

vector <double> McmcMachinery::calcTmpTitre(){
    vector <double> tmpTitre;
    for ( size_t k = 0; k < this->kStrain_; k++){
        double dt = (*this->deltaX_normal_distribution_)((*std_generator_));
        tmpTitre.push_back( currentTitre_[k] + dt );
    }
    return tmpTitre;
}




void McmcMachinery::updateSingleHap(){
    dout << "update Single Hap "<<endl;
}


void McmcMachinery::updatePairHaps(){
    dout << "update Pair Hap "<<endl;
}
