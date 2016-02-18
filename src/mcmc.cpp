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

#include "mcmc.hpp"
#include "utility.hpp"
#include <math.h>       /* ceil */
#include <random>
#include <iomanip>      // std::setw
#include "updateHap.hpp"

McmcMachinery::McmcMachinery(Input* input, Panel *panel, McmcSample *mcmcSample,
                       size_t nSample, size_t McmcMachineryRate, size_t seed ){ // initialiseMCMCmachinery
    this->seed_ = (unsigned)(time(0));
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
    this->panel_ = panel;
    this->mcmcSample_ = mcmcSample;
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

    this->initializeExpectedWsaf(); // This requires currentHap_ and currentProp_
    this->currentLLks_ = calcLLKs( this->input_->refCount, this->input_->altCount, this->currentExpectedWsaf_ );
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
    return sumOfVec(tmp);
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
    double tmpSum = sumOfVec(tmpExpTitre);

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

    ofstream hap_file( "tmp.hap", ios::out | ios::app | ios::binary );
    for ( size_t i = 0; i < this->currentHap_.size(); i++ ){
        for ( size_t ii = 0; ii < this->kStrain_; ii++){
            hap_file << this->currentHap_[i][ii];
            hap_file << ((ii < (this->kStrain_-1)) ? "\t" : "\n") ;
        }
    }
    hap_file.close();
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
        for ( size_t k = 0; k < kStrain_; k++){
            expectedWsaf[i] += currentHap_[i][k] * proportion[k];
        }
    }
    return expectedWsaf;
}


void McmcMachinery::recordMcmcMachinery(){
    dout << "Record mcmc sample " <<endl;
    this->mcmcSample_->proportion.push_back(this->currentProp_);

    ofstream prop_file( "tmp.prop", ios::out | ios::app | ios::binary );
    for ( size_t i = 0; i < this->kStrain_; i++){
        prop_file << (10) << mcmcSample_->proportion.back()[i];
        prop_file << ((i < (this->kStrain_-1)) ? "\t" : "\n") ;
    }
    prop_file.close();

    this->mcmcSample_->sumLLKs.push_back(sumOfVec(this->currentLLks_));
    ofstream llk_file( "tmp.llk", ios::out | ios::app | ios::binary );
    llk_file << this->mcmcSample_->sumLLKs.back() << endl;
    llk_file.close();
}


void McmcMachinery::updateProportion(){
    if ( this->kStrain_ < 2 ) return;

    // calculate dt
    vector <double> tmpTitre = calcTmpTitre();
    vector <double> tmpProp = titre2prop(tmpTitre);
    if ( minOfVec(tmpProp) < 0 || maxOfVec(tmpProp) > 1 ) return;

    vector <double> tmpExpecedWsaf = calcExpectedWsaf(tmpProp);
    vector <double> tmpLLKs = calcLLKs (this->input_->refCount, this->input_->altCount, tmpExpecedWsaf);
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
    this->currentProp_ = tmpProp;
}

double McmcMachinery::deltaLLKs ( vector <double> &newLLKs ){
    //double tmp = 0.0;
    //for ( size_t i = 0; i < newLLKs.size(); i++){
        //tmp += ( newLLKs[i] - this->currentLLks_[i] );
    //}
    //return tmp;
    vector <double> tmpdiff = vecDiff ( newLLKs,  this->currentLLks_);
    return sumOfVec(tmpdiff);
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
    UpdateSingleHap updating( this->input_->refCount,
                              this->input_->altCount,
                              this->currentExpectedWsaf_,
                              this->currentProp_, this->currentHap_, this->rg_, this->panel_);

    dout << "update Single Hap "<<endl;
    for ( size_t i = 0 ; i < this->nLoci_; i++ ){
        this->currentHap_[i][updating.strainIndex_] = updating.hap_[i];
    }
    this->currentLLks_ = updating.newLLK;
}



void McmcMachinery::updatePairHaps(){
    UpdatePairHap updating( this->input_->refCount,
                            this->input_->altCount,
                            this->currentExpectedWsaf_,
                            this->currentProp_, this->currentHap_, this->rg_, this->panel_);
    dout << "update Pair Hap "<<endl;

    //for ( size_t i = 0 ; i < this->nLoci_; i++ ){
        //this->currentHap_[i][updating.strainIndex1_] = updating.hap1_[i];
        //this->currentHap_[i][updating.strainIndex2_] = updating.hap2_[i];
    //}
    //this->currentLLks_ = updating.newLLK;
}
