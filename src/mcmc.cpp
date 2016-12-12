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

#include "mcmc.hpp"
#include "utility.hpp"
#include <math.h>       /* ceil */
#include <random>
#include "updateHap.hpp"
#include <stdio.h>

McmcSample::McmcSample(){};
McmcSample::~McmcSample(){};

McmcMachinery::McmcMachinery(DEploidIO* dEploidIO, Panel *panel, McmcSample *mcmcSample, RandomGenerator* rg_ ){ // initialiseMCMCmachinery

    this->dEploidIO_ = dEploidIO;
    this->panel_ = panel;
    this->mcmcSample_ = mcmcSample;
    this->seed_ = rg_->seed();

    //this->hapRg_ = new MersenneTwister(this->seed_);
    this->hapRg_ = rg_ ;
    this->mcmcEventRg_ = this->hapRg_;
    this->propRg_ = this->hapRg_;
    this->initialHapRg_ = this->hapRg_;
    //this->mcmcEventRg_ = new MersenneTwister(this->seed_);
    //this->propRg_  = new MersenneTwister(this->seed_);
    //this->initialHapRg_ = new MersenneTwister(this->seed_);

    this->calcMaxIteration( dEploidIO_->nMcmcSample_ , dEploidIO_->mcmcMachineryRate_, dEploidIO_->mcmcBurn_ );

    this->MN_LOG_TITRE = 0.0;
    this->SD_LOG_TITRE = 3.0;
    this->PROP_SCALE = 40.0;

    stdNorm_ = new StandNormalRandomSample(this->seed_);

    this->kStrain_ = this->dEploidIO_->kStrain_;
    this->nLoci_ = this->dEploidIO_->plaf_.size();
    this->initializeMcmcChain( );

}


McmcMachinery::~McmcMachinery(){
    if ( this->stdNorm_ ){
        delete stdNorm_;
    }
}


void McmcMachinery::calcMaxIteration( size_t nSample, size_t McmcMachineryRate, double burnIn ){
    this->burnIn_ = burnIn;
    this->McmcMachineryRate_ = McmcMachineryRate;
    this->maxIteration_ = (size_t)ceil((double)nSample*(double)McmcMachineryRate/(1.0-this->burnIn_))+1;
    this->mcmcThresh_ = (size_t)ceil((double)nSample*(double)McmcMachineryRate*this->burnIn_/(1.0-this->burnIn_));
}


void McmcMachinery::initializeMcmcChain( ){
    // Initialization
    dout << "###########################################"<< endl;
    dout << "#            Initialization               #"<< endl;
    dout << "###########################################"<< endl;

    this->initializeTitre();
    this->currentLogPriorTitre_ = this->calcLogPriorTitre(this->currentTitre_);
    this->initializeHap();
    this->initializeProp();
    this->initializeExpectedWsaf(); // This requires currentHap_ and currentProp_
    this->currentLLks_ = calcLLKs( this->dEploidIO_->refCount_, this->dEploidIO_->altCount_, this->currentExpectedWsaf_ , 0, this->currentExpectedWsaf_.size());

    assert (doutProp());
    assert (doutLLK());
    dout << "###########################################"<< endl;
    dout << "#        Initialization finished          #"<< endl;
    dout << "###########################################"<< endl;
}


void McmcMachinery::initializeHap(){
    assert( currentHap_.size() == 0);
    for ( size_t i = 0; i < this->dEploidIO_->plaf_.size(); i++ ){
        double currentPlaf = this->dEploidIO_->plaf_[i];
        vector <double> tmpVec;
        for ( size_t k = 0; k < this->kStrain_; k++){
            tmpVec.push_back( this->rBernoulli(currentPlaf) );
        }
        this->currentHap_.push_back(tmpVec);
    }
}


double McmcMachinery::rBernoulli(double p){
    double u = this->initialHapRg_->sample();
    return ( u < p ) ? 1.0 : 0.0;
}


void McmcMachinery::initializeExpectedWsaf(){
    assert( this->currentExpectedWsaf_.size() == 0);
    this->currentExpectedWsaf_ = this->calcExpectedWsaf( this->currentProp_ );
    assert( this->currentExpectedWsaf_.size() == this->nLoci_ );
}


void McmcMachinery::initializellk(){
    assert( this->currentLLks_.size() == (size_t)0);
    this->currentLLks_ = vector <double> (this->nLoci_, 0.0);
    assert( this->currentLLks_.size() == this->nLoci_);
}


void McmcMachinery::initializeProp( ){
    assert( this->currentProp_.size() == (size_t)0 );
    this->currentProp_ = ( this->dEploidIO_ -> initialPropWasGiven()) ?
                          this->dEploidIO_ ->initialProp:
                          this->titre2prop( this->currentTitre_ );
    if ( this->dEploidIO_ -> initialPropWasGiven() ){
        this->currentTitre_.clear();
        for ( size_t i = 0; i <  this->dEploidIO_ ->initialProp.size(); i++ ) {
            this->currentTitre_.push_back( log(this->dEploidIO_ ->initialProp[i]));
        }
    }
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

    if ( this->dEploidIO_->doUpdateProp() ){
        for ( size_t k = 0; k < this->kStrain_; k++){
            this->currentTitre_[k] = this->initialTitreNormalVariable();
        }
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
        assert (tmpProp.back() > 0);
        assert (tmpProp.back() < 1);
    }
    return tmpProp;
}


void McmcMachinery::runMcmcChain( bool showProgress ){
    for ( this->currentMcmcIteration_ = 0 ; currentMcmcIteration_ < this->maxIteration_ ; currentMcmcIteration_++){
        dout << endl;
        dout << "MCMC iteration: " << this->currentMcmcIteration_ << endl;
        if ( this->currentMcmcIteration_ > 0 && this->currentMcmcIteration_%100 == 0 && showProgress ){
            clog << "\r" << " MCMC step" << setw(4) << int(currentMcmcIteration_ * 100 / this->maxIteration_) << "% completed."<<flush;
        }
        this->sampleMcmcEvent();
    }
    clog << "\r" << " MCMC step" << setw(4) << 100 << "% completed."<<endl;
    this->mcmcSample_->hap = this->currentHap_;

    this->writeLastFwdProb();
}


void McmcMachinery::sampleMcmcEvent( ){
    this->recordingMcmcBool_ = ( currentMcmcIteration_ > this->mcmcThresh_ && currentMcmcIteration_ % this->McmcMachineryRate_ == 0 );

    this->eventInt_ = this->mcmcEventRg_->sampleInt(3);
    if ( (this->eventInt_ == 0) && (this->dEploidIO_->doUpdateProp() == true) ){
        this->updateProportion();
    } else if ( (this->eventInt_ == 1) && (this->dEploidIO_->doUpdateSingle() == true) ){
        this->updateSingleHap();
    } else if ( (this->eventInt_ == 2) && (this->dEploidIO_->doUpdatePair() == true) ){
        this->updatePairHaps();
    }

    assert(doutLLK());

    if ( this->recordingMcmcBool_ ){
        this->recordMcmcMachinery();
    }
}


vector <double> McmcMachinery::calcExpectedWsaf( vector <double> &proportion ){
    //assert ( sumOfVec(proportion) == 1.0); // this fails ...
    vector <double> expectedWsaf (this->nLoci_, 0.0);
    for ( size_t i = 0; i < currentHap_.size(); i++ ){
        assert( kStrain_ == currentHap_[i].size() );
        for ( size_t k = 0; k < kStrain_; k++){
            expectedWsaf[i] += currentHap_[i][k] * proportion[k];
        }
        assert ( expectedWsaf[i] >= 0 );
        //assert ( expectedWsaf[i] <= 1.0 );
    }
    return expectedWsaf;
}


void McmcMachinery::recordMcmcMachinery(){
    dout << "***Record mcmc sample " <<endl;
    this->mcmcSample_->proportion.push_back(this->currentProp_);
    this->mcmcSample_->sumLLKs.push_back(sumOfVec(this->currentLLks_));
    this->mcmcSample_->moves.push_back(this->eventInt_);
}


void McmcMachinery::updateProportion(){
    dout << " Attempt of update proportion";

    if ( this->kStrain_ < 2 ) {
        dout << "(failed)" << endl;
        return;
    }

    // calculate dt
    vector <double> tmpTitre = calcTmpTitre();
    vector <double> tmpProp = titre2prop(tmpTitre);

    if ( min_value(tmpProp) < 0 || max_value(tmpProp) > 1 ) {
        dout << "(failed)" << endl;
        return;
    }

    vector <double> tmpExpecedWsaf = calcExpectedWsaf(tmpProp);
    vector <double> tmpLLKs = calcLLKs (this->dEploidIO_->refCount_, this->dEploidIO_->altCount_, tmpExpecedWsaf, 0, tmpExpecedWsaf.size());
    double diffLLKs = this->deltaLLKs(tmpLLKs);
    double tmpLogPriorTitre = calcLogPriorTitre( tmpTitre );
    double priorPropRatio = exp(tmpLogPriorTitre - this->currentLogPriorTitre_ );
    double hastingsRatio = 1.0;

    //runif(1)<prior.prop.ratio*hastings.ratio*exp(del.llk))
    if ( this->propRg_->sample() > priorPropRatio*hastingsRatio*exp(diffLLKs) ) {
        dout << "(failed)" << endl;
        return;
    }

    dout << "(successed) " << endl;

    this->currentExpectedWsaf_ = tmpExpecedWsaf;
    this->currentLLks_ = tmpLLKs;
    this->currentLogPriorTitre_ = tmpLogPriorTitre;
    this->currentTitre_ = tmpTitre;
    this->currentProp_ = tmpProp;

    assert (doutProp());
}


double McmcMachinery::deltaLLKs ( vector <double> &newLLKs ){
    vector <double> tmpdiff = vecDiff ( newLLKs,  this->currentLLks_);
    return sumOfVec(tmpdiff);
}


vector <double> McmcMachinery::calcTmpTitre(){
    vector <double> tmpTitre;
    for ( size_t k = 0; k < this->kStrain_; k++){
        double dt = this->deltaXnormalVariable();
        tmpTitre.push_back( currentTitre_[k] + dt );
    }
    return tmpTitre;
}


void McmcMachinery::updateSingleHap(){
    dout << " Update Single Hap "<<endl;
    this->findUpdatingStrainSingle();

    for ( size_t chromi = 0 ; chromi < this->dEploidIO_->indexOfChromStarts_.size(); chromi++ ){
        size_t start = this->dEploidIO_->indexOfChromStarts_[chromi];
        size_t length = this->dEploidIO_->position_[chromi].size();
        dout << "   Update Chrom with index " << chromi << ", starts at "<< start << ", with " << length << " sites" << endl;
        UpdateSingleHap updating( this->dEploidIO_->refCount_,
                                  this->dEploidIO_->altCount_,
                                  this->dEploidIO_->plaf_,
                                  this->currentExpectedWsaf_,
                                  this->currentProp_, this->currentHap_, this->hapRg_,
                                  start, length,
                                  this->panel_, this->dEploidIO_->missCopyProb_,
                                  this->strainIndex_);
        updating.core ( this->dEploidIO_->refCount_, this->dEploidIO_->altCount_, this->dEploidIO_->plaf_, this->currentExpectedWsaf_, this->currentProp_, this->currentHap_);

        size_t updateIndex = 0;
        for ( size_t ii = start ; ii < (start+length); ii++ ){
            this->currentHap_[ii][this->strainIndex_] = updating.hap_[updateIndex];
            this->currentLLks_[ii] = updating.newLLK[updateIndex];
            updateIndex++;
        }

        if ( this->dEploidIO_->doExportSwitchMissCopy() ){
            this->dEploidIO_->ofstreamExportTmp.open( this->dEploidIO_->strExportOneSwitchOne.c_str(), ios::out | ios::app | ios::binary );
            for ( size_t i = 0; i < updating.siteOfOneSwitchOne.size(); i++ ){
                this->dEploidIO_->ofstreamExportTmp << this->dEploidIO_->chrom_[chromi] << "\t" << (size_t)this->dEploidIO_->position_[chromi][updating.siteOfOneSwitchOne[i]] << endl;
            }
            this->dEploidIO_->ofstreamExportTmp.close();

            this->dEploidIO_->ofstreamExportTmp.open( this->dEploidIO_->strExportOneMissCopyOne.c_str(), ios::out | ios::app | ios::binary );
            for ( size_t i = 0; i < updating.siteOfOneMissCopyOne.size(); i++ ){
                this->dEploidIO_->ofstreamExportTmp << this->dEploidIO_->chrom_[chromi] << "\t" << (size_t)this->dEploidIO_->position_[chromi][updating.siteOfOneMissCopyOne[i]] << endl;
            }
            this->dEploidIO_->ofstreamExportTmp.close();
        }
    }

    this->currentExpectedWsaf_ = this->calcExpectedWsaf( this->currentProp_ );
}


void McmcMachinery::updatePairHaps(){
    dout << " Update Pair Hap "<<endl;
    this->findUpdatingStrainPair();

    for ( size_t chromi = 0 ; chromi < this->dEploidIO_->indexOfChromStarts_.size(); chromi++ ){
        size_t start = this->dEploidIO_->indexOfChromStarts_[chromi];
        size_t length = this->dEploidIO_->position_[chromi].size();
        dout << "   Update Chrom with index " << chromi << ", starts at "<< start << ", with " << length << " sites" << endl;

        UpdatePairHap updating( this->dEploidIO_->refCount_,
                                this->dEploidIO_->altCount_,
                                this->dEploidIO_->plaf_,
                                this->currentExpectedWsaf_,
                                this->currentProp_, this->currentHap_, this->hapRg_,
                                start, length,
                                this->panel_, this->dEploidIO_->missCopyProb_, this->dEploidIO_->forbidCopyFromSame(),
                                this->strainIndex1_,
                                this->strainIndex2_);
        updating.core ( this->dEploidIO_->refCount_, this->dEploidIO_->altCount_, this->dEploidIO_->plaf_, this->currentExpectedWsaf_, this->currentProp_, this->currentHap_);

        size_t updateIndex = 0;
        for ( size_t ii = start ; ii < (start+length); ii++ ){
            this->currentHap_[ii][this->strainIndex1_] = updating.hap1_[updateIndex];
            this->currentHap_[ii][this->strainIndex2_] = updating.hap2_[updateIndex];
            this->currentLLks_[ii] = updating.newLLK[updateIndex];
            updateIndex++;
        }

        if ( this->dEploidIO_->doExportSwitchMissCopy() ){
            this->dEploidIO_->ofstreamExportTmp.open( this->dEploidIO_->strExportTwoSwitchOne.c_str(), ios::out | ios::app | ios::binary );
            for ( size_t i = 0; i < updating.siteOfTwoSwitchOne.size(); i++ ){
                this->dEploidIO_->ofstreamExportTmp << this->dEploidIO_->chrom_[chromi] << "\t" << (size_t)this->dEploidIO_->position_[chromi][updating.siteOfTwoSwitchOne[i]] << endl;
            }
            this->dEploidIO_->ofstreamExportTmp.close();

            this->dEploidIO_->ofstreamExportTmp.open( this->dEploidIO_->strExportTwoMissCopyOne.c_str(), ios::out | ios::app | ios::binary );
            for ( size_t i = 0; i < updating.siteOfTwoMissCopyOne.size(); i++ ){
                this->dEploidIO_->ofstreamExportTmp << this->dEploidIO_->chrom_[chromi] << "\t" << (size_t)this->dEploidIO_->position_[chromi][updating.siteOfTwoMissCopyOne[i]] << endl;
            }
            this->dEploidIO_->ofstreamExportTmp.close();

            this->dEploidIO_->ofstreamExportTmp.open( this->dEploidIO_->strExportTwoSwitchTwo.c_str(), ios::out | ios::app | ios::binary );
            for ( size_t i = 0; i < updating.siteOfTwoSwitchTwo.size(); i++ ){
                this->dEploidIO_->ofstreamExportTmp << this->dEploidIO_->chrom_[chromi] << "\t" << (size_t)this->dEploidIO_->position_[chromi][updating.siteOfTwoSwitchTwo[i]] << endl;
            }
            this->dEploidIO_->ofstreamExportTmp.close();

            this->dEploidIO_->ofstreamExportTmp.open( this->dEploidIO_->strExportTwoMissCopyTwo.c_str(), ios::out | ios::app | ios::binary );
            for ( size_t i = 0; i < updating.siteOfTwoMissCopyTwo.size(); i++ ){
                this->dEploidIO_->ofstreamExportTmp << this->dEploidIO_->chrom_[chromi] << "\t" << (size_t)this->dEploidIO_->position_[chromi][updating.siteOfTwoMissCopyTwo[i]] << endl;
            }
            this->dEploidIO_->ofstreamExportTmp.close();
        }
    }

    this->currentExpectedWsaf_ = this->calcExpectedWsaf( this->currentProp_ );
}


void McmcMachinery::findUpdatingStrainSingle( ){
    vector <double> eventProb (this->kStrain_, 1);
    (void)normalizeBySum(eventProb);
    this->strainIndex_ = sampleIndexGivenProp ( this->mcmcEventRg_, eventProb );
    dout << "  Updating hap: "<< this->strainIndex_ <<endl;
}


void McmcMachinery::findUpdatingStrainPair( ){
    vector <size_t> strainIndex (2, 0);
    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    double u;

    while (m < 2) {
        u = this->mcmcEventRg_->sample(); // call a uniform(0,1) random number generator
        if ( ( this->kStrain_ - t)*u < 2 - m ) {
            strainIndex[m] = t;
            m++;
        }
        t++;
    }

    this->strainIndex1_ = strainIndex[0];
    this->strainIndex2_ = strainIndex[1];
    assert( strainIndex1_ != strainIndex2_ );
    dout << "  Updating hap: "<< this->strainIndex1_ << " and " << strainIndex2_ <<endl;
}


