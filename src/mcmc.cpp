/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
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

#include <stdio.h>
#include <math.h>       // ceil
#include <random>
#include <limits>       // std::numeric_limits< double >::min()
#include <numeric>      // std::accumulate, std::inner_product
#include "global.hpp"     // dout
#include "updateHap.hpp"
#include "mcmc.hpp"
#include "utility.hpp"

McmcSample::McmcSample() {}
McmcSample::~McmcSample() {}

McmcMachinery::McmcMachinery( vector <double> * plaf,
                              vector <double> * refCount,
                              vector <double> * altCount,
                              Panel *panel_ptr,
                              DEploidIO* dEploidIO,
                              string mcmcJob,
                              string jobbrief,
                              McmcSample *mcmcSample,
                              RandomGenerator* rg_,
                              bool useIBD ) {  // initialiseMCMCmachinery
    this->mcmcJob = mcmcJob;
    this->jobbrief = jobbrief;
    this->plaf_ptr_ = plaf;
    this->refCount_ptr_ = refCount;
    this->altCount_ptr_ = altCount;
    this->panel_ = panel_ptr;
    this->dEploidIO_ = dEploidIO;
    //this->panel_ = dEploidIO->panel;
    this->mcmcSample_ = mcmcSample;
    this->seed_ = rg_->seed();

    //this->hapRg_ = new MersenneTwister(this->seed_);
    this->hapRg_ = rg_ ;
    this->mcmcEventRg_ = this->hapRg_;
    this->propRg_ = this->hapRg_;
    this->initialHapRg_ = this->hapRg_;
    //this->ibdRg_ = this->hapRg_;
    //this->mcmcEventRg_ = new MersenneTwister(this->seed_);
    //this->propRg_  = new MersenneTwister(this->seed_);
    //this->initialHapRg_ = new MersenneTwister(this->seed_);
    if (useIBD == true) {
        this->calcMaxIteration(100, 10, 0.5);
        //this->calcMaxIteration(10, 10, 0.5);
    } else {
        this->calcMaxIteration( dEploidIO_->nMcmcSample_ , dEploidIO_->mcmcMachineryRate_, dEploidIO_->mcmcBurn_ );
    }
    this->MN_LOG_TITRE = 0.0;
    this->SD_LOG_TITRE = (useIBD == true) ? this->dEploidIO_->ibdSigma() : this->dEploidIO_->parameterSigma();
    //this->SD_LOG_TITRE = this->dEploidIO_->parameterSigma();
    this->PROP_SCALE = 40.0;

    stdNorm_ = new StandNormalRandomSample(this->seed_);

    this->setKstrain(this->dEploidIO_->kStrain());
    this->setNLoci(this->plaf_ptr_->size());
    this->initializeMcmcChain(useIBD);
}


McmcMachinery::~McmcMachinery() {
    if ( this->stdNorm_ ) {
        delete stdNorm_;
    }
}


void McmcMachinery::calcMaxIteration( size_t nSample, size_t McmcMachineryRate, double burnIn ) {
    this->burnIn_ = burnIn;
    this->McmcMachineryRate_ = McmcMachineryRate;
    this->maxIteration_ = (size_t)ceil((double)nSample*(double)McmcMachineryRate/(1.0-this->burnIn_))+1;
    this->mcmcThresh_ = (size_t)ceil((double)nSample*(double)McmcMachineryRate*this->burnIn_/(1.0-this->burnIn_));
}


void McmcMachinery::initializeMcmcChain(bool useIBD) {
    // Initialization
    dout << "###########################################"<< endl;
    dout << "#            Initialization               #"<< endl;
    dout << "###########################################"<< endl;
    dout << "plaf.size() " << plaf_ptr_->size() <<"   ";
    if (panel_ != NULL) {
        dout << "Panel size = " << this->panel_->truePanelSize() << endl;
    }
    this->initializeTitre();
    this->currentLogPriorTitre_ = this->calcLogPriorTitre(this->currentTitre_);
    this->initializeHap();
    this->initializeProp();
    this->initializeExpectedWsaf(); // This requires currentHap_ and currentProp_
    this->currentLLks_ = calcLLKs( *this->refCount_ptr_, *this->altCount_ptr_, this->currentExpectedWsaf_ , 0, this->currentExpectedWsaf_.size(), this->dEploidIO_->scalingFactor());
    this->acceptUpdate = 0;

    if ( this->dEploidIO_->doAllowInbreeding() == true ) {
        this->initializeUpdateReferencePanel(this->panel_->truePanelSize()+kStrain_-1);
    }

    if ( useIBD == true ) {
        this->ibdInitializeEssentials();
    }

    this->mcmcSample_->siteOfTwoSwitchOne = vector <double> (this->nLoci());
    this->mcmcSample_->siteOfTwoMissCopyOne = vector <double> (this->nLoci());
    this->mcmcSample_->siteOfTwoSwitchTwo = vector <double> (this->nLoci());
    this->mcmcSample_->siteOfTwoMissCopyTwo = vector <double> (this->nLoci());
    this->mcmcSample_->siteOfOneSwitchOne = vector <double> (this->nLoci());
    this->mcmcSample_->siteOfOneMissCopyOne = vector <double> (this->nLoci());

    this->mcmcSample_->currentsiteOfTwoSwitchOne = vector <double> (this->nLoci());
    this->mcmcSample_->currentsiteOfTwoMissCopyOne = vector <double> (this->nLoci());
    this->mcmcSample_->currentsiteOfTwoSwitchTwo = vector <double> (this->nLoci());
    this->mcmcSample_->currentsiteOfTwoMissCopyTwo = vector <double> (this->nLoci());
    this->mcmcSample_->currentsiteOfOneSwitchOne = vector <double> (this->nLoci());
    this->mcmcSample_->currentsiteOfOneMissCopyOne = vector <double> (this->nLoci());

    assert (doutProp());
    assert (doutLLK());
    dout << "###########################################"<< endl;
    dout << "#        Initialization finished          #"<< endl;
    dout << "###########################################"<< endl;
}


void McmcMachinery::initializeHap() {
    assert( currentHap_.size() == 0);
    if ( this->dEploidIO_ -> initialHapWasGiven() ) {
        this->currentHap_ = this->dEploidIO_->initialHap;
        dout << "given initial hap ?" << endl;
    } else {
        for ( size_t i = 0; i < this->plaf_ptr_->size(); i++ ) {
            double currentPlaf = this->plaf_ptr_->at(i);
            vector <double> tmpVec;
            for ( size_t k = 0; k < this->kStrain_; k++) {
                tmpVec.push_back( this->rBernoulli(currentPlaf) );
            }
            this->currentHap_.push_back(tmpVec);
        }
    }
    assert(this->currentHap_.size() == this->plaf_ptr_->size());
}


void McmcMachinery::initializeUpdateReferencePanel(size_t inbreedingPanelSizeSetTo) {
    if ( this->dEploidIO_->doAllowInbreeding() != true ) {
        return;
    }

    this->panel_->initializeUpdatePanel(inbreedingPanelSizeSetTo);
}


void McmcMachinery::updateReferencePanel(size_t inbreedingPanelSizeSetTo, size_t excludedStrain) {
    if ( this->burnIn_ > this->currentMcmcIteration_ ) {
        return;
    }

    //if ( this->dEploidIO_->doAllowInbreeding() != true ) {
        //return;
    //}
    this->panel_->updatePanelWithHaps( inbreedingPanelSizeSetTo, excludedStrain, this->currentHap_);
}


double McmcMachinery::rBernoulli(double p) {
    double u = this->initialHapRg_->sample();
    return ( u < p ) ? 1.0 : 0.0;
}


void McmcMachinery::initializeExpectedWsaf() {
    assert( this->currentExpectedWsaf_.size() == 0);
    this->currentExpectedWsaf_ = this->calcExpectedWsaf( this->currentProp_ );
    assert( this->currentExpectedWsaf_.size() == this->nLoci_ );
    this->cumExpectedWsaf_ = this->currentExpectedWsaf_;
}


void McmcMachinery::initializellk() {
    assert( this->currentLLks_.size() == (size_t)0);
    this->currentLLks_ = vector <double> (this->nLoci_, 0.0);
    assert( this->currentLLks_.size() == this->nLoci_);
}


void McmcMachinery::initializeProp() {
    assert( this->currentProp_.size() == (size_t)0 );
    this->currentProp_ = ( this->dEploidIO_ -> initialPropWasGiven()) ?
                          this->dEploidIO_ ->initialProp:
                          this->titre2prop( this->currentTitre_ );
    if ( this->dEploidIO_ -> initialPropWasGiven() ) {
        this->currentTitre_.clear();
        for ( size_t i = 0; i <  this->dEploidIO_ ->initialProp.size(); i++ ) {
            this->currentTitre_.push_back( log(this->dEploidIO_ ->initialProp[i]));
        }
    }
}


double McmcMachinery::calcLogPriorTitre(const vector <double> &tmpTitre) {
    //sum(dnorm(titre, MN_LOG_TITRE, SD_LOG_TITRE, log=TRUE));
    vector <double> tmp;
    for ( auto const& value: tmpTitre ) {
        tmp.push_back(log(normal_pdf(value, MN_LOG_TITRE, SD_LOG_TITRE)));
    }
    return sumOfVec(tmp);
}


void McmcMachinery::initializeTitre() {
    /*   titre<-rnorm(initial.k, MN_LOG_TITRE, SD_LOG_TITRE); */
    assert( currentTitre_.size() == 0);
    currentTitre_ = vector <double> (this->kStrain_, 0.0);

    if ( this->dEploidIO_->doUpdateProp() ) {
        for ( size_t k = 0; k < this->kStrain_; k++) {
            this->currentTitre_[k] = this->initialTitreNormalVariable();
            //this->currentTitre_[k] = 10.0;
        }
    }
    assert( currentTitre_.size() == this->kStrain_);
}


vector <double> McmcMachinery::titre2prop(const vector <double> & tmpTitre) {
    vector <double> tmpExpTitre;
    for ( auto const& value: tmpTitre ) {
        tmpExpTitre.push_back( exp(value) );
    }
    double tmpSum = sumOfVec(tmpExpTitre);

    vector <double> tmpProp;
    for ( auto const& value: tmpExpTitre ) {
        tmpProp.push_back( value/tmpSum );
        assert (tmpProp.back() > 0);
        assert (tmpProp.back() <= 1);
    }
    return tmpProp;
}




void McmcMachinery::initializePropIBD() {
    //#Initialise titres and convert to proportions
    //this->initializeTitre();
    this->currentProp_ = ( this->dEploidIO_ -> initialPropWasGiven()) ?
                          this->dEploidIO_ ->initialProp:
                          this->titre2prop( this->currentTitre_ );

    //this->currentTitre_ = vector <double> (this->kStrain(), 0.0);
    //this->currentProp_ = ( this->dEploidIO_ -> initialPropWasGiven()) ?
                          //this->dEploidIO_ ->initialProp :
                          //vector <double> (this->kStrain(), 1.0/(double)kStrain());
}


void McmcMachinery::runMcmcChain( bool showProgress, bool useIBD, bool notInR, bool averageP) {
    for ( this->currentMcmcIteration_ = 0 ; currentMcmcIteration_ < this->maxIteration_ ; currentMcmcIteration_++) {
        dout << endl;
        dout << "MCMC iteration: " << this->currentMcmcIteration_ << endl;
        if ( this->currentMcmcIteration_ > 0 && this->currentMcmcIteration_%30 == 0 && showProgress ) {
            #ifndef RBUILD
                clog << "\r" << " MCMC step" << setw(4) << int(currentMcmcIteration_ * 100 / this->maxIteration_) << "% completed ("<<this->mcmcJob<<")"<<flush;
            #endif
        }
        this->sampleMcmcEvent(useIBD);

        //printArray(this->currentProp_);
    }

    #ifndef RBUILD
        clog << "\r" << " MCMC step" << setw(4) << 100 << "% completed ("<<this->mcmcJob<<")"<<endl;
    #endif
    printArray(this->currentProp_);

    this->mcmcSample_->hap = this->currentHap_;

    this->writeLastFwdProb(useIBD);
    this->dEploidIO_->finalProp.clear();
    this->dEploidIO_->finalProp = this->mcmcSample_->proportion.back();

    if (averageP){
      double remainingP = 1.0;
      this->dEploidIO_->finalProp.clear();
      for (size_t i = 0; i < (this->mcmcSample_->proportion.back().size()-1); i ++){
        double strainp = 0.0;
        size_t np = this->mcmcSample_->proportion.size();
        for (size_t j = 0; j < np; j++){
          strainp += this->mcmcSample_->proportion[j][i];
        }
        strainp = strainp/(np*1.0);
        remainingP -= strainp;
        this->dEploidIO_->finalProp.push_back(strainp);
      }
      this->dEploidIO_->finalProp.push_back(remainingP);
      // this->dEploidIO_->finalProp = this->mcmcSample_->proportion.back();
    }

    for (size_t atSiteI = 0; atSiteI < nLoci(); atSiteI++ ) {
        this->mcmcSample_->siteOfTwoSwitchOne[atSiteI] /= (double)this->maxIteration_;
        this->mcmcSample_->siteOfTwoMissCopyOne[atSiteI] /= (double)this->maxIteration_;
        this->mcmcSample_->siteOfTwoSwitchTwo[atSiteI] /= (double)this->maxIteration_;
        this->mcmcSample_->siteOfTwoMissCopyTwo[atSiteI] /= (double)this->maxIteration_;
        this->mcmcSample_->siteOfOneSwitchOne[atSiteI] /= (double)this->maxIteration_;
        this->mcmcSample_->siteOfOneMissCopyOne[atSiteI] /= (double)this->maxIteration_;
    }

    if ( notInR & ((jobbrief == "lassoK") | (jobbrief == "ibd") | (jobbrief == "classic")) ) { // notInPython
        this->dEploidIO_->writeMcmcRelated(this->mcmcSample_, jobbrief, useIBD);
    }

    if ( useIBD == true ) {
        for (size_t atSiteI = 0; atSiteI < nLoci(); atSiteI++ ) {
            this->ibdPath.IBDpathChangeAt[atSiteI] /= (double)this->maxIteration_;
        }
        //vector < vector <double> > reshapedProbs = this->reshapeFm(hprior.stateIdx);
        //this->dEploidIO_->ibdProbsHeader = getIBDprobsHeader();
        //this->dEploidIO_->ibdProbsIntegrated = getIBDprobsIntegrated(reshapedProbs);
        //this->dEploidIO_->writeIBDpostProb(reshapedProbs, this->dEploidIO_->ibdProbsHeader);
        //clog << "Proportion update acceptance rate: "<<acceptUpdate / (this->kStrain()*1.0*this->maxIteration_)<<endl;
        this->dEploidIO_->initialProp = averageProportion();
        this->dEploidIO_->setInitialPropWasGiven(true);
        this->dEploidIO_->setDoUpdateProp(false);
        //this->dEploidIO_->initialHap = this->mcmcSample_->hap;
        //this->dEploidIO_->setInitialHapWasGiven(true);
    }
    clog << "Proportion update acceptance rate: "<<acceptUpdate / (this->kStrain()*1.0*this->maxIteration_)<<endl;

    this->computeDiagnostics();
    dout << "###########################################"<< endl;
    dout << "#            MCMC RUN finished            #"<< endl;
    dout << "###########################################"<< endl;
}




void McmcMachinery::computeDiagnostics() {
    // clog << "Proportion update acceptance rate: "<<acceptUpdate / (this->kStrain()*1.0*this->maxIteration_)<<endl;
    this->dEploidIO_->setacceptRatio(acceptUpdate / (1.0*this->maxIteration_));

    // average cumulate expectedWSAF
    for ( size_t i = 0; i < this->cumExpectedWsaf_.size(); i++) {
        //cout << "cumExpectedWsaf_ i = "<<i <<" " << this->cumExpectedWsaf_[i] << " " << this->dEploidIO_->nMcmcSample_<<endl;
        this->cumExpectedWsaf_[i] /= static_cast<double>(this->dEploidIO_->nMcmcSample_);
        if (this->cumExpectedWsaf_[i]>1){
            this->cumExpectedWsaf_[i] = 1;
        }
    }
    vector <double> tmpLLKs1 = calcLLKs (*this->refCount_ptr_, *this->altCount_ptr_, this->cumExpectedWsaf_, 0, this->cumExpectedWsaf_.size(), this->dEploidIO_->scalingFactor());
    this->dEploidIO_->setmeanThetallks( sumOfVec(tmpLLKs1) );

    vector <double> wsaf_vec;
    for ( size_t i = 0; i < nLoci(); i++) {
        double wsaf = this->altCount_ptr_->at(i) / (this->refCount_ptr_->at(i) + this->altCount_ptr_->at(i) + 0.00000000000001);
        double adjustedWsaf = wsaf*(1-0.01) + (1-wsaf)*0.01;
        wsaf_vec.push_back(adjustedWsaf);
        //llkOfData.push_back( logBetaPdf(adjustedWsaf, this->llkSurf[i][0], this->llkSurf[i][1]));
    }
    vector <double> tmpLLKs = calcLLKs (*this->refCount_ptr_, *this->altCount_ptr_, wsaf_vec, 0, wsaf_vec.size(), this->dEploidIO_->scalingFactor());
    this->dEploidIO_->setmaxLLKs( sumOfVec(tmpLLKs) );

    double sum = std::accumulate(this->mcmcSample_->sumLLKs.begin(), this->mcmcSample_->sumLLKs.end(), 0.0);
    double mean = sum / this->mcmcSample_->sumLLKs.size();
    double sq_sum = std::inner_product(this->mcmcSample_->sumLLKs.begin(), this->mcmcSample_->sumLLKs.end(), this->mcmcSample_->sumLLKs.begin(), 0.0);
    double varLLKs = sq_sum / this->mcmcSample_->sumLLKs.size() - mean * mean;
    double stdev = std::sqrt(varLLKs);
    this->dEploidIO_->setmeanllks(mean);
    this->dEploidIO_->setstdvllks(stdev);

    double dicByVar = (-2*mean) + 4*varLLKs/2;
    this->dEploidIO_->setdicByVar(dicByVar);
     //return (  mean(-2*tmpllk) + var(-2*tmpllk)/2 )# D_bar + 1/2 var (D_theta), where D_theta = -2*tmpllk, and D_bar = mean(D_theta)

    double dicWSAFBar = -2 * sumOfVec(tmpLLKs1);
    double dicByTheta = (-2*mean) + (-2*mean) - dicWSAFBar;
    this->dEploidIO_->setdicByTheta(dicByTheta);
    //DIC.WSAF.bar = -2 * sum(thetallk)
    //return (  mean(-2*tmpllk) + (mean(-2*tmpllk) - DIC.WSAF.bar) ) # D_bar + pD, where pD = D_bar - D_theta, and D_bar = mean(D_theta)
}


vector <double> McmcMachinery::averageProportion() {
    assert(this->mcmcSample_->proportion.size()>0);
    vector <double> ret(this->kStrain());
    for ( size_t i = 0; i < kStrain(); i++ ) {
        for (vector <double> p : this->mcmcSample_->proportion) {
            ret[i] += p[i];
        }
        ret[i] /= (1.0*this->mcmcSample_->proportion.size());
    }
    (void)normalizeBySum(ret);
    return ret;
}


void McmcMachinery::sampleMcmcEvent( bool useIBD ) {
    this->recordingMcmcBool_ = ( currentMcmcIteration_ > this->mcmcThresh_ && currentMcmcIteration_ % this->McmcMachineryRate_ == 0 );
    if ( useIBD == true ) {
        ibdSampleMcmcEventStep();
        assert(doutProp());
    } else if (this->jobbrief == "classic"){
        this->eventInt_ = this->mcmcEventRg_->sampleInt(3);
        if ( (this->eventInt_ == 0) && (this->dEploidIO_->doUpdateProp() == true) ) {
            this->updateProportion();
        } else if ( (this->eventInt_ == 1) && (this->dEploidIO_->doUpdateSingle() == true) ) {
            this->updateSingleHap(this->panel_);
        } else if ( (this->eventInt_ == 2) && (this->dEploidIO_->doUpdatePair() == true) ) {
            this->updatePairHaps(this->panel_);
        }
    }else {  // Best practice
        this->eventInt_ = this->mcmcEventRg_->sampleInt(4);
        if ( (this->eventInt_ == 0) && (this->dEploidIO_->doUpdateProp() == true) ) {
            this->updateProportion();
        } else if ( (this->eventInt_ == 1) && (this->dEploidIO_->doUpdateSingle() == true) ) {
            this->updateSingleHap(this->panel_);
        } else if ( (this->eventInt_ == 2) && (this->dEploidIO_->doUpdatePair() == true) ) {
            this->updatePairHaps(this->panel_);
        } else if ( (this->eventInt_ == 3) && (this->dEploidIO_->doUpdateSingle() == true) ) {
            this->updateSingleHap(NULL);
            this->updateSingleHap(NULL);
            this->updateSingleHap(NULL);
            this->updateSingleHap(NULL);
        }
    }

    assert(doutLLK());

    if ( this->recordingMcmcBool_ ) {
        this->recordMcmcMachinery();
    }
}


void McmcMachinery::ibdInitializeEssentials() {

    this->initializePropIBD();
    this->ibdPath.init(*this->dEploidIO_, this->hapRg_);

    vector <double> llkOfData;
    for ( size_t i = 0; i < nLoci(); i++) {
        double wsaf = this->altCount_ptr_->at(i) / (this->refCount_ptr_->at(i) + this->altCount_ptr_->at(i) + 0.00000000000001);
        double adjustedWsaf = wsaf*(1-0.01) + (1-wsaf)*0.01;
        llkOfData.push_back( logBetaPdf(adjustedWsaf, this->ibdPath.llkSurf[i][0], this->ibdPath.llkSurf[i][1]));
    }
    dout << "LLK of data = " << sumOfVec(llkOfData) << endl;

}



void McmcMachinery::ibdSampleMcmcEventStep() {
    vector <double> effectiveKPrior = this->ibdPath.computeEffectiveKPrior(this->ibdPath.theta());
    vector <double> statePrior = this->ibdPath.computeStatePrior(effectiveKPrior);
    // First building the path likelihood
    this->ibdPath.computeIbdPathFwdProb(this->currentProp_, statePrior);

    ////#Now sample path given matrix
    this->ibdPath.ibdSamplePath(statePrior);

    //#Get haplotypes and update LLK for each site
    this->ibdUpdateHaplotypesFromPrior();
    vector <double> llkAtAllSites = computeLlkAtAllSites();
    ////#Given current haplotypes, sample titres 1 by 1 using MH
    vector <double> updatedllkAtAllSites = this->ibdUpdateProportionGivenHap(llkAtAllSites);
    // Compute new theta after all proportion and haplotypes are up to date.
    this->ibdPath.computeAndUpdateTheta();

    this->currentLLks_ = updatedllkAtAllSites;
    this->currentExpectedWsaf_ = this->calcExpectedWsaf( this->currentProp_ );
}


void McmcMachinery::ibdUpdateHaplotypesFromPrior() {
    for (size_t i = 0; i < this->nLoci(); i++) {
        for ( size_t j = 0; j < kStrain(); j++) {
            this->currentHap_[i][j] = (double)this->ibdPath.hprior.hSet[ibdPath.ibdConfigurePath[i]][j];
        }
    }
}


vector <double> McmcMachinery::ibdUpdateProportionGivenHap(
                        const vector <double> &llkAtAllSites) {
    vector <double> ret(llkAtAllSites.begin(), llkAtAllSites.end());
    for (size_t i = 0; i < kStrain(); i++) {
        double v0 = this->currentTitre_[i];
        vector <double> oldProp = this->currentProp_;
        //this->currentTitre_[i] += (this->stdNorm_->genReal() * 0.1 + 0.0); // tit.0[i]+rnorm(1, 0, scale.t.prop);
        this->currentTitre_[i] += (this->stdNorm_->genReal() * SD_LOG_TITRE* 1.0/PROP_SCALE + 0.0); // tit.0[i]+rnorm(1, 0, scale.t.prop);
        this->currentProp_ = this->titre2prop(this->currentTitre_);
        vector <double> vv = computeLlkAtAllSites();
        double rr = normal_pdf( this->currentTitre_[i], 0, 1) /
                    normal_pdf( v0, 0, 1) * exp( sumOfVec(vv) - sumOfVec(ret));

        if ( this->propRg_->sample() < rr) {
            //llkAtAllSites = vv;
            ret.clear();
            ret = vv;
            acceptUpdate++;
        } else {
            this->currentTitre_[i] = v0;
            this->currentProp_ = oldProp;
        }
    }
    return ret;
}


vector <double> McmcMachinery::computeLlkAtAllSites(double err) {
    vector <double > ret;
    for ( size_t site = 0; site < this->nLoci(); site++ ) {
        double qs = 0;
        for ( size_t j = 0; j < this->kStrain() ; j++ ) {
            qs += (double)this->currentHap_[site][j] * this->currentProp_[j];
        }
        double qs2 = qs*(1-err) + (1-qs)*err ;
        ret.push_back(logBetaPdf(qs2, this->ibdPath.llkSurf[site][0], this->ibdPath.llkSurf[site][1]));
    }
    return ret;
}


vector <double> McmcMachinery::calcExpectedWsaf(const vector <double> &proportion ) {
    //assert ( sumOfVec(proportion) == 1.0); // this fails ...
    vector <double> expectedWsaf (this->nLoci_, 0.0);
    for ( size_t i = 0; i < currentHap_.size(); i++ ) {
        assert( kStrain_ == currentHap_[i].size() );
        for ( size_t k = 0; k < kStrain_; k++) {
            expectedWsaf[i] += currentHap_[i][k] * proportion[k];
        }
        assert ( expectedWsaf[i] >= 0 );
        //assert ( expectedWsaf[i] <= 1.0 );
    }
    return expectedWsaf;
}


void McmcMachinery::recordMcmcMachinery() {
    dout << "***Record mcmc sample " <<endl;
    this->mcmcSample_->proportion.push_back(this->currentProp_);
    this->mcmcSample_->sumLLKs.push_back(sumOfVec(this->currentLLks_));
    this->mcmcSample_->moves.push_back(this->eventInt_);

    // Cumulate expectedWSAF for computing the mean expectedWSAF
    for ( size_t i = 0; i < this->cumExpectedWsaf_.size(); i++) {
        this->cumExpectedWsaf_[i] += this->currentExpectedWsaf_[i];
    }
}


void McmcMachinery::updateProportion() {
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
    vector <double> tmpLLKs = calcLLKs (*this->refCount_ptr_, *this->altCount_ptr_, tmpExpecedWsaf, 0, tmpExpecedWsaf.size(), this->dEploidIO_->scalingFactor());
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
    this->acceptUpdate++;

    this->currentExpectedWsaf_ = tmpExpecedWsaf;
    this->currentLLks_ = tmpLLKs;
    this->currentLogPriorTitre_ = tmpLogPriorTitre;
    this->currentTitre_ = tmpTitre;
    this->currentProp_ = tmpProp;

    assert (doutProp());
}


double McmcMachinery::deltaLLKs (const vector <double> &newLLKs ) {
    vector <double> tmpdiff = vecDiff ( newLLKs,  this->currentLLks_);
    return sumOfVec(tmpdiff);
}


vector <double> McmcMachinery::calcTmpTitre() {
    vector <double> tmpTitre;

    //size_t tmpEvet = this->mcmcEventRg_->sampleInt(2);
/*
    if (tmpEvet == 0){
*/
        for ( size_t k = 0; k < this->kStrain_; k++) {
            double dt = this->deltaXnormalVariable();
            tmpTitre.push_back( currentTitre_[k] + dt );
        }
/*
    } else {
        double dt;
        for ( size_t k = 0; k < this->kStrain_-1; k++) {
            dt = this->deltaXnormalVariable();
            tmpTitre.push_back( currentTitre_[k] + dt );
        }
        tmpTitre.push_back( currentTitre_.back() + dt );
    }
*/
    return tmpTitre;
}


void McmcMachinery::updateSingleHap(Panel *useThisPanel) {
    dout << " Update Single Hap "<<endl;
    this->findUpdatingStrainSingle();

    if ( this->dEploidIO_->doAllowInbreeding() == true ) {
        this->updateReferencePanel(this->panel_->truePanelSize()+kStrain_-1, this->strainIndex_);
    }

    for ( size_t chromi = 0 ; chromi < this->dEploidIO_->indexOfChromStarts_.size(); chromi++ ) {
        size_t start = this->dEploidIO_->indexOfChromStarts_[chromi];
        size_t length = this->dEploidIO_->position_[chromi].size();
        dout << "   Update Chrom with index " << chromi << ", starts at "<< start << ", with " << length << " sites" << endl;
        UpdateSingleHap updating( *this->refCount_ptr_,
                                  *this->altCount_ptr_,
                                  *this->plaf_ptr_,
                                  this->currentExpectedWsaf_,
                                  this->currentProp_, this->currentHap_, this->hapRg_,
                                  start, length,
                                  useThisPanel, this->dEploidIO_->missCopyProb_, this->dEploidIO_->scalingFactor(),
                                  this->strainIndex_);

        if ( this->dEploidIO_->doAllowInbreeding() == true ) {
            updating.setPanelSize(this->panel_->inbreedingPanelSize());
        }

        updating.core ( *this->refCount_ptr_, *this->altCount_ptr_, *this->plaf_ptr_, this->currentExpectedWsaf_, this->currentProp_, this->currentHap_);

        size_t updateIndex = 0;
        for ( size_t ii = start ; ii < (start+length); ii++ ) {
            this->currentHap_[ii][this->strainIndex_] = updating.hap_[updateIndex];
            this->currentLLks_[ii] = updating.newLLK[updateIndex];
            updateIndex++;
        }

        for (size_t siteI = 0; siteI< length; siteI++) {
            this->mcmcSample_->siteOfOneSwitchOne[start+siteI] += updating.siteOfOneSwitchOne[siteI];
            this->mcmcSample_->siteOfOneMissCopyOne[start+siteI] += updating.siteOfOneMissCopyOne[siteI];
            this->mcmcSample_->siteOfOneSwitchOne[start+siteI] = updating.siteOfOneSwitchOne[siteI];
            this->mcmcSample_->siteOfOneMissCopyOne[start+siteI] = updating.siteOfOneMissCopyOne[siteI];
        }
    }
    this->currentExpectedWsaf_ = this->calcExpectedWsaf( this->currentProp_ );
}


void McmcMachinery::updatePairHaps(Panel *useThisPanel) {
    if ( this->kStrain() == 1 ) {
        return;
    }

    dout << " Update Pair Hap "<<endl;
    this->findUpdatingStrainPair();

    for ( size_t chromi = 0 ; chromi < this->dEploidIO_->indexOfChromStarts_.size(); chromi++ ) {
        size_t start = this->dEploidIO_->indexOfChromStarts_[chromi];
        size_t length = this->dEploidIO_->position_[chromi].size();
        dout << "   Update Chrom with index " << chromi << ", starts at "<< start << ", with " << length << " sites" << endl;

        UpdatePairHap updating( *this->refCount_ptr_,
                                *this->altCount_ptr_,
                                *this->plaf_ptr_,
                                this->currentExpectedWsaf_,
                                this->currentProp_, this->currentHap_, this->hapRg_,
                                start, length,
                                useThisPanel, this->dEploidIO_->missCopyProb_, this->dEploidIO_->scalingFactor(),
                                this->dEploidIO_->forbidCopyFromSame(),
                                this->strainIndex1_,
                                this->strainIndex2_);

        updating.core(*this->refCount_ptr_, *this->altCount_ptr_, *this->plaf_ptr_, this->currentExpectedWsaf_, this->currentProp_, this->currentHap_);

        size_t updateIndex = 0;
        for ( size_t ii = start ; ii < (start+length); ii++ ) {
            this->currentHap_[ii][this->strainIndex1_] = updating.hap1_[updateIndex];
            this->currentHap_[ii][this->strainIndex2_] = updating.hap2_[updateIndex];
            this->currentLLks_[ii] = updating.newLLK[updateIndex];
            updateIndex++;
        }

        for (size_t siteI = 0; siteI< length; siteI++) {
            this->mcmcSample_->siteOfTwoSwitchOne[start+siteI] += updating.siteOfTwoSwitchOne[siteI];
            this->mcmcSample_->siteOfTwoMissCopyOne[start+siteI] += updating.siteOfTwoMissCopyOne[siteI];
            this->mcmcSample_->siteOfTwoSwitchTwo[start+siteI] += updating.siteOfTwoSwitchTwo[siteI];
            this->mcmcSample_->siteOfTwoMissCopyTwo[start+siteI] += updating.siteOfTwoMissCopyTwo[siteI];
            this->mcmcSample_->currentsiteOfTwoSwitchOne[start+siteI] = updating.siteOfTwoSwitchOne[siteI];
            this->mcmcSample_->currentsiteOfTwoMissCopyOne[start+siteI] = updating.siteOfTwoMissCopyOne[siteI];
            this->mcmcSample_->currentsiteOfTwoSwitchTwo[start+siteI] = updating.siteOfTwoSwitchTwo[siteI];
            this->mcmcSample_->currentsiteOfTwoMissCopyTwo[start+siteI] = updating.siteOfTwoMissCopyTwo[siteI];
        }
    }

    this->currentExpectedWsaf_ = this->calcExpectedWsaf( this->currentProp_ );
}


void McmcMachinery::findUpdatingStrainSingle( ) {
    vector <double> eventProb (this->kStrain_, 1);
    (void)normalizeBySum(eventProb);
    this->strainIndex_ = sampleIndexGivenProp ( this->mcmcEventRg_, eventProb );
    dout << "  Updating hap: "<< this->strainIndex_ <<endl;
}


void McmcMachinery::findUpdatingStrainPair( ) {
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
