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

McmcMachinery::McmcMachinery(DEploidIO* dEploidIO, McmcSample *mcmcSample, RandomGenerator* rg_ ){ // initialiseMCMCmachinery

    this->dEploidIO_ = dEploidIO;
    this->panel_ = dEploidIO->panel;
    this->mcmcSample_ = mcmcSample;
    this->seed_ = rg_->seed();

    //this->hapRg_ = new MersenneTwister(this->seed_);
    this->hapRg_ = rg_ ;
    this->mcmcEventRg_ = this->hapRg_;
    this->propRg_ = this->hapRg_;
    this->initialHapRg_ = this->hapRg_;
    this->ibdRg_ = this->hapRg_;
    //this->mcmcEventRg_ = new MersenneTwister(this->seed_);
    //this->propRg_  = new MersenneTwister(this->seed_);
    //this->initialHapRg_ = new MersenneTwister(this->seed_);

    this->calcMaxIteration( dEploidIO_->nMcmcSample_ , dEploidIO_->mcmcMachineryRate_, dEploidIO_->mcmcBurn_ );

    this->MN_LOG_TITRE = 0.0;
    this->SD_LOG_TITRE = 3.0;
    this->PROP_SCALE = 40.0;

    stdNorm_ = new StandNormalRandomSample(this->seed_);

    this->setKstrain(this->dEploidIO_->kStrain());
    this->setNLoci(this->dEploidIO_->plaf_.size());
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

    if ( this->dEploidIO_->doAllowInbreeding() == true ){
        this->initializeUpdateReferencePanel(this->panel_->truePanelSize()+kStrain_-1);
    }

    assert (doutProp());
    assert (doutLLK());
    dout << "###########################################"<< endl;
    dout << "#        Initialization finished          #"<< endl;
    dout << "###########################################"<< endl;
}


void McmcMachinery::initializeHap(){
    assert( currentHap_.size() == 0);
    if ( this->dEploidIO_ -> initialHapWasGiven() ){
        this->currentHap_ = this->dEploidIO_->initialHap;
    } else {
        for ( size_t i = 0; i < this->dEploidIO_->plaf_.size(); i++ ){
            double currentPlaf = this->dEploidIO_->plaf_[i];
            vector <double> tmpVec;
            for ( size_t k = 0; k < this->kStrain_; k++){
                tmpVec.push_back( this->rBernoulli(currentPlaf) );
            }
            this->currentHap_.push_back(tmpVec);
        }
    }
    assert(this->currentHap_.size() == this->dEploidIO_->plaf_.size());
}


void McmcMachinery::initializeUpdateReferencePanel(size_t inbreedingPanelSizeSetTo){
    if ( this->dEploidIO_->doAllowInbreeding() != true ){
        return;
    }

    this->panel_->initializeUpdatePanel(inbreedingPanelSizeSetTo);
}


void McmcMachinery::updateReferencePanel(size_t inbreedingPanelSizeSetTo, size_t excludedStrain){
    if ( this->burnIn_ > this->currentMcmcIteration_ ){
        return;
    }

    //if ( this->dEploidIO_->doAllowInbreeding() != true ){
        //return;
    //}
    this->panel_->updatePanelWithHaps( inbreedingPanelSizeSetTo, excludedStrain, this->currentHap_);
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


void McmcMachinery::makeLlkSurf(vector <double> altCount, vector <double> refCount, double scalingConst, double err, size_t gridSize){
    //n.site<-length(alt.cts);
    double pGridSpacing = 1.0 / (double)(gridSize+1);
    vector <double> pGrid;
    pGrid.push_back(pGridSpacing);
    for (size_t i = 1; i < gridSize; i++){
        pGrid.push_back(pGrid.back() + pGridSpacing);
    }
    assert(pGrid.size() == gridSize);

    assert(llkSurf.size() == 0);

    for ( size_t i = 0 ; i < altCount.size(); i++) {
        double alt = altCount[i];
        double ref = refCount[i];

        vector <double> ll;
        for ( double unadjustedP : pGrid ){
            ll.push_back(calcLLK( ref, alt, unadjustedP, err, scalingConst));
        }

        double llmax = max_value(ll);
        vector <double> ln;
        for ( double lltmp : ll ){
            ln.push_back(exp(lltmp-llmax));
        }

        double lnSum = sumOfVec(ln);
        for (size_t i = 0; i < ln.size(); i++){
            ln[i] = ln[i]/lnSum;
        }

        //double mn;
        //mn<-sum(ln*p.grid); # Compute the posterior, then normalize by the mean?
        vector <double> tmpVec1 = vecProd(ln, pGrid);
        double mn = sumOfVec(tmpVec1);

        //double vr;
        //vr<-sum(ln*p.grid^2)-mn^2;
        vector <double> pGridSq = vecProd(pGrid, pGrid);
        vector <double> tmpVec2 = vecProd(ln, pGridSq);
        double vr = sumOfVec(tmpVec2) - mn*mn;

        double comm = (mn*(1.0-mn)/vr-1.0);
        llkSurf.push_back(vector <double> {mn*comm, (1-mn)*comm});
        //llk.apx[site,]<-(mn*(1-mn)/vr-1)*c(mn, 1-mn);

    }
    assert(llkSurf.size() == this->nLoci());
}


void McmcMachinery::initializePropIBD(){
    //#Initialise titres and convert to proportions
    this->currentTitre_ = vector <double> (this->kStrain(), 0.0);
    this->currentProp_ = vector <double> (this->kStrain(), 1.0/(double)kStrain());
}


void McmcMachinery::runMcmcChain( bool showProgress, bool useIBD ){
    if ( useIBD ){
        this->initializeIbdEssentials();
    }

    for ( this->currentMcmcIteration_ = 0 ; currentMcmcIteration_ < this->maxIteration_ ; currentMcmcIteration_++){
        dout << endl;
        dout << "MCMC iteration: " << this->currentMcmcIteration_ << endl;
        if ( this->currentMcmcIteration_ > 0 && this->currentMcmcIteration_%100 == 0 && showProgress ){
            clog << "\r" << " MCMC step" << setw(4) << int(currentMcmcIteration_ * 100 / this->maxIteration_) << "% completed."<<flush;
        }
        this->sampleMcmcEvent(useIBD);
    }
    clog << "\r" << " MCMC step" << setw(4) << 100 << "% completed."<<endl;
    this->mcmcSample_->hap = this->currentHap_;

    this->writeLastFwdProb();

    this->dEploidIO_->filnalProp = this->mcmcSample_->proportion.back();
    this->dEploidIO_->writeMcmcRelated(this->mcmcSample_, useIBD);

    if ( useIBD ){
        this->dEploidIO_->initialProp = this->mcmcSample_->proportion.back();
        this->dEploidIO_->initialHap = this->mcmcSample_->hap;
        this->dEploidIO_->setInitialHapWasGiven(true);
    }
}


void McmcMachinery::sampleMcmcEvent( bool useIBD ){
    this->recordingMcmcBool_ = ( currentMcmcIteration_ > this->mcmcThresh_ && currentMcmcIteration_ % this->McmcMachineryRate_ == 0 );
    if ( useIBD == true ){
        sampleMcmcEventIbdStep();
    } else {
        this->eventInt_ = this->mcmcEventRg_->sampleInt(3);
        if ( (this->eventInt_ == 0) && (this->dEploidIO_->doUpdateProp() == true) ){
            this->updateProportion();
        } else if ( (this->eventInt_ == 1) && (this->dEploidIO_->doUpdateSingle() == true) ){
            this->updateSingleHap();
        } else if ( (this->eventInt_ == 2) && (this->dEploidIO_->doUpdatePair() == true) ){
            this->updatePairHaps();
        }
    }

    assert(doutLLK());

    if ( this->recordingMcmcBool_ ){
        this->recordMcmcMachinery();
    }
}


vector <size_t> McmcMachinery::findWhichIsSomething(vector <size_t> tmpOp, size_t something){
    vector <size_t> ret;
    for ( size_t i = 0; i < tmpOp.size(); i++){
        if ( tmpOp[i] == something){
            ret.push_back(i);
        }
    }
    return ret;
}


void McmcMachinery::makeIbdTransProbs(){
    size_t nPattern = hprior.nPattern();
    size_t nState = hprior.nState();
    //#Transition probs
    assert(ibdTransProbs.size() == 0);

    //tij<-array(0, c(n.pattern, n.state));

    for ( size_t i = 0; i < nPattern; i++ ){
        vector <double> transProbRow(nState);
        //wi<-which(state.idx==i);
        vector <size_t> wi = findWhichIsSomething(hprior.stateIdx, i);
        //tij[i,wi]<-1;
        for (size_t wii : wi){
            transProbRow[wii] = 1;
        }
        ibdTransProbs.push_back(transProbRow);
    }
}


void McmcMachinery::computeUniqueEffectiveKCount(){
    this->uniqueEffectiveKCount = vector <int> (this->kStrain());
    for (size_t effectiveKtmp : this->hprior.effectiveK) {
        int effectiveKidx = effectiveKtmp-1;
        assert(effectiveKidx>=0);
        this->uniqueEffectiveKCount[effectiveKidx]++;
    }
}


vector <double> McmcMachinery::computeStatePrior( double theta ){
    //#Calculate state prior given theta (theta is prob IBD)
    //pr.0<-dbinom(0:(k.max-1), k.max-1, theta.0);
    vector <double> pr0(this->kStrain());
    for (int i = 0; i < (int)pr0.size(); i++){
        pr0[i] = binomialPdf(i, (int)(this->kStrain()-1), theta);
    }
    //a.prior<-pr.0[ibd.st$k.eff]/a.ct[ibd.st$k.eff];
    //st.prior<-a.prior[state.idx];

    vector <double> effectiveKPrior;
    for ( size_t effectiveKtmp : this->hprior.effectiveK){
        int effectiveKidx = effectiveKtmp-1;
        assert(effectiveKidx >= 0);
        assert(effectiveKidx < (int)this->kStrain());
        effectiveKPrior.push_back(pr0[effectiveKidx]/uniqueEffectiveKCount[effectiveKidx]);
    }

    vector <double> ret;
    for (size_t stateIdxTmp : this->hprior.stateIdx){
        ret.push_back(effectiveKPrior[stateIdxTmp]);
    }
    return ret;
}


void McmcMachinery::initializeIbdEssentials(){
    // initialize haplotype prior
    this->hprior.buildHprior(this->kStrain(), this->dEploidIO_->plaf_);
    this->hprior.transposePriorProbs();
    // compute likelihood surface
    this->makeLlkSurf(this->dEploidIO_->altCount_, this->dEploidIO_->refCount_);
    this->initializePropIBD();
    this->setTheta(1.0 / (double)kStrain());
    this->makeIbdTransProbs();
    this->computeUniqueEffectiveKCount();


    // initialize fm

    this->fSumState = vector <double> (this->hprior.nPattern());

    // initialize ibdPath
    this->ibdPath = vector <size_t> (this->nLoci());
}


vector <double> McmcMachinery::computeLlkOfStatesAtSiteI( size_t siteI, double err ){

    //qs<-h.set %*% t(prop.0); # estimated proportion
    //qs2<-qs*(1-err)+(1-qs)*err; # adjusted estimated proportion
    //lk.data<-dbeta(qs2, llk.apx[1,1], llk.apx[1,2], log=T); # likelihood of the proportion
    //lk.norm<-exp(lk.data-max(lk.data));

    vector <double> qs2;
    for ( size_t i = 0; i < this->hprior.hSet.size(); i++ ){
        double qs = 0;
        for ( size_t j = 0; j < this->hprior.hSet[i].size() ; j++ ){
            qs += (double)this->hprior.hSet[i][j] * this->currentProp_[j];
        }
        qs2.push_back( qs*(1-err) + (1-qs)*err );
    }

    vector <double> llks;
    for (double qs2Tmp : qs2){
        llks.push_back(logBetaPdf(qs2Tmp, this->llkSurf[siteI][0], this->llkSurf[siteI][1]));
    }

    double maxllk = max_value(llks);
    vector <double> ret;
    for ( double llk : llks ){
        ret.push_back( exp(llk-maxllk));
    }
    return ret;
}


void McmcMachinery::updateFmAtSiteI(vector <double> & prior, vector <double> & llk){
    //cout << "prior = ";
    //for ( double p : prior ){
        //cout << p <<" ";
    //}
    //cout << "llk = ";
    //for ( double p : llk ){
        //cout << p <<" ";
    //}

    vector <double> postAtSiteI = vecProd(prior, llk);
    //cout << "postAtSiteI = ";
    //for ( double p : postAtSiteI ){
        //cout << p <<" ";
    //}
    //cout<<endl;
    normalizeByMax(postAtSiteI);
    //cout << "postAtSiteI = ";
    //for ( double p : postAtSiteI ){
        //cout << p <<" ";
    //}
    //cout<<endl;
    this->fm.push_back(postAtSiteI);
    this->fSum = sumOfVec(postAtSiteI);
//cout << "fSum =" << fSum<<endl;
    for ( size_t i = 0; i < fSumState.size(); i++){
        this->fSumState[i] = 0;
        for ( size_t j = 0; j < hprior.nState(); j++ ){
            this->fSumState[i] += ibdTransProbs[i][j]*postAtSiteI[j];
        }
    }

}


void McmcMachinery::sampleMcmcEventIbdStep(){
    double pRecomb = 0.01;
    double pNoRecomb = 0.99;
    vector <double> statePrior = this->computeStatePrior(this->theta());

    ////#Given current proportions sample allelic and inbreeding states at each site
    ////#Initialise:
    //v.prior<-st.prior * h.prior[,1];
    vector <double> vPrior = vecProd(statePrior, this->hprior.priorProbTrans[0]);
    vector <double> lk = computeLlkOfStatesAtSiteI(0);
    this->updateFmAtSiteI(vPrior, lk);
    for ( size_t siteI = 1; siteI < this->nLoci(); siteI++ ){
        vector <double> vNoRec;
        for ( size_t stateIdxTmp : hprior.stateIdx ){
            vNoRec.push_back(this->fSumState[stateIdxTmp]);
        }

        //v.norec<-f.sum.state[state.idx];
        for ( size_t i = 0; i < hprior.nState(); i++ ){
            vPrior[i] = (vNoRec[i] * pNoRecomb + fSum * pRecomb * statePrior[i]) * hprior.priorProbTrans[siteI][i];
        }
        //fm[,site]<-(v.norec*(1-p.rec)+f.sum*p.rec*st.prior)*h.prior[,site];
        lk = computeLlkOfStatesAtSiteI(siteI);
        this->updateFmAtSiteI(vPrior, lk);
    }

    ////#Now sample path given matrix
    int lociIdx = this->nLoci()-1;
    vector <double> tmpProp = fm[lociIdx];
    normalizeBySum(tmpProp);
    ibdPath[lociIdx] = sampleIndexGivenProp(this->ibdRg_, tmpProp);
    //cout << "tmpProp = ";
    //for ( double p : tmpProp ){
        //cout << p <<" ";
    //}
    //cout<<endl;
    //cout << "sum(tmpProp) = " << sumOfVec(tmpProp) <<endl;
    //cout << "ibdPath[lociIdx] = "<<ibdPath[lociIdx]<<endl;

    while ( lociIdx > 0 ){
        lociIdx--;
        //cout << "lociIdx = "<<lociIdx <<", ";
        vector <double> vNoRecomb = vecProd(this->ibdTransProbs[this->hprior.stateIdx[ibdPath[lociIdx+1]]], fm[lociIdx]);
        assert(vNoRecomb.size() == this->hprior.nState());
        vector <double> vRecomb = fm[lociIdx];
        assert(vRecomb.size() == this->hprior.nState());

        vector <double> prop (this->hprior.nState());
        //cout << "prop = ";
        for ( size_t i = 0; i < prop.size(); i++){
            prop[i] = vNoRecomb[i]*pNoRecomb + vRecomb[i]*pRecomb*statePrior[ibdPath[lociIdx+1]];
            //cout << prop[i]<<" ";
        }
        tmpProp = prop;
        normalizeBySum(tmpProp);
        //cout <<endl;
        //v.norec<-tij[state.idx[a.path[site+1]],] * fm[,site] * (1-p.rec);
        //v.rec<-fm[,site] * p.rec * st.prior[a.path[site+1]];
        //a.path[site]<-sample(n.state, 1, p=v.norec+v.rec);
        ibdPath[lociIdx] = sampleIndexGivenProp(this->ibdRg_, tmpProp);
        //cout << "tmpProp = ";
        //for ( double p : tmpProp ){
            //cout << p <<" ";
        //}
        //cout<<endl;
        //cout << "sum(tmpProp) = " << sumOfVec(tmpProp) <<endl;
        //cout << "ibdPath["<<lociIdx<<"] = "<<ibdPath[lociIdx]<<endl;


        assert( ibdPath[lociIdx] < this->hprior.nState() );
        assert( ibdPath[lociIdx] >= 0 );
    }

    //a.path[n.loci]<-sample(n.state, 1, p=fm[,n.loci]);
    //for (site in (n.loci-1):1) {
        //v.norec<-tij[state.idx[a.path[site+1]],] * fm[,site] * (1-p.rec);
        //v.rec<-fm[,site] * p.rec * st.prior[a.path[site+1]];
        //a.path[site]<-sample(n.state, 1, p=v.norec+v.rec);
    //}

    //#Get haplotypes and update LLK for each site
    for (size_t i = 0; i < this->nLoci(); i++){
        //cout << "ibdPath["<<i<<"] = "<<ibdPath[i]<<endl;
        for ( size_t j = 0; j < kStrain(); j++){
            //cout << " "<<this->hprior.hSet[ibdPath[i]][j];
            this->currentHap_[i][j] = (double)this->hprior.hSet[ibdPath[i]][j];
        }
        //cout<<endl;
        //h.0[,site]<-h.set[a.path[site],];
    }
    //cout<<"hea"<<endl;
    //qs<-prop.0 %*% h.0;
    //qs2<-qs*(1-err)+(1-qs)*err;
    //llk.site<-dbeta(qs2, llk.apx[,1], llk.apx[,2], log=T);

    ////#Given current haplotypes, sample titres 1 by 1 using MH
    //for (i in 1:k.max) {
        //v0<-tit.0[i];
        //p.old<-prop.0;
        //tit.0[i]<-tit.0[i]+rnorm(1, 0, scale.t.prop);
        //prop.0[1,]<-exp(tit.0)/sum(exp(tit.0));
        //qs<-prop.0 %*% h.0;
        //qs2<-qs*(1-err)+(1-qs)*err;
        //vv<-dbeta(qs2, llk.apx[,1], llk.apx[,2], log=T);
        //rr<-dnorm(tit.0[i],0,1)/dnorm(v0, 0, 1) * exp(sum(vv)-sum(llk.site));
        //if (runif(1)<rr) {
            //llk.site[]<-vv;
            //prop.accpt<-prop.accpt+1;
        //} else {
            //tit.0[i]<-v0;
            //prop.0<-p.old;
        //}
    //}

    ////#Now update inbreeding parameter (theta);
    //s.path<-c(0,state.idx[a.path]);
    //obs.states<-s.path[which(diff(s.path)!=0)+1];
    //k.eff.states<-ibd.st$k.eff[obs.states];
    //sccs<-sum(k.max-k.eff.states);
    //theta.0<-rbeta(1, sccs+1, sum(k.eff.states-1)+1);

    //llk.hist = c(llk.hist, sum(llk.site));
    //prop.hist = rbind(prop.hist, prop.0)

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

    if ( this->dEploidIO_->doAllowInbreeding() == true ){
        this->updateReferencePanel(this->panel_->truePanelSize()+kStrain_-1, this->strainIndex_);
    }

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

        if ( this->dEploidIO_->doAllowInbreeding() == true ){
            updating.setPanelSize(this->panel_->inbreedingPanelSize());
        }

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


