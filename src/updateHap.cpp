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

#include "updateHap.hpp"
#include <algorithm>    // std::reverse


UpdateHap::UpdateHap( vector <double> &refCount,
                      vector <double> &altCount,
                      vector <double> &expectedWsaf,
                      vector <double> &proportion,
                      vector < vector <double> > &haplotypes, MersenneTwister* rg, Panel* panel){
    this->panel_ = panel;
    this->nPanel_ = this->panel_->nPanel_;

    this->missCopyProb = 0.01;
    this->kStrain_ = proportion.size();
    this->nLoci_ = expectedWsaf.size();
    this->rg_ = rg;
}


UpdateSingleHap::UpdateSingleHap( vector <double> &refCount,
                                  vector <double> &altCount,
                                  vector <double> &expectedWsaf,
                                  vector <double> &proportion,
                                  vector < vector <double> > &haplotypes, MersenneTwister* rg, Panel* panel):
                UpdateHap(refCount, altCount, expectedWsaf, proportion, haplotypes, rg, panel){
    this->findUpdatingStrain( proportion );
    this->calcExpectedWsaf( expectedWsaf, proportion, haplotypes);
    this->calcHapLLKs(refCount, altCount);
    this->buildEmission();
    this->calcFwdProbs();
    this->samplePaths();
    this->addMissCopying();
}


void UpdateSingleHap::findUpdatingStrain( vector <double> proportion ){
    this->strainIndex_ = sampleIndexGivenProp ( proportion, this->rg_ );
    ////vector <size_t> strainIndex = sampleNoReplace( this->currentProp_, this->rg_ );
    ////assert( strainIndex.size() == 1);
    ////size_t ws = strainIndex[0];
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
        expectedWsaf1_[i] += proportion[strainIndex_] ;
    }
    cout << "expectedWsaf[0] = "<<expectedWsaf[0]<< " expectedWsaf0_[0] = "<<expectedWsaf0_[0]<<" expectedWsaf1_[0] = "<<expectedWsaf1_[0]<<endl;
}


void UpdateSingleHap::buildEmission(){
    vector <double> noMissProb (this->nLoci_, log(1.0 - this->missCopyProb));
    vector <double> t1omu = vecSum(llk0_, noMissProb);
    vector <double> t2omu = vecSum(llk1_, noMissProb);

    vector <double> missProb (this->nLoci_, log(this->missCopyProb));
    vector <double> t1u = vecSum(llk0_, missProb);
    vector <double> t2u = vecSum(llk1_, missProb);

    assert(emission_.size() == 0 );
    for ( size_t i = 0; i < this->nLoci_; i++){
        vector <double> tmp ({t1omu[i], t2omu[i], t1u[i], t2u[i]});
        double tmaxTmp = maxOfVec(tmp);
        vector <double> emissRow ({exp(t1omu[i] - tmaxTmp) + exp(t2u[i] - tmaxTmp),
                                   exp(t2omu[i] - tmaxTmp) + exp(t1u[i] - tmaxTmp)});
        this->emission_.push_back(emissRow);
    }
}


void UpdateSingleHap::calcFwdProbs(){
    assert ( this->fwdProbs_.size() == 0 );
    vector <double> fwd1st (this->nPanel_, 0.0);
    for ( size_t i = 0 ; i < this->nPanel_; i++){
        fwd1st[i] = this->emission_[0][this->panel_->content_[0][i]];
    }
    (void)normalizeBySum(fwd1st);
    this->fwdProbs_.push_back(fwd1st);

    for ( size_t j = 1; j < this->nLoci_; j++ ){
        double pRec = this->panel_->recombProbs_[j-1];
        double pRecEachHap = pRec / nPanel_;
        double pNoRec = 1.0 - pRec;

        double massFromRec = sumOfVec(fwdProbs_[j-1]) * pRecEachHap;
        vector <double> fwdTmp (this->nPanel_, 0.0);
        for ( size_t i = 0 ; i < this->nPanel_; i++){
            fwdTmp[i] = this->emission_[j][this->panel_->content_[j][i]] * (fwdProbs_[j-1][i] * pNoRec + massFromRec);
        }
        (void)normalizeBySum(fwdTmp);
        this->fwdProbs_.push_back(fwdTmp);
    }
}


void UpdateSingleHap::calcHapLLKs( vector <double> &refCount,
                                   vector <double> &altCount){
    llk0_ = calcLLKs( refCount, altCount, expectedWsaf0_ );
    llk1_ = calcLLKs( refCount, altCount, expectedWsaf1_ );
}


void UpdateSingleHap::samplePaths(){
    assert ( this->path_.size() == 0 );
    size_t pathTmp = sampleIndexGivenProp ( fwdProbs_[nLoci_-1], this->rg_ );
    //cout << pathTmp << endl;
    this->path_.push_back( this->panel_->content_.back()[pathTmp]);
    //cout << path_.back()<<endl;
    for ( size_t j = (this->nLoci_ - 1); j > 0; j--){
        dout << "j="<<j<<", ";
        double pRec = this->panel_->recombProbs_[j-1];
        double pRecEachHap = pRec / nPanel_;
        double pNoRec = 1.0 - pRec;

        vector <double> previousDist = fwdProbs_[j-1];
        vector <double> weightOfNoRecAndRec ({ previousDist[pathTmp]*pNoRec,
                                               sumOfVec(previousDist)*pRecEachHap});
        (void)normalizeBySum(weightOfNoRecAndRec);
        size_t tmpState = sampleIndexGivenProp(weightOfNoRecAndRec, this->rg_ );
        dout <<" recomb state = " << tmpState<< ", ";
        if ( tmpState == (size_t)1){
            pathTmp = sampleIndexGivenProp( previousDist, this->rg_ );
        }
        dout <<"path = "<<pathTmp<<endl;
        this->path_.push_back(this->panel_->content_[j][pathTmp]);
    }
    reverse(path_.begin(), path_.end());
    assert(path_.size() == nLoci_);
}


void UpdateSingleHap::addMissCopying(){
    newLLK = vector <double> (this->nLoci_, 0.0);
    for ( size_t i = 0; i < this->nLoci_; i++){
        double tmpLLKmax = (this->llk0_[i] > this->llk1_[i] ? this->llk0_[i] : this->llk1_[i]) +0.00000001; // Add 0.00000001, in order to prevent that tmpLLKmax = 0
        dout <<"site "<<i<<" tmpLLKmax = "<<tmpLLKmax <<" "<<this->llk0_[i]<<" "<<this->llk1_[i]<<endl;
        dout <<"site "<<i<<" expectedWsaf0_[i] = "<<expectedWsaf0_[i] <<" expectedWsaf1_[i] = "<<expectedWsaf1_[i]<<endl;
        vector <double> emissionTmp ({exp(this->llk0_[i]/tmpLLKmax), exp(this->llk1_[i]/tmpLLKmax)});
        dout << "emissionTmp[0] = "<<emissionTmp[0]<<" emissionTmp[1] = "<<emissionTmp[1]<<endl;
        vector <double> sameDiffDist ({emissionTmp[path_[i]]*(1.0 - this->missCopyProb),
                                       emissionTmp[ 1 -path_[i]] * this->missCopyProb });
        (void)normalizeBySum(sameDiffDist);
        if ( sampleIndexGivenProp(sameDiffDist, this->rg_ ) == 1 ){
            this->hap_.push_back( 1 - this->path_[i] ); // differ
        } else {
            this->hap_.push_back( this->path_[i] ); // same
        }

        if ( this->hap_[i] == 0){
            newLLK[i] = llk0_[i];
        } else if (this->hap_[i] == 1){
            newLLK[i] = llk1_[i];
        } else {
            throw("should never get here!");
        }
        dout <<"site "<<i<<" "<<this->hap_[i] <<endl;
    }
}


