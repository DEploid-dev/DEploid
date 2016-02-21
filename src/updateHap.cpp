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
    dout << "Update single hap: "<< this->strainIndex_ << " strain"<<endl;
}

void UpdateSingleHap::calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes ){
    //expected.WSAF.0 <- bundle$expected.WSAF - (bundle$prop[ws] * bundle$h[,ws]);
    this->expectedWsaf0_ = expectedWsaf;
    for ( size_t i = 0; i < expectedWsaf0_.size(); i++ ){
        expectedWsaf0_[i] -= proportion[strainIndex_] * haplotypes[i][strainIndex_];
        //dout << expectedWsaf[i] << " " << expectedWsaf0_[i] << endl;
        assert (expectedWsaf0_[i] >= 0 );
        assert (expectedWsaf0_[i] < 1 );
    }

    //expected.WSAF.1 <- expected.WSAF.0 + bundle$prop[ws] ;
    this->expectedWsaf1_ = expectedWsaf0_;
    for ( size_t i = 0; i < expectedWsaf1_.size(); i++ ){
        expectedWsaf1_[i] += proportion[strainIndex_] ;
    }
    //cout << "expectedWsaf[0] = "<<expectedWsaf[0]<< " expectedWsaf0_[0] = "<<expectedWsaf0_[0]<<" expectedWsaf1_[0] = "<<expectedWsaf1_[0]<<endl;
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
    this->llk0_ = calcLLKs( refCount, altCount, expectedWsaf0_ );
    this->llk1_ = calcLLKs( refCount, altCount, expectedWsaf1_ );
}


void UpdateSingleHap::samplePaths(){
    assert ( this->path_.size() == 0 );
    size_t pathTmp = sampleIndexGivenProp ( fwdProbs_[nLoci_-1], this->rg_ );
    //cout << pathTmp << endl;
    this->path_.push_back( this->panel_->content_.back()[pathTmp]);
    //cout << path_.back()<<endl;
    for ( size_t j = (this->nLoci_ - 1); j > 0; j--){
        //dout << "j="<<j<<", ";
        double pRec = this->panel_->recombProbs_[j-1];
        double pRecEachHap = pRec / nPanel_;
        double pNoRec = 1.0 - pRec;

        vector <double> previousDist = fwdProbs_[j-1];
        vector <double> weightOfNoRecAndRec ({ previousDist[pathTmp]*pNoRec,
                                               sumOfVec(previousDist)*pRecEachHap});
        (void)normalizeBySum(weightOfNoRecAndRec);
        size_t tmpState = sampleIndexGivenProp(weightOfNoRecAndRec, this->rg_ );
        //dout <<" recomb state = " << tmpState<< ", ";
        if ( tmpState == (size_t)1){
            pathTmp = sampleIndexGivenProp( previousDist, this->rg_ );
        }
        //dout <<"path = "<<pathTmp<<endl;
        this->path_.push_back(this->panel_->content_[j][pathTmp]);
    }
    reverse(path_.begin(), path_.end());
    assert(path_.size() == nLoci_);
}


void UpdateSingleHap::addMissCopying(){
    assert( this->hap_.size() == 0 );
    newLLK = vector <double> (this->nLoci_, 0.0);
    for ( size_t i = 0; i < this->nLoci_; i++){
          //double tmpLLKmax = (this->llk0_[i] > this->llk1_[i] ? this->llk0_[i] : this->llk1_[i]) +0.00000001; // Add 0.00000001, in order to prevent that tmpLLKmax = 0
        //dout <<"site "<<i<<" tmpLLKmax = "<<tmpLLKmax <<" "<<this->llk0_[i]<<" "<<this->llk1_[i]<<endl;
        //dout <<"site "<<i<<" expectedWsaf0_[i] = "<<expectedWsaf0_[i] <<" expectedWsaf1_[i] = "<<expectedWsaf1_[i]<<endl;
        //vector <double> emissionTmp ({exp(this->llk0_[i]/tmpLLKmax), exp(this->llk1_[i]/tmpLLKmax)});
        vector <double> emissionTmp ({exp(this->llk0_[i]), exp(this->llk1_[i])});
        vector <double> sameDiffDist ({emissionTmp[path_[i]]*(1.0 - this->missCopyProb), // probability of the same
                                       emissionTmp[(size_t)(1 -path_[i])] * this->missCopyProb }); // probability of differ
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
        //dout <<"site "<<i<<" "<<this->hap_[i] <<endl;
    }
}


UpdatePairHap::UpdatePairHap( vector <double> &refCount,
                              vector <double> &altCount,
                              vector <double> &expectedWsaf,
                              vector <double> &proportion,
                              vector < vector <double> > &haplotypes, MersenneTwister* rg, Panel* panel):
                UpdateHap(refCount, altCount, expectedWsaf, proportion, haplotypes, rg, panel){
    this->findUpdatingStrain( proportion );
    this->calcExpectedWsaf( expectedWsaf, proportion, haplotypes);
    this->calcHapLLKs(refCount, altCount);
    this->buildEmission();
    //this->calcFwdProbs();
    //this->samplePaths();
    //this->addMissCopying();
}

void UpdatePairHap:: findUpdatingStrain( vector <double> proportion ){
    vector <size_t> strainIndex = sampleNoReplace( proportion, this->rg_ , 2);
    assert( strainIndex.size() == 2);
    this->strainIndex1_ = strainIndex[0];
    this->strainIndex2_ = strainIndex[1];
    assert( strainIndex1_ != strainIndex2_ );
}


void UpdatePairHap:: calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes){
  //expected.WSAF.00 <- expected.WSAF-(prop[ws[1]]*h[,ws[1]] + prop[ws[2]]*h[,ws[2]]);
  //expected.WSAF.10 <- expected.WSAF.00 + prop[ws[1]];
  //expected.WSAF.01 <- expected.WSAF.00 + prop[ws[2]];
  //expected.WSAF.11 <- expected.WSAF.00 + prop[ws[1]] + prop[ws[2]];    //expected.WSAF.0 <- bundle$expected.WSAF - (bundle$prop[ws] * bundle$h[,ws]);
    this->expectedWsaf00_ = expectedWsaf;
    for ( size_t i = 0; i < expectedWsaf00_.size(); i++ ){
        //expectedWsaf00_[i] -= proportion[strainIndex1_] * haplotypes[i][strainIndex1_];
        //expectedWsaf00_[i] -= proportion[strainIndex2_] * haplotypes[i][strainIndex2_];
        expectedWsaf00_[i] -= (proportion[strainIndex1_] * haplotypes[i][strainIndex1_] + proportion[strainIndex2_] * haplotypes[i][strainIndex2_]);
        dout << expectedWsaf[i] << " " << expectedWsaf00_[i] << endl;
        assert (expectedWsaf00_[i] >= 0 );
        assert (expectedWsaf00_[i] < 1 );
    }

    this->expectedWsaf10_ = expectedWsaf00_;
    for ( size_t i = 0; i < expectedWsaf10_.size(); i++ ){
        expectedWsaf10_[i] += proportion[strainIndex1_] ;
    }

    this->expectedWsaf01_ = expectedWsaf00_;
    for ( size_t i = 0; i < expectedWsaf01_.size(); i++ ){
        expectedWsaf01_[i] += proportion[strainIndex2_] ;
    }

    this->expectedWsaf11_ = expectedWsaf00_;
    for ( size_t i = 0; i < expectedWsaf11_.size(); i++ ){
        expectedWsaf11_[i] += (proportion[strainIndex1_] + proportion[strainIndex2_]);
    }
}


void UpdatePairHap:: calcHapLLKs( vector <double> &refCount, vector <double> &altCount){
    this->llk00_ = calcLLKs( refCount, altCount, expectedWsaf00_ );
    this->llk10_ = calcLLKs( refCount, altCount, expectedWsaf10_ );
    this->llk01_ = calcLLKs( refCount, altCount, expectedWsaf01_ );
    this->llk11_ = calcLLKs( refCount, altCount, expectedWsaf11_ );
}


void UpdatePairHap:: buildEmission(){
    //llk.00 = logemiss[,1]
    //llk.10 = logemiss[,2]
    //llk.01 = logemiss[,3]
    //llk.11 = logemiss[,4]

    //log.omu = log(1-miss.copy.rate)
    //log.u = log(miss.copy.rate)
    vector <double> noMissProb (this->nLoci_, log(1.0 - this->missCopyProb));
    vector <double> missProb (this->nLoci_, log(this->missCopyProb));
    vector <double> noNo = vecSum(noMissProb, noMissProb);
    vector <double> misMis = vecSum(missProb, missProb);
    vector <double> misNo = vecSum(noMissProb, missProb);

    vector <double> tmp_00_1 = vecSum(llk00_, noNo);
    vector <double> tmp_00_2 = vecSum(llk10_, misNo);
    vector <double> tmp_00_3 = vecSum(llk01_, misNo);
    vector <double> tmp_00_4 = vecSum(llk11_, misMis);

    vector <double> tmp_01_1 = vecSum(llk01_, noNo);
    vector <double> tmp_01_2 = vecSum(llk00_, misNo);
    vector <double> tmp_01_3 = vecSum(llk11_, misNo);
    vector <double> tmp_01_4 = vecSum(llk10_, misMis);

    vector <double> tmp_10_1 = vecSum(llk10_, noNo);
    vector <double> tmp_10_2 = vecSum(llk00_, misNo);
    vector <double> tmp_10_3 = vecSum(llk11_, misNo);
    vector <double> tmp_10_4 = vecSum(llk01_, misMis);

    vector <double> tmp_11_1 = vecSum(llk11_, noNo);
    vector <double> tmp_11_2 = vecSum(llk10_, misNo);
    vector <double> tmp_11_3 = vecSum(llk01_, misNo);
    vector <double> tmp_11_4 = vecSum(llk00_, misMis);

    //tmp.max = apply(cbind(tmp.00.1, tmp.00.2, tmp.00.3, tmp.00.4,
                          //tmp.10.1, tmp.10.2, tmp.10.3, tmp.10.4,
                          //tmp.01.1, tmp.01.2, tmp.01.3, tmp.01.4,
                          //tmp.11.1, tmp.11.2, tmp.11.3, tmp.11.4), 1, max)
    //emiss = cbind ( exp( tmp.00.1-tmp.max ) + exp( tmp.00.2-tmp.max ) + exp( tmp.00.3-tmp.max ) + exp( tmp.00.4-tmp.max ),
                    //exp( tmp.10.1-tmp.max ) + exp( tmp.10.2-tmp.max ) + exp( tmp.10.3-tmp.max ) + exp( tmp.10.4-tmp.max ),
                    //exp( tmp.01.1-tmp.max ) + exp( tmp.01.2-tmp.max ) + exp( tmp.01.3-tmp.max ) + exp( tmp.01.4-tmp.max ),
                    //exp( tmp.11.1-tmp.max ) + exp( tmp.11.2-tmp.max ) + exp( tmp.11.3-tmp.max ) + exp( tmp.11.4-tmp.max ))

    assert(emission_.size() == 0 );
    for ( size_t i = 0; i < this->nLoci_; i++){
        vector <double> tmp ({tmp_00_1[i], tmp_00_2[i], tmp_00_3[i], tmp_00_4[i],
                              tmp_01_1[i], tmp_01_2[i], tmp_01_3[i], tmp_01_4[i],
                              tmp_10_1[i], tmp_10_2[i], tmp_10_3[i], tmp_10_4[i],
                              tmp_11_1[i], tmp_11_2[i], tmp_11_3[i], tmp_11_4[i]});
        double tmaxTmp = maxOfVec(tmp);
        vector <double> emissRow ({exp(tmp_00_1[i] - tmaxTmp) + exp(tmp_00_2[i] - tmaxTmp) + exp(tmp_00_3[i] - tmaxTmp) + exp(tmp_00_4[i] - tmaxTmp),
                                   exp(tmp_01_1[i] - tmaxTmp) + exp(tmp_01_2[i] - tmaxTmp) + exp(tmp_01_3[i] - tmaxTmp) + exp(tmp_01_4[i] - tmaxTmp),
                                   exp(tmp_10_1[i] - tmaxTmp) + exp(tmp_10_2[i] - tmaxTmp) + exp(tmp_10_3[i] - tmaxTmp) + exp(tmp_10_4[i] - tmaxTmp),
                                   exp(tmp_11_1[i] - tmaxTmp) + exp(tmp_11_2[i] - tmaxTmp) + exp(tmp_11_3[i] - tmaxTmp) + exp(tmp_11_4[i] - tmaxTmp)});
        this->emission_.push_back(emissRow);
    }
}


vector <double> UpdatePairHap::computeRowMarginalDist( vector < vector < double > > & probDist ){
    vector <double> marginalDist (probDist.size(), 0.0);
    for ( size_t i = 0; i < probDist.size(); i++ ){
        marginalDist[i] = sumOfVec(probDist[i]);
    }
    return marginalDist;
}


vector <double> UpdatePairHap::computeColMarginalDist( vector < vector < double > > & probDist ){
    vector <double> marginalDist (probDist.size(), 0.0);
    for ( size_t coli = 0; coli < probDist[0].size(); coli++ ){
        for ( size_t rowi = 0; rowi < probDist.size(); rowi++ ){
            marginalDist[coli] += probDist[rowi][coli];
        }
    }
    assert ( sumOfVec ( marginalDist ) == 1 );
    return marginalDist;
}



void UpdatePairHap:: calcFwdProbs(){
    assert ( this->fwdProbs_.size() == 0 );
    vector < vector < double > > fwd1st;
    for ( size_t i = 0 ; i < this->nPanel_; i++){
        size_t rowObs = (size_t)this->panel_->content_[0][i];
        vector <double> fwd1stRow (this->nPanel_, 0.0);

        for ( size_t ii = 0 ; ii < this->nPanel_; ii++){
            if ( i == ii ) continue;

            size_t colObs = (size_t)this->panel_->content_[0][ii];
            size_t obs = rowObs*2 + colObs;
            fwd1stRow[ii] = this->emission_[0][obs];
        }
        fwd1st.push_back(fwd1stRow);
    }

    (void)normalizeBySumMat(fwd1st);
    this->fwdProbs_.push_back(fwd1st);

    for ( size_t j = 1; j < this->nLoci_; j++ ){
        double pRec = this->panel_->recombProbs_[j-1];
        double pRecEachHap = pRec / nPanel_;
        double pNoRec = 1.0 - pRec;

        double recRec = pRecEachHap * pRecEachHap;
        double recNorec = pRecEachHap * pNoRec;
        double norecNorec = pNoRec * pNoRec;

        //fwdList[[j]] =  ( rec.rec * sum(fwdList[[j-1]]) +
                           //norec.norec * fwdList[[j-1]] +
                           //rec.norec * marginalOfRows * marginalOfCols ) *
                         //fun.buildEmission( ref.panel[,j], emiss[j,] )
        //diag( fwdList[[j]] ) <- 0
        //fwdList[[j]] = fun.normalize.bymax(fwdList[[j]])

        //vector <double> marginalOfRows =


        vector < vector < double > > fwdTmp;
        for ( size_t i = 0 ; i < this->nPanel_; i++){
            size_t rowObs = (size_t)this->panel_->content_[0][i];
            vector <double> fwdTmpRow (this->nPanel_, 0.0);
            for ( size_t ii = 0 ; ii < this->nPanel_; ii++){
                if ( i == ii ) continue;

                size_t colObs = (size_t)this->panel_->content_[0][ii];
                size_t obs = rowObs*2 + colObs;
                // TODO, get row sum and column sum
                //fwdTmpRow[ii] = this->emission_[0][obs] * (sumOfMat(this->fwdProbs_.back())*recRec +
                                                           //fwdProbs_.back()[i][ii]*norecNorec+
                                                           //);
            }
            fwdTmp.push_back(fwdTmpRow);
        }
        //(void)normalizeBySumMat(fwdTmp);
        //this->fwdProbs_.push_back(fwdTmp);

            //double massFromRec = sumOfVec(fwdProbs_[j-1]) * pRecEachHap;
        //vector <double> fwdTmp (this->nPanel_, 0.0);
        //for ( size_t i = 0 ; i < this->nPanel_; i++){
            //fwdTmp[i] = this->emission_[j][this->panel_->content_[j][i]] * (fwdProbs_[j-1][i] * pNoRec + massFromRec);
        //}
        //(void)normalizeBySum(fwdTmp);
        //this->fwdProbs_.push_back(fwdTmp);
    }

}
void UpdatePairHap:: samplePaths(){}
void UpdatePairHap:: addMissCopying(){}

