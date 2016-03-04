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
#include <cstdlib> // div
#include <valarray>

UpdateHap::UpdateHap( vector <double> &refCount,
                      vector <double> &altCount,
                      vector <double> &expectedWsaf,
                      vector <double> &proportion,
                      vector < vector <double> > &haplotypes,
                      MersenneTwister* rg,
                      size_t segmentStartIndex,
                      size_t nLoci,
                      Panel* panel){
    this->panel_ = panel;

    if ( this->panel_ != NULL ){
        this->nPanel_ = this->panel_->nPanel_;
    } else {
        this->nPanel_ = 0;
    }

    this->missCopyProb = 0.01;
    this->kStrain_ = proportion.size();

    this->rg_ = rg;

    this->segmentStartIndex_ = segmentStartIndex;
    this->nLoci_ = nLoci;

}


size_t UpdateHap::sampleIndexGivenProp ( vector <double> proportion ){
    vector <size_t> strainIndex = this->sampleNoReplace( proportion );
    assert( strainIndex.size() == 1);
    return strainIndex[0];
}


vector <size_t> UpdateHap::sampleNoReplace( vector <double> proportion, size_t nSample ){
    vector <size_t> indexReturn;
    assert( indexReturn.size() == 0 );
    vector <double> tmpDist(proportion) ;
    vector <size_t> tmpIndex;
    for ( size_t i = 0; i < proportion.size(); i++ ){
        tmpIndex.push_back(i);
    }
    for ( size_t nSampleRemaining = nSample; nSampleRemaining > 0; nSampleRemaining-- ){
        // Compute cdf of tmpDist
        vector <double> tmpCdf = computeCdf(tmpDist);
        double u = this->rg_->sample();
        size_t i = 0;
        for ( ; i < tmpCdf.size() ; i++){
            if ( u < tmpCdf[i] ){
                indexReturn.push_back(tmpIndex[i]);
                break;
            }
        }
        // Reduce tmpDist and tmpIndex
        tmpDist.erase(tmpDist.begin()+i);
        (void)normalizeBySum(tmpDist);
        tmpIndex.erase(tmpIndex.begin()+i);
    }
    return indexReturn;
}


UpdateSingleHap::UpdateSingleHap( vector <double> &refCount,
                                  vector <double> &altCount,
                                  vector <double> &expectedWsaf,
                                  vector <double> &proportion,
                                  vector < vector <double> > &haplotypes,
                                  MersenneTwister* rg,
                                  size_t segmentStartIndex,
                                  size_t nLoci,
                                  Panel* panel):
                UpdateHap(refCount, altCount, expectedWsaf, proportion, haplotypes, rg, segmentStartIndex, nLoci, panel){
    this->findUpdatingStrain( proportion );
    this->calcExpectedWsaf( expectedWsaf, proportion, haplotypes);
    this->calcHapLLKs(refCount, altCount);
    if ( this->panel_ != NULL ){
        this->buildEmission();
        this->calcFwdProbs();
        this->samplePaths();
        this->addMissCopying();
    } else {
        this->sampleHapIndependently();
    }
    this->updateLLK();
}


void UpdateSingleHap::findUpdatingStrain( vector <double> proportion ){
    this->strainIndex_ = sampleIndexGivenProp ( proportion );
    dout << "  Updating hap: "<< this->strainIndex_ <<endl;
}


void UpdateSingleHap::calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes ){
    //expected.WSAF.0 <- bundle$expected.WSAF - (bundle$prop[ws] * bundle$h[,ws]);
    this->expectedWsaf0_ = vector <double> (expectedWsaf.begin()+this->segmentStartIndex_, expectedWsaf.begin()+(this->segmentStartIndex_+this->nLoci_));
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
        size_t previous_site = j-1;
        double pRecEachHap = this->panel_->pRecEachHap_[previous_site];
        double pNoRec = this->panel_->pNoRec_[previous_site];

        double massFromRec = sumOfVec(fwdProbs_.back()) * pRecEachHap;
        vector <double> fwdTmp (this->nPanel_, 0.0);
        for ( size_t i = 0 ; i < this->nPanel_; i++){
            fwdTmp[i] = this->emission_[j][this->panel_->content_[j][i]] * (fwdProbs_.back()[i] * pNoRec + massFromRec);
        }
        (void)normalizeBySum(fwdTmp);
        this->fwdProbs_.push_back(fwdTmp);
    }
}


void UpdateSingleHap::calcHapLLKs( vector <double> &refCount,
                                   vector <double> &altCount){
    this->llk0_ = calcLLKs( refCount, altCount, expectedWsaf0_, this->segmentStartIndex_, this->nLoci_ );
    this->llk1_ = calcLLKs( refCount, altCount, expectedWsaf1_, this->segmentStartIndex_, this->nLoci_ );
}


void UpdateSingleHap::samplePaths(){
    assert ( this->path_.size() == 0 );
    size_t pathTmp = sampleIndexGivenProp ( fwdProbs_[nLoci_-1] );
    this->path_.push_back( this->panel_->content_.back()[pathTmp]);
    for ( size_t j = (this->nLoci_ - 1); j > 0; j--){
        size_t previous_site = j-1;
        double pRecEachHap = this->panel_->pRecEachHap_[previous_site];
        double pNoRec = this->panel_->pNoRec_[previous_site];

        vector <double> previousDist = fwdProbs_[previous_site];
        vector <double> weightOfNoRecAndRec ({ previousDist[pathTmp]*pNoRec,
                                               sumOfVec(previousDist)*pRecEachHap});
        (void)normalizeBySum(weightOfNoRecAndRec);
        size_t tmpState = sampleIndexGivenProp(weightOfNoRecAndRec);
        if ( tmpState == (size_t)1){
            pathTmp = sampleIndexGivenProp( previousDist );
        }

        this->path_.push_back(this->panel_->content_[j][pathTmp]);
    }
    reverse(path_.begin(), path_.end());
    assert(path_.size() == nLoci_);
}


void UpdateSingleHap::addMissCopying(){
    assert( this->hap_.size() == 0 );
    for ( size_t i = 0; i < this->nLoci_; i++){
        vector <double> emissionTmp ({exp(this->llk0_[i]), exp(this->llk1_[i])});
        vector <double> sameDiffDist ({emissionTmp[path_[i]]*(1.0 - this->missCopyProb), // probability of the same
                                       emissionTmp[(size_t)(1 -path_[i])] * this->missCopyProb }); // probability of differ
        (void)normalizeBySum(sameDiffDist);
        if ( sampleIndexGivenProp(sameDiffDist) == 1 ){
            this->hap_.push_back( 1 - this->path_[i] ); // differ
        } else {
            this->hap_.push_back( this->path_[i] ); // same
        }
    }
    assert ( this->hap_.size() == this->nLoci_ );
}


void UpdateSingleHap::sampleHapIndependently(){
    assert( this->hap_.size() == 0 );
    for ( size_t i = 0; i < this->nLoci_; i++){
        valarray <double> llkArray( {llk0_[i], llk1_[i]} );
        double tmpMax = llkArray.max();
        vector <double> tmpDist ( {exp(llk0_[i] - tmpMax),
                                   exp(llk1_[i] - tmpMax)} );
        (void)normalizeBySum(tmpDist);
        this->hap_.push_back ( (double)sampleIndexGivenProp(tmpDist) );
    }
    assert ( this->hap_.size() == this->nLoci_ );
}


void UpdateSingleHap::updateLLK(){
    newLLK = vector <double> (this->nLoci_, 0.0);
    for ( size_t i = 0; i < this->nLoci_; i++){
        if ( this->hap_[i] == 0){
            newLLK[i] = llk0_[i];
        } else if (this->hap_[i] == 1){
            newLLK[i] = llk1_[i];
        } else {
            throw("should never get here!");
        }
    }
}


UpdatePairHap::UpdatePairHap( vector <double> &refCount,
                              vector <double> &altCount,
                              vector <double> &expectedWsaf,
                              vector <double> &proportion,
                              vector < vector <double> > &haplotypes,
                              MersenneTwister* rg,
                              size_t segmentStartIndex,
                              size_t nLoci,
                              Panel* panel ):
                UpdateHap(refCount, altCount, expectedWsaf, proportion, haplotypes, rg, segmentStartIndex, nLoci, panel){
    this->findUpdatingStrain( proportion );
    this->calcExpectedWsaf( expectedWsaf, proportion, haplotypes);
    this->calcHapLLKs(refCount, altCount);
    if ( this->panel_ != NULL ){
        this->buildEmission();
        this->calcFwdProbs();
        this->samplePaths();
        this->addMissCopying();
    } else {
        this->sampleHapIndependently();
    }
    this->updateLLK();
}


void UpdatePairHap:: findUpdatingStrain( vector <double> proportion ){
    vector <size_t> strainIndex = sampleNoReplace( proportion, 2);
    assert( strainIndex.size() == 2);
    this->strainIndex1_ = strainIndex[0];
    this->strainIndex2_ = strainIndex[1];
    assert( strainIndex1_ != strainIndex2_ );
    dout << "  Updating hap: "<< this->strainIndex1_ << " and " << strainIndex2_ <<endl;

}


void UpdatePairHap:: calcExpectedWsaf( vector <double> & expectedWsaf, vector <double> &proportion, vector < vector <double> > &haplotypes){
  //expected.WSAF.00 <- expected.WSAF-(prop[ws[1]]*h[,ws[1]] + prop[ws[2]]*h[,ws[2]]);
  //expected.WSAF.10 <- expected.WSAF.00 + prop[ws[1]];
  //expected.WSAF.01 <- expected.WSAF.00 + prop[ws[2]];
  //expected.WSAF.11 <- expected.WSAF.00 + prop[ws[1]] + prop[ws[2]];    //expected.WSAF.0 <- bundle$expected.WSAF - (bundle$prop[ws] * bundle$h[,ws]);
    this->expectedWsaf00_ = vector <double> (expectedWsaf.begin()+this->segmentStartIndex_, expectedWsaf.begin()+(this->segmentStartIndex_+this->nLoci_));
    for ( size_t i = 0; i < expectedWsaf00_.size(); i++ ){
        expectedWsaf00_[i] -= (proportion[strainIndex1_] * haplotypes[i][strainIndex1_] + proportion[strainIndex2_] * haplotypes[i][strainIndex2_]);
        //dout << expectedWsaf[i] << " " << expectedWsaf00_[i] << endl;
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
    this->llk00_ = calcLLKs( refCount, altCount, expectedWsaf00_, this->segmentStartIndex_, this->nLoci_ );
    this->llk10_ = calcLLKs( refCount, altCount, expectedWsaf10_, this->segmentStartIndex_, this->nLoci_ );
    this->llk01_ = calcLLKs( refCount, altCount, expectedWsaf01_, this->segmentStartIndex_, this->nLoci_ );
    this->llk11_ = calcLLKs( refCount, altCount, expectedWsaf11_, this->segmentStartIndex_, this->nLoci_ );
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

    assert(this->emission_.size() == 0 );
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
    assert(this->emission_.size() == this->nLoci_ );
}


vector <double> UpdatePairHap::computeRowMarginalDist( vector < vector < double > > & probDist ){
    vector <double> marginalDist (probDist.size(), 0.0);
    for ( size_t i = 0; i < probDist.size(); i++ ){
        marginalDist[i] = sumOfVec(probDist[i]);
    }
    //assert ( sumOfVec (marginalDist) == sumOfMat(probDist));
    return marginalDist;
}


vector <double> UpdatePairHap::computeColMarginalDist( vector < vector < double > > & probDist ){
    vector <double> marginalDist (probDist.size(), 0.0);
    for ( size_t coli = 0; coli < probDist[0].size(); coli++ ){
        for ( size_t rowi = 0; rowi < probDist.size(); rowi++ ){
            marginalDist[coli] += probDist[rowi][coli];
        }
    }
    //assert ( sumOfVec ( marginalDist ) == 1.0 );
    //assert ( sumOfVec (marginalDist) == sumOfMat(probDist));
    return marginalDist;
}



void UpdatePairHap:: calcFwdProbs(){
    assert ( this->fwdProbs_.size() == 0 );
    vector < vector < double > > fwd1st;
    for ( size_t i = 0 ; i < this->nPanel_; i++){ // Row of the matrix
        size_t rowObs = (size_t)this->panel_->content_[0][i];
        vector <double> fwd1stRow (this->nPanel_, 0.0);

        for ( size_t ii = 0 ; ii < this->nPanel_; ii++){ // Column of the matrix
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
        size_t previous_site = j-1;
        double recRec = this->panel_->pRecRec_[previous_site];
        double recNorec = this->panel_->pRecNoRec_[previous_site];
        double norecNorec = this->panel_->pNoRecNoRec_[previous_site];

        vector <double> marginalOfRows = this->computeRowMarginalDist( this->fwdProbs_.back() );
        vector <double> marginalOfCols = this->computeColMarginalDist( this->fwdProbs_.back() );

        vector < vector < double > > fwdTmp;
        for ( size_t i = 0 ; i < this->nPanel_; i++){
            size_t rowObs = (size_t)this->panel_->content_[j][i];
            vector <double> fwdTmpRow (this->nPanel_, 0.0);
            for ( size_t ii = 0 ; ii < this->nPanel_; ii++){
                if ( i == ii ) continue;

                size_t colObs = (size_t)this->panel_->content_[j][ii];
                size_t obs = rowObs*2 + colObs;
                fwdTmpRow[ii] = this->emission_[j][obs] * (sumOfMat(this->fwdProbs_.back())*recRec +
                                                           fwdProbs_.back()[i][ii]*norecNorec+
                                                           recNorec *marginalOfRows[ii]*marginalOfCols[i]);
            }
            fwdTmp.push_back(fwdTmpRow);
        }
        (void)normalizeBySumMat(fwdTmp);
        this->fwdProbs_.push_back(fwdTmp);
    }
}


vector <double> UpdatePairHap::reshapeMatToVec ( vector < vector <double> > &Mat ){
    vector <double> tmp;
    for (auto const& array: Mat){
        for (auto const& value: array){
            tmp.push_back(value);
        }
    }
    return tmp;
}


vector <size_t> UpdatePairHap::sampleMatrixIndex( vector < vector < double > > &probDist ){
    size_t tmp = sampleIndexGivenProp ( reshapeMatToVec(probDist));
    div_t divresult;
    divresult = div((int)tmp, (int)this->nPanel_);
    //vector <size_t> returnVec({(size_t)divresult.quot, (size_t)divresult.rem});
    //return returnVec;
    return vector <size_t> ({(size_t)divresult.quot, (size_t)divresult.rem});
}

void UpdatePairHap::samplePaths(){
    assert ( this->path1_.size() == 0 );
    assert ( this->path2_.size() == 0 );

    vector <size_t> tmpPath = sampleMatrixIndex(fwdProbs_[nLoci_-1]);
    size_t rowI = tmpPath[0];
    size_t colJ = tmpPath[1];
    this->path1_.push_back(this->panel_->content_.back()[rowI]);
    this->path2_.push_back(this->panel_->content_.back()[colJ]);

    for ( size_t j = (this->nLoci_ - 1); j > 0; j--){
        size_t previous_site = j-1;
        double recRec = this->panel_->pRecRec_[previous_site];
        double recNorec = this->panel_->pRecNoRec_[previous_site];
        double norecNorec = this->panel_->pNoRecNoRec_[previous_site];

        vector < vector < double > > previousDist = fwdProbs_[previous_site];
        double previousProbij =previousDist[rowI][colJ];

        double tmpRowSum = sumOfVec(previousDist[rowI]) - previousProbij;
        double tmpColSum = 0.0;
        for ( size_t i = 0; i < previousDist.size(); i++){
            tmpColSum += previousDist[i][colJ];
        }
        tmpColSum -= previousProbij;

        vector <double> weightOfFourCases ({ recRec     * sumOfMat(previousDist),           // recombination happened on both strains
                                             recNorec   * tmpRowSum,  // first strain no recombine, second strain recombine
                                             recNorec   * tmpColSum,  // first strain recombine, second strain no recombine
                                             norecNorec * previousProbij }); // no recombine on either strain
        (void)normalizeBySum(weightOfFourCases);
        size_t tmpCase = sampleIndexGivenProp( weightOfFourCases );

        if ( tmpCase == (size_t)0 ){ // switching both strains
            tmpPath = sampleMatrixIndex(previousDist);
            rowI = tmpPath[0];
            colJ = tmpPath[1];
            assert (rowI != colJ);
            //switch.two = switch.two + 1
            //switch.table = rbind(switch.table, c("twoSwitchTwo", j ))
        } else if ( tmpCase == (size_t)1 ){ // switching second strain
            rowI = rowI;
            vector <double> rowIdist = previousDist[rowI];
            (void)normalizeBySum(rowIdist);
            colJ = sampleIndexGivenProp(rowIdist);
            assert (rowI != colJ);
            //switch.one = switch.one + 1
            //switch.table = rbind(switch.table, c("twoSwitchOne", j ))
        } else if ( tmpCase == (size_t)2 ){ // switching first strain
            vector <double> colJdist;
            for ( auto const& array: previousDist ){
                colJdist.push_back( array[colJ] );
            }
            assert(this->nPanel_ == colJdist.size());
            (void)normalizeBySum(colJdist);
            rowI = sampleIndexGivenProp(colJdist);
            colJ = colJ;
            assert (rowI != colJ);
            //switch.one = switch.one + 1
            //switch.table = rbind(switch.table, c("twoSwitchOne", j ))
        } else if ( tmpCase == (size_t)3 ) { // no switching
            rowI = rowI;
            colJ = colJ;
            assert (rowI != colJ);
        } else {
            throw ("Unknow case ... Should never reach here!");
        }
        this->path1_.push_back(this->panel_->content_[j][rowI]);
        this->path2_.push_back(this->panel_->content_[j][colJ]);

    }
    reverse(path1_.begin(), path1_.end());
    reverse(path2_.begin(), path2_.end());
    assert(path1_.size() == nLoci_);
    assert(path2_.size() == nLoci_);
}


void UpdatePairHap::addMissCopying(){
    assert( this->hap1_.size() == 0 );
    assert( this->hap2_.size() == 0 );

    for ( size_t i = 0; i < this->nLoci_; i++){
        vector <double> emissionTmp ({exp(this->llk00_[i]), exp(this->llk01_[i]), exp(this->llk10_[i]), exp(this->llk11_[i])});
        vector <double> casesDist ( { emissionTmp[(size_t)(2*path1_[i]     +path2_[i])]     * (1.0 - this->missCopyProb) * (1.0 - this->missCopyProb), // probability of both same
                                      emissionTmp[(size_t)(2*path1_[i]     +(1-path2_[i]))] * (1.0 - this->missCopyProb) * this->missCopyProb,         // probability of same1diff2
                                      emissionTmp[(size_t)(2*(1 -path1_[i])+path2_[i])]     * this->missCopyProb * (1.0 - this->missCopyProb),         // probability of same2diff1
                                      emissionTmp[(size_t)(2*(1 -path1_[i])+(1-path2_[i]))] * this->missCopyProb * this->missCopyProb });              // probability of both differ
        (void)normalizeBySum(casesDist);
        size_t tmpCase = sampleIndexGivenProp(casesDist);

        if ( tmpCase == 0 ){
            this->hap1_.push_back( this->path1_[i] );
            this->hap2_.push_back( this->path2_[i] );
        } else if ( tmpCase == 1 ){
            this->hap1_.push_back( this->path1_[i] );
            this->hap2_.push_back( 1.0 - this->path2_[i] );
        } else if ( tmpCase == 2 ){
            this->hap1_.push_back( 1.0 - this->path1_[i] );
            this->hap2_.push_back( this->path2_[i] );
        } else if ( tmpCase == 3 ){
            this->hap1_.push_back( 1.0 - this->path1_[i] );
            this->hap2_.push_back( 1.0 - this->path2_[i] );
        } else {
            throw ("add missing copy should never reach here" );
        }
    }

    assert ( this->hap1_.size() == this->nLoci_ );
    assert ( this->hap2_.size() == this->nLoci_ );
}


void UpdatePairHap::sampleHapIndependently(){
    assert( this->hap1_.size() == 0 );
    assert( this->hap2_.size() == 0 );

    for ( size_t i = 0; i < this->nLoci_; i++){
        valarray <double> llkArray( {llk00_[i], llk01_[i], llk10_[i], llk11_[i]} );
        double tmpMax = llkArray.max();
        vector <double> tmpDist ( {exp(llk00_[i] - tmpMax),
                                   exp(llk01_[i] - tmpMax),
                                   exp(llk10_[i] - tmpMax),
                                   exp(llk11_[i] - tmpMax) } );
        (void)normalizeBySum(tmpDist);

        size_t tmpCase = sampleIndexGivenProp(tmpDist);

        if ( tmpCase == 0 ){
            this->hap1_.push_back( 0.0 );
            this->hap2_.push_back( 0.0 );
        } else if ( tmpCase == 1 ){
            this->hap1_.push_back( 0.0 );
            this->hap2_.push_back( 1.0 );
        } else if ( tmpCase == 2 ){
            this->hap1_.push_back( 1.0 );
            this->hap2_.push_back( 0.0 );
        } else if ( tmpCase == 3 ){
            this->hap1_.push_back( 1.0 );
            this->hap2_.push_back( 1.0 );
        } else {
            throw ("add missing copy should never reach here" );
        }
    }

    assert ( this->hap1_.size() == this->nLoci_ );
    assert ( this->hap2_.size() == this->nLoci_ );
}



void UpdatePairHap::updateLLK(){
    newLLK = vector <double> (this->nLoci_, 0.0);
    for ( size_t i = 0; i < this->nLoci_; i++){
        if ( this->hap1_[i] == 0 && this->hap2_[i] == 0 ){
            newLLK[i] = llk00_[i];
        } else if (this->hap1_[i] == 0 && this->hap2_[i] == 1){
            newLLK[i] = llk01_[i];
        } else if (this->hap1_[i] == 1 && this->hap2_[i] == 0){
            newLLK[i] = llk10_[i];
        } else if (this->hap1_[i] == 1 && this->hap2_[i] == 1){
            newLLK[i] = llk11_[i];
        } else {
            throw("add missing copy, update llk should never reach here");
        }
    }
}
