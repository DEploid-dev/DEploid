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

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <limits>
#include "ibd.hpp"

IBDconfiguration::IBDconfiguration() {}


void IBDconfiguration::buildIBDconfiguration(size_t k) {
    this->setKstrain(k);
    this->enumerateOp();
    this->makePairList();
    this->makePairToEmission();
    this->findUniqueState();
    this->findEffectiveK();
}


IBDconfiguration::~IBDconfiguration() {}


void IBDconfiguration::enumerateOp() {
    // #For each configuration, identify which pairs are IBD
    this->op = enumerateBinaryMatrixOfK(nchoose2(this->kStrain()));
}


void IBDconfiguration::makePairList() {
    // #Make a map of pairs to pair value
    assert(pairList.size() == 0);
    for ( size_t i = 0; i < this->kStrain(); i++ ) {  // 0-indexed
        for ( size_t j = i+1; j < this->kStrain(); j++ ) {
            pairList.push_back(vector <size_t> ({i, j}));
        }
    }
    assert(pairList.size() == (size_t)nchoose2(this->kStrain()));
}


void IBDconfiguration::makePairToEmission() {
    assert(pairToEmission.size() == 0);
    for ( vector<int> tmpOp : op ) {
        vector <int> tmpRow = makeTmpRow();
        // for ( size_t i = 0 ; i <(*opIt).size(); i++) {
            // cout << (*opIt)[i] << " ";
        // }
        // cout<<endl;

        vector <size_t> ii = findWhichIsOne(tmpOp);
        if (ii.size() > 0) {
            vector < vector <size_t> > tmpIBDPairs;
            for (size_t j = 0; j < ii.size(); j++) {
                tmpIBDPairs.push_back(pairList[ii[j]]);
            }

            int tmpIndex = (tmpIBDPairs.size()-1);
            while (tmpIndex >= 0) {
                tmpRow[tmpIBDPairs[tmpIndex][0]] =
                        tmpRow[tmpIBDPairs[tmpIndex][1]];
                tmpIndex--;
            }
        }
        pairToEmission.push_back(tmpRow);
    }
}


vector <int> IBDconfiguration::makeTmpRow() {
    vector <int> ret(this->kStrain());
    for (size_t i =0; i < ret.size(); i++) {
        ret[i] = static_cast<int>(i);
    }
    return ret;
}


vector <size_t> IBDconfiguration::findWhichIsOne(vector <int> tmpOp) {
    vector <size_t> ret;
    for (size_t i = 0; i < tmpOp.size(); i++) {
        if (tmpOp[i] == 1) {
            ret.push_back(i);
        }
    }
    return ret;
}


void IBDconfiguration::findUniqueState() {
    assert(states.size() == 0);
    states = unique(this->pairToEmission);
}


void IBDconfiguration::findEffectiveK() {
    assert(effectiveK.size() == 0);
    for (vector<int> state : states) {
        set <int> tmpSet(state.begin(), state.end());
        effectiveK.push_back(tmpSet.size());
    }
    assert(effectiveK.size() == states.size());
}


vector <string> IBDconfiguration::getIBDconfigureHeader() {
    vector <string> ret;
    for (size_t i = 0; i < this->states.size(); i++) {
        string tmp;
        for (size_t j = 0; j < this->states[i].size(); j++) {
            stringstream tmp_ss;
            tmp_ss << this->states[i][j];
            tmp += tmp_ss.str() + ((j < (this->states[i].size()-1)) ? "-":"");
        }
        ret.push_back(tmp);
    }
    return ret;
}


Hprior::Hprior() {}


void Hprior::buildHprior(size_t kStrain, const vector <double> &plaf) {
    ibdConfig.buildIBDconfiguration(kStrain);
    this->effectiveK = ibdConfig.effectiveK;
    this->nState_ = 0;
    this->plaf_ = plaf;
    this->setKstrain(kStrain);
    this->setnLoci(this->plaf_.size());
    vector < vector<int> > hSetBase = enumerateBinaryMatrixOfK(this->kStrain());
    size_t stateI = 0;
    for ( vector<int> state : ibdConfig.states ) {
        set <int> stateUnique(state.begin(), state.end());
        assert(stateUnique.size() == effectiveK[stateI]);
        vector < vector<int> > hSetBaseTmp = hSetBase;
        for (size_t j = 0; j < (size_t)this->kStrain(); j++) {
            for ( size_t i = 0; i < hSetBase.size(); i++ ) {
                hSetBaseTmp[i][j] = hSetBase[i][state[j]];
            }
        }

        vector < vector<int> > hSetBaseTmpUnique = unique(hSetBaseTmp);  // uu
        size_t sizeOfhSetBaseTmpUnique = hSetBaseTmpUnique.size();
        stateIdxFreq.push_back(sizeOfhSetBaseTmpUnique);

        for (size_t i = 0; i < sizeOfhSetBaseTmpUnique; i++) {
            int tmpSum = 0;
            for (int uniqSt : stateUnique) {
                tmpSum += hSetBaseTmpUnique[i][uniqSt];
            }
            // sumOfVec(hSetBaseTmpUnique[i]);
            int tmpDiff = stateUnique.size()-tmpSum;
            vector <double> hPriorTmp(nLoci());
            for (size_t site = 0; site < nLoci(); site++) {
                hPriorTmp[site] =
                        pow(plaf_[site], static_cast<double>(tmpSum))*
                        pow((1.0-plaf_[site]), static_cast<double>(tmpDiff));
            }
            priorProb.push_back(hPriorTmp);
            hSet.push_back(hSetBaseTmpUnique[i]);

            nState_++;
            stateIdx.push_back(stateI);
        }
        stateI++;
    }
}


void Hprior::transposePriorProbs() {
    assert(priorProbTrans.size() == 0);
    for ( size_t i = 0; i < nLoci(); i++ ) {
        vector <double> priorProbTransTmp(nState());
        for ( size_t j = 0; j < nState(); j++ ) {
            priorProbTransTmp[j] = priorProb[j][i];
        }
        priorProbTrans.push_back(priorProbTransTmp);
    }
}


vector <string> Hprior::getIBDconfigureHeader() {
    return this->ibdConfig.getIBDconfigureHeader();
}


Hprior::~Hprior() {}


vector < vector<int> > unique(const vector < vector<int> > &mat ) {
    vector < vector<int> > ret;
    ret.push_back(mat[0]);
    for (size_t i = 1; i < mat.size(); i++) {
        bool aNewState = true;
        for (vector<int> state : ret) {
            if ( twoVectorsAreSame(state, mat[i]) ) {
                aNewState = false;
                break;
            }
        }
        if ( aNewState ) {
            ret.push_back(mat[i]);
        }
    }

    return ret;
}


vector < vector<int> > enumerateBinaryMatrixOfK(size_t k) {
    // This function enumerate all possible binary combinations of k elements
    int ksq = pow(2, k);
    vector < vector<int> > ret;
    for (int i = 0; i < ksq; i++) {
        ret.push_back(convertIntToBinary(i, k));
    }
    return ret;
}


vector<int> convertIntToBinary(int x, size_t len) {
    vector<int> ret(len);
    size_t idx = 0;
    while (x) {
        ret[idx] = (x&1) ? 1:0;
        idx++;
        if ( idx > len ) {
            throw OutOfVectorSize();
        }
        x >>= 1;
    }
    reverse(ret.begin(), ret.end());
    return ret;
}


int nchoose2(int n) {
    if ( n < 2 ) {
        throw InvalidInput("Input must be at least 2!");
    }
    int ret = n*(n-1)/2;
    return ret;
}


bool twoVectorsAreSame(vector<int> vec1, vector<int> vec2) {
    if (vec1.size() != vec2.size()) {
        throw InvalidInput("Input vectors have different length!");
    }

    bool ret = true;
    for (size_t i = 0; i < vec1.size(); i++) {
        if (vec1[i] != vec2[i]) {
            ret = false;
            break;
        }
    }
    return ret;
}


IBDpath::IBDpath() {}


void IBDpath::init(const DEploidIO &dEploidIO, RandomGenerator* rg) {
    this->ibdRg_ = rg;
    this->setNLoci(dEploidIO.nLoci());
    this->setKstrain(dEploidIO.kStrain());
    this->setTheta(1.0 / static_cast<double>(kStrain()));

    this->IBDpathChangeAt = vector <double> (this->nLoci());

    // compute likelihood surface
    this->makeLlkSurf(dEploidIO.altCount_, dEploidIO.refCount_);

    // initialize haplotype prior
    this->hprior.buildHprior(kStrain(), dEploidIO.plaf_);
    this->hprior.transposePriorProbs();

    this->makeIbdTransProbs();

    // initialize fm
    this->fSumState = vector <double> (this->hprior.nPattern());

    // initialize ibdConfigurePath
    this->ibdConfigurePath = vector <size_t> (this->nLoci());

    // initialize recombination probabilities;
    this->ibdRecombProbs = IBDrecombProbs(dEploidIO.position_,
                                          dEploidIO.nLoci());
    this->ibdRecombProbs.computeRecombProbs(
                                        dEploidIO.averageCentimorganDistance(),
                                        dEploidIO.parameterG(),
                                        dEploidIO.useConstRecomb(),
                                        dEploidIO.constRecombProb());
    this->currentIBDpathChangeAt = vector <double> (this->nLoci());

    this->computeUniqueEffectiveKCount();
}


void IBDpath::ibdSamplePath(vector <double> statePrior) {
    int lociIdx = this->nLoci()-1;
    vector <double> tmpProp = fm[lociIdx];
    (void)normalizeBySum(tmpProp);
    ibdConfigurePath[lociIdx] = sampleIndexGivenProp(this->ibdRg_, tmpProp);

    assert(this->fm.size() == nLoci());
    while (lociIdx > 0) {
        lociIdx--;
        vector <double> vNoRecomb = vecProd(
        this->ibdTransProbs[this->hprior.stateIdx[ibdConfigurePath[lociIdx+1]]],
        fm[lociIdx]);
        assert(vNoRecomb.size() == this->hprior.nState());
        vector <double> vRecomb = fm[lociIdx];
        assert(vRecomb.size() == this->hprior.nState());
        vector <double> prop(this->hprior.nState());
        for (size_t i = 0; i < prop.size(); i++) {
            prop[i] = vNoRecomb[i]*this->ibdRecombProbs.pNoRec_[lociIdx] +
                        vRecomb[i]*this->ibdRecombProbs.pRec_[lociIdx] *
                        statePrior[ibdConfigurePath[lociIdx+1]];
        }
        tmpProp = prop;
        (void)normalizeBySum(tmpProp);
        ibdConfigurePath[lociIdx] = sampleIndexGivenProp(this->ibdRg_, tmpProp);
        assert(ibdConfigurePath[lociIdx] < this->hprior.nState());
        assert(ibdConfigurePath[lociIdx] >= 0);
    }
}


vector <size_t> IBDpath::findWhichIsSomething(vector <size_t> tmpOp,
                                              size_t something) {
    vector <size_t> ret;
    for (size_t i = 0; i < tmpOp.size(); i++) {
        if (tmpOp[i] == something) {
            ret.push_back(i);
        }
    }
    return ret;
}


void IBDpath::buildPathProbabilityForPainting(vector <double> proportion) {
    // vector <double> effectiveKPrior =
    //                    this->computeEffectiveKPrior(this->theta());
    vector <double> effectiveKPrior = vector <double> (this->hprior.nPattern(),
                                                1.0/this->hprior.nPattern());
    vector <double> statePrior = this->computeStatePrior(effectiveKPrior);
    // First building the path likelihood
    this->computeIbdPathFwdProb(proportion, statePrior);
    // Reshape Fwd
    vector < vector <double>> reshapedFwd = reshapeProbs(this->fm);

    this->computeIbdPathBwdProb(proportion, effectiveKPrior, statePrior);
    // Reshape Bwd
    vector < vector <double>> reshapedBwd = reshapeProbs(this->bwd);

    // Combine Fwd Bwd
    this->combineFwdBwd(reshapedFwd, reshapedBwd);
}



void IBDpath::computeIbdPathBwdProb(vector <double> proportion,
                                    vector <double> effectiveKPrior,
                                    vector <double> statePrior) {
    // # assuming each ibd state has equal probabilities,
    // # transform it into ibd configurations
    // dout << " start building ibd bwd "<< endl;
    vector <double> tmp = vector <double> (hprior.stateIdxFreq.size());
    assert(effectiveKPrior.size() == hprior.stateIdxFreq.size());
    for (size_t i = 0; i < tmp.size(); i++) {
        tmp[i] = effectiveKPrior[i] /
                   static_cast<double>(hprior.stateIdxFreq[i]);
    }

    vector <double> tmpBw = vector <double> (hprior.nState());
    for (size_t j = 0; j < tmpBw.size(); j++) {
        for (size_t i = 0; i < tmp.size(); i++) {
            tmpBw[j] += tmp[i] * ibdTransProbs[i][j];
        }
    }

    this->bwd.push_back(tmpBw);
    for ( size_t rev_siteI = 1; rev_siteI < this->nLoci(); rev_siteI++ ) {
        size_t siteI = this->nLoci()-rev_siteI;

        vector<double> lk = computeLlkOfStatesAtSiteI(proportion, siteI);
        vector <double> bSumState = vector <double> (hprior.nPattern());
        for (size_t i = 0; i < bSumState.size(); i++) {
            for (size_t j = 0; j < hprior.nState(); j++) {
                bSumState[i] += ibdTransProbs[i][j]*this->bwd.back()[j];
            }
        }
        vector <double> vNoRecomb(hprior.nState());
        for (size_t i = 0; i < hprior.stateIdx.size(); i++) {
            vNoRecomb[i] = bSumState[hprior.stateIdx[i]];
        }

        for (size_t i = 0; i < hprior.nState(); i++) {
            tmpBw[i] = 0;
            for (size_t j = 0; j < lk.size(); j++) {
                tmpBw[i] += (lk[j] * bwd.back()[j]) *
                            this->ibdRecombProbs.pRec_[siteI-1];
            }
            tmpBw[i] *= statePrior[i];
            tmpBw[i] += lk[i] * (this->ibdRecombProbs.pNoRec_[siteI-1]) *
                        vNoRecomb[i];
            tmpBw[i] *= hprior.priorProb[i][siteI];
        }
        normalizeBySum(tmpBw);
        this->bwd.push_back(tmpBw);
    }
    reverse(bwd.begin(), bwd.end());
}


void IBDpath::computeIbdPathFwdProb(vector <double> proportion,
                                    vector <double> statePrior) {
    this->fm.clear();
    vector <double> vPrior = vecProd(statePrior,
                                     this->hprior.priorProbTrans[0]);

    vector <double> lk = computeLlkOfStatesAtSiteI(proportion, 0);
    this->updateFmAtSiteI(vPrior, lk);
    for ( size_t siteI = 1; siteI < this->nLoci(); siteI++ ) {
        vector <double> vNoRec;
        for ( size_t stateIdxTmp : hprior.stateIdx ) {
            vNoRec.push_back(this->fSumState[stateIdxTmp]);
        }
        for ( size_t i = 0; i < hprior.nState(); i++ ) {
            vPrior[i] = (vNoRec[i] * this->ibdRecombProbs.pNoRec_[siteI] +
                         fSum * this->ibdRecombProbs.pRec_[siteI] *
                         statePrior[i]) *
                         hprior.priorProbTrans[siteI][i];
        }

        lk = computeLlkOfStatesAtSiteI(proportion, siteI);
        this->updateFmAtSiteI(vPrior, lk);
    }
}


void IBDpath::updateFmAtSiteI(const vector <double> & prior,
                              const vector <double> & llk) {
    vector <double> postAtSiteI = vecProd(prior, llk);
    normalizeBySum(postAtSiteI);
    this->fm.push_back(postAtSiteI);
    this->fSum = sumOfVec(postAtSiteI);
    for (size_t i = 0; i < fSumState.size(); i++) {
        this->fSumState[i] = 0;
        for ( size_t j = 0; j < hprior.nState(); j++ ) {
            this->fSumState[i] += ibdTransProbs[i][j]*postAtSiteI[j];
        }
    }
}


double IBDpath::bestPath(vector <double> proportion, double err) {
    double sumLLK = 0.0;
    for (size_t i = 0; i < nLoci(); i++) {
        vector <double> tmp;
        for (size_t j = 0; j < fm[i].size(); j++) {
            tmp.push_back(exp(log(fm[i][j])+log(bwd[i][j])));
        }
        normalizeBySum(tmp);
        size_t indx = distance(tmp.begin(),
                               max_element(tmp.begin(), tmp.end()));

        vector <int> hSetI = this->hprior.hSet[indx];
        double qs = 0;
        for (size_t j = 0; j < this->kStrain(); j++) {
            qs += static_cast<double>(hSetI[j]) * proportion[j];
        }
        double qs2 = qs*(1-err) + (1-qs)*err;

        if ( (qs > 0) & (qs < 1) ) {
            sumLLK += logBetaPdf(qs2, this->llkSurf[i][0], this->llkSurf[i][1]);
        }
    }
    return sumLLK;
}


void IBDpath::combineFwdBwd(const vector < vector <double>> &reshapedFwd,
                            const vector < vector <double>> &reshapedBwd) {
    for (size_t i = 0; i < nLoci(); i++) {
        vector <double> tmp;
        for (size_t j = 0; j < reshapedFwd[i].size(); j++) {
            tmp.push_back(exp(log(reshapedFwd[i][j])+log(reshapedBwd[i][j])));
        }
        normalizeBySum(tmp);
        fwdbwd.push_back(tmp);
    }
}


void IBDpath::makeIbdTransProbs() {
    assert(this->ibdTransProbs.size() == 0);
    for ( size_t i = 0; i < hprior.nPattern(); i++ ) {
        vector <double> transProbRow(hprior.nState());
        vector <size_t> wi = findWhichIsSomething(hprior.stateIdx, i);
        for (size_t wii : wi) {
            transProbRow[wii] = 1;
        }
        ibdTransProbs.push_back(transProbRow);
    }
}


vector <string> IBDpath::getIBDprobsHeader() {
    return this->hprior.getIBDconfigureHeader();
}


vector < vector <double> > IBDpath::reshapeProbs(
        const vector < vector <double> >& probs) {
    assert(this->nLoci() == probs.size());
    vector < vector <double> > ret;
    for (size_t siteIndex = 0; siteIndex < this->nLoci(); siteIndex++) {
        size_t previousStateIdx = 0;
        vector <double> tmpRow;
        double cumProb = 0;
        for (size_t prob_ij = 0; prob_ij < probs[siteIndex].size(); prob_ij++) {
            cumProb += probs[siteIndex][prob_ij];
            if (previousStateIdx != this->hprior.stateIdx[prob_ij]) {
                cumProb -= probs[siteIndex][prob_ij];
                previousStateIdx++;
                tmpRow.push_back(cumProb);
                cumProb = probs[siteIndex][prob_ij];
            }
        }
        tmpRow.push_back(cumProb);
        normalizeBySum(tmpRow);
        ret.push_back(tmpRow);
    }
    return ret;
}


vector <double> IBDpath::computeEffectiveKPrior(double theta) {
    // #Calculate state prior given theta (theta is prob IBD)
    vector <double> pr0(this->kStrain());
    for (int i = 0; i < static_cast<int>(pr0.size()); i++) {
        pr0[i] = binomialPdf(i, static_cast<int>(this->kStrain()-1), theta);
    }
    vector <double> effectiveKPrior;
    for (size_t effectiveKtmp : this->hprior.effectiveK) {
        int effectiveKidx = effectiveKtmp-1;
        assert(effectiveKidx >= 0);
        assert(effectiveKidx < static_cast<int>(this->kStrain()));
        effectiveKPrior.push_back(
            pr0[effectiveKidx]/uniqueEffectiveKCount[effectiveKidx]);
    }
    return effectiveKPrior;
}


vector <double> IBDpath::computeStatePrior(vector <double> effectiveKPrior) {
    vector <double> ret;
    for (size_t stateIdxTmp : this->hprior.stateIdx) {
        ret.push_back(effectiveKPrior[stateIdxTmp]);
    }
    return ret;
}


void IBDpath::computeAndUpdateTheta() {
    vector <size_t> obsState;
    size_t previousState = 0;
    size_t atSiteI = 0;
    for (size_t a : this->ibdConfigurePath) {
        if (a != previousState) {
            obsState.push_back(a);
        }
        if (this->hprior.stateIdx[a] != this->hprior.stateIdx[previousState]) {
            this->IBDpathChangeAt[atSiteI] += 1.0;
            this->currentIBDpathChangeAt[atSiteI] = 1.0;
        } else {
            this->currentIBDpathChangeAt[atSiteI] = 0.0;
        }
        previousState = a;
        atSiteI++;
    }

    size_t sumOfKeffStates = 0;
    size_t sccs = 0;
    for (size_t obs : obsState) {
        sumOfKeffStates += this->hprior.effectiveK[obs] - 1;
        sccs += this->kStrain() - this->hprior.effectiveK[obs];
    }
    this->setTheta(rBeta(sccs+1.0, sumOfKeffStates+1.0, this->ibdRg_));
}


void IBDpath::computeUniqueEffectiveKCount() {
    this->uniqueEffectiveKCount = vector <int> (this->kStrain());
    for (size_t effectiveKtmp : this->hprior.effectiveK) {
        int effectiveKidx = effectiveKtmp-1;
        assert(effectiveKidx >= 0);
        this->uniqueEffectiveKCount[effectiveKidx]++;
    }
}


void IBDpath::makeLlkSurf(vector <double> altCount, vector <double> refCount,
                          double scalingConst, double err, size_t gridSize) {
    double pGridSpacing = 1.0 / static_cast<double>(gridSize+1);
    vector <double> pGrid;
    pGrid.push_back(pGridSpacing);
    for (size_t i = 1; i < gridSize; i++) {
        pGrid.push_back(pGrid.back() + pGridSpacing);
    }
    assert(pGrid.size() == gridSize);

    assert(llkSurf.size() == 0);

    for (size_t i = 0 ; i < altCount.size(); i++) {
        double alt = altCount[i];
        double ref = refCount[i];

        vector <double> ll;
        for ( double unadjustedP : pGrid ) {
            ll.push_back(calcLLK(ref, alt, unadjustedP, err, scalingConst));
        }

        double llmax = max_value(ll);
        vector <double> ln;
        for ( double lltmp : ll ) {
            ln.push_back(exp(lltmp-llmax));
        }

        double lnSum = sumOfVec(ln);
        for (size_t i = 0; i < ln.size(); i++) {
            ln[i] = ln[i]/lnSum;
        }

        vector <double> tmpVec1 = vecProd(ln, pGrid);
        double mn = sumOfVec(tmpVec1);
        vector <double> pGridSq = vecProd(pGrid, pGrid);
        vector <double> tmpVec2 = vecProd(ln, pGridSq);
        double vr = sumOfVec(tmpVec2) - mn*mn;

        double comm = (mn*(1.0-mn)/vr-1.0);
        llkSurf.push_back(vector <double> {mn*comm, (1-mn)*comm});
    }
    assert(llkSurf.size() == this->nLoci());
}


vector <double> IBDpath::computeLlkOfStatesAtSiteI(vector<double> proportion,
                                                   size_t siteI, double err ) {
    vector <double> llks;
    for ( vector <int> hSetI : this->hprior.hSet ) {
        double qs = 0;
        for ( size_t j = 0; j < this->kStrain() ; j++ ) {
            qs += static_cast<double>(hSetI[j]) * proportion[j];
        }
        double qs2 = qs*(1-err) + (1-qs)*err;
        llks.push_back(
            logBetaPdf(qs2, this->llkSurf[siteI][0], this->llkSurf[siteI][1]));
    }

    double maxllk = max_value(llks);
    vector <double> ret;
    for ( double llk : llks ) {
        double normalized = exp(llk-maxllk);
        if ( normalized == 0 ) {
            // normalized = std::numeric_limits< double >::min();
            normalized = 2.22507e-308;
        }
        ret.push_back(normalized);
    }
    return ret;
}


IBDpath::~IBDpath() {}
