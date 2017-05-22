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
#include "ibd.hpp"
#include <algorithm>
#include <cassert>

IBDconfiguration::IBDconfiguration(){}


void IBDconfiguration::buildIBDconfiguration(size_t k){
    this->setKstrain(k);
    this->enumerateOp();
    this->makePairList();
    this->makePairToEmission();
    this->findUniqueState();
    this->findEffectiveK();
}


IBDconfiguration::~IBDconfiguration(){}


void IBDconfiguration::enumerateOp(){
    //#For each configuration, identify which pairs are IBD
    this->op = enumerateBinaryMatrixOfK(nchoose2(this->kStrain()));
}


void IBDconfiguration::makePairList(){
    //#Make a map of pairs to pair value
    assert(pairList.size() == 0);
    for ( size_t i = 0; i < this->kStrain(); i++ ){ // 0-indexed
        for ( size_t j = i+1; j < this->kStrain(); j++ ){
            pairList.push_back(vector <size_t> ({i, j}));
        }
    }
    assert(pairList.size() == (size_t)nchoose2(this->kStrain()));
}


void IBDconfiguration::makePairToEmission(){
    assert(pairToEmission.size() == 0);
    //prs2ems<-array(0, c(nrow(op), k));
    //for ( size_t i = 0; i < this->op.size(); i++ ){
    for ( vector<int> tmpOp : op ){
        vector <int> tmpRow = makeTmpRow();
        //for ( size_t i = 0 ; i <(*opIt).size(); i++){
            //cout << (*opIt)[i] << " ";
        //}
        //cout<<endl;

        vector <size_t> ii = findWhichIsOne(tmpOp);
        //ii <- which(op[rowI,]==1);
        //cout << ii.size()<<endl;
        //cout << "##############" <<endl;
        if (ii.size() > 0){
            vector < vector <size_t> > tmpIBDPairs;
            for ( size_t j = 0; j < ii.size(); j++){
                //cout << "j = " << j <<" ii[j] = "<<ii[j]<<endl;
                tmpIBDPairs.push_back(pairList[ii[j]]);

                //cout << tmpIBDPairs.back()[0]<< " "<<tmpIBDPairs.back()[1]<<endl;
            }

            int tmpIndex = (tmpIBDPairs.size()-1);
            while (tmpIndex >= 0){
                //cout << "replacing element "<< tmpIBDPairs[tmpIndex][0] << " by element " << tmpIBDPairs[tmpIndex][1] <<endl;
                tmpRow[tmpIBDPairs[tmpIndex][0]] = tmpRow[tmpIBDPairs[tmpIndex][1]];
                tmpIndex--;
            }
        }
        pairToEmission.push_back(tmpRow);
        //for (size_t i = 0 ; i < pairToEmission.back().size(); i++){
            //cout << pairToEmission.back()[i]<<" ";
        //}
        //cout <<endl;
    }
    //for (size_t i = 0; i < pairToEmission.size(); i++){
        //for (size_t ii = 0 ; ii < pairToEmission.back().size(); ii++){
            //cout << pairToEmission[i][ii]<<" ";
        //}
        //cout <<endl;
    //}
}


vector <int> IBDconfiguration::makeTmpRow(){
    vector <int> ret(this->kStrain());
    for (size_t i =0; i < ret.size(); i++){
        ret[i] = (int)i;
    }
    return ret;
}


vector <size_t> IBDconfiguration::findWhichIsOne(vector <int> tmpOp){
    vector <size_t> ret;
    for ( size_t i = 0; i < tmpOp.size(); i++){
        if ( tmpOp[i] == 1){
            ret.push_back(i);
        }
    }
    return ret;
}


void IBDconfiguration::findUniqueState(){
    assert (states.size() == 0);
    //states.push_back(pairToEmission[0]);
    //for (size_t i = 1; i < this->pairToEmission.size(); i++){
        //bool aNewState = true;
        //for ( vector<int> state : states){
            //if ( twoVectorsAreSame(state, this->pairToEmission[i]) ){
                //aNewState = false;
                //break;
            //}
        //}
        //if ( aNewState ){
            //states.push_back(this->pairToEmission[i]);
        //}
    //}
    states = unique(this->pairToEmission);
    //for ( vector<int> state : states){
        //for (int i : state){
            //cout << i <<" ";
        //}
        //cout <<endl;
    //}
}


void IBDconfiguration::findEffectiveK(){
    assert(effectiveK.size() == 0);
    for ( vector<int> state : states){
        set <int> tmpSet (state.begin(), state.end());
        //cout << tmpSet.size() <<endl;
        effectiveK.push_back(tmpSet.size());
    }
    assert(effectiveK.size() == states.size());
}


Hprior::Hprior(){}


void Hprior::buildHprior(size_t kStrain, vector <double> &plaf){
    ibdConfig.buildIBDconfiguration(kStrain);
    this->effectiveK = ibdConfig.effectiveK;
    this->nState_ = 0;
    this->plaf_ = plaf;
    this->setKstrain(kStrain);
    this->setnLoci(this->plaf_.size());
    vector < vector<int> > hSetBase = enumerateBinaryMatrixOfK(this->kStrain());

    size_t stateI = 0;
    for ( vector<int> state : ibdConfig.states ) {
        set <int> stateUnique (state.begin(), state.end());
        assert(stateUnique.size() == effectiveK[stateI]);
        vector < vector<int> > hSetBaseTmp = hSetBase;
        for (size_t j = 0; j < (size_t)this->kStrain(); j++){
            for ( size_t i = 0; i < hSetBase.size(); i++ ){
                hSetBaseTmp[i][j] = hSetBase[i][state[j]];
            }
        }

        vector < vector<int> > hSetBaseTmpUnique = unique(hSetBaseTmp); // uu
        size_t sizeOfhSetBaseTmpUnique = hSetBaseTmpUnique.size();
        //h.prior.i<-array(0, c(size.h.set.i, n.loci));

        for ( size_t i = 0; i < sizeOfhSetBaseTmpUnique; i++){
            //vector<int> hSetBaseTmpUniqueSubSet(); // uu[i,a.u,drop=F]
            int tmpSum = 0;
            for (int uniqSt : stateUnique){
                tmpSum += hSetBaseTmpUnique[i][uniqSt];
            }
            // sumOfVec(hSetBaseTmpUnique[i]);
            int tmpDiff = stateUnique.size()-tmpSum;
            //cout << stateUnique.size() << " " << tmpSum << " " << tmpDiff<< endl;
            vector <double> hPriorTmp (nLoci());
            for (size_t site = 0; site < nLoci(); site++ ){
                //cout << (1.0-plaf_[site]) << " " << tmpDiff << " " <<pow((1.0-plaf_[site]),(double)tmpDiff) << endl;
                hPriorTmp[site] = pow(plaf_[site], (double)tmpSum) * pow((1.0-plaf_[site]),(double)tmpDiff);
            }
            priorProb.push_back(hPriorTmp);
            hSet.push_back(hSetBaseTmpUnique[i]);

            nState_++;
            stateIdx.push_back(stateI);
        }
        stateI++;
    }
}


void Hprior::transposePriorProbs(){
    assert(priorProbTrans.size() == 0);
    for ( size_t i = 0; i < nLoci(); i++ ){
        vector <double> priorProbTransTmp(nState());
        for ( size_t j = 0; j < nState(); j++ ){
            priorProbTransTmp[j] = priorProb[j][i];
            //cout << priorProb[j][i] << " ";
        }
        priorProbTrans.push_back(priorProbTransTmp);
        //cout << endl;
    }
}



Hprior::~Hprior(){}


vector < vector<int> > unique( vector < vector<int> > &mat ){
    vector < vector<int> > ret;
    ret.push_back(mat[0]);
    for (size_t i = 1; i < mat.size(); i++){
        bool aNewState = true;
        for ( vector<int> state : ret){
            if ( twoVectorsAreSame(state, mat[i]) ){
                aNewState = false;
                break;
            }
        }
        if ( aNewState ){
            ret.push_back(mat[i]);
        }
    }

    return ret;
}


vector < vector<int> > enumerateBinaryMatrixOfK( size_t k ){
    // This function enumerate all possible binary combinations of k elements
    int ksq = pow(2,k);
    vector < vector<int> > ret;
    for ( int i = 0; i < ksq; i++){
        ret.push_back(convertIntToBinary(i, k));
    }
    return ret;
}


vector<int> convertIntToBinary(int x, size_t len) {
    vector<int> ret(len);
    size_t idx = 0;
    while(x) {
        ret[idx] = (x&1) ? 1:0;
        idx++;
        //cout << "x " <<x<< " idx "<<idx<<" len "<< len<<endl;
        if ( idx > len ){
            throw OutOfVectorSize();
        }
        x >>= 1;
    }
    reverse(ret.begin(),ret.end());
    //for (size_t i = 0; i < ret.size(); i++){
            //cout << ret[i] << " ";
        //}
        //cout<<endl;
    return ret;
}


int nchoose2(int n){
    if ( n < 2 ){
        throw InvalidInput("Input must be at least 2!");
    }
    int ret = n*(n-1)/2;
    return ret;
}


bool twoVectorsAreSame(vector<int> vec1, vector<int> vec2){
    if (vec1.size() != vec2.size()){
        throw InvalidInput("Input vectors have different length!");
    }

    bool ret = true;
    for (size_t i = 0; i < vec1.size(); i++){
        if (vec1[i] != vec2[i]){
            ret = false;
            break;
        }
    }
    return ret;
}

