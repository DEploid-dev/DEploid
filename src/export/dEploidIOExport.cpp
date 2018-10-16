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

#include "dEploidIO.hpp"
#include "mcmc.hpp"


void DEploidIO::wrapUp(){
    this->writeRecombProb( panel );

    // Get End time before writing the log
    this->getTime(false);

    this->writeLog (&std::cout);

    ofstreamExportTmp.open( strExportLog.c_str(), ios::out | ios::app | ios::binary );
    this->writeLog (&ofstreamExportTmp);
    ofstreamExportTmp.close();
}


void DEploidIO::writeRecombProb ( Panel * panel ){
    if ( !doExportRecombProb() ) return;

    if ( panel != NULL ){
        ofstreamExportTmp.open( strExportRecombProb.c_str(), ios::out | ios::app | ios::binary );
        ofstreamExportTmp << "p.recomb"       << "\t"
                                 << "p.each"         << "\t"
                                 << "p.no.recomb"    << "\t"
                                 << "p.rec.rec"      << "\t"
                                 << "p.rec.norec"    << "\t"
                                 << "p.norec.norec"  << "\n";
        for ( size_t i = 0; i < panel->pRec_.size(); i++ ){
            ofstreamExportTmp << panel->pRec_[i]           << "\t"
                                     << panel->pRecEachHap_[i]    << "\t"
                                     << panel->pNoRec_[i]         << "\t"
                                     << panel->pRecRec_[i]        << "\t"
                                     << panel->pRecNoRec_[i]      << "\t"
                                     << panel->pNoRecNoRec_[i]    << "\n";
        }
        ofstreamExportTmp.close();
    }
}



void DEploidIO::writeLog ( ostream * writeTo ){
    size_t nHash = 30 + string(VERSION).size();
    for ( size_t i = 0; i < nHash; i++){
        (*writeTo) << "#";
    }
    (*writeTo) << "\n";
    (*writeTo) << "#        dEploid "<< setw(10) << VERSION << " log        #\n";
    for ( size_t i = 0; i < nHash; i++){
        (*writeTo) << "#";
    }
    (*writeTo) << "\n";
    (*writeTo) << "Program was compiled on: " << compileTime_ << endl;
    (*writeTo) << "dEploid version: " << dEploidGitVersion_ << endl;
    (*writeTo) << "lasso version: " << lassoGitVersion_ << endl;
    (*writeTo) << "\n";
    (*writeTo) << "Input data: \n";
    if (panelFileName_.size() > 0){
        (*writeTo) << setw(12) << "Panel: "     << panelFileName_  << "\n";
    }
    (*writeTo) << setw(12) << "PLAF: "      << plafFileName_   << "\n";
    if ( useVcf() ) (*writeTo) << setw(12) << "VCF: " << vcfFileName_    << "\n";
    if ( refFileName_.size()>0) (*writeTo) << setw(12) << "REF count: " << refFileName_    << "\n";
    if ( altFileName_.size()>0) (*writeTo) << setw(12) << "ALT count: " << altFileName_    << "\n";
    if ( excludeSites() ){ (*writeTo) << setw(12) << "Exclude: " << excludeFileName_    << "\n"; }
    (*writeTo) << "\n";
    if ( (this->doLsPainting() == false) & (this->doIbdPainting() == false) ) {
        (*writeTo) << "MCMC parameters: "<< "\n";
        (*writeTo) << setw(19) << " MCMC burn: " << mcmcBurn_ << "\n";
        (*writeTo) << setw(19) << " MCMC sample: " << nMcmcSample_ << "\n";
        (*writeTo) << setw(19) << " MCMC sample rate: " << mcmcMachineryRate_ <<"\n";
        (*writeTo) << setw(19) << " Random seed: " << this->randomSeed() << "\n";
        if (this->useIBD()){
            (*writeTo) << setw(19) << "  IBD Method used: YES" << "\n";
        }
        (*writeTo) << setw(19) << " Update Prop: "   << (this->doUpdateProp()  ? "YES":"NO") << "\n";
        (*writeTo) << setw(19) << " Update Single: " << (this->doUpdateSingle()? "YES":"NO") << "\n";
        (*writeTo) << setw(19) << " Update Pair: "   << (this->doUpdatePair()  ? "YES":"NO") << "\n";
        (*writeTo) << "\n";
    }
    (*writeTo) << "Other parameters:"<< "\n";
    if ( forbidCopyFromSame_ ){ (*writeTo) << " Update pair haplotypes move forbid copying from the same strain!!! \n"; }
    (*writeTo) << setw(20) << " Miss copy prob: "   << this->missCopyProb_ << "\n";
    (*writeTo) << setw(20) << " Avrg Cent Morgan: " << this->averageCentimorganDistance_ << "\n";
    (*writeTo) << setw(20) << " G: "               << this->parameterG() << "\n";
    if (this->useIBD()){
    (*writeTo) << setw(20) << " IBD sigma: "               << this->ibdSigma() << "\n";
    } else {
    (*writeTo) << setw(20) << " sigma: "               << this->parameterSigma() << "\n";
    }
    (*writeTo) << setw(20) << " ScalingFactor: "    << this->scalingFactor() << "\n";
    if ( this->initialPropWasGiven() ){
        (*writeTo) << setw(20) << " Initial prob: " ;
        for ( size_t i = 0; i < this->initialProp.size(); i++ ){
            (*writeTo) << this->initialProp[i]
                       << ( ( i != (this->kStrain_-1) ) ? " " : "\n" );
        }
    }
    (*writeTo) << "\n";
    if ((this->useLasso() == false) & (this->doLsPainting() == false) & (this->doIbdPainting() == false) & (this->doComputeLLK() == false) ) {
        (*writeTo) << "MCMC diagnostic:"<< "\n";
        (*writeTo) << setw(19) << " Accept_ratio: " << acceptRatio_ << "\n";
        (*writeTo) << setw(19) << " Max_llks: " << maxLLKs_ << "\n";
        (*writeTo) << setw(19) << " Final_theta_llks: " << meanThetallks_ << "\n";
        (*writeTo) << setw(19) << " Mean_llks: " << meanllks_ << "\n";
        (*writeTo) << setw(19) << " Stdv_llks: " << stdvllks_ << "\n";
        (*writeTo) << setw(19) << " DIC_by_Dtheta: " << dicByTheta_ << "\n";
        (*writeTo) << setw(19) << " DIC_by_varD: " << dicByVar_ << "\n";
        (*writeTo) << "\n";
    }
    (*writeTo) << "Run time:\n";
    (*writeTo) << setw(14) << "Start at: "  << startingTime_  ;
    (*writeTo) << setw(14) << "End at: "    << endTime_  ;
    (*writeTo) << "\n";
    if ( this->doComputeLLK() ){
        (*writeTo) << "Input likelihood: " << llkFromInitialHap_;
        (*writeTo) << "\n";
    } else {
        (*writeTo) << "Output saved to:\n";
        if ( this->doLsPainting() ){
            for ( size_t i = 0; i < kStrain(); i++ ){
                (*writeTo) << "Posterior probability of strain " << i << ": "<< strExportSingleFwdProbPrefix << i <<endl;
            }
        } else if (this->doIbdPainting()){
            if (this->ibdProbsIntegrated.size()>1){
                (*writeTo) << setw(14) << "IBD probs: "  << strIbdExportProbs  << "\n\n";
                (*writeTo) << " IBD probabilities:\n";
                for ( size_t stateI = 0; stateI < this->ibdProbsHeader.size(); stateI++ ){
                    (*writeTo) << setw(14) << this->ibdProbsHeader[stateI] << ": " << this->ibdProbsIntegrated[stateI] << "\n";
                }
            }
        } else {
            if (this->useLasso() == false) {
                (*writeTo) << setw(14) << "Likelihood: "  << strExportLLK  << "\n";
                (*writeTo) << setw(14) << "Proportions: " << strExportProp << "\n";
            }
            (*writeTo) << setw(14) << "Haplotypes: "  << strExportHap  << "\n";
            if ( doExportVcf() ) { (*writeTo) << setw(14) << "Vcf: "  << strExportVcf  << "\n"; }
            if (this->useIBD()){
                (*writeTo) << " IBD method output saved to:\n";
                (*writeTo) << setw(14) << "Likelihood: "  << strIbdExportLLK  << "\n";
                (*writeTo) << setw(14) << "Proportions: " << strIbdExportProp << "\n";
                (*writeTo) << setw(14) << "Haplotypes: "  << strIbdExportHap  << "\n";
            }
            if (this->ibdProbsIntegrated.size()>1){
                (*writeTo) << setw(14) << "IBD probs: "  << strIbdExportProbs  << "\n\n";
                (*writeTo) << " IBD probabilities:\n";
                for ( size_t stateI = 0; stateI < this->ibdProbsHeader.size(); stateI++ ){
                    (*writeTo) << setw(14) << this->ibdProbsHeader[stateI] << ": " << this->ibdProbsIntegrated[stateI] << "\n";
                }
            }
        }
        (*writeTo) << "\n";
        if (this->useLasso() == false) {
            (*writeTo) << " IBD best path llk: " << ibdLLK_ << "\n\n";
        }
        this->computeEffectiveKstrain(this->finalProp);
        (*writeTo) << "         Effective_K: " << this->effectiveKstrain_ <<"\n";
        this->computeInferredKstrain(this->finalProp);
        (*writeTo) << "          Inferred_K: " << this->inferredKstrain_ <<"\n";
        this->computeAdjustedEffectiveKstrain();
        (*writeTo) << "Adjusted_effective_K: " << this->adjustedEffectiveKstrain_ <<"\n";
    }
    (*writeTo) << "\n";
    (*writeTo) << "Proportions:\n";
    for ( size_t ii = 0; ii < this->finalProp.size(); ii++){
        (*writeTo) << setw(10) << this->finalProp[ii];
        (*writeTo) << ((ii < (this->finalProp.size()-1)) ? "\t" : "\n") ;
    }
}


void DEploidIO::writeEventCount(){
    ofstreamExportTmp.open( strExportExtra.c_str(), ios::out | ios::app | ios::binary );

    // HEADER
    ofstreamExportTmp << "CHROM" << "\t"
                      << "POS" << "\t"
                      << "IBDpathChangeAt" << "\t"
                      << "finalIBDpathChangeAt" << "\t"

                      << "siteOfTwoSwitchOne" << "\t"
                      << "finalSiteOfTwoSwitchOne" << "\t"

                      << "siteOfTwoMissCopyOne" << "\t"
                      << "finalSiteOfTwoMissCopyOne" << "\t"

                      << "siteOfTwoSwitchTwo" << "\t"
                      << "finalSiteOfTwoSwitchTwo" << "\t"

                      << "siteOfTwoMissCopyTwo" << "\t"
                      << "finalSiteOfTwoMissCopyTwo" << "\t"

                      << "siteOfOneSwitchOne" << "\t"
                      << "finalSiteOfOneSwitchOne" << "\t"

                      << "siteOfOneMissCopyOne" << "\t"
                      << "finalSiteOfOneMissCopyOne" << endl;

    size_t siteIndex = 0;
    for ( size_t chromI = 0; chromI < chrom_.size(); chromI++ ){
        for ( size_t posI = 0; posI < position_[chromI].size(); posI++){
            ofstreamExportTmp << chrom_[chromI] << "\t"
                              << (int)position_[chromI][posI] << "\t"

                              << this->IBDpathChangeAt[siteIndex] << "\t"
                              << this->finalIBDpathChangeAt[siteIndex] << "\t"

                              << this->siteOfTwoSwitchOne[siteIndex] << "\t"
                              << this->finalSiteOfTwoSwitchOne[siteIndex] << "\t"

                              << this->siteOfTwoMissCopyOne[siteIndex] << "\t"
                              << this->finalSiteOfTwoMissCopyOne[siteIndex] << "\t"

                              << this->siteOfTwoSwitchTwo[siteIndex] << "\t"
                              << this->finalSiteOfTwoSwitchTwo[siteIndex] << "\t"

                              << this->siteOfTwoMissCopyTwo[siteIndex] << "\t"
                              << this->finalSiteOfTwoMissCopyTwo[siteIndex] << "\t"

                              << this->siteOfOneSwitchOne[siteIndex] << "\t"
                              << this->finalSiteOfOneSwitchOne[siteIndex] << "\t"

                              << this->siteOfOneMissCopyOne[siteIndex] << "\t"
                              << this->finalSiteOfOneMissCopyOne[siteIndex] << endl;
            siteIndex++;
        }
    }

    assert(siteIndex == this->IBDpathChangeAt.size());
    ofstreamExportTmp.close();
}


void DEploidIO::writeIBDpostProb(vector < vector <double> > & reshapedProbs, vector <string> header){
    ostream * writeTo;
    #ifdef UNITTEST
        writeTo = &std::cout;
    #endif

    #ifndef UNITTEST
        ofstreamExportTmp.open( strIbdExportProbs.c_str(), ios::out | ios::app | ios::binary );
        writeTo = &ofstreamExportTmp;
    #endif

    (*writeTo) << "CHROM" << "\t" << "POS" << "\t";
    for (string tmp : header){
        (*writeTo) << tmp << ((tmp!=header[header.size()-1])?"\t":"\n");
    }

    size_t siteIndex = 0;
    for ( size_t chromIndex = 0; chromIndex < position_.size(); chromIndex++){
        for ( size_t posI = 0; posI < position_[chromIndex].size(); posI++){
            (*writeTo) << chrom_[chromIndex] << "\t" << (int)position_[chromIndex][posI] << "\t";
            for (size_t ij = 0; ij < reshapedProbs[siteIndex].size(); ij++){
                (*writeTo) << reshapedProbs[siteIndex][ij] << "\t";
            }
            (*writeTo) << endl;
            siteIndex++;
        }
    }
    assert(siteIndex == nLoci());

    #ifndef UNITTEST
        ofstreamExportTmp.close();
    #endif
}


void DEploidIO::paintIBD(){
    vector <double> goodProp;
    vector <size_t> goodStrainIdx;

    if ( this->doIbdPainting() ){
        this->finalProp = this->initialProp;
    }

    for ( size_t i = 0; i < this->finalProp.size(); i++){
        if (this->finalProp[i] > 0.01){
            goodProp.push_back(this->finalProp[i]);
            goodStrainIdx.push_back(i);
        }
    }

    if (goodProp.size() == 1){
        return;
    }

    DEploidIO tmpDEploidIO; // (*this);
    tmpDEploidIO.setKstrain(goodProp.size());
    tmpDEploidIO.setInitialPropWasGiven(true);
    tmpDEploidIO.initialProp = goodProp;
    tmpDEploidIO.finalProp = goodProp;
    tmpDEploidIO.refCount_ = this->refCount_;
    tmpDEploidIO.altCount_ = this->altCount_;
    tmpDEploidIO.plaf_ = this->plaf_;
    tmpDEploidIO.nLoci_= this->nLoci();
    tmpDEploidIO.position_ = this->position_;
    tmpDEploidIO.chrom_ = this->chrom_;
    //tmpDEploidIO.useConstRecomb_ = true;
    //tmpDEploidIO.constRecombProb_ = 0.000001;

    //tmpDEploidIO.writeLog (&std::cout);

    MersenneTwister tmpRg(this->randomSeed());
    IBDpath tmpIBDpath;
    tmpIBDpath.init(tmpDEploidIO, &tmpRg);
    tmpIBDpath.buildPathProbabilityForPainting(goodProp);
    this->ibdLLK_ = tmpIBDpath.bestPath(goodProp);
    this->ibdProbsHeader = tmpIBDpath.getIBDprobsHeader();
    this->getIBDprobsIntegrated(tmpIBDpath.fwdbwd);
    this->writeIBDpostProb(tmpIBDpath.fwdbwd, this->ibdProbsHeader);
}

