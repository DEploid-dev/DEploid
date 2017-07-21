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
    (*writeTo) << "\n";
    (*writeTo) << "Input data: \n";
    (*writeTo) << setw(12) << "Panel: "     << panelFileName_  << "\n";
    (*writeTo) << setw(12) << "PLAF: "      << plafFileName_   << "\n";
    if ( useVcf() ) (*writeTo) << setw(12) << "VCF: " << vcfFileName_    << "\n";
    if ( refFileName_.size()>0) (*writeTo) << setw(12) << "REF count: " << refFileName_    << "\n";
    if ( altFileName_.size()>0) (*writeTo) << setw(12) << "ALT count: " << altFileName_    << "\n";
    if ( excludeSites() ){ (*writeTo) << setw(12) << "Exclude: " << excludeFileName_    << "\n"; }
    (*writeTo) << "\n";
    if ( this->doPainting() == false ) {
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
    (*writeTo) << setw(20) << " sigma: "               << this->parameterSigma() << "\n";
    (*writeTo) << setw(20) << " ScalingFactor: "    << this->scalingFactor() << "\n";
    if ( this->initialPropWasGiven() ){
        (*writeTo) << setw(20) << " Initial prob: " ;
        for ( size_t i = 0; i < this->initialProp.size(); i++ ){
            (*writeTo) << this->initialProp[i]
                       << ( ( i != (this->kStrain_-1) ) ? " " : "\n" );
        }
    }
    (*writeTo) << "\n";
    if ( this->doPainting() == false ) {
        (*writeTo) << "MCMC diagnostic:"<< "\n";
        (*writeTo) << setw(19) << " Accept_ratio: " << acceptRatio_ << "\n";
        (*writeTo) << setw(19) << " Max_llks: " << maxLLKs_ << "\n";
        (*writeTo) << setw(19) << " Mean_theta_llks: " << meanThetallks_ << "\n";
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
    (*writeTo) << "Output saved to:\n";
    if ( this->doPainting() ){
        for ( size_t i = 0; i < kStrain(); i++ ){
            (*writeTo) << "Posterior probability of strain " << i << ": "<< strExportSingleFwdProbPrefix << i <<endl;
        }
    } else {
        (*writeTo) << setw(14) << "Likelihood: "  << strExportLLK  << "\n";
        (*writeTo) << setw(14) << "Proportions: " << strExportProp << "\n";
        (*writeTo) << setw(14) << "Haplotypes: "  << strExportHap  << "\n";
        if ( doExportVcf() ) { (*writeTo) << setw(14) << "Vcf: "  << strExportVcf  << "\n"; }
        if (this->useIBD()){
            (*writeTo) << " IBD method output saved to:\n";
            (*writeTo) << setw(14) << "Likelihood: "  << strIbdExportProp  << "\n";
            (*writeTo) << setw(14) << "Proportions: " << strIbdExportProp << "\n";
            (*writeTo) << setw(14) << "Haplotypes: "  << strIbdExportHap  << "\n";
        }
    }
    (*writeTo) << "\n";
    (*writeTo) << "Proportions:\n";
    for ( size_t ii = 0; ii < this->filnalProp.size(); ii++){
        (*writeTo) << setw(10) << this->filnalProp[ii];
        (*writeTo) << ((ii < (this->filnalProp.size()-1)) ? "\t" : "\n") ;
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

