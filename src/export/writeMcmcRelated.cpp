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

void DEploidIO::writeMcmcRelated (McmcSample * mcmcSample, bool useIBD){
    this->writeProp( mcmcSample, useIBD );
    this->writeLLK( mcmcSample, useIBD );
    this->writeHap( mcmcSample, useIBD);

    if ( useIBD == false ){
        this->writeVcf( mcmcSample );
        this->siteOfTwoSwitchOne = mcmcSample->siteOfTwoSwitchOne;
        this->siteOfTwoMissCopyOne = mcmcSample->siteOfTwoMissCopyOne;
        this->siteOfTwoSwitchTwo = mcmcSample->siteOfTwoSwitchTwo;
        this->siteOfTwoMissCopyTwo = mcmcSample->siteOfTwoMissCopyTwo;
        this->siteOfOneSwitchOne = mcmcSample->siteOfOneSwitchOne;
        this->siteOfOneMissCopyOne = mcmcSample->siteOfOneMissCopyOne;

        this->finalSiteOfTwoSwitchOne = mcmcSample->currentsiteOfTwoSwitchOne;
        this->finalSiteOfTwoMissCopyOne = mcmcSample->currentsiteOfTwoMissCopyOne;
        this->finalSiteOfTwoSwitchTwo = mcmcSample->currentsiteOfTwoSwitchTwo;
        this->finalSiteOfTwoMissCopyTwo = mcmcSample->currentsiteOfTwoMissCopyTwo;
        this->finalSiteOfOneSwitchOne = mcmcSample->currentsiteOfOneSwitchOne;
        this->finalSiteOfOneMissCopyOne = mcmcSample->currentsiteOfOneMissCopyOne;

        //this->writeEventCount( );
    } else {
        this->IBDpathChangeAt = mcmcSample->IBDpathChangeAt;
        this->finalIBDpathChangeAt = mcmcSample->currentIBDpathChangeAt;
    }
}


void DEploidIO::writeProp( McmcSample * mcmcSample, bool useIBD){
    if ( useIBD ){
        ofstreamExportTmp.open( strIbdExportProp.c_str(), ios::out | ios::app | ios::binary );
    } else {
        ofstreamExportTmp.open( strExportProp.c_str(), ios::out | ios::app | ios::binary );
    }
    for ( size_t i = 0; i < mcmcSample->proportion.size(); i++){
        for ( size_t ii = 0; ii < mcmcSample->proportion[i].size(); ii++){
            ofstreamExportTmp << setw(10) << mcmcSample->proportion[i][ii];
            ofstreamExportTmp << ((ii < (mcmcSample->proportion[i].size()-1)) ? "\t" : "\n") ;
        }
    }
    ofstreamExportTmp.close();
}


void DEploidIO::writeLLK( McmcSample * mcmcSample, bool useIBD){
    if ( useIBD ){
        ofstreamExportTmp.open( strIbdExportLLK.c_str(), ios::out | ios::app | ios::binary );
    } else {
        ofstreamExportTmp.open( strExportLLK.c_str(), ios::out | ios::app | ios::binary );
    }
    for ( size_t i = 0; i < mcmcSample->sumLLKs.size(); i++){
        ofstreamExportTmp << mcmcSample->moves[i] << "\t" << mcmcSample->sumLLKs[i] << endl;
    }
    ofstreamExportTmp.close();
}


void DEploidIO::writeHap( McmcSample * mcmcSample, bool useIBD){
    if ( useIBD ){
        ofstreamExportTmp.open( strIbdExportHap.c_str(), ios::out | ios::app | ios::binary );
    } else {
        ofstreamExportTmp.open( strExportHap.c_str(), ios::out | ios::app | ios::binary );
    }
    // HEADER
    ofstreamExportTmp << "CHROM" << "\t" << "POS" << "\t";;
    for ( size_t ii = 0; ii < kStrain_; ii++){
        ofstreamExportTmp << "h" << (ii+1) ;
        ofstreamExportTmp << ((ii < (kStrain_-1)) ? "\t" : "\n") ;
    }

    size_t siteIndex = 0;
    for ( size_t chromI = 0; chromI < chrom_.size(); chromI++ ){
        for ( size_t posI = 0; posI < position_[chromI].size(); posI++){
            ofstreamExportTmp << chrom_[chromI] << "\t" << (int)position_[chromI][posI] << "\t";
            for ( size_t ii = 0; ii < mcmcSample->hap[siteIndex].size(); ii++){
                ofstreamExportTmp << mcmcSample->hap[siteIndex][ii];
                ofstreamExportTmp << ((ii < (mcmcSample->hap[siteIndex].size()-1)) ? "\t" : "\n") ;
            }
            siteIndex++;
        }
    }

    assert ( siteIndex == mcmcSample->hap.size());
    ofstreamExportTmp.close();
}


void DEploidIO::writeVcf( McmcSample * mcmcSample ){
    if ( !doExportVcf() ) return;

    ogzstream ogstreamExport;
    ostream * writeTo;
    if ( compressVcf() ){
        ogstreamExport.open( strExportVcf.c_str(), ios::out );
        writeTo = &ogstreamExport;
    } else {
        ofstreamExportTmp.open( strExportVcf.c_str(), ios::out | ios::app | ios::binary );
        writeTo = &ofstreamExportTmp;
    }

    // VCF HEADER
    if ( this->useVcf() ){
        for ( auto const& headerLine: this->vcfReaderPtr_->headerLines){
            (*writeTo) << headerLine << endl;
        }
    } else {
        (*writeTo) << "##fileformat=VCFv4.2" << endl;
    }
    // Include proportions
    for ( size_t ii = 0; ii < kStrain_; ii++){
        (*writeTo) << "##Proportion of strain "
                   << ( this->useVcf() ? this->vcfReaderPtr_->sampleName : "h" )
                   << "." << (ii+1)
                   << "=" << mcmcSample->proportion.back()[ii] << endl;
    }

    // HEADER
    (*writeTo) << "#CHROM" << "\t"
               << "POS"    << "\t"
               << "ID"     << "\t"
               << "REF"    << "\t"
               << "ALT"    << "\t"
               << "QUAL"   << "\t"
               << "FILTER" << "\t"
               << "INFO"   << "\t"
               << "FORMAT" << "\t";
    for ( size_t ii = 0; ii < kStrain_; ii++){
        (*writeTo) << ( this->useVcf() ? this->vcfReaderPtr_->sampleName : "h" )
                          << "." << (ii+1) ;
        (*writeTo) << ((ii < (kStrain_-1)) ? "\t" : "\n") ;
    }

    size_t siteIndex = 0;
    for ( size_t chromI = 0; chromI < chrom_.size(); chromI++ ){
        for ( size_t posI = 0; posI < position_[chromI].size(); posI++){
            if ( useVcf() ) {
                (*writeTo) << this->vcfReaderPtr_->variants[siteIndex].chromStr  << "\t"
                           << this->vcfReaderPtr_->variants[siteIndex].posStr    << "\t"
                           << this->vcfReaderPtr_->variants[siteIndex].idStr     << "\t"
                           << this->vcfReaderPtr_->variants[siteIndex].refStr    << "\t"
                           << this->vcfReaderPtr_->variants[siteIndex].altStr    << "\t"
                           << this->vcfReaderPtr_->variants[siteIndex].qualStr   << "\t"
                           << this->vcfReaderPtr_->variants[siteIndex].filterStr << "\t"
                           << this->vcfReaderPtr_->variants[siteIndex].infoStr   << "\t"
                           << "GT"                                      << "\t";
            } else {
                (*writeTo) << chrom_[chromI]               << "\t"
                           << (int)position_[chromI][posI] << "\t"
                           << "."                          << "\t"
                           << "."                          << "\t"
                           << "."                          << "\t"
                           << "."                          << "\t"
                           << "."                          << "\t"
                           << "."                          << "\t"
                           << "GT"                         << "\t";
            }

            for ( size_t ii = 0; ii < mcmcSample->hap[siteIndex].size(); ii++){
                (*writeTo) << mcmcSample->hap[siteIndex][ii];
                (*writeTo) << ((ii < (mcmcSample->hap[siteIndex].size()-1)) ? "\t" : "\n") ;
            }
            siteIndex++;
        }
    }

    assert ( siteIndex == mcmcSample->hap.size());
    if ( compressVcf() ){
        ogstreamExport.close();
    } else {
        ofstreamExportTmp.close();
    }
}

