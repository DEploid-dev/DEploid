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

#include "pfDeconvIO.hpp"
#include "updateHap.hpp"
#include "mcmc.hpp"

void McmcMachinery::writeLastFwdProb(){
    if ( this->pfDeconvIO_ ->doExportPostProb() != true ){
        return;
    }

    for ( size_t chromi = 0 ; chromi < this->pfDeconvIO_->indexOfChromStarts_.size(); chromi++ ){
        size_t start = this->pfDeconvIO_->indexOfChromStarts_[chromi];
        size_t length = this->pfDeconvIO_->position_[chromi].size();

        UpdateSingleHap updating0( this->pfDeconvIO_->refCount_,
                                  this->pfDeconvIO_->altCount_,
                                  this->pfDeconvIO_->plaf_,
                                  this->currentExpectedWsaf_,
                                  this->currentProp_, this->currentHap_, this->hapRg_,
                                  start, length,
                                  this->panel_, this->pfDeconvIO_->missCopyProb_,
                                  (size_t)0);
        updating0.core ( this->pfDeconvIO_->refCount_, this->pfDeconvIO_->altCount_, this->pfDeconvIO_->plaf_, this->currentExpectedWsaf_, this->currentProp_, this->currentHap_);
        this->pfDeconvIO_->writeLastSingleFwdProb( updating0, chromi, (size_t)0 );

        UpdateSingleHap updating1( this->pfDeconvIO_->refCount_,
                                  this->pfDeconvIO_->altCount_,
                                  this->pfDeconvIO_->plaf_,
                                  this->currentExpectedWsaf_,
                                  this->currentProp_, this->currentHap_, this->hapRg_,
                                  start, length,
                                  this->panel_, this->pfDeconvIO_->missCopyProb_,
                                  (size_t)1);
        updating1.core ( this->pfDeconvIO_->refCount_, this->pfDeconvIO_->altCount_, this->pfDeconvIO_->plaf_, this->currentExpectedWsaf_, this->currentProp_, this->currentHap_);
        this->pfDeconvIO_->writeLastSingleFwdProb( updating1, chromi, (size_t)1 );

        UpdatePairHap updating( this->pfDeconvIO_->refCount_,
                                this->pfDeconvIO_->altCount_,
                                this->pfDeconvIO_->plaf_,
                                this->currentExpectedWsaf_,
                                this->currentProp_, this->currentHap_, this->hapRg_,
                                start, length,
                                this->panel_, this->pfDeconvIO_->missCopyProb_, this->pfDeconvIO_->forbidCopyFromSame(),
                                (size_t)0,
                                (size_t)1);
        updating.core ( this->pfDeconvIO_->refCount_, this->pfDeconvIO_->altCount_, this->pfDeconvIO_->plaf_, this->currentExpectedWsaf_, this->currentProp_, this->currentHap_);
        this->pfDeconvIO_->writeLastPairFwdProb( updating, chromi );
    }
}


void PfDeconvIO::writeLastSingleFwdProb( UpdateSingleHap & updateSingle, size_t chromIndex, size_t strainIndex ){
    string strExportFwdProb = (strainIndex == 0) ? strExportSingleFwdProb0 : strExportSingleFwdProb1;
    ofstreamExportFwdProb.open( strExportFwdProb.c_str(), ios::out | ios::app | ios::binary );

    if ( chromIndex == 0 ){ // Print header
        ofstreamExportFwdProb << "CHROM" << "\t" << "POS" << "\t";;
        for ( size_t ii = 0; ii < updateSingle.fwdProbs_[0].size(); ii++){
            ofstreamExportFwdProb << (ii+1) ;
            ofstreamExportFwdProb << ((ii < (updateSingle.fwdProbs_[0].size()-1)) ? "\t" : "\n") ;
        }
    }

    size_t siteIndex = 0;
    for ( size_t posI = 0; posI < position_[chromIndex].size(); posI++){
        ofstreamExportFwdProb << chrom_[chromIndex] << "\t" << (int)position_[chromIndex][posI] << "\t";
        for ( size_t ii = 0; ii < updateSingle.fwdProbs_[siteIndex].size(); ii++){
            ofstreamExportFwdProb << updateSingle.fwdProbs_[siteIndex][ii];
            ofstreamExportFwdProb << ((ii < (updateSingle.fwdProbs_[siteIndex].size()-1)) ? "\t" : "\n") ;
        }
        siteIndex++;
    }

    ofstreamExportFwdProb.close();
}


void PfDeconvIO::writeLastPairFwdProb( UpdatePairHap & updatePair, size_t chromIndex ){
    cout << "starts printing writeLastPairFwdProb " << endl;
    ofstreamExportFwdProb.open( strExportPairFwdProb.c_str(), ios::out | ios::app | ios::binary );
cout<<"prod = "<<updatePair.fwdProbs_[0].size()*updatePair.fwdProbs_[0][3].size()<<endl;
    if ( chromIndex == 0 ){ // Print header
        ofstreamExportFwdProb << "CHROM" << "\t" << "POS" << "\t";;
        for ( size_t ii = 0; ii < updatePair.fwdProbs_[0].size(); ii++){
            for ( size_t ij = 0; ij < updatePair.fwdProbs_[0][ii].size(); ij++){
                ofstreamExportFwdProb << (ii+1) << "X" << (ij+1);
                ofstreamExportFwdProb << ((((ii+1) * (ij+1)) < (updatePair.fwdProbs_[0].size()*updatePair.fwdProbs_[0][ii].size()))  ? "\t" : "\n") ;
            }
        }
    }

    size_t siteIndex = 0;
    for ( size_t posI = 0; posI < position_[chromIndex].size(); posI++){
        ofstreamExportFwdProb << chrom_[chromIndex] << "\t" << (int)position_[chromIndex][posI] << "\t";
        for ( size_t ii = 0; ii < updatePair.fwdProbs_[siteIndex].size(); ii++){
            for ( size_t ij = 0; ij < updatePair.fwdProbs_[siteIndex][ii].size(); ij++){
                ofstreamExportFwdProb << updatePair.fwdProbs_[siteIndex][ii][ij];
                ofstreamExportFwdProb << ((((ii+1) * (ij+1)) < (updatePair.fwdProbs_[0].size()*updatePair.fwdProbs_[0][ii].size()))  ? "\t" : "\n") ;
            }
        }
        siteIndex++;
    }

    ofstreamExportFwdProb.close();
}

