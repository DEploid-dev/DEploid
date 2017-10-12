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
#include "updateHap.hpp"
#include "mcmc.hpp"

void McmcMachinery::writeLastFwdProb(bool useIBD){
    if ( this->dEploidIO_ ->doExportPostProb() != true ){
        return;
    }

    for ( size_t tmpk = 0; tmpk < this->kStrain_; tmpk++ ){
        if ( this->dEploidIO_->doAllowInbreeding() == true ){
            this->updateReferencePanel(this->panel_->truePanelSize()+kStrain_-1, tmpk);
        }

        for ( size_t chromi = 0 ; chromi < this->dEploidIO_->indexOfChromStarts_.size(); chromi++ ){
            size_t start = this->dEploidIO_->indexOfChromStarts_[chromi];
            size_t length = this->dEploidIO_->position_[chromi].size();

            UpdateSingleHap updatingSingle( this->dEploidIO_->refCount_,
                                      this->dEploidIO_->altCount_,
                                      this->dEploidIO_->plaf_,
                                      this->currentExpectedWsaf_,
                                      this->currentProp_, this->currentHap_, this->hapRg_,
                                      start, length,
                                      this->panel_, this->dEploidIO_->missCopyProb_, this->dEploidIO_->scalingFactor(),
                                      tmpk);
            if ( this->dEploidIO_->doAllowInbreeding() == true ){
                updatingSingle.setPanelSize(this->panel_->inbreedingPanelSize());
            }

            updatingSingle.core ( this->dEploidIO_->refCount_, this->dEploidIO_->altCount_, this->dEploidIO_->plaf_, this->currentExpectedWsaf_, this->currentProp_, this->currentHap_);
            this->dEploidIO_->writeLastSingleFwdProb( updatingSingle.fwdProbs_, chromi, tmpk, useIBD );
        }
        //UpdatePairHap updating( this->dEploidIO_->refCount_,
                                //this->dEploidIO_->altCount_,
                                //this->dEploidIO_->plaf_,
                                //this->currentExpectedWsaf_,
                                //this->currentProp_, this->currentHap_, this->hapRg_,
                                //start, length,
                                //this->panel_, this->dEploidIO_->missCopyProb_, this->dEploidIO_->forbidCopyFromSame(),
                                //(size_t)0,
                                //(size_t)1);
        //updating.core ( this->dEploidIO_->refCount_, this->dEploidIO_->altCount_, this->dEploidIO_->plaf_, this->currentExpectedWsaf_, this->currentProp_, this->currentHap_);
        //this->dEploidIO_->writeLastPairFwdProb( updating, chromi );
    }
}


void DEploidIO::writeLastSingleFwdProb( vector < vector <double> >& probabilities, size_t chromIndex, size_t strainIndex, bool useIBD ){
    if ( probabilities.size() == 0 ){
        return;
    }

    size_t panelSize = probabilities[0].size();

    string strExportFwdProb = (( useIBD == true ) ? strIbdExportSingleFwdProbPrefix :
                                                    strExportSingleFwdProbPrefix) + to_string(strainIndex);
    ofstreamExportFwdProb.open( strExportFwdProb.c_str(), ios::out | ios::app | ios::binary );

    if ( chromIndex == 0 ){ // Print header
        ofstreamExportFwdProb << "CHROM" << "\t" << "POS" << "\t";
        for ( size_t ii = 0; ii < probabilities[0].size(); ii++){
            if (this->doAllowInbreeding() == true){
                if ( ii <= (panelSize - this->kStrain()) ){
                    ofstreamExportFwdProb << "P" << (ii+1) ;
                } else {
                    ofstreamExportFwdProb << "I" << (ii)-(panelSize - this->kStrain()) ;
                }
            } else {
                ofstreamExportFwdProb << (ii+1) ;
            }
            ofstreamExportFwdProb << ((ii < (panelSize-1)) ? "\t" : "\n") ;
        }
    }

    size_t siteIndex = 0;
    for ( size_t posI = 0; posI < position_[chromIndex].size(); posI++){
        ofstreamExportFwdProb << chrom_[chromIndex] << "\t" << (int)position_[chromIndex][posI] << "\t";
        for ( size_t ii = 0; ii < probabilities[siteIndex].size(); ii++){
            ofstreamExportFwdProb << probabilities[siteIndex][ii];
            ofstreamExportFwdProb << ((ii < (probabilities[siteIndex].size()-1)) ? "\t" : "\n") ;
        }
        siteIndex++;
    }

    ofstreamExportFwdProb.close();
}


//void DEploidIO::writeLastPairFwdProb( UpdatePairHap & updatePair, size_t chromIndex ){
    //ofstreamExportFwdProb.open( strExportPairFwdProb.c_str(), ios::out | ios::app | ios::binary );
    //if ( chromIndex == 0 ){ // Print header
        //ofstreamExportFwdProb << "CHROM" << "\t" << "POS" << "\t";;
        //for ( size_t ii = 0; ii < updatePair.fwdProbs_[0].size(); ii++){
            //for ( size_t ij = 0; ij < updatePair.fwdProbs_[0][ii].size(); ij++){
                //ofstreamExportFwdProb << (ii+1) << "X" << (ij+1);
                //ofstreamExportFwdProb << ((((ii+1) * (ij+1)) < (updatePair.fwdProbs_[0].size()*updatePair.fwdProbs_[0][ii].size()))  ? "\t" : "\n") ;
            //}
        //}
    //}

    //size_t siteIndex = 0;
    //for ( size_t posI = 0; posI < position_[chromIndex].size(); posI++){
        //ofstreamExportFwdProb << chrom_[chromIndex] << "\t" << (int)position_[chromIndex][posI] << "\t";
        //for ( size_t ii = 0; ii < updatePair.fwdProbs_[siteIndex].size(); ii++){
            //for ( size_t ij = 0; ij < updatePair.fwdProbs_[siteIndex][ii].size(); ij++){
                //ofstreamExportFwdProb << updatePair.fwdProbs_[siteIndex][ii][ij];
                //ofstreamExportFwdProb << ((((ii+1) * (ij+1)) < (updatePair.fwdProbs_[0].size()*updatePair.fwdProbs_[0][ii].size()))  ? "\t" : "\n") ;
            //}
        //}
        //siteIndex++;
    //}

    //ofstreamExportFwdProb.close();
//}

