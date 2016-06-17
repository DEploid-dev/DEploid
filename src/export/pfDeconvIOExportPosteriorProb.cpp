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


void PfDeconvIO::writeLastSingFwdProb( UpdateSingleHap & updateSingle ){

}


void PfDeconvIO::writeLastPairFwdProb( UpdatePairHap & updatePair ){
    //ofstreamExportHap.open( strExportHap.c_str(), ios::out | ios::app | ios::binary );

    //// HEADER
    //ofstreamExportHap << "CHROM" << "\t" << "POS" << "\t";;
    //for ( size_t ii = 0; ii < kStrain_; ii++){
        //ofstreamExportHap << "h" << (ii+1) ;
        //ofstreamExportHap << ((ii < (kStrain_-1)) ? "\t" : "\n") ;
    //}

    //size_t siteIndex = 0;
    //for ( size_t chromI = 0; chromI < chrom_.size(); chromI++ ){
        //for ( size_t posI = 0; posI < position_[chromI].size(); posI++){
            //ofstreamExportHap << chrom_[chromI] << "\t" << (int)position_[chromI][posI] << "\t";
            //for ( size_t ii = 0; ii < mcmcSample->hap[siteIndex].size(); ii++){
                //ofstreamExportHap << mcmcSample->hap[siteIndex][ii];
                //ofstreamExportHap << ((ii < (mcmcSample->hap[siteIndex].size()-1)) ? "\t" : "\n") ;
            //}
            //siteIndex++;
        //}
    //}

    //assert ( siteIndex == mcmcSample->hap.size());
    //ofstreamExportHap.close();
}

