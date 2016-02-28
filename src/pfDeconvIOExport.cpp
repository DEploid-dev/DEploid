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
#include "mcmc.hpp"

void PfDeconvIO::write( McmcSample * mcmcSample){
    this->writeProp( mcmcSample );
    this->writeLLK( mcmcSample );
    this->writeHap( mcmcSample );

    this->writeLog ( &std::cout );
    ofstreamExportLog.open( strExportLog.c_str(), ios::out | ios::app | ios::binary );
    this->writeLog ( &ofstreamExportLog );
    ofstreamExportLog.close();
}


void PfDeconvIO::writeLog ( ostream * writeTo ){
    (*writeTo) << "#########################################\n";
    (*writeTo) << "#        pfDeconv "<< setw(10) << VERSION << " log        #\n";
    (*writeTo) << "#########################################\n";
    (*writeTo) << "Program was compiled on: " << compileTime_ << endl;
    (*writeTo) << "pfDeconv version: " << pfDeconvVersion_ << endl;
    (*writeTo) << "input data: \n";
    (*writeTo) << "PLAF: " << plafFileName_ << "\n";
    (*writeTo) << "REF count: " << refFileName_  << "\n";
    (*writeTo) << "ALT count: " << altFileName_  << "\n";
    (*writeTo) << "\n";
    (*writeTo) << "Mcmc chain parameters:"<< "\n";
    (*writeTo) << "MCMC sample: " << nMcmcSample_ << ", sample rate: " << mcmcMachineryRate_ <<"\n";
    if ( seed_set_ ) { (*writeTo) << "Random seed: "<< random_seed_ << "\n";}
    (*writeTo) << "\n";
    (*writeTo) << "\n";
}


void PfDeconvIO::writeProp( McmcSample * mcmcSample){
    ofstreamExportProp.open( strExportProp.c_str(), ios::out | ios::app | ios::binary );
    for ( size_t i = 0; i < mcmcSample->proportion.size(); i++){
        for ( size_t ii = 0; ii < mcmcSample->proportion[i].size(); ii++){
            ofstreamExportProp << setw(10) << mcmcSample->proportion[i][ii];
            ofstreamExportProp << ((ii < (mcmcSample->proportion[i].size()-1)) ? "\t" : "\n") ;
        }
    }
    ofstreamExportProp.close();
}


void PfDeconvIO::writeLLK( McmcSample * mcmcSample){
    ofstreamExportLLK.open( strExportLLK.c_str(), ios::out | ios::app | ios::binary );
    for ( size_t i = 0; i < mcmcSample->sumLLKs.size(); i++){
        //if ( this->moves[i] != 2)
        ofstreamExportLLK << mcmcSample->moves[i] << "\t" << mcmcSample->sumLLKs[i] << endl;
    }
    ofstreamExportLLK.close();
}


void PfDeconvIO::writeHap( McmcSample * mcmcSample){
    ofstreamExportHap.open( strExportHap.c_str(), ios::out | ios::app | ios::binary );
    for ( size_t i = 0; i < mcmcSample->hap.size(); i++ ){
        for ( size_t ii = 0; ii < mcmcSample->hap[i].size(); ii++){
            ofstreamExportHap << mcmcSample->hap[i][ii];
            ofstreamExportHap << ((ii < (mcmcSample->hap[i].size()-1)) ? "\t" : "\n") ;
        }
    }
    ofstreamExportHap.close();
}
