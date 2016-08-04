/*
 * dEploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016, Sha (Joe) Zhu, Jacob Almagro and Prof. Gil McVean
 *
 * This file is part of dEploid.
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

#include "dEploidIO.hpp"
#include "mcmc.hpp"

void DEploidIO::write( McmcSample * mcmcSample, Panel * panel ){
    this->writeProp( mcmcSample );
    this->writeLLK( mcmcSample );
    this->writeHap( mcmcSample );
    this->writeVcf( mcmcSample );
    this->writeRecombProb( panel );

    // Get End time before writing the log
    this->getTime(false);

    this->writeLog ( mcmcSample, &std::cout );

    ofstreamExportLog.open( strExportLog.c_str(), ios::out | ios::app | ios::binary );
    this->writeLog ( mcmcSample, &ofstreamExportLog );
    ofstreamExportLog.close();
}


void DEploidIO::writeRecombProb ( Panel * panel ){
    if ( !doExportRecombProb() ) return;

    if ( panel != NULL ){
        ofstreamExportRecombProb.open( strExportRecombProb.c_str(), ios::out | ios::app | ios::binary );
        ofstreamExportRecombProb << "p.recomb"       << "\t"
                                 << "p.each"         << "\t"
                                 << "p.no.recomb"    << "\t"
                                 << "p.rec.rec"      << "\t"
                                 << "p.rec.norec"    << "\t"
                                 << "p.norec.norec"  << "\n";
        for ( size_t i = 0; i < panel->pRec_.size(); i++ ){
            ofstreamExportRecombProb << panel->pRec_[i]           << "\t"
                                     << panel->pRecEachHap_[i]    << "\t"
                                     << panel->pNoRec_[i]         << "\t"
                                     << panel->pRecRec_[i]        << "\t"
                                     << panel->pRecNoRec_[i]      << "\t"
                                     << panel->pNoRecNoRec_[i]    << "\n";
        }
        ofstreamExportProp.close();
    }
}



void DEploidIO::writeLog ( McmcSample * mcmcSample, ostream * writeTo ){
    (*writeTo) << "#########################################\n";
    (*writeTo) << "#        dEploid "<< setw(10) << VERSION << " log        #\n";
    (*writeTo) << "#########################################\n";
    (*writeTo) << "Program was compiled on: " << compileTime_ << endl;
    (*writeTo) << "dEploid version: " << dEploidVersion_ << endl;
    (*writeTo) << "\n";
    (*writeTo) << "Input data: \n";
    (*writeTo) << setw(12) << "Panel: "     << panelFileName_  << "\n";
    (*writeTo) << setw(12) << "PLAF: "      << plafFileName_   << "\n";
    if ( useVcf() ) (*writeTo) << setw(12) << "VCF: " << vcfFileName_    << "\n";
    if ( refFileName_.size()>0) (*writeTo) << setw(12) << "REF count: " << refFileName_    << "\n";
    if ( altFileName_.size()>0) (*writeTo) << setw(12) << "ALT count: " << altFileName_    << "\n";
    if ( excludeSites() ){ (*writeTo) << setw(12) << "Exclude: " << excludeFileName_    << "\n"; }
    (*writeTo) << "\n";
    (*writeTo) << "MCMC parameters: "<< "\n";
    (*writeTo) << setw(19) << " MCMC burn: " << mcmcBurn_ << "\n";
    (*writeTo) << setw(19) << " MCMC sample: " << nMcmcSample_ << "\n";
    (*writeTo) << setw(19) << " MCMC sample rate: " << mcmcMachineryRate_ <<"\n";
    (*writeTo) << setw(19) << " Random seed: " << this->randomSeed() << "\n";
    (*writeTo) << setw(19) << " Update Prop: "   << (this->doUpdateProp()  ? "YES":"NO") << "\n";
    (*writeTo) << setw(19) << " Update Single: " << (this->doUpdateSingle()? "YES":"NO") << "\n";
    (*writeTo) << setw(19) << " Update Pair: "   << (this->doUpdatePair()  ? "YES":"NO") << "\n";
    (*writeTo) << "\n";
    (*writeTo) << "Other parameters:"<< "\n";
    if ( forbidCopyFromSame_ ){ (*writeTo) << " Update pair haplotypes move forbid copying from the same strain!!! \n"; }
    (*writeTo) << setw(20) << " Miss copy prob: "   << this->missCopyProb_ << "\n";
    (*writeTo) << setw(20) << " Avrg Cent Morgan: " << this->averageCentimorganDistance_ << "\n";
    (*writeTo) << setw(20) << " Ne: "               << this->Ne_ << "\n";
    if ( this->initialPropWasGiven() ){
        (*writeTo) << setw(20) << " Initial prob: " ;
        for ( size_t i = 0; i < this->initialProp.size(); i++ ){
            (*writeTo) << this->initialProp[i]
                       << ( ( i != (this->kStrain_-1) ) ? " " : "\n" );
        }
    }
    (*writeTo) << "\n";
    (*writeTo) << "Run time:\n";
    (*writeTo) << setw(14) << "Start at: "  << startingTime_  ;
    (*writeTo) << setw(14) << "End at: "    << endTime_  ;
    (*writeTo) << "\n";
    (*writeTo) << "Output saved to:\n";
    (*writeTo) << setw(14) << "Likelihood: "  << strExportLLK  << "\n";
    (*writeTo) << setw(14) << "Proportions: " << strExportProp << "\n";
    (*writeTo) << setw(14) << "Haplotypes: "  << strExportHap  << "\n";
    if ( doExportVcf() ) { (*writeTo) << setw(14) << "Vcf: "  << strExportVcf  << "\n"; }
    (*writeTo) << "\n";
    (*writeTo) << "Proportions at the last iteration:\n";
    for ( size_t ii = 0; ii < mcmcSample->proportion.back().size(); ii++){
        (*writeTo) << setw(10) << mcmcSample->proportion.back()[ii];
        (*writeTo) << ((ii < (mcmcSample->proportion.back().size()-1)) ? "\t" : "\n") ;
    }

}


void DEploidIO::writeProp( McmcSample * mcmcSample){
    ofstreamExportProp.open( strExportProp.c_str(), ios::out | ios::app | ios::binary );
    for ( size_t i = 0; i < mcmcSample->proportion.size(); i++){
        for ( size_t ii = 0; ii < mcmcSample->proportion[i].size(); ii++){
            ofstreamExportProp << setw(10) << mcmcSample->proportion[i][ii];
            ofstreamExportProp << ((ii < (mcmcSample->proportion[i].size()-1)) ? "\t" : "\n") ;
        }
    }
    ofstreamExportProp.close();
}


void DEploidIO::writeLLK( McmcSample * mcmcSample){
    ofstreamExportLLK.open( strExportLLK.c_str(), ios::out | ios::app | ios::binary );
    for ( size_t i = 0; i < mcmcSample->sumLLKs.size(); i++){
        ofstreamExportLLK << mcmcSample->moves[i] << "\t" << mcmcSample->sumLLKs[i] << endl;
    }
    ofstreamExportLLK.close();
}


void DEploidIO::writeHap( McmcSample * mcmcSample ){
    ofstreamExportHap.open( strExportHap.c_str(), ios::out | ios::app | ios::binary );

    // HEADER
    ofstreamExportHap << "CHROM" << "\t" << "POS" << "\t";;
    for ( size_t ii = 0; ii < kStrain_; ii++){
        ofstreamExportHap << "h" << (ii+1) ;
        ofstreamExportHap << ((ii < (kStrain_-1)) ? "\t" : "\n") ;
    }

    size_t siteIndex = 0;
    for ( size_t chromI = 0; chromI < chrom_.size(); chromI++ ){
        for ( size_t posI = 0; posI < position_[chromI].size(); posI++){
            ofstreamExportHap << chrom_[chromI] << "\t" << (int)position_[chromI][posI] << "\t";
            for ( size_t ii = 0; ii < mcmcSample->hap[siteIndex].size(); ii++){
                ofstreamExportHap << mcmcSample->hap[siteIndex][ii];
                ofstreamExportHap << ((ii < (mcmcSample->hap[siteIndex].size()-1)) ? "\t" : "\n") ;
            }
            siteIndex++;
        }
    }

    assert ( siteIndex == mcmcSample->hap.size());
    ofstreamExportHap.close();
}


void DEploidIO::writeVcf( McmcSample * mcmcSample ){
    if ( !doExportVcf() ) return;

    ofstreamExportTmp.open( strExportVcf.c_str(), ios::out | ios::app | ios::binary );
    // VCF HEADER
    if ( this->useVcf() ){
        for ( auto const& headerLine: this->vcfReaderPtr_->headerLines){
            ofstreamExportTmp << headerLine << endl;
        }
    } else {
        ofstreamExportTmp << "##fileformat=VCFv4.2" << endl;
    }
    // Include proportions
    for ( size_t ii = 0; ii < kStrain_; ii++){
        ofstreamExportTmp << "##Proportion of strain "
                          << ( this->useVcf() ? this->vcfReaderPtr_->sampleName : "h" )
                          << "." << (ii+1)
                          << "=" << mcmcSample->proportion.back()[ii] << endl;
    }

    // HEADER
    ofstreamExportTmp << "#CHROM" << "\t"
                      << "POS"    << "\t"
                      << "ID"     << "\t"
                      << "REF"    << "\t"
                      << "ALT"    << "\t"
                      << "QUAL"   << "\t"
                      << "FILTER" << "\t"
                      << "INFO"   << "\t"
                      << "FORMAT" << "\t";
    for ( size_t ii = 0; ii < kStrain_; ii++){
        ofstreamExportTmp << ( this->useVcf() ? this->vcfReaderPtr_->sampleName : "h" )
                          << "." << (ii+1) ;
        ofstreamExportTmp << ((ii < (kStrain_-1)) ? "\t" : "\n") ;
    }

    size_t siteIndex = 0;
    for ( size_t chromI = 0; chromI < chrom_.size(); chromI++ ){
        for ( size_t posI = 0; posI < position_[chromI].size(); posI++){
            if ( useVcf() ) {
                ofstreamExportTmp << this->vcfReaderPtr_->variants[siteIndex].chromStr  << "\t"
                                  << this->vcfReaderPtr_->variants[siteIndex].posStr    << "\t"
                                  << this->vcfReaderPtr_->variants[siteIndex].idStr     << "\t"
                                  << this->vcfReaderPtr_->variants[siteIndex].refStr    << "\t"
                                  << this->vcfReaderPtr_->variants[siteIndex].altStr    << "\t"
                                  << this->vcfReaderPtr_->variants[siteIndex].qualStr   << "\t"
                                  << this->vcfReaderPtr_->variants[siteIndex].filterStr << "\t"
                                  << this->vcfReaderPtr_->variants[siteIndex].infoStr   << "\t"
                                  << "GT"                                      << "\t";
            } else {
                ofstreamExportTmp << chrom_[chromI]               << "\t"
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
                ofstreamExportTmp << mcmcSample->hap[siteIndex][ii];
                ofstreamExportTmp << ((ii < (mcmcSample->hap[siteIndex].size()-1)) ? "\t" : "\n") ;
            }
            siteIndex++;
        }
    }

    assert ( siteIndex == mcmcSample->hap.size());
    ofstreamExportTmp.close();
}
