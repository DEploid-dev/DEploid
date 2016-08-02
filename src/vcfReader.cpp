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

#include "global.h"
#include <cassert>       // assert
#include <iostream>      // std::cout
#include "vcfReader.hpp"

using namespace std;

/*! Initialize vcf file, search for the end of the vcf header.
 *  Extract the first block of data ( "buffer_length" lines ) into buff
 */
VcfReader::VcfReader(string fileName){

    /*! Initialize by read in the vcf header file */
    this->init( fileName );
    this->readHeader();
    this->readVariants();
    this->finalize();
}


void VcfReader::init( string fileName ){
    /*! Initialize other VcfReader class members
     */
    this->fileName_ = fileName;
    this->inFile.open( this->fileName_.c_str());
}


void VcfReader::finalize(){
    for ( size_t i = 0; i < this->variants.size(); i++ ){
        this->refCount.push_back( (double)this->variants[i].ref );
        this->altCount.push_back( (double)this->variants[i].alt );
    }
    this->inFile.close();
}


void VcfReader::readHeader(){
    if ( inFile.good() ){
        getline ( inFile, this->tmpLine_ );
        while ( this->tmpLine_.size()>0 ){
            if ( this->tmpLine_[0]=='#' ){
                if ( this->tmpLine_[1]=='#' ){
                    this->headerLines.push_back( this->tmpLine_ );
                    getline ( inFile, this->tmpLine_ );
                } else {
                    this->checkFeilds();
                    break; //end of the header
                }
            }
        }
    }
    else {
        throw InvalidInputFile( this->fileName_ );
    }

    dout << " there are "<< this->headerLines.size() << " lines" <<endl;
}


void VcfReader::checkFeilds(){
    size_t feild_start = 0;
    size_t field_end = 0;
    size_t field_index = 0;

    while( field_end < this->tmpLine_.size() ) {
        field_end = min( this->tmpLine_.find('\t',feild_start), this->tmpLine_.find('\n',feild_start) );
        this->tmpStr_ = this->tmpLine_.substr( feild_start, field_end-feild_start );
        string correctFieldValue;
        switch ( field_index ){
            case 0: correctFieldValue = "#CHROM";  break;
            case 1: correctFieldValue =  "POS";    break;
            case 2: correctFieldValue =  "ID";     break;
            case 3: correctFieldValue =  "REF";    break;
            case 4: correctFieldValue =  "ALT";    break;
            case 5: correctFieldValue =  "QUAL";   break;
            case 6: correctFieldValue =  "FILTER"; break;
            case 7: correctFieldValue =  "INFO";   break;
            case 8: correctFieldValue =  "FORMAT"; break;
            case 9: sampleName = this->tmpStr_;    break;
        }

        if ( this->tmpStr_ != correctFieldValue && field_index < 9){
            throw VcfInvalidFieldHeader( correctFieldValue, this->tmpStr_ );
        }

        if ( field_index == 9 ){
            break;
        }

        feild_start = field_end+1;
        field_index++;
    } // End of while loop: field_end < line.size()

    assert( field_index == 9 );
    assert( printSampleName() );
}


bool VcfReader::printSampleName(){
    dout << "Sample name is " << this->sampleName << endl;
    return true;
}


void VcfReader::readVariants(){
    getline( inFile, this->tmpLine_ );
    while ( inFile.good() && this->tmpLine_.size()>0 ){
        VariantLine newVariant ( this->tmpLine_ );
        this->variants.push_back(newVariant);
        getline(inFile, this->tmpLine_);
    }
}


VariantLine::VariantLine ( string tmpLine ){
    this->init( tmpLine );

    while ( fieldEnd_ < this->tmpLine_.size() ){
        fieldEnd_ = min ( this->tmpLine_.find('\t',feildStart_), this->tmpLine_.find('\n', feildStart_) );
        this->tmpStr_ = this->tmpLine_.substr( feildStart_, fieldEnd_ - feildStart_ );
        switch ( fieldIndex_ ){
            case 0: this->extract_field_CHROM();   break;
            case 1: this->extract_field_POS ();    break;
            case 2: this->extract_field_ID ();     break;
            case 3: this->extract_field_REF () ;   break;
            case 4: this->extract_field_ALT () ;   break;
            case 5: this->extract_field_QUAL() ;   break;
            case 6: this->extract_field_FILTER();  break;
            case 7: this->extract_field_INFO();    break;
            case 8: this->extract_field_FORMAT();  break;
            case 9: this->extract_field_VARIANT(); break;
        }

        feildStart_ = fieldEnd_+1;
        fieldIndex_++;
    }
}


void VariantLine::init( string tmpLine ){
    this->tmpLine_ = tmpLine;
    this->feildStart_ = 0;
    this->fieldEnd_ = 0;
    this->fieldIndex_  = 0;
    this->adFieldIndex_ = -1;
}


void VariantLine::extract_field_CHROM () {
    this->chromStr = this->tmpStr_;
}


void VariantLine::extract_field_POS ( ){
    this->posStr = this->tmpStr_;
}


void VariantLine::extract_field_ID ( ){
    this->idStr = this->tmpStr_;
}


void VariantLine::extract_field_REF ( ){
    this->refStr = this->tmpStr_;
}


void VariantLine::extract_field_ALT ( ){
    this->altStr = this->tmpStr_;
}


void VariantLine::extract_field_QUAL ( ){
    this->qualStr = this->tmpStr_;
}


void VariantLine::extract_field_FILTER ( ){
    this->filterStr = this->tmpStr_;
}

void VariantLine::extract_field_INFO ( ){
    this->infoStr = this->tmpStr_;
}


void VariantLine::extract_field_FORMAT ( ){
    this->formatStr = this->tmpStr_;
    size_t feild_start = 0;
    size_t field_end = 0;
    size_t field_index = 0;

    while( field_end < this->formatStr.size() ) {
        field_end = min( this->formatStr.find(':',feild_start), this->formatStr.find('\n',feild_start) );
        if ( "AD" == this->formatStr.substr( feild_start, field_end-feild_start ) ){
            break;
        }
        feild_start = field_end+1;
        field_index++;
    }
    adFieldIndex_ = field_index;
    assert ( adFieldIndex_ > -1 );
    //cout << formatStr << "  " << adFieldIndex_ << endl;
}


void VariantLine::extract_field_VARIANT ( ){
    size_t feild_start = 0;
    size_t field_end = 0;
    int field_index = 0;

    while( field_end < this->tmpStr_.size() ) {
        field_end = min( this->tmpStr_.find(':',feild_start), this->tmpStr_.find('\n',feild_start) );
        if ( field_index == adFieldIndex_ ){
            string adStr = this->tmpStr_.substr( feild_start, field_end-feild_start );
            size_t commaIndex = adStr.find(',',0);
            ref = stoi(adStr.substr(0, commaIndex) );
            alt = stoi(adStr.substr(commaIndex+1, adStr.size()) );
            break;
        }
        feild_start = field_end+1;
        field_index++;
    }

}
