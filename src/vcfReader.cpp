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
#include <iostream> // std::cout
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
    //this->reset_chrom_site(); // reinitilize VariantPosition
    //this->filter_window_             = 1;
    //this->current_line_index_        = 0;
    //this->empty_file_line_counter_   = 0;
    //this->nsam_                      = 0;
    ////this->nfield_                    = 0;
    ////this->current_block_line_        = 0;
    //this->missing_data_threshold_ = INT_MAX;
    //this->file_length_                = 0;
    //this->even_interval_             = 0.0;
    //this->ghost_num_mut            = 0;
    //this->eof_                       = false;
    //this->end_data_                  = false;
    //this->empty_block();
}


void VcfReader::finalize(){
    this->inFile.close();

}


void VcfReader::readHeader(){
    //ifstream in_file;
    //in_file.open( this->fileName_.c_str());
    //in_file.seekg (0, in_file.end);
    //this->file_length_ = in_file.tellg();

    //in_file.seekg (0, in_file.beg);
    //header_end_pos_=0;
    //header_end_line=0;
    //in_file.seekg (0, in_file.end);
    //in_file.seekg (0, in_file.beg);

    if ( inFile.good() ){
        getline ( inFile, this->tmpLine_ );
        //header_end_pos_ += this->tmpLine_.size() + 1;
        while ( this->tmpLine_.size()>0 ){
            if ( this->tmpLine_[0]=='#' ){
                if ( this->tmpLine_[1]=='#' ){
                    this->headerLines.push_back( this->tmpLine_ );
                    getline ( inFile, this->tmpLine_ );

                    //header_end_pos_ += this->tmpLine_.size() + 1;
                    //header_end_line++;

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

    //this->end_pos_ = this->header_end_pos_;
    cout << " there are "<< this->headerLines.size() << " lines" <<endl;
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
    //int line_counter = 0;
    getline( inFile, this->tmpLine_ );
    while ( inFile.good() && this->tmpLine_.size()>0 ){
        VariantLine newVariant ( this->tmpLine_ );
        getline(inFile, this->tmpLine_);
    }
}




VariantLine::VariantLine ( string tmpLine ){
    cout << tmpLine << endl;

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
}


void VariantLine::extract_field_CHROM () {

    //this->chrom_ = strtol( tmpStr.c_str(), NULL, 0);
    //if ( pervious_chrom_ != chrom_ && pervious_chrom_!= 0 ) { throw std::invalid_argument ( "Two different chroms" ); }
    }


void VariantLine::extract_field_POS ( ){

    //this->site_ = strtol( tmpStr.c_str(), NULL, 0);

    //if (  this->pervious_chrom_ > 0) {assert ( this->pervious_chrom_ == this->chrom_ );}
    //assert ( this->previous_site_at_ >= 0);

    //if ( ( this->site_ - this->previous_site_at_ ) < this->filter_window_   ) { // Making 100, as a filtering process, to screen mutations that are too close
        //cout << "Skip reads at chrom " << this->chrom_<<" at position " <<  this->site_ << ", due to it's too close to the previous variant (at " << this->previous_site_at_ <<")." << endl;
        //this->skip_tmpLine = true;
        //}
    }


void VariantLine::extract_field_ID ( ){}


void VariantLine::extract_field_REF ( ){
    /*! Check the REF field, if there REF length is greater than 1, it is indel.
     *  if it is gvcf file, and REF is "." or "N", it is INVARIANT,
     */

    //this->ref = tmpStr;

    //if ( this->ref.size() > 1 ){ // BAD LINE, SKIP EXACTING INFORMATION FOR FIELD
        //cout << "Skip reads at chrom " << this->chrom_<<" at position " <<  this->site_<<", due to deletion or replacement" << endl;
        ////this->current_variant_state = OTHER_VARIANT; // IGNORE OTHER_VARIANT FOR NOW, JUST SKIP
        //this->skip_tmpLine = true;
        //return ;
    //}

    //if ( this->ref == "." || this->ref == "N" ) {
        //if ( this->FileType == GVCF  ) {
             //this->current_variant_state = INVARIANT;
             //this->current_seg_state = SEQ_INVARIANT;
             //}
        //else if ( this->FileType == RGVCF ) {
             //this->current_variant_state = INVARIANT;
             //this->current_seg_state = MISSING;
             //}
        //}
    //else { // assume besides ".", "N", it can only be "A","T", "C", "G"
        //this->current_variant_state = SNP;
        //this->current_seg_state = ZERO_SEG;
        //}
}


void VariantLine::extract_field_ALT ( ){
    //if ( this->tmpStr == "." || this->tmpStr == "N" ) {
        //if ( this->FileType == GVCF  ) {
             //this->current_variant_state = INVARIANT;
             //this->current_seg_state = SEQ_INVARIANT;
             //}
        //if ( this->FileType == RGVCF ) {
             //this->current_variant_state = INVARIANT;
             //this->current_seg_state = MISSING;
             //}
        //return;
    //}
    //size_t alt_start = 0;
    //size_t alt_end = 0;
    //string alt_tmpStr;
    //bool alt_not_valid = false; // Depends on the number of samples, need to check each ALT string for each sample
    //while ( alt_end < this->tmpStr.size()){
        //alt_end = min( this->tmpStr.find(',', alt_start), this->tmpStr.size() );
        //alt_tmpStr = this->tmpStr.substr(alt_start, alt_end);
        //if ( alt_tmpStr.size() > 1 ){ // BAD LINE, SKIP EXACTING INFORMATION FOR FIELD
            //cout << "Skip reads at chrom " << this->chrom_<<" at position " <<  this->site_ <<", due to insertion" << endl;
            //this->current_variant_state = OTHER_VARIANT;
            //alt_not_valid = true;
            //break;
            //}
        //alt.push_back(alt_tmpStr);
        //alt_start=alt_end+1;
        //}
    //if ( alt_not_valid ){
        //this->skip_tmpLine = true;
        //}
}

void VariantLine::extract_field_QUAL ( ){}

void VariantLine::extract_field_FILTER ( ){
    //this->skip_tmpLine = false;
    //if ( this->tmpStr.find( "PASS") == std::string::npos && this->tmpStr.find( "REFCALL") == std::string::npos ){ // BAD LINE, SKIP EXACTING INFORMATION FOR FIELD
        //cout << "Skip reads at chrom " << chrom_<<" at position " <<  site_<<", due to low qualitiy." << endl;
        ////cout << "skipping: "<< line << endl; // DEBUG
        //this->skip_tmpLine = true;
        //}
}

void VariantLine::extract_field_INFO ( ){
    //if ( this->current_variant_state == SNP ){
        //assert( this->current_seg_state == ZERO_SEG );
        //this->seg_end_site_ = this->site_;
        //}
    //else {
        //assert( this->tmpStr.find( "END=", 0 ) != std::string::npos );
        //this->seg_end_site_ = strtol( tmpStr.substr( (size_t)4 ).c_str(), NULL, 0);

        //assert ( (this->FileType == GVCF) || ( this->FileType == RGVCF ) );
        //this->current_seg_state = this->FileType == GVCF ? SEQ_INVARIANT : MISSING ;
        ////if      ( this->FileType == VCF )  this->current_seg_state = SEQ_INVARIANT;
        ////else if ( this->FileType == GVCF )  this->current_seg_state = SEQ_INVARIANT;
        ////else    ( this->FileType == RGVCF )  this->current_seg_state = MISSING;

        ////if ( this->FileType == GVCF ) { this->current_seg_state = SEQ_INVARIANT; }
        ////else if ( this->FileType == RGVCF ) { this->current_seg_state = MISSING; }
        //}
}


void VariantLine::extract_field_FORMAT ( ){ }

void VariantLine::extract_field_VARIANT ( ){}
    ////assert ( this->skip_tmpLine = false );
    //if ( this->current_variant_state == INVARIANT ){
        //this->int_vec_of_sample_alt.push_back( 0 );
        //this->int_vec_of_sample_alt.push_back( 0 );
        //return;
    //}

    //size_t bar_index   = this->tmpStr.find('|',0);
    //size_t slash_index = this->tmpStr.find('/',0);
    //size_t colon_index = this->tmpStr.find(':',0);
    //size_t break_index = min(bar_index, slash_index);
    //assert( break_index < colon_index );
    //this->vec_of_sample_alt.push_back( extract_field_ALT_str( 0, break_index));
    //this->vec_of_sample_alt.push_back( extract_field_ALT_str( break_index+1, colon_index));

    //size_t alt_index_0 = strtol (tmpStr.substr(0,1).c_str(), NULL, 0);;
    //this->int_vec_of_sample_alt.push_back( (alt_index_0 == (size_t)0) ? 0 : 1);

    //size_t alt_index_2 = strtol (tmpStr.substr(2,1).c_str(), NULL, 0);;
    //this->int_vec_of_sample_alt.push_back( (alt_index_2 == (size_t)0) ? 0 : 1);
    //}


//string VariantLine::extract_field_ALT_str( size_t start, size_t end ){
    ///*! Extract haplotype */
    ////size_t alt_index = strtol ( this->tmpStr.substr(start,end-start).c_str(), NULL, 0);
    ////return ( alt_index==0 ) ? ref : alt[alt_index-1];
//}



