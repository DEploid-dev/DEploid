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

#include <string>  /* string */
#include <vector>  /* vector */
#include <fstream>
#include <stdlib.h>     /* strtol, strtod */
#include "exceptions.hpp"


#ifndef VCF
#define VCF

using namespace std;

struct InvalidVcf : public InvalidInput{
    InvalidVcf( string str ):InvalidInput( str ){
    }
    virtual ~InvalidVcf() throw() {}
    virtual const char* what () const noexcept {
    return throwMsg.c_str();
    }
};


struct VcfInvalidFieldHeader : public InvalidVcf{
    VcfInvalidFieldHeader( string str1, string str2 ):InvalidVcf( str1 ){
        this->reason = " VCF field header should be ";
        throwMsg = this->reason + this->src + ", " + str2 + " was found!" ;
    }
    ~VcfInvalidFieldHeader() throw() {}
};

// More informative exceptions for vcf related errors


class VariantLine{
 friend class VcfReader;
  public:
    VariantLine ( string tmpLine );
    ~VariantLine(){}


  private:
    string tmpLine_;
    string tmpStr_;

    void init( string tmpLine );

    void extract_field_CHROM   ( );
    void extract_field_POS     ( );
    void extract_field_ID      ( );
    void extract_field_REF     ( );
    void extract_field_ALT     ( );
    void extract_field_QUAL    ( );
    void extract_field_FILTER  ( );
    void extract_field_INFO    ( );
    void extract_field_FORMAT  ( );
    void extract_field_VARIANT ( );

    size_t feildStart_;
    size_t fieldEnd_;
    size_t fieldIndex_;

    string chromStr;
    string posStr;
    string idStr;
    string refStr;
    string altStr;
    string qualStr;
    string filterStr;
    string infoStr;
    string formatStr;

    int ref;
    int alt;

};



/*! \brief VCF file reader @ingroup group_data */
class VcfReader{

#ifdef UNITTEST
 friend class TestVCF;
#endif

  public:
    // Constructors and Destructors
    VcfReader(string fileName);
    ~VcfReader(){};

    void init( string fileName );
    ifstream inFile;

    void finalize();

    // Methods
    //void read_new_line();
    //void reset_data_to_first_entry();
    //void force_to_end_data() { this->end_data_ = true; }

    //// DEBUG
    //void print_vcf_line();
    //void print();

    bool end_data() const { return this->end_data_; }
    //INPUT_FILETYPE FileType;

  private:

    vector <VariantLine> variants;

    vector <string> headerLines ;
    string sampleName;
    // Methods
    //void empty_block();
    void readVariants();

    //  Setters and getters:
    //size_t nfield() const { return this->nfield_; }

    //bool eof() const { return this->eof_; }
    //int even_interval() const { return this-> even_interval_; }
    //void set_even_interval( int interval ) { this->even_interval_ = interval; }


    string tmpLine_;
    string tmpStr_;
    void readHeader( );
    void checkFeilds( );

    //void set_empty_line_entry();
    void initialize_read_newLine();
    void check_and_update_block();
    void check_and_update_newLine();
    void finalize_read_new_line();


    string extract_field_ALT_str( size_t start, size_t end );




    // Members

    // Line related


    // FILE related
    bool eof_;
    string fileName_;
    size_t current_line_index_; // line counter in the entire vcf file
    size_t file_length_;
    bool end_data_;
    size_t end_pos_;
    int buffer_max_number_of_lines;

    // Header related
    size_t header_end_pos_;
    //size_t nfield_;
    size_t header_end_line;

    // Block related
    size_t current_block_line_;  // line counter in the current data block
    size_t empty_file_line_counter_;

    size_t ghost_num_mut;
    int filter_window_;  // If two snps are too close, i.e. difference between the site is less than this window, should skip to the next read.
    int missing_data_threshold_; // if two snps are too far away apart, treat as missing data
    int even_interval_;

    //size_t vcf_file_length;

    vector <string> buffer_lines;

    bool skip_tmpLine;


    // Debug tools
    bool printSampleName();


};

#endif
