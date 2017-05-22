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

#include <string>  /* string */
#include <vector>  /* vector */
#include <fstream>
#include <stdlib.h>     /* strtol, strtod */
#include "exceptions.hpp"
#include "variantIndex.hpp"
#include "gzstream.h"

#ifndef VCF
#define VCF

using namespace std;

struct InvalidVcf : public InvalidInput{
    InvalidVcf( string str ):InvalidInput( str ){
    }
    virtual ~InvalidVcf() throw() {}
    //virtual const char* what () const noexcept {
        //return throwMsg.c_str();
    //}
};


struct VcfInvalidHeaderFieldNames : public InvalidVcf{
    VcfInvalidHeaderFieldNames( string str1, string str2 ):InvalidVcf( str1 ){
        this->reason = " VCF field header expects: ";
        throwMsg = this->reason + this->src + ", " + str2 + " was found!" ;
    }
    ~VcfInvalidHeaderFieldNames() throw() {}
};


struct VcfInvalidVariantEntry : public InvalidVcf{
    VcfInvalidVariantEntry( string str ):InvalidVcf( str ){}
    virtual ~VcfInvalidVariantEntry() throw() {}
    //virtual const char* what () const noexcept {
        //return throwMsg.c_str();
    //}
};


struct VcfCoverageFieldNotFound : public VcfInvalidVariantEntry{
    VcfCoverageFieldNotFound( string str ):VcfInvalidVariantEntry( str ){
        this->reason = "Coverage field AD was not found in the FORMAT, found: ";
        throwMsg = this->reason + this->src ;
    }
    ~VcfCoverageFieldNotFound() throw() {}
};

// This requires more thinking, check for data type?


// More informative exceptions for vcf related errors


class VariantLine{
 friend class VcfReader;
 friend class DEploidIO;
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
    int adFieldIndex_;

    int ref;
    int alt;
};



/*! \brief VCF file reader @ingroup group_data */
class VcfReader : public VariantIndex {

#ifdef UNITTEST
 friend class TestVCF;
#endif
 friend class DEploidIO;
  public:
    // Constructors and Destructors
    VcfReader(string fileName); // parse in exclude sites
    ~VcfReader(){};

  private:
    vector <VariantLine> variants;
    vector <VariantLine> keptVariants;
    vector <double> refCount;
    vector <double> altCount;
    vector <string> headerLines ;
    string fileName_;
    ifstream inFile;
    igzstream inFileGz;
    bool isCompressed_;
    bool isCompressed() const { return this->isCompressed_; }
    void setIsCompressed ( const bool compressed ){ this->isCompressed_ = compressed; }
    void checkFileCompressed();
    string sampleName;
    string tmpLine_;
    string tmpStr_;

    // Methods
    void init( string fileName );
    void finalize();
    void readVariants();
    void readHeader( );
    void checkFeilds( );

    void getChromList();
    void removeMarkers ( );

    // Debug tools
    bool printSampleName();
};

#endif
